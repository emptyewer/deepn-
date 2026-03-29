#include "gcworker.h"

#include <Xlsx/Workbook.h>

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QStandardPaths>
#include <QStringList>
#include <QTextStream>
#include <QUuid>

GCWorker::GCWorker(GCStat *gcstat) {
  stat = gcstat;
  QFileInfo fileInfo(stat->input);
  QString basename = fileInfo.completeBaseName();
  mappedOuputDBName =
      fileInfo.absolutePath() + "/../analyzed_files/" + basename + ".sqlite";
  // Unique connection names so multiple workers don't collide
  QString uid = QUuid::createUuid().toString(QUuid::Id128);
  writeDbConn = "writeDB_" + uid;
  fileDbConn = "fileDb_" + uid;
  memDbConn = "memoryDb_" + uid;
  qDebug() << mappedOuputDBName;
}

void GCWorker::run() {
  setupDB();
  readMapOutput();
  createInMemoryDB();
  writeGeneCount();
  emit sig->gc_finished_sig();
  emit finished();
}

void GCWorker::writeGeneCount() {
  qDebug() << "Writing Gene Count";
  QSqlQuery q(mem_db);

  // Totals from the maps table
  q.exec("SELECT COUNT(DISTINCT read) FROM maps");
  int totalReads = q.next() ? q.value(0).toInt() : 0;
  q.exec("SELECT COUNT(*) FROM maps");
  int totalHits = q.next() ? q.value(0).toInt() : 0;
  q.clear();

  // Detach the previously attached fileDb if still attached
  QSqlQuery detachQ(mem_db);
  detachQ.exec("DETACH DATABASE IF EXISTS fileDb");
  detachQ.clear();

  // Attach the output file from the in-memory DB and write directly
  QSqlQuery w(mem_db);
  if (!w.exec(QString("ATTACH DATABASE '%1' AS outFile").arg(mappedOuputDBName))) {
    qDebug() << "Error attaching output database:" << w.lastError().text();
    return;
  }

  w.exec("DROP TABLE IF EXISTS outFile.gene_counts");
  w.exec("DROP TABLE IF EXISTS outFile.summary");

  // Gene counts: aggregate distinct reads per gene with PPM
  w.exec("CREATE TABLE outFile.gene_counts ("
         "gene TEXT PRIMARY KEY, count INTEGER, ppm REAL)");
  w.exec(QString(
      "INSERT INTO outFile.gene_counts (gene, count, ppm) "
      "SELECT gene, COUNT(DISTINCT read), "
      "COUNT(DISTINCT read) * 1000000.0 / %1 "
      "FROM maps GROUP BY gene ORDER BY COUNT(DISTINCT read) DESC")
      .arg(qMax(totalReads, 1)));

  // Summary metadata
  w.exec("CREATE TABLE outFile.summary (key TEXT PRIMARY KEY, value TEXT)");
  QFileInfo fileInfo(stat->input);
  w.exec(QString("INSERT INTO outFile.summary VALUES ('file', '%1')")
      .arg(fileInfo.fileName()));
  w.exec(QString("INSERT INTO outFile.summary VALUES ('total_reads', '%1')")
      .arg(totalReads));
  w.exec(QString("INSERT INTO outFile.summary VALUES ('total_hits', '%1')")
      .arg(totalHits));

  w.exec("DETACH DATABASE outFile");
  w.clear();

  // Clean up memory DB
  mem_db.close();
  QSqlDatabase::removeDatabase(memDbConn);

  qDebug() << "Gene counts written to:" << mappedOuputDBName;
}

void GCWorker::createInMemoryDB() {
  QSqlDatabase fileDb = QSqlDatabase::addDatabase("QSQLITE", fileDbConn);
  fileDb.setDatabaseName(mappedOuputDBName);
  if (!fileDb.open()) {
    qDebug() << "Error opening file-based database:"
             << fileDb.lastError().text();
  }
  mem_db = QSqlDatabase::addDatabase("QSQLITE", memDbConn);
  mem_db.setDatabaseName(":memory:");
  if (!mem_db.open()) {
    qDebug() << "Error opening in-memory database:"
             << mem_db.lastError().text();
  }
  // Start a transaction to speed up the copying process
  mem_db.transaction();

  // Attach the file-based database to the in-memory database
  QSqlQuery attachQuery(mem_db);
  if (!attachQuery.exec(QString("ATTACH DATABASE '%1' AS fileDb")
                            .arg(fileDb.databaseName()))) {
    qDebug() << "Error attaching file-based database to in-memory database:"
             << attachQuery.lastError().text();
  }

  // Get the table schema
  QSqlQuery fileQuery(fileDb);
  fileQuery.exec(QString(
      "SELECT sql FROM sqlite_master WHERE type='table' AND name='maps'"));
  if (fileQuery.next()) {
    QString createTableSql = fileQuery.value(0).toString();
    // Create the table in the in-memory database
    QSqlQuery memoryQuery(mem_db);
    if (!memoryQuery.exec(createTableSql)) {
      qDebug() << "Error creating table in in-memory database:"
               << memoryQuery.lastError().text();
    }
    // Copy data from file-based database to in-memory database
    if (!memoryQuery.exec(
            QString("INSERT INTO maps SELECT * FROM fileDb.maps"))) {
      qDebug() << "Error copying data to in-memory database:"
               << memoryQuery.lastError().text();
    }
  }
  // Commit the transaction
  mem_db.commit();
  // Close the file-based database connection
  fileQuery.clear();
  attachQuery.clear();
  fileDb.close();
  fileDb = QSqlDatabase();
  QSqlDatabase::removeDatabase(fileDbConn);

  //  // Now you can query the in-memory database (memoryDb) for faster
  //  performance
  //  // Example:
  //  QSqlQuery m_query(mem_db);
  //  if (m_query.exec("SELECT * FROM maps")) {
  //    while (m_query.next()) {
  //      qDebug() << m_query.value(0).toString()
  //               << m_query.value(1)
  //                      .toString();  // Adjust according to your table schema
  //    }
  //  } else {
  //    qDebug() << "Error executing query on in-memory database:"
  //             << m_query.lastError().text();
  //  }

  // Close the in-memory database connection when you're done
  //  mem_db.close();
  //  QSqlDatabase::removeDatabase("memoryDb");
}

void GCWorker::readMapOutput() {
  elapsedTimer.start();
  qDebug() << stat->input;
  QFile f(stat->input);
  if (!f.open(QIODevice::ReadOnly)) {
    qDebug() << "Unable to open file";
    return;
  }
  stat->readCount = 1;
  QList<ReadHits> rCollect;
  ReadHits rHits;
  QTextStream in(&f);
  while (!in.atEnd()) {
    QString line = in.readLine();
    MapHit hit = MapHit(line);
    if (hit.identity_ > 97 && hit.bitscore_ > 50) {
      if (rHits.size() > 0) {
        if (rHits == hit) {
          if (hit.bitscore_ > 0.98 * rHits.lastBitscore()) {
            rHits.add(hit);
          }
        } else {
          rCollect.append(ReadHits(rHits));
          if (stat->readCount % 10000 == 0) {
            stat->elapsedTime = elapsedTimer.elapsed() / 1000;
            stat->readName = hit.read_;
            emit sig->gc_update_progress_sig();
          }
          rHits.clear();
          stat->readCount += 1;
          rHits.add(hit);
        }
      } else {
        rHits.add(hit);
      }
    }
  }
  if (rCollect.length() > 0) {
    writeToDatabase(&rCollect);
  }
  query.exec("CREATE INDEX maps_idx ON maps (read, gene, refseq)");
  f.close();
  stat->running = false;
  elapsedTimer.invalidate();
  query.clear();
  db.commit();
  db.close();
  db = QSqlDatabase();
  QSqlDatabase::removeDatabase(writeDbConn);
  emit sig->gc_update_progress_sig();
}

void GCWorker::setupDB() {
  // Ensure output directory exists
  QFileInfo dbInfo(mappedOuputDBName);
  QDir().mkpath(dbInfo.absolutePath());

  db = QSqlDatabase::addDatabase("QSQLITE", writeDbConn);
  db.setDatabaseName(mappedOuputDBName);
  if (!db.open()) {
    qDebug() << "Error opening database:" << db.lastError().text();
    return;
  }
  query = QSqlQuery(db);
  // Create tables (or clear existing)
  query.exec("CREATE TABLE IF NOT EXISTS maps ("
             "id INTEGER PRIMARY KEY AUTOINCREMENT, "
             "read TEXT, gene TEXT, qstart INTEGER, qend INTEGER, "
             "refseq TEXT, frame TEXT, location TEXT)");
  query.exec("DELETE FROM maps");
  query.prepare("INSERT INTO maps (read, gene, qstart, qend, refseq, frame, location) "
                "VALUES (:read, :gene, :qstart, :qend, :refseq, :frame, :location)");
  db.transaction();
}

void GCWorker::writeToDatabase(QList<ReadHits> *collectedReads) {
  QVariantList reads;
  QVariantList genes;
  QVariantList qstarts;
  QVariantList qends;
  QVariantList nms;
  QVariantList frames;
  QVariantList locations;
  for (ReadHits &hits : *collectedReads) {
    for (const MapHit &hit : hits) {
      reads << hit.read_;
      genes << hit.gene_;
      qstarts << hit.qstart_;
      qends << hit.qend_;
      nms << hit.nm_number_;
      frames << hit.frame();
      locations << hit.location();
    }
  }

  if (reads.length() != genes.length() || reads.length() != qstarts.length() ||
      reads.length() != qends.length() || reads.length() != nms.length() ||
      reads.length() != frames.length() ||
      reads.length() != locations.length()) {
    qDebug() << "Error: Length of QVariantList objects is not equal.";
    return;
  }

  query.prepare(
      "INSERT INTO maps (read, gene, qstart, qend, refseq, frame, location) "
      "VALUES (:read, :gene, :qstart, :qend, :refseq, :frame, :location)");
  query.bindValue(":read", reads);
  query.bindValue(":gene", genes);
  query.bindValue(":qstart", qstarts);
  query.bindValue(":qend", qends);
  query.bindValue(":refseq", nms);
  query.bindValue(":frame", frames);
  query.bindValue(":location", locations);
  if (!query.execBatch()) {
    qDebug() << "Error:" << query.lastError().text();
    return;
  }
  collectedReads->clear();
}

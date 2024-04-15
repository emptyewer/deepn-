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
  qDebug() << mappedOuputDBName;
  setupDB();
}

void GCWorker::run() {
  readMapOutput();
  createInMemoryDB();
  writeGeneCount();
}

void GCWorker::writeGeneCount() { qDebug() << "Writing Gene Count"; }

void GCWorker::createInMemoryDB() {
  QSqlDatabase fileDb = QSqlDatabase::addDatabase("QSQLITE", "fileDb");
  fileDb.setDatabaseName(mappedOuputDBName);
  if (!fileDb.open()) {
    qDebug() << "Error opening file-based database:"
             << fileDb.lastError().text();
  }
  mem_db = QSqlDatabase::addDatabase("QSQLITE", "memoryDb");
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
  fileDb.close();
  QSqlDatabase::removeDatabase("fileDb");

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
  QSqlDatabase::removeDatabase("writeDB");
  emit sig->gc_update_progress_sig();
}

void GCWorker::setupDB() {
  db = QSqlDatabase::addDatabase("QSQLITE", "writeDB");
  db.setDatabaseName(mappedOuputDBName);
  QFileInfo dbInfo(mappedOuputDBName);
  if (!dbInfo.exists()) {
    return;
  }
  if (!db.open()) {
    qDebug() << "Error opening database";
  }
  query = QSqlQuery(db);
  query.exec("DELETE FROM maps");
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

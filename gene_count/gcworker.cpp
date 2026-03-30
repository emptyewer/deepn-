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
  writeDbConn = "writeDB_" + QUuid::createUuid().toString(QUuid::Id128);
  qDebug() << mappedOuputDBName;
}

void GCWorker::run() {
  setupDB();
  readMapOutput();
  writeGeneCount();
  emit sig->gc_finished_sig();
  emit finished();
}

void GCWorker::writeGeneCount() {
  qDebug() << "Writing Gene Count";
  // Reopen the on-disk database for the aggregate queries.
  // The index on (read, gene, refseq) created in readMapOutput() makes
  // these GROUP BY queries fast without needing an in-memory copy.
  QString readConn = "gcRead_" + QUuid::createUuid().toString(QUuid::Id128);
  {
    QSqlDatabase fileDb = QSqlDatabase::addDatabase("QSQLITE", readConn);
    fileDb.setDatabaseName(mappedOuputDBName);
    if (!fileDb.open()) {
      qDebug() << "Error opening database for gene count:" << fileDb.lastError().text();
      return;
    }

    QSqlQuery q(fileDb);
    q.exec("PRAGMA cache_size = -2048");
    q.exec("PRAGMA mmap_size = 0");

    // Totals from the maps table
    q.exec("SELECT COUNT(DISTINCT read) FROM maps");
    int totalReads = q.next() ? q.value(0).toInt() : 0;
    q.exec("SELECT COUNT(*) FROM maps");
    int totalHits = q.next() ? q.value(0).toInt() : 0;
    q.clear();

    QSqlQuery w(fileDb);
    w.exec("DROP TABLE IF EXISTS gene_counts");
    w.exec("DROP TABLE IF EXISTS summary");

    // Gene counts: aggregate distinct reads per gene with PPM
    w.exec("CREATE TABLE gene_counts ("
           "gene TEXT PRIMARY KEY, count INTEGER, ppm REAL)");
    w.exec(QString(
        "INSERT INTO gene_counts (gene, count, ppm) "
        "SELECT gene, COUNT(DISTINCT read), "
        "COUNT(DISTINCT read) * 1000000.0 / %1 "
        "FROM maps GROUP BY gene ORDER BY COUNT(DISTINCT read) DESC")
        .arg(qMax(totalReads, 1)));

    // Summary metadata
    w.exec("CREATE TABLE summary (key TEXT PRIMARY KEY, value TEXT)");
    QFileInfo fileInfo(stat->input);
    w.exec(QString("INSERT INTO summary VALUES ('file', '%1')")
        .arg(fileInfo.fileName()));
    w.exec(QString("INSERT INTO summary VALUES ('total_reads', '%1')")
        .arg(totalReads));
    w.exec(QString("INSERT INTO summary VALUES ('total_hits', '%1')")
        .arg(totalHits));

    // Checkpoint and switch to DELETE mode so WAL/SHM files are removed
    w.exec("PRAGMA wal_checkpoint(TRUNCATE)");
    w.exec("PRAGMA journal_mode = DELETE");
    w.clear();
    fileDb.close();
  }
  QSqlDatabase::removeDatabase(readConn);

  qDebug() << "Gene counts written to:" << mappedOuputDBName;
}

void GCWorker::readMapOutput() {
  static const int kCommitInterval = 50000;
  elapsedTimer.start();
  qDebug() << stat->input;
  QFile f(stat->input);
  if (!f.open(QIODevice::ReadOnly)) {
    qDebug() << "Unable to open file";
    return;
  }
  stat->readCount = 1;
  int insertsSinceCommit = 0;
  ReadHits rHits;

  while (!f.atEnd()) {
    QByteArray rawLine = f.readLine();
    if (rawLine.size() < 20) continue;

    // Pre-filter on identity/bitscore from raw tab fields to avoid
    // constructing MapHit + breakDownMap for rejected lines.
    int tabCount = 0;
    int tab2 = -1, tab3 = -1, tab11 = -1, tab12 = -1;
    for (int i = 0; i < rawLine.size(); ++i) {
      if (rawLine[i] == '\t') {
        tabCount++;
        if (tabCount == 2)  tab2 = i + 1;
        if (tabCount == 3)  tab3 = i;
        if (tabCount == 11) tab11 = i + 1;
        if (tabCount == 12) { tab12 = i; break; }
      }
    }
    if (tabCount < 11) continue;
    if (tab12 < 0) tab12 = rawLine.size();
    // Trim trailing newline for the last field
    while (tab12 > tab11 && (rawLine[tab12 - 1] == '\n' || rawLine[tab12 - 1] == '\r'))
      tab12--;

    double identity = QByteArray::fromRawData(rawLine.constData() + tab2, tab3 - tab2).toDouble();
    double bitscore = QByteArray::fromRawData(rawLine.constData() + tab11, tab12 - tab11).toDouble();
    if (identity <= 97 || bitscore <= 50) continue;

    QString line = QString::fromUtf8(rawLine).trimmed();
    MapHit hit(line);

    if (rHits.size() > 0) {
      if (rHits == hit) {
        if (hit.bitscore_ > 0.98 * rHits.lastBitscore()) {
          rHits.add(hit);
        }
      } else {
        // Write the completed read group directly to DB — no accumulation
        writeReadHitsToDB(rHits);
        insertsSinceCommit += rHits.size();

        // Commit periodically so SQLite doesn't hold a huge WAL
        if (insertsSinceCommit >= kCommitInterval) {
          query.clear();
          db.commit();
          db.transaction();
          query.prepare(
              "INSERT INTO maps (read, gene, qstart, qend, refseq, frame, location, rstart, rend) "
              "VALUES (:read, :gene, :qstart, :qend, :refseq, :frame, :location, :rstart, :rend)");
          insertsSinceCommit = 0;
        }

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
  // Write final read group
  if (rHits.size() > 0) {
    writeReadHitsToDB(rHits);
  }
  query.clear();
  db.commit();
  // Create indexes optimized for downstream queries:
  //   MultiQuery++: WHERE gene = ? OR refseq = ?
  //   ReadDepth++:  WHERE (gene = ? OR refseq = ?) AND rstart <= ? AND rend >= ?
  //   Gene list:    SELECT DISTINCT gene FROM maps
  QSqlQuery idxQ(db);
  idxQ.exec("DROP INDEX IF EXISTS maps_idx");
  idxQ.exec("CREATE INDEX IF NOT EXISTS idx_maps_gene ON maps (gene)");
  idxQ.exec("CREATE INDEX IF NOT EXISTS idx_maps_refseq ON maps (refseq)");
  idxQ.exec("CREATE INDEX IF NOT EXISTS idx_maps_depth ON maps (gene, rstart, rend)");
  idxQ.exec("CREATE INDEX IF NOT EXISTS idx_maps_depth_ref ON maps (refseq, rstart, rend)");
  // Checkpoint WAL and switch back to DELETE journal so WAL/SHM files are removed
  idxQ.exec("PRAGMA wal_checkpoint(TRUNCATE)");
  idxQ.exec("PRAGMA journal_mode = DELETE");
  idxQ.clear();
  f.close();
  stat->running = false;
  elapsedTimer.invalidate();
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
  // Constrain SQLite memory usage: small page cache, no mmap, WAL mode
  query.exec("PRAGMA journal_mode = WAL");
  query.exec("PRAGMA cache_size = -2048");    // 2 MB page cache max
  query.exec("PRAGMA mmap_size = 0");         // disable memory-mapped I/O
  query.exec("PRAGMA wal_autocheckpoint = 1000");  // checkpoint frequently
  // Drop and recreate maps table (preserves reads table from JunctionDice++)
  query.exec("DROP TABLE IF EXISTS maps");
  query.exec("DROP INDEX IF EXISTS maps_idx");
  query.exec("DROP TABLE IF EXISTS gene_counts");
  query.exec("DROP TABLE IF EXISTS summary");
  query.exec("CREATE TABLE maps ("
             "id INTEGER PRIMARY KEY AUTOINCREMENT, "
             "read TEXT, gene TEXT, qstart INTEGER, qend INTEGER, "
             "refseq TEXT, frame TEXT, location TEXT, "
             "rstart INTEGER, rend INTEGER)");
  // Schema metadata for version tracking
  query.exec("CREATE TABLE IF NOT EXISTS schema_meta ("
             "key TEXT PRIMARY KEY, value TEXT)");
  query.exec("INSERT OR REPLACE INTO schema_meta VALUES ('schema_version', '2')");
  query.prepare("INSERT INTO maps (read, gene, qstart, qend, refseq, frame, location, rstart, rend) "
                "VALUES (:read, :gene, :qstart, :qend, :refseq, :frame, :location, :rstart, :rend)");
  db.transaction();
}

void GCWorker::writeReadHitsToDB(ReadHits& hits) {
  for (const MapHit &hit : hits) {
    query.bindValue(":read", hit.read_);
    query.bindValue(":gene", hit.gene_);
    query.bindValue(":qstart", hit.qstart_);
    query.bindValue(":qend", hit.qend_);
    query.bindValue(":refseq", hit.nm_number_);
    query.bindValue(":frame", hit.frame());
    query.bindValue(":location", hit.location());
    query.bindValue(":rstart", hit.rstart_);
    query.bindValue(":rend", hit.rend_);
    if (!query.exec()) {
      qDebug() << "Insert error:" << query.lastError().text();
    }
  }
}

void GCWorker::writeToDatabase(QList<ReadHits> *collectedReads) {
  // Estimate capacity to avoid realloc
  int totalHits = 0;
  for (const ReadHits &hits : *collectedReads) {
    totalHits += hits.size();
  }

  QVariantList reads, genes, qstarts, qends, nms, frames, locations, rstarts, rends;
  reads.reserve(totalHits);
  genes.reserve(totalHits);
  qstarts.reserve(totalHits);
  qends.reserve(totalHits);
  nms.reserve(totalHits);
  frames.reserve(totalHits);
  locations.reserve(totalHits);
  rstarts.reserve(totalHits);
  rends.reserve(totalHits);

  // Move data out of ReadHits into flat lists, then free each ReadHits
  while (!collectedReads->isEmpty()) {
    ReadHits hits = std::move(collectedReads->takeFirst());
    for (const MapHit &hit : hits) {
      reads << hit.read_;
      genes << hit.gene_;
      qstarts << hit.qstart_;
      qends << hit.qend_;
      nms << hit.nm_number_;
      frames << hit.frame();
      locations << hit.location();
      rstarts << hit.rstart_;
      rends << hit.rend_;
    }
    // hits goes out of scope here, freeing its m_hits list
  }
  // collectedReads is now empty

  if (reads.isEmpty()) return;

  query.prepare(
      "INSERT INTO maps (read, gene, qstart, qend, refseq, frame, location, rstart, rend) "
      "VALUES (:read, :gene, :qstart, :qend, :refseq, :frame, :location, :rstart, :rend)");
  query.bindValue(":read", reads);
  query.bindValue(":gene", genes);
  query.bindValue(":qstart", qstarts);
  query.bindValue(":qend", qends);
  query.bindValue(":refseq", nms);
  query.bindValue(":frame", frames);
  query.bindValue(":location", locations);
  query.bindValue(":rstart", rstarts);
  query.bindValue(":rend", rends);
  if (!query.execBatch()) {
    qDebug() << "Error:" << query.lastError().text();
  }
}

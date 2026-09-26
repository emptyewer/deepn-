#include "jdworker.h"

#include <QCoreApplication>
#include <QDir>
#include <QJsonDocument>
#include <QJsonObject>
#include <QStandardPaths>
#include <QRegularExpression>
#include <QtSql>

#include "qglobal.h"
#include "zstr/zstr.hpp"

JDWorker::JDWorker(JDStat* jdstat, int fcount) {
  fileCount = fcount;
  stat = jdstat;
  dbConnectionName = "jdDb_" + QUuid::createUuid().toString(QUuid::Id128);
  QFileInfo fileInfo(stat->dstat.input);
  stat->dstat.fileSize = fileInfo.size();
  QString basename = fileInfo.baseName();
  stat->dstat.summary = fileInfo.absolutePath() + "/../junction_diced_fasta/" +
                        basename + ".json";
  stat->dstat.output = fileInfo.absolutePath() + "/../junction_diced_fasta/" +
                       basename + ".jdice.fasta";

  QString suffix = ".blastn";
  if (stat->mstat.algo == blat) {
    suffix = ".blat";
  } else if (stat->mstat.algo == megablast) {
    suffix = ".megablast";
  }
  readDepthFileName = fileInfo.absolutePath() + "/../analyzed_files/" +
                      basename + suffix + ".sqlite";
  stat->mstat.output = fileInfo.absolutePath() + "/../mapped_files/" +
                       basename + suffix + ".txt";
  // Build tier-1 junction pattern (full matchLength bases).
  // First (matchLength - 3) positions: allow expected nucleotide OR N (sequencer
  // ambiguity code). Last 3 positions: any nucleotide (loose 3' boundary).
  int char_count = 0;
  stat->dstat.match = stat->dstat.junction.right(stat->dstat.matchLength);
  foreach (QChar c, stat->dstat.match) {
    if (char_count < stat->dstat.matchLength - 3) {
      jseq_pattern += QString("[%1N]").arg(c);
    } else {
      jseq_pattern += QString("[ATGCN]");
    }
    char_count++;
  }

  // Build tier-2 pattern (rightmost 2/3 of matchLength).
  // All positions allow the expected nucleotide OR N — recovers reads where
  // sequencing errors cluster in the 5' (outer) part of the junction.
  {
    int matchLen2 = qMax(5, stat->dstat.matchLength * 2 / 3);
    if (matchLen2 < stat->dstat.matchLength) {
      QString match2 = stat->dstat.junction.right(matchLen2);
      foreach (QChar c, match2)
        jseq_pattern2 += QString("[%1N]").arg(c);
    }
  }

  // Build tier-3 pattern (rightmost 1/3 of matchLength, minimum 5 bases).
  // Each position allows expected base OR N — catches reads with N-calls
  // across much of the junction but clear sequence at the 3' boundary.
  {
    int matchLen3 = qMax(5, stat->dstat.matchLength / 3);
    if (matchLen3 < stat->dstat.matchLength * 2 / 3) {
      QString match3 = stat->dstat.junction.right(matchLen3);
      foreach (QChar c, match3)
        jseq_pattern3 += QString("[%1N]").arg(c);
    }
  }

  repeats_sequence = QString("(?:([ATGCN])(\\1{20,}))");
}

void JDWorker::createDepthDatabase() {
  db = QSqlDatabase::addDatabase("QSQLITE", dbConnectionName);
  db.setDatabaseName(readDepthFileName);
  QFileInfo dbInfo(readDepthFileName);
  QDir().mkpath(dbInfo.absolutePath());
  if (dbInfo.exists()) {
    QFile::remove(dbInfo.absoluteFilePath());
  }
  if (!db.open()) {
    stat->errorMessage = QString("Cannot open depth database: %1\n%2")
        .arg(readDepthFileName, db.lastError().text());
    stat->dicing = false;
    emit sig->jd_update_progress_sig();
    emit sig->jd_finished_sig();
    emit finished();
    return;
  }
  query = QSqlQuery(QSqlDatabase::database(dbConnectionName));
  query.exec("PRAGMA auto_vacuum = FULL;");
  query.exec("PRAGMA journal_mode = WAL;");

  query.exec(
      "CREATE TABLE reads (id INTEGER PRIMARY KEY AUTOINCREMENT, read TEXT, "
      "sequence TEXT)");
  query.exec(
      "CREATE TABLE maps (id INTEGER PRIMARY KEY AUTOINCREMENT, read TEXT, "
      "gene "
      "TEXT, qstart INTEGER, qend INTEGER, refseq TEXT, frame TEXT, location "
      "TEXT)");
  db.transaction();
  query.prepare("INSERT INTO reads (read, sequence) VALUES (:read, :sequence)");
}

JDWorker::~JDWorker() { process.close(); }

void JDWorker::mapCallBack() {
  if (mapFinished) return;
  QFile file(stat->mstat.output);
  stat->mstat.elapsedTime = elapsedTimer.elapsed() / 1000;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QTextStream in(&file);
    if (file.size() > 4096) {
      in.seek(file.size() - 4096);
      in.readLine();  // skip partial first line
    }
    // Read the last complete line with a read name
    QString lastRead;
    while (!in.atEnd()) {
      QString line = in.readLine();
      if (!line.isEmpty()) {
        lastRead = line.split("\t").at(0);
      }
    }
    file.close();
    if (!lastRead.isEmpty()) {
      stat->mstat.currentRead = lastRead;
      int idx = readIndex.value(lastRead, -1);
      if (idx >= 0 && readIndex.size() > 0) {
        stat->mstat.percentComplete = static_cast<float>(idx) * 100.0f / readIndex.size();
      }
    }
  }

  // Stall detection: abort if output file never appears (60s) or stays
  // empty (2 min) — covers both "aligner never started" and "aligner crashed"
  mapPollCount++;
  QFileInfo mappedFileInfo(stat->mstat.output);
  bool stallDetected = (mapPollCount >= 20 && !mappedFileInfo.exists()) ||
                       (mapPollCount >= 40 && mappedFileInfo.size() == 0);
  if (stallDetected) {
    stat->errorMessage = QString(
        "Mapping produced no output after 60 seconds. "
        "The reference database may be inaccessible or the aligner failed.\n"
        "Command: %1").arg(stat->mstat.map_command);
    stat->blasting = false;
    elapsedTimer.invalidate();
    mapFinished = true;
    mapTimer->stop();
    emit sig->jd_update_progress_sig();
    emit sig->jd_finished_sig();
    emit finished();
    return;
  }

  if (stat->mstat.percentComplete < 98) {
    emit sig->jd_update_progress_sig();
  } else {
    stat->mstat.percentComplete = 100;
    stat->blasting = false;
    elapsedTimer.invalidate();
    mapFinished = true;
    mapTimer->stop();
    emit sig->jd_update_progress_sig();
    emit sig->jd_finished_sig();
    emit finished();
  }
}

void JDWorker::run() {
  mapTimer = new QTimer(this);
  QObject::connect(mapTimer, &QTimer::timeout, this, &JDWorker::mapCallBack);
  elapsedTimer.start();
  if (!stat->readdice) {
    createDepthDatabase();
    doDice();
  } else {
    readDiceSummary();
    readDice();
  }
  elapsedTimer.restart();
  doMapping();
}

// A function to compute the reverse complement of a DNA sequence
QString JDWorker::reverseComplement(QString dna_sequence) {
  // Reverse the DNA sequence using the standard library `std::reverse` function
  std::reverse(dna_sequence.begin(), dna_sequence.end());

  // Create a vector of the complement nucleotides, in the order 'T', 'G', 'C',
  // 'A', 'N'
  QVector<QChar> complement_map = {'T', 'G', 'C', 'A', 'N'};

  // Complement each nucleotide in the reversed DNA sequence and store in a new
  // string
  QString complement_sequence;
  for (QChar &base : dna_sequence) {
    if (base == 'A') {
      base = complement_map[0];  // 'A' complements to 'T'
    } else if (base == 'C') {
      base = complement_map[1];  // 'C' complements to 'G'
    } else if (base == 'G') {
      base = complement_map[2];  // 'G' complements to 'C'
    } else if (base == 'T') {
      base = complement_map[3];  // 'T' complements to 'A'
    } else if (base == 'N') {
      base = complement_map[4];  // 'N' (unknown base) complements to 'N'
    }
    complement_sequence.push_back(base);
  }
  // Return the reverse complement sequence
  return complement_sequence;
}

QString JDWorker::translate(QString dna) {
  // Initialize the genetic code table
  QMap<QString, QString> geneticCode = {
      {"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "L"}, {"CTT", "L"},
      {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"}, {"ATT", "I"}, {"ATC", "I"},
      {"ATA", "I"}, {"ATG", "M"}, {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"},
      {"GTG", "V"}, {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},
      {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"}, {"ACT", "T"},
      {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"}, {"GCT", "A"}, {"GCC", "A"},
      {"GCA", "A"}, {"GCG", "A"}, {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "*"},
      {"TAG", "*"}, {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
      {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"}, {"GAT", "D"},
      {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"}, {"TGT", "C"}, {"TGC", "C"},
      {"TGA", "*"}, {"TGG", "W"}, {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"},
      {"CGG", "R"}, {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
      {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}};
  QString protein = "";

  for (int i = 0; i < dna.length() - 2; i += 3) {
    QString codon = dna.mid(i, 3);
    QString aminoAcid = geneticCode.value(codon, "X");
    protein += aminoAcid;
  }
  return protein;
}

void JDWorker::readDice() {
  QFile file(stat->dstat.output);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    stat->errorMessage = QString(
        "Cannot re-read diced FASTA: %1\n"
        "The file may have been moved or deleted.")
        .arg(stat->dstat.output);
    emit finished();
    return;
  }

  QTextStream in(&file);
  while (!in.atEnd()) {
    QString line = in.readLine();
    if (line.startsWith(">")) {
      line.replace(">", "");
      readIndex.insert(line, readIndex.size());
    }
    if (readIndex.size() % 50000 == 0) {
      stat->dstat.elapsedTime = elapsedTimer.elapsed() / 1000;
      emit sig->jd_update_progress_sig();
    }
  }
  file.close();
}

void JDWorker::doDice() {
  QString line;
  QRegularExpression pattern(jseq_pattern);
  QRegularExpression pattern2(!jseq_pattern2.isEmpty() ? jseq_pattern2 : jseq_pattern);
  QRegularExpression pattern3(!jseq_pattern3.isEmpty() ? jseq_pattern3 : jseq_pattern);
  const bool hasTier2 = !jseq_pattern2.isEmpty();
  const bool hasTier3 = !jseq_pattern3.isEmpty();
  QRegularExpression repeat_pattern(repeats_sequence);
  QDir().mkpath(QFileInfo(stat->dstat.output).absolutePath());
  QFile of(stat->dstat.output);
  if (!of.open(QIODevice::WriteOnly | QIODevice::Text)) {
    stat->errorMessage = QString(
        "Cannot create output file: %1\n"
        "Check that the output directory exists and is writable.")
        .arg(stat->dstat.output);
    stat->dicing = false;
    emit sig->jd_update_progress_sig();
    emit finished();
    return;
  }
  QTextStream out(&of);

  std::unique_ptr<std::istream> is_p = std::unique_ptr<std::istream>(
      new zstr::ifstream(stat->dstat.input.toStdString()));
  const std::streamsize buff_size = 65536;
  char* buff = new char[buff_size];
  QVariantList sql_ids;
  QVariantList sql_seqs;
  while (true) {
    is_p->read(buff, buff_size);
    std::streamsize cnt = is_p->gcount();
    int line_slice_len = 0;
    for (int i = 0; i < cnt; i++) {
      char c = buff[i];
      if (c == '\n') {
        if (line.startsWith("@")) {
          stat->dstat.totalReads += 1;
          stat->dstat.currentRead = line.split(QRegularExpression("\\s+")).at(0);
        } else if (line.contains(QRegularExpression("^[ATGCN]+$"))) {
          QRegularExpressionMatch repeat_match = repeat_pattern.match(line);
          if (!repeat_match.hasMatch()) {
            QString reverseSeq = reverseComplement(line);
            // Three-tier junction search (Piper 2016): try progressively shorter
            // patterns to recover reads with sequencing errors in the 5' part of the
            // junction sequence.
            QRegularExpressionMatch match = pattern.match(line);
            QRegularExpressionMatch rev_match = pattern.match(reverseSeq);
            if (!match.hasMatch() && !rev_match.hasMatch() && hasTier2) {
              match = pattern2.match(line);
              rev_match = pattern2.match(reverseSeq);
            }
            if (!match.hasMatch() && !rev_match.hasMatch() && hasTier3) {
              match = pattern3.match(line);
              rev_match = pattern3.match(reverseSeq);
            }
            if (match.hasMatch()) {
              line_slice_len = line.length() - (match.capturedStart(0) +
                                                match.capturedLength(0));
              if (line_slice_len > 25) {
                QString seq = line.right(line_slice_len);
                out << ">" << stat->dstat.currentRead << Qt::endl;
                out << seq << Qt::endl;
                sql_ids << stat->dstat.currentRead;
                sql_seqs << seq;
                readIndex.insert(stat->dstat.currentRead, readIndex.size());
                stat->dstat.forwardMatches += 1;
                stat->dstat.writtenReads += 1;
              } else {
                stat->dstat.removedReads += 1;
              }
            } else if (rev_match.hasMatch()) {
              line_slice_len = line.size() - (rev_match.capturedStart(0) +
                                              rev_match.capturedLength(0));
              if (line_slice_len > 25) {
                QString seq = reverseSeq.right(line_slice_len);
                out << ">" << stat->dstat.currentRead << Qt::endl;
                out << seq << Qt::endl;
                sql_ids << stat->dstat.currentRead;
                sql_seqs << seq;
                readIndex.insert(stat->dstat.currentRead, readIndex.size());
                stat->dstat.reverseMatches += 1;
                stat->dstat.writtenReads += 1;
              } else {
                stat->dstat.removedReads += 1;
              }
            } else {
              // No junction match in forward or reverse — discard read.
              stat->dstat.removedReads += 1;
            }
            if (stat->dstat.totalReads % 50000 == 0) {
              query.bindValue(":read", sql_ids);
              query.bindValue(":sequence", sql_seqs);
              query.execBatch();
              sql_ids.clear();
              sql_seqs.clear();
              stat->dstat.elapsedTime = elapsedTimer.elapsed() / 1000;
              emit sig->jd_update_progress_sig();
            }
          }
        } else if (line.startsWith("+")) {
          //          out << line << Qt::endl;
        } else {
          //          out << line.left(line_slice_len) << Qt::endl;
        }
        line.clear();
      } else {
        line += c;
      }
    }
    if (cnt == 0) break;
  }
  if (sql_ids.length() > 0) {
    query.bindValue(":read", sql_ids);
    query.bindValue(":sequence", sql_seqs);
    query.execBatch();
  }
  query.exec("CREATE INDEX reads_idx ON reads (read)");
  // Commit transaction before WAL checkpoint — checkpoint is a no-op inside
  // an open write transaction.
  db.commit();
  // Checkpoint WAL and switch to DELETE journal so WAL/SHM files are removed
  query.exec("PRAGMA wal_checkpoint(TRUNCATE)");
  query.exec("PRAGMA journal_mode = DELETE");
  query.clear();
  stat->dicing = false;
  writeDiceSummary();
  delete[] buff;
  out.flush();
  of.close();
  // Release query reference before removing the connection
  query = QSqlQuery();
  db.close();
  QSqlDatabase::removeDatabase(dbConnectionName);
  emit sig->jd_update_progress_sig();
}

void JDWorker::writeDiceSummary() {
  QJsonObject junction;
  junction["forward_match_count"] = stat->dstat.forwardMatches;
  junction["reverse_match_count"] = stat->dstat.reverseMatches;
  junction["sequence"] = stat->dstat.junction;
  junction["match_length"] = stat->dstat.matchLength;
  junction["match_sequence"] = stat->dstat.match;
  QJsonObject dicesum;
  dicesum["input"] = stat->dstat.input;
  dicesum["output"] = stat->dstat.output;
  dicesum["time"] = stat->dstat.elapsedTime;
  dicesum["total_reads"] = stat->dstat.totalReads;
  dicesum["diced_reads"] = stat->dstat.writtenReads;
  dicesum["junction"] = junction;

  QJsonDocument doc(dicesum);
  // Write the JSON document to a file
  QFile file(stat->dstat.summary);
  if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    file.write(doc.toJson());
    file.close();
  }
}

void JDWorker::readDiceSummary() {
  QFile file(stat->dstat.summary);
  if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    QString contents = file.readAll();
    QJsonDocument doc = QJsonDocument::fromJson(contents.toUtf8());
    QJsonObject root = doc.object();
    stat->dstat.totalReads = root["total_reads"].toInt();
    stat->dstat.writtenReads = root["diced_reads"].toInt();
    stat->dstat.elapsedTime = root["time"].toInt();
    QJsonObject junction = root["junction"].toObject();
    stat->dstat.junction = junction["sequence"].toString();
    stat->dstat.forwardMatches = junction["forward_match_count"].toInt();
    stat->dstat.reverseMatches = junction["reverse_match_count"].toInt();
  }
}

void JDWorker::doMapping() {
  stat->blasting = true;
  QString exec_path;
  QStringList args;
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    exec_path = QDir::cleanPath(application_directory.path() +
                                QDir::separator() + "Tools");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    application_directory.cdUp();
    exec_path = application_directory.path();
  } else {
    application_directory.cdUp();
    exec_path = application_directory.path();
  }
  if (stat->mstat.algo == blat) {
    exec_path = QDir::cleanPath(exec_path + QDir::separator() + "blat");
    args << stat->mstat.db << stat->dstat.output << stat->mstat.output
         << "-minIdentity=98"
         << "-out=blast8";
  } else if (stat->mstat.algo == blastn || stat->mstat.algo == megablast) {
    exec_path = QDir::cleanPath(exec_path + QDir::separator() + "blastn");
    QString options =
        "-ungapped -dust no -outfmt 6 -evalue 0.2 "
        "-use_index true -perc_identity 98";
    if (stat->mstat.algo == blastn) {
      int threads = QThread::idealThreadCount() / fileCount;
      if (threads <= 0) {
        threads = 1;
      }
      args << "-task"
           << "blastn"
           << "-num_threads" << QString("%1").arg(threads);
    } else {
      args << "-task"
           << "megablast"
           << "-num_threads" << QString("%1").arg(QThread::idealThreadCount());
    }
    args << options.split(QRegularExpression("\\s+")) << "-query" << stat->dstat.output
         << "-db" << stat->mstat.db << "-out" << stat->mstat.output;
  }
  QDir().mkpath(QFileInfo(stat->mstat.output).absolutePath());
  stat->mstat.map_command = exec_path + " " + args.join(" ");
  qDebug() << "Launching:" << stat->mstat.map_command;
  mapTimer->start(3000);
  // Launch as background process via system() since QProcess::startDetached
  // silently fails from worker threads on macOS.
  // Double-quotes inside path components are escaped so the shell does not
  // interpret them as argument delimiters.
  auto shellQuote = [](const QString &s) -> QString {
    return "\"" + QString(s).replace("\"", "\\\"") + "\"";
  };
  QString cmd = shellQuote(QDir::toNativeSeparators(exec_path));
  for (const QString &arg : args) {
    cmd += " " + shellQuote(arg);
  }
  cmd += " &";
  qDebug() << "system():" << cmd;
  system(cmd.toUtf8().constData());
}

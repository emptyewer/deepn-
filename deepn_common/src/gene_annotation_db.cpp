#include "gene_annotation_db.h"

#include <QCryptographicHash>
#include <QDateTime>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QTextStream>
#include <QUuid>

namespace deepn {

bool GeneAnnotationDB::loadFromFasta(const QString& fastaPath)
{
    m_error.clear();
    m_byRefseq.clear();
    m_byGeneName.clear();

    QFileInfo fastaInfo(fastaPath);
    if (!fastaInfo.exists()) {
        m_error = QStringLiteral("FASTA file not found: %1").arg(fastaPath);
        qDebug() << m_error;
        return false;
    }

    QString cachePath = fastaPath + ".annotations.sqlite";
    QFileInfo cacheInfo(cachePath);

    // Check if cache exists and is valid
    if (cacheInfo.exists()) {
        QString currentChecksum = computeChecksum(fastaPath);
        if (currentChecksum.isEmpty()) {
            m_error = QStringLiteral("Failed to compute checksum for: %1").arg(fastaPath);
            qDebug() << m_error;
            return false;
        }

        // Open cache and verify checksum
        QString connName = "annot_check_" + QUuid::createUuid().toString(QUuid::Id128);
        {
            QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
            db.setDatabaseName(cachePath);
            if (db.open()) {
                QSqlQuery q(db);
                q.prepare("SELECT value FROM metadata WHERE key = 'source_checksum'");
                if (q.exec() && q.next()) {
                    QString storedChecksum = q.value(0).toString();
                    q.clear();
                    db.close();
                    if (storedChecksum == currentChecksum) {
                        QSqlDatabase::removeDatabase(connName);
                        return loadCache(cachePath);
                    }
                    qDebug() << "FASTA checksum changed, rebuilding cache";
                } else {
                    q.clear();
                    db.close();
                }
            } else {
                db.close();
            }
        }
        QSqlDatabase::removeDatabase(connName);
    }

    // Parse FASTA and create cache
    return parseFasta(fastaPath, cachePath);
}

bool GeneAnnotationDB::parseFasta(const QString& fastaPath, const QString& cachePath)
{
    QFile fastaFile(fastaPath);
    if (!fastaFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        m_error = QStringLiteral("Cannot open FASTA file: %1").arg(fastaPath);
        qDebug() << m_error;
        return false;
    }

    // Compute checksum
    QString checksum = computeChecksum(fastaPath);
    if (checksum.isEmpty()) {
        m_error = QStringLiteral("Failed to compute checksum for: %1").arg(fastaPath);
        qDebug() << m_error;
        return false;
    }

    // Remove old cache if exists
    if (QFileInfo::exists(cachePath)) {
        QFile::remove(cachePath);
    }

    // Create SQLite cache
    QString connName = "annot_write_" + QUuid::createUuid().toString(QUuid::Id128);
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(cachePath);
        if (!db.open()) {
            m_error = QStringLiteral("Cannot create cache database: %1").arg(db.lastError().text());
            qDebug() << m_error;
            QSqlDatabase::removeDatabase(connName);
            return false;
        }

        QSqlQuery q(db);
        q.exec("PRAGMA journal_mode = WAL");
        q.exec("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)");
        q.exec("CREATE TABLE genes ("
               "refseq TEXT PRIMARY KEY, "
               "gene_name TEXT, "
               "orf_start INTEGER, "
               "orf_end INTEGER, "
               "mrna_length INTEGER, "
               "sequence TEXT)");
        q.exec("CREATE INDEX idx_gene_name ON genes(gene_name)");

        db.transaction();

        // Parse the FASTA file
        QTextStream in(&fastaFile);
        QString currentRefseq;
        QString currentGeneName;
        int currentOrfStart = 0;
        int currentOrfEnd = 0;
        QString currentSequence;
        bool inEntry = false;
        int geneCount = 0;

        q.prepare("INSERT INTO genes (refseq, gene_name, orf_start, orf_end, mrna_length, sequence) "
                  "VALUES (:refseq, :gene_name, :orf_start, :orf_end, :mrna_length, :sequence)");

        auto flushEntry = [&]() {
            if (!inEntry || currentRefseq.isEmpty())
                return;

            int seqLen = currentSequence.length();

            // Insert into database
            q.bindValue(":refseq", currentRefseq);
            q.bindValue(":gene_name", currentGeneName);
            q.bindValue(":orf_start", currentOrfStart);
            q.bindValue(":orf_end", currentOrfEnd);
            q.bindValue(":mrna_length", seqLen);
            q.bindValue(":sequence", currentSequence);
            if (!q.exec()) {
                qDebug() << "Failed to insert gene:" << currentGeneName << q.lastError().text();
            }

            // Store in memory
            GeneAnnotation ann;
            ann.refseq = currentRefseq;
            ann.geneName = currentGeneName;
            ann.orfStart = currentOrfStart;
            ann.orfEnd = currentOrfEnd;
            ann.mRNALength = seqLen;
            ann.sequence = currentSequence;

            m_byRefseq.insert(currentRefseq, ann);
            m_byGeneName.insert(currentGeneName, ann);
            geneCount++;
        };

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty())
                continue;

            if (line.startsWith('>')) {
                // Flush previous entry
                flushEntry();

                // Reset for new entry
                currentRefseq.clear();
                currentGeneName.clear();
                currentOrfStart = 0;
                currentOrfEnd = 0;
                currentSequence.clear();
                inEntry = false;

                // Parse header: >GeneSymbol:TBC1D15|RefSeqNM:NM_146001|ORFstart:542|ORFstop:2180
                QString header = line.mid(1);  // strip '>'
                QStringList fields = header.split('|');

                for (const QString& field : fields) {
                    int colonPos = field.indexOf(':');
                    if (colonPos < 0)
                        continue;
                    QString key = field.left(colonPos).trimmed();
                    QString value = field.mid(colonPos + 1).trimmed();

                    if (key == "GeneSymbol") {
                        currentGeneName = value;
                    } else if (key == "RefSeqNM") {
                        currentRefseq = value;
                    } else if (key == "ORFstart") {
                        currentOrfStart = value.toInt();
                    } else if (key == "ORFstop") {
                        currentOrfEnd = value.toInt();
                    }
                }

                // Validation: reject headers missing GeneSymbol or RefSeqNM
                if (currentGeneName.isEmpty() || currentRefseq.isEmpty()) {
                    qDebug() << "Skipping FASTA entry with missing GeneSymbol or RefSeqNM:" << header;
                    currentRefseq.clear();
                    currentGeneName.clear();
                    inEntry = false;
                    continue;
                }

                inEntry = true;
            } else if (inEntry) {
                // Sequence line -- append (strip whitespace)
                currentSequence += line.toUpper();
            }
        }

        // Flush last entry
        flushEntry();

        // Write metadata
        q.prepare("INSERT INTO metadata (key, value) VALUES (:key, :value)");

        q.bindValue(":key", "source_checksum");
        q.bindValue(":value", checksum);
        q.exec();

        q.bindValue(":key", "source_path");
        q.bindValue(":value", fastaPath);
        q.exec();

        q.bindValue(":key", "created_at");
        q.bindValue(":value", QDateTime::currentDateTime().toString(Qt::ISODate));
        q.exec();

        q.bindValue(":key", "gene_count");
        q.bindValue(":value", QString::number(geneCount));
        q.exec();

        db.commit();
        q.clear();
        db.close();
    }
    QSqlDatabase::removeDatabase(connName);

    fastaFile.close();
    qDebug() << "Gene annotation cache created:" << cachePath << "with" << m_byRefseq.size() << "genes";
    return true;
}

bool GeneAnnotationDB::loadCache(const QString& cachePath)
{
    m_byRefseq.clear();
    m_byGeneName.clear();

    QString connName = "annot_load_" + QUuid::createUuid().toString(QUuid::Id128);
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(cachePath);
        if (!db.open()) {
            m_error = QStringLiteral("Cannot open cache database: %1").arg(db.lastError().text());
            qDebug() << m_error;
            QSqlDatabase::removeDatabase(connName);
            return false;
        }

        QSqlQuery q(db);
        if (!q.exec("SELECT refseq, gene_name, orf_start, orf_end, mrna_length, sequence FROM genes")) {
            m_error = QStringLiteral("Cannot query cache: %1").arg(q.lastError().text());
            qDebug() << m_error;
            q.clear();
            db.close();
            QSqlDatabase::removeDatabase(connName);
            return false;
        }

        while (q.next()) {
            GeneAnnotation ann;
            ann.refseq = q.value(0).toString();
            ann.geneName = q.value(1).toString();
            ann.orfStart = q.value(2).toInt();
            ann.orfEnd = q.value(3).toInt();
            ann.mRNALength = q.value(4).toInt();
            ann.sequence = q.value(5).toString();

            m_byRefseq.insert(ann.refseq, ann);
            m_byGeneName.insert(ann.geneName, ann);
        }

        q.clear();
        db.close();
    }
    QSqlDatabase::removeDatabase(connName);

    qDebug() << "Loaded" << m_byRefseq.size() << "gene annotations from cache:" << cachePath;
    return true;
}

QString GeneAnnotationDB::computeChecksum(const QString& filePath) const
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Cannot open file for checksum:" << filePath;
        return {};
    }

    QCryptographicHash hash(QCryptographicHash::Sha256);
    // Read first 1MB only for speed on large reference files
    constexpr qint64 maxBytes = 1024 * 1024;
    QByteArray data = file.read(maxBytes);
    hash.addData(data);
    file.close();

    return hash.result().toHex();
}

GeneAnnotation GeneAnnotationDB::findByRefseq(const QString& refseq) const
{
    auto it = m_byRefseq.find(refseq);
    if (it != m_byRefseq.end())
        return it.value();
    return {};
}

GeneAnnotation GeneAnnotationDB::findByGeneName(const QString& geneName) const
{
    auto it = m_byGeneName.find(geneName);
    if (it != m_byGeneName.end())
        return it.value();
    return {};
}

QVector<GeneAnnotation> GeneAnnotationDB::search(const QString& query, int maxResults) const
{
    QVector<GeneAnnotation> results;
    if (query.isEmpty())
        return results;

    QString upper = query.toUpper();

    // Exact match first
    auto exactIt = m_byGeneName.find(query);
    if (exactIt != m_byGeneName.end()) {
        results.append(exactIt.value());
    }

    // Prefix matches
    for (auto it = m_byGeneName.constBegin(); it != m_byGeneName.constEnd(); ++it) {
        if (results.size() >= maxResults)
            break;
        if (it.key().toUpper().startsWith(upper) && it.key() != query) {
            results.append(it.value());
        }
    }

    // Substring matches (if still room)
    if (results.size() < maxResults) {
        for (auto it = m_byGeneName.constBegin(); it != m_byGeneName.constEnd(); ++it) {
            if (results.size() >= maxResults)
                break;
            QString nameUpper = it.key().toUpper();
            if (nameUpper.contains(upper) && !nameUpper.startsWith(upper)) {
                results.append(it.value());
            }
        }
    }

    // Also search refseq if still room
    if (results.size() < maxResults) {
        for (auto it = m_byRefseq.constBegin(); it != m_byRefseq.constEnd(); ++it) {
            if (results.size() >= maxResults)
                break;
            if (it.key().toUpper().contains(upper)) {
                // Avoid duplicates -- check by refseq since we might already have it via gene name
                bool found = false;
                for (const auto& existing : results) {
                    if (existing.refseq == it.key()) {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    results.append(it.value());
            }
        }
    }

    return results;
}

QStringList GeneAnnotationDB::allGeneNames() const
{
    return m_byGeneName.keys();
}

QStringList GeneAnnotationDB::allRefseqs() const
{
    return m_byRefseq.keys();
}

int GeneAnnotationDB::count() const
{
    return m_byRefseq.size();
}

QString GeneAnnotationDB::lastError() const
{
    return m_error;
}

bool GeneAnnotationDB::loadFromSqlite(const QString& sqlitePath)
{
    m_error.clear();
    m_byRefseq.clear();
    m_byGeneName.clear();
    m_sqlitePath.clear();

    QFileInfo fi(sqlitePath);
    if (!fi.exists()) {
        m_error = QStringLiteral("SQLite annotation file not found: %1").arg(sqlitePath);
        qDebug() << m_error;
        return false;
    }

    QString connName = "annot_" + QUuid::createUuid().toString(QUuid::Id128);
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(sqlitePath);
        db.setConnectOptions("QSQLITE_OPEN_READONLY");
        if (!db.open()) {
            m_error = QStringLiteral("Cannot open annotation DB: %1").arg(db.lastError().text());
            qDebug() << m_error;
            return false;
        }

        QSqlQuery q(db);
        q.exec("PRAGMA cache_size = -2048");
        q.exec("PRAGMA mmap_size = 0");

        // Load only small metadata columns — skip sequence and length(sequence)
        // which force SQLite to read the entire 169MB+ of sequence data.
        // mRNALength and sequence are fetched per-gene on demand.
        if (!q.exec("SELECT name, nm_number, start, stop FROM gene")) {
            m_error = QStringLiteral("Query failed: %1").arg(q.lastError().text());
            qDebug() << m_error;
            db.close();
            QSqlDatabase::removeDatabase(connName);
            return false;
        }

        while (q.next()) {
            GeneAnnotation ann;
            ann.geneName = q.value(0).toString().trimmed();
            ann.refseq = q.value(1).toString().trimmed();
            ann.orfStart = q.value(2).toInt();
            ann.orfEnd = q.value(3).toInt();
            ann.mRNALength = 0;  // loaded on demand via populateGeneDetails()

            if (!ann.refseq.isEmpty()) {
                m_byRefseq.insert(ann.refseq, ann);
            }
            if (!ann.geneName.isEmpty()) {
                m_byGeneName.insert(ann.geneName, ann);
            }
        }

        q.clear();
        db.close();
    }
    QSqlDatabase::removeDatabase(connName);

    m_sqlitePath = sqlitePath;
    qDebug() << "Loaded" << m_byRefseq.size() << "gene annotations from" << sqlitePath;
    return m_byRefseq.size() > 0;
}

QString GeneAnnotationDB::sequenceForRefseq(const QString& refseq) const
{
    if (m_sqlitePath.isEmpty()) {
        // FASTA-loaded path: sequence is already in memory
        auto it = m_byRefseq.find(refseq);
        if (it != m_byRefseq.end()) {
            return it->sequence;
        }
        return {};
    }

    // SQLite path: query on demand
    QString result;
    QString connName = "seq_" + QUuid::createUuid().toString(QUuid::Id128);
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(m_sqlitePath);
        db.setConnectOptions("QSQLITE_OPEN_READONLY");
        if (db.open()) {
            QSqlQuery q(db);
            q.prepare("SELECT sequence FROM gene WHERE nm_number = :nm LIMIT 1");
            q.bindValue(":nm", refseq);
            if (q.exec() && q.next()) {
                result = q.value(0).toString();
            }
            q.clear();
            db.close();
        }
    }
    QSqlDatabase::removeDatabase(connName);
    return result;
}

void GeneAnnotationDB::populateGeneDetails(GeneAnnotation& annotation) const
{
    if (annotation.refseq.isEmpty()) return;

    // Already populated
    if (annotation.mRNALength > 0 && !annotation.sequence.isEmpty()) return;

    if (m_sqlitePath.isEmpty()) {
        // FASTA path: details already in memory
        auto it = m_byRefseq.find(annotation.refseq);
        if (it != m_byRefseq.end()) {
            if (annotation.mRNALength <= 0) annotation.mRNALength = it->mRNALength;
            if (annotation.sequence.isEmpty()) annotation.sequence = it->sequence;
        }
        return;
    }

    // SQLite path: single query for this gene
    QString connName = "det_" + QUuid::createUuid().toString(QUuid::Id128);
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(m_sqlitePath);
        db.setConnectOptions("QSQLITE_OPEN_READONLY");
        if (db.open()) {
            QSqlQuery q(db);
            q.prepare("SELECT sequence FROM gene WHERE nm_number = :nm LIMIT 1");
            q.bindValue(":nm", annotation.refseq);
            if (q.exec() && q.next()) {
                QString seq = q.value(0).toString();
                annotation.sequence = seq;
                annotation.mRNALength = seq.length();
            }
            q.clear();
            db.close();
        }
    }
    QSqlDatabase::removeDatabase(connName);
}

}  // namespace deepn

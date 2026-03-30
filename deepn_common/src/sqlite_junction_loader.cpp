#include "sqlite_junction_loader.h"

#include <QDebug>
#include <QFileInfo>
#include <QMap>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QSqlRecord>
#include <QUuid>

namespace deepn {

bool SqliteJunctionLoader::open(const QString& dbPath)
{
    m_error.clear();
    m_totalReads = -1;

    QFileInfo fi(dbPath);
    if (!fi.exists()) {
        m_error = QStringLiteral("Database file not found: %1").arg(dbPath);
        qDebug() << m_error;
        return false;
    }

    // Close previous connection if any
    close();

    m_dbPath = dbPath;
    m_connName = "jloader_" + QUuid::createUuid().toString(QUuid::Id128);

    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", m_connName);
    db.setDatabaseName(m_dbPath);
    db.setConnectOptions("QSQLITE_OPEN_READONLY");
    if (!db.open()) {
        m_error = QStringLiteral("Cannot open database: %1").arg(db.lastError().text());
        qDebug() << m_error;
        QSqlDatabase::removeDatabase(m_connName);
        m_connName.clear();
        return false;
    }

    // Verify the maps table exists
    QSqlQuery q(db);
    if (!q.exec("SELECT name FROM sqlite_master WHERE type='table' AND name='maps'")) {
        m_error = QStringLiteral("Cannot query schema: %1").arg(q.lastError().text());
        qDebug() << m_error;
        q.clear();
        db.close();
        QSqlDatabase::removeDatabase(m_connName);
        m_connName.clear();
        return false;
    }
    if (!q.next()) {
        m_error = QStringLiteral("No 'maps' table found in database: %1").arg(dbPath);
        qDebug() << m_error;
        q.clear();
        db.close();
        QSqlDatabase::removeDatabase(m_connName);
        m_connName.clear();
        return false;
    }
    q.clear();

    qDebug() << "Opened junction database:" << dbPath << "schema v" << schemaVersion();
    return true;
}

void SqliteJunctionLoader::close()
{
    if (m_connName.isEmpty())
        return;

    {
        QSqlDatabase db = QSqlDatabase::database(m_connName, false);
        if (db.isOpen())
            db.close();
    }
    QSqlDatabase::removeDatabase(m_connName);
    m_connName.clear();
    m_dbPath.clear();
    m_totalReads = -1;
}

bool SqliteJunctionLoader::hasPositionColumns() const
{
    return schemaVersion() >= 2;
}

int SqliteJunctionLoader::schemaVersion() const
{
    if (m_connName.isEmpty())
        return 0;

    QSqlDatabase db = QSqlDatabase::database(m_connName, false);
    if (!db.isOpen())
        return 0;

    QSqlQuery q(db);
    // Check for rstart/rend columns in maps table
    if (!q.exec("PRAGMA table_info(maps)")) {
        q.clear();
        return 0;
    }

    bool hasRstart = false;
    bool hasRend = false;
    while (q.next()) {
        QString colName = q.value(1).toString();
        if (colName == "rstart") hasRstart = true;
        if (colName == "rend") hasRend = true;
    }
    q.clear();

    if (hasRstart && hasRend)
        return 2;

    // Check for metadata table presence
    QSqlQuery q2(db);
    q2.exec("SELECT name FROM sqlite_master WHERE type='table' AND name='metadata'");
    bool hasMeta = q2.next();
    q2.clear();

    return hasMeta ? 1 : 1;  // v1: has maps but no rstart/rend
}

GeneJunctionProfile SqliteJunctionLoader::loadGeneJunctions(
    const QString& gene, const GeneAnnotation& annotation) const
{
    GeneJunctionProfile profile;
    profile.annotation = annotation;
    profile.datasetLabel = gene;
    profile.sourceFile = m_dbPath;

    if (m_connName.isEmpty()) {
        m_error = "No database open";
        qDebug() << m_error;
        return profile;
    }

    QSqlDatabase db = QSqlDatabase::database(m_connName, false);
    if (!db.isOpen()) {
        m_error = "Database not open";
        qDebug() << m_error;
        return profile;
    }

    int total = totalDistinctReads();
    profile.totalReads = total;

    bool v2 = hasPositionColumns();

    // Query junction sites for this gene
    QString sql;
    if (v2) {
        sql = "SELECT read, gene, qstart, qend, refseq, frame, location, rstart, rend "
              "FROM maps WHERE gene = :gene";
    } else {
        sql = "SELECT read, gene, qstart, qend, refseq, frame, location "
              "FROM maps WHERE gene = :gene";
    }

    QSqlQuery q(db);
    q.prepare(sql);
    q.bindValue(":gene", gene);
    if (!q.exec()) {
        m_error = QStringLiteral("Query failed: %1").arg(q.lastError().text());
        qDebug() << m_error;
        q.clear();
        return profile;
    }

    // Collect raw rows, grouping by position to count distinct reads
    // For PPM: we need distinct reads per position
    // position -> set of read names
    QMap<int, QMap<QString, JunctionSite>> positionReads;

    while (q.next()) {
        JunctionSite site;
        QString readName = q.value(0).toString();
        site.geneName = q.value(1).toString();
        site.queryStart = q.value(2).toInt();
        site.queryEnd = q.value(3).toInt();
        site.refseq = q.value(4).toString();
        site.frame = q.value(5).toString();
        site.cdsClass = q.value(6).toString();

        if (v2) {
            site.position = q.value(7).toInt();
            site.positionEnd = q.value(8).toInt();
        } else {
            // Fallback: use qstart as position proxy
            site.position = site.queryStart;
            site.positionEnd = site.queryEnd;
        }

        // Use position as grouping key; track unique reads per position
        if (!positionReads.contains(site.position)) {
            positionReads[site.position] = {};
        }
        // Store per-read (use first occurrence)
        if (!positionReads[site.position].contains(readName)) {
            positionReads[site.position].insert(readName, site);
        }
    }
    q.clear();

    // Build JunctionSite list with PPM
    double totalForPpm = qMax(total, 1);

    for (auto posIt = positionReads.constBegin(); posIt != positionReads.constEnd(); ++posIt) {
        const auto& readMap = posIt.value();
        int distinctCount = readMap.size();
        double ppm = (distinctCount * 1000000.0) / totalForPpm;

        for (auto readIt = readMap.constBegin(); readIt != readMap.constEnd(); ++readIt) {
            JunctionSite site = readIt.value();
            site.rawCount = distinctCount;
            site.ppm = ppm;
            profile.sites.append(site);
        }
    }

    return profile;
}

QVector<CollapsedJunction> SqliteJunctionLoader::collapseByPosition(
    const QVector<JunctionSite>& sites)
{
    // Group by position
    QMap<int, QVector<JunctionSite>> groups;
    for (const auto& site : sites) {
        groups[site.position].append(site);
    }

    QVector<CollapsedJunction> collapsed;
    collapsed.reserve(groups.size());

    for (auto it = groups.constBegin(); it != groups.constEnd(); ++it) {
        const auto& variants = it.value();
        CollapsedJunction cj;
        cj.position = it.key();
        cj.variants = variants;

        // Count distinct queryStart values
        QMap<int, int> qstartCounts;
        // Count frame occurrences
        QMap<QString, int> frameCounts;
        double totalPpm = 0.0;
        QString cdsClass;

        for (const auto& v : variants) {
            qstartCounts[v.queryStart]++;
            frameCounts[v.frame]++;
            totalPpm += v.ppm;
            if (cdsClass.isEmpty())
                cdsClass = v.cdsClass;
        }

        cj.variantCount = qstartCounts.size();

        // PPM: sum across all variants at this position (they share the same position PPM,
        // so just use the first one's PPM since they are all the same)
        // Actually, PPM was set per-position already, so just take from first variant
        if (!variants.isEmpty()) {
            cj.totalPpm = variants.first().ppm;
        } else {
            cj.totalPpm = totalPpm;
        }

        cj.cdsClass = cdsClass;

        // Dominant frame: most frequent
        QString bestFrame;
        int bestCount = 0;
        for (auto fi = frameCounts.constBegin(); fi != frameCounts.constEnd(); ++fi) {
            if (fi.value() > bestCount) {
                bestCount = fi.value();
                bestFrame = fi.key();
            }
        }
        cj.dominantFrame = bestFrame;

        collapsed.append(cj);
    }

    return collapsed;
}

QStringList SqliteJunctionLoader::availableGenes() const
{
    QStringList genes;
    if (m_connName.isEmpty()) {
        m_error = "No database open";
        return genes;
    }

    QSqlDatabase db = QSqlDatabase::database(m_connName, false);
    if (!db.isOpen()) {
        m_error = "Database not open";
        return genes;
    }

    QSqlQuery q(db);
    if (!q.exec("SELECT DISTINCT gene FROM maps ORDER BY gene")) {
        m_error = QStringLiteral("Query failed: %1").arg(q.lastError().text());
        qDebug() << m_error;
        q.clear();
        return genes;
    }

    while (q.next()) {
        genes.append(q.value(0).toString());
    }
    q.clear();
    return genes;
}

int SqliteJunctionLoader::totalDistinctReads() const
{
    if (m_totalReads >= 0)
        return m_totalReads;

    if (m_connName.isEmpty()) {
        m_error = "No database open";
        return 0;
    }

    QSqlDatabase db = QSqlDatabase::database(m_connName, false);
    if (!db.isOpen()) {
        m_error = "Database not open";
        return 0;
    }

    QSqlQuery q(db);
    if (!q.exec("SELECT COUNT(DISTINCT read) FROM maps")) {
        m_error = QStringLiteral("Query failed: %1").arg(q.lastError().text());
        qDebug() << m_error;
        q.clear();
        return 0;
    }

    if (q.next()) {
        m_totalReads = q.value(0).toInt();
    } else {
        m_totalReads = 0;
    }
    q.clear();
    return m_totalReads;
}

int SqliteJunctionLoader::geneReadCount(const QString& gene) const
{
    if (m_connName.isEmpty()) {
        m_error = "No database open";
        return 0;
    }

    QSqlDatabase db = QSqlDatabase::database(m_connName, false);
    if (!db.isOpen()) {
        m_error = "Database not open";
        return 0;
    }

    QSqlQuery q(db);
    q.prepare("SELECT COUNT(DISTINCT read) FROM maps WHERE gene = :gene");
    q.bindValue(":gene", gene);
    if (!q.exec()) {
        m_error = QStringLiteral("Query failed: %1").arg(q.lastError().text());
        qDebug() << m_error;
        q.clear();
        return 0;
    }

    int count = 0;
    if (q.next()) {
        count = q.value(0).toInt();
    }
    q.clear();
    return count;
}

QString SqliteJunctionLoader::lastError() const
{
    return m_error;
}

QString SqliteJunctionLoader::databasePath() const
{
    return m_dbPath;
}

}  // namespace deepn

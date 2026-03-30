#include "depth_calculator.h"

#include <QDebug>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QUuid>

#include <algorithm>
#include <cmath>

namespace deepn {

DepthProfile DepthCalculator::calculate(const QString& dbPath,
                                         const QString& gene,
                                         const GeneAnnotation& annotation,
                                         int intervalWidth,
                                         int intervalSpacing) const
{
    DepthProfile profile;
    profile.annotation = annotation;
    profile.intervalWidth = intervalWidth;
    profile.intervalSpacing = intervalSpacing;
    profile.datasetLabel = gene;
    profile.sourceFile = dbPath;

    m_error.clear();

    if (!annotation.isValid()) {
        m_error = QStringLiteral("Invalid gene annotation for: %1").arg(gene);
        qDebug() << m_error;
        return profile;
    }

    QString connName = "depth_" + QUuid::createUuid().toString(QUuid::Id128);
    const QString geneKey = annotation.geneName.isEmpty() ? gene : annotation.geneName;
    const QString refseqKey = annotation.refseq.isEmpty() ? gene : annotation.refseq;
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(dbPath);
        db.setConnectOptions("QSQLITE_OPEN_READONLY");
        if (!db.open()) {
            m_error = QStringLiteral("Cannot open database: %1").arg(db.lastError().text());
            qDebug() << m_error;
            QSqlDatabase::removeDatabase(connName);
            return profile;
        }

        // Get total distinct reads for normalization
        QSqlQuery totalQ(db);
        if (!totalQ.exec("SELECT COUNT(DISTINCT read) FROM maps")) {
            m_error = QStringLiteral("Total reads query failed: %1").arg(totalQ.lastError().text());
            qDebug() << m_error;
            totalQ.clear();
            db.close();
            QSqlDatabase::removeDatabase(connName);
            return profile;
        }
        int totalReads = 0;
        if (totalQ.next())
            totalReads = totalQ.value(0).toInt();
        totalQ.clear();
        profile.totalReads = totalReads;

        // Check for v2 schema (rstart/rend columns)
        QSqlQuery pragmaQ(db);
        pragmaQ.exec("PRAGMA table_info(maps)");
        bool hasRstart = false;
        bool hasRend = false;
        while (pragmaQ.next()) {
            QString colName = pragmaQ.value(1).toString();
            if (colName == "rstart") hasRstart = true;
            if (colName == "rend") hasRend = true;
        }
        pragmaQ.clear();

        if (!hasRstart || !hasRend) {
            m_error = "Database missing rstart/rend columns in maps table";
            qDebug() << m_error;
            db.close();
            QSqlDatabase::removeDatabase(connName);
            return profile;
        }

        int mRNALen = annotation.mRNALength;
        double totalForNorm = qMax(totalReads, 1);

        // Fetch all read spans for this gene in a single query, then
        // count interval coverage in memory. Much faster than per-interval queries.
        QSqlQuery q(db);
        q.prepare("SELECT DISTINCT read, rstart, rend FROM maps "
                  "WHERE gene = :gene "
                  "UNION "
                  "SELECT DISTINCT read, rstart, rend FROM maps "
                  "WHERE refseq = :refseq");
        q.bindValue(":gene", geneKey);
        q.bindValue(":refseq", refseqKey);

        struct ReadSpan { int rstart; int rend; };
        QVector<ReadSpan> spans;
        if (q.exec()) {
            while (q.next()) {
                spans.append({q.value(1).toInt(), q.value(2).toInt()});
            }
        }
        q.clear();

        for (int p = 0; p + intervalWidth <= mRNALen; p += intervalSpacing) {
            int pEnd = p + intervalWidth;
            int count = 0;
            for (const auto& s : spans) {
                if (s.rstart <= p && s.rend >= pEnd) {
                    count++;
                }
            }

            DepthPoint pt;
            pt.position = p;
            pt.count = count;
            pt.normalized = (count * 1000000.0) / totalForNorm;
            profile.points.append(pt);
        }

        db.close();
    }
    QSqlDatabase::removeDatabase(connName);

    return profile;
}

BoundaryResult DepthCalculator::detect3PrimeBoundary(const DepthProfile& profile)
{
    BoundaryResult result;
    result.confidenceLabel = "LOW";

    if (profile.points.size() < 3)
        return result;

    // Step 1: Smooth the depth profile
    QVector<DepthPoint> smoothed = smooth(profile.points, 3);

    // Step 2: Find max depth
    int maxDepth = 0;
    for (const auto& pt : smoothed) {
        if (pt.count > maxDepth)
            maxDepth = pt.count;
    }

    if (maxDepth == 0)
        return result;

    // Step 3: Compute gradient (difference between consecutive points)
    QVector<double> gradient;
    gradient.reserve(smoothed.size() - 1);
    for (int i = 1; i < smoothed.size(); i++) {
        gradient.append(smoothed[i].count - smoothed[i - 1].count);
    }

    // Step 4: Find largest negative gradient
    int largestDropIdx = -1;
    double largestDrop = 0.0;
    for (int i = 0; i < gradient.size(); i++) {
        if (gradient[i] < largestDrop) {
            largestDrop = gradient[i];
            largestDropIdx = i;
        }
    }

    if (largestDropIdx < 0)
        return result;

    // Step 5: Walk backward from largest drop to find last interval with depth > 10% of max
    double threshold = maxDepth * 0.10;
    int boundaryIdx = largestDropIdx;

    // The drop is between index largestDropIdx-1 and largestDropIdx in smoothed
    // Walk backward to last point above threshold
    for (int i = largestDropIdx; i >= 0; i--) {
        if (smoothed[i].count > threshold) {
            boundaryIdx = i;
            break;
        }
    }

    result.position = smoothed[boundaryIdx].position;
    result.depthBefore = smoothed[boundaryIdx].count;
    result.depthAfter = (boundaryIdx + 1 < smoothed.size()) ? smoothed[boundaryIdx + 1].count : 0;

    // Step 6: Score confidence
    double dropPercent = 0.0;
    if (result.depthBefore > 0) {
        dropPercent = 1.0 - (static_cast<double>(result.depthAfter) / result.depthBefore);
    }

    // Count intervals over which the drop occurs
    int dropIntervals = 0;
    if (largestDropIdx >= 0 && boundaryIdx < smoothed.size()) {
        // Walk forward from boundary to find where depth stays below threshold
        for (int i = boundaryIdx + 1; i < smoothed.size(); i++) {
            dropIntervals++;
            if (smoothed[i].count <= threshold)
                break;
        }
    }
    if (dropIntervals == 0) dropIntervals = 1;

    if (dropPercent > 0.80 && dropIntervals <= 2) {
        result.confidence = 0.9;
        result.confidenceLabel = "HIGH";
    } else if (dropPercent >= 0.50 && dropIntervals <= 5) {
        result.confidence = 0.6;
        result.confidenceLabel = "MEDIUM";
    } else {
        result.confidence = 0.3;
        result.confidenceLabel = "LOW";
    }

    // Look for secondary boundaries (other significant drops)
    for (int i = 0; i < gradient.size(); i++) {
        if (i == largestDropIdx)
            continue;
        // Secondary if drop is at least 30% of max drop
        if (gradient[i] < largestDrop * 0.3 && gradient[i] < -maxDepth * 0.1) {
            BoundaryResult secondary;
            secondary.position = smoothed[i].position;
            secondary.depthBefore = smoothed[i].count;
            secondary.depthAfter = (i + 1 < smoothed.size()) ? smoothed[i + 1].count : 0;
            secondary.confidence = 0.3;
            secondary.confidenceLabel = "LOW";
            result.secondary.append(secondary);
        }
    }

    return result;
}

QVector<DepthPoint> DepthCalculator::smooth(const QVector<DepthPoint>& points, int windowSize)
{
    if (points.isEmpty() || windowSize <= 1)
        return points;

    QVector<DepthPoint> result;
    result.reserve(points.size());

    int halfWindow = windowSize / 2;

    for (int i = 0; i < points.size(); i++) {
        DepthPoint pt;
        pt.position = points[i].position;

        double sumCount = 0.0;
        double sumNorm = 0.0;
        int n = 0;

        for (int j = i - halfWindow; j <= i + halfWindow; j++) {
            if (j >= 0 && j < points.size()) {
                sumCount += points[j].count;
                sumNorm += points[j].normalized;
                n++;
            }
        }

        pt.count = static_cast<int>(std::round(sumCount / n));
        pt.normalized = sumNorm / n;
        result.append(pt);
    }

    return result;
}

QVector<DepthPoint> DepthCalculator::downsample(const QVector<DepthPoint>& points, int targetCount)
{
    int n = points.size();
    if (n <= targetCount || targetCount < 3)
        return points;

    // Largest-Triangle-Three-Buckets (LTTB) algorithm
    QVector<DepthPoint> sampled;
    sampled.reserve(targetCount);

    // Always include first point
    sampled.append(points.first());

    // Bucket size (excluding first and last points)
    double bucketSize = static_cast<double>(n - 2) / (targetCount - 2);

    int prevSelected = 0;

    for (int i = 0; i < targetCount - 2; i++) {
        // Current bucket range
        int bucketStart = static_cast<int>(std::floor((i) * bucketSize)) + 1;
        int bucketEnd = static_cast<int>(std::floor((i + 1) * bucketSize)) + 1;
        if (bucketEnd > n - 1) bucketEnd = n - 1;

        // Next bucket range (for computing average point)
        int nextBucketStart = static_cast<int>(std::floor((i + 1) * bucketSize)) + 1;
        int nextBucketEnd = static_cast<int>(std::floor((i + 2) * bucketSize)) + 1;
        if (nextBucketEnd > n - 1) nextBucketEnd = n - 1;
        if (i == targetCount - 3) {
            // Last bucket before the final point
            nextBucketStart = n - 1;
            nextBucketEnd = n;
        }

        // Compute average point in next bucket
        double avgX = 0.0, avgY = 0.0;
        int nextCount = 0;
        for (int j = nextBucketStart; j < nextBucketEnd && j < n; j++) {
            avgX += points[j].position;
            avgY += points[j].count;
            nextCount++;
        }
        if (nextCount > 0) {
            avgX /= nextCount;
            avgY /= nextCount;
        }

        // Previous selected point
        double prevX = points[prevSelected].position;
        double prevY = points[prevSelected].count;

        // Find point in current bucket with largest triangle area
        double maxArea = -1.0;
        int bestIdx = bucketStart;

        for (int j = bucketStart; j < bucketEnd && j < n; j++) {
            double curX = points[j].position;
            double curY = points[j].count;

            // Triangle area (using cross product / 2)
            double area = std::abs((prevX - avgX) * (curY - prevY) -
                                   (prevX - curX) * (avgY - prevY)) * 0.5;

            if (area > maxArea) {
                maxArea = area;
                bestIdx = j;
            }
        }

        sampled.append(points[bestIdx]);
        prevSelected = bestIdx;
    }

    // Always include last point
    sampled.append(points.last());

    return sampled;
}

QString DepthCalculator::lastError() const
{
    return m_error;
}

}  // namespace deepn

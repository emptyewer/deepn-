#include "main_window.h"
#include "ui_mainwindow.h"
#include <QApplication>
#include <QFileDialog>
#include <QMessageBox>
#include <QHeaderView>
#include <QStandardPaths>
#include <QDateTime>
#include <QElapsedTimer>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QSettings>
#include <QDir>
#include <QProgressBar>
#include <QThread>
#include <QMutex>
#include <QMutexLocker>
#include <QPainter>
#include <QPrinter>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlError>
#include <QSqlRecord>
#include <QClipboard>
#include <QProcess>
#include <QInputDialog>
#include <QRegularExpression>
#include <QSvgGenerator>
#include <QUuid>
#include <stdexcept>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <QLegendMarker>

// Y2H-SCORES library
#include "y2h_scores.h"

// Qt Charts (replaces QCustomPlot)
#include <QScatterSeries>
#include <QLineSeries>
#include <QValueAxis>
#include <QLegend>

namespace deseq2
{

    namespace
    {
        constexpr int kDeseqColumnCount = 6;
        constexpr int kResultsColumnGene = 0;
        constexpr int kResultsColumnBaseMean = 1;
        constexpr int kResultsColumnLog2FoldChange = 2;
        constexpr int kResultsColumnPValue = 5;
        constexpr int kResultsColumnPadj = 6;
        constexpr int kResultsColumnEnrichmentCall = 7;
        constexpr int kResultsColumnEnrichmentScore = 8;
        constexpr int kResultsColumnSpecificityScore = 9;
        constexpr int kResultsColumnInFrameScore = 10;
        constexpr int kResultsColumnBordaScore = 11;
        constexpr int kResultsColumnInFrameTranscripts = 12;
        constexpr int kInteractivePointLimit = 5000;

        QString joinPaths(const QStringList &paths)
        {
            return paths.join("; ");
        }

        QString sanitizeGeneName(const std::string &geneName)
        {
            QString sanitized;
            sanitized.reserve(static_cast<int>(geneName.size()));
            for (unsigned char ch : geneName)
            {
                if (ch == '\0')
                {
                    sanitized += QChar(0xFFFD);
                }
                else if (ch < 32 || ch == 127)
                {
                    sanitized += '?';
                }
                else
                {
                    sanitized += QChar(ch);
                }
            }
            return sanitized.trimmed();
        }

        QString csvEscape(const QString &value)
        {
            QString escaped = value;
            escaped.replace('"', "\"\"");
            if (escaped.contains(',') || escaped.contains('"') || escaped.contains('\n'))
                return QString("\"%1\"").arg(escaped);
            return escaped;
        }

        bool geneNameHasInvalidBytes(const std::string &geneName)
        {
            for (unsigned char ch : geneName)
            {
                if (ch == '\0' || ch < 32 || ch == 127)
                    return true;
            }
            return false;
        }

        QStringList splitPaths(const QString &raw)
        {
            QStringList parts;
            for (const QString &path : raw.split(';', Qt::SkipEmptyParts))
            {
                const QString trimmed = path.trimmed();
                if (!trimmed.isEmpty())
                    parts << QFileInfo(trimmed).absoluteFilePath();
            }
            parts.removeDuplicates();
            return parts;
        }

        QString inferGroupNameForSample(const QString &sampleName)
        {
            const QString lower = sampleName.toLower();
            if (lower.contains("non-selected") || lower.contains("nonselected") ||
                lower.contains("non_selection") || lower.contains("_non_") ||
                lower.contains(" non ") || lower.endsWith("_non") ||
                lower.startsWith("non_") || lower.contains("vector"))
            {
                return QStringLiteral("Non-Selected");
            }
            if (lower.contains("selected") || lower.contains("_sel_") ||
                lower.endsWith("_sel") || lower.startsWith("sel_"))
            {
                return QStringLiteral("Selected");
            }
            if (lower.contains("control"))
                return QStringLiteral("Vector Control");
            if (lower.contains("bait"))
                return QStringLiteral("Bait");
            return {};
        }

        QStringList discoverFiles(const QString &directoryPath, const QStringList &filters)
        {
            QDir dir(directoryPath);
            if (!dir.exists())
                return {};

            QStringList files;
            const QFileInfoList infoList = dir.entryInfoList(filters, QDir::Files, QDir::Name);
            for (const QFileInfo &info : infoList)
                files << info.absoluteFilePath();
            return files;
        }

        QString databaseBaitName(const QStringList &groupNames)
        {
            for (const QString &groupName : groupNames)
            {
                if (groupName.compare("Vector Control", Qt::CaseInsensitive) != 0 &&
                    groupName.compare("Control", Qt::CaseInsensitive) != 0 &&
                    groupName.compare("Non-Selected", Qt::CaseInsensitive) != 0)
                {
                    return groupName;
                }
            }
            return groupNames.isEmpty() ? QStringLiteral("analysis") : groupNames.last();
        }

        QList<QPointF> downsamplePointsByXBucket(const QList<QPointF> &points,
                                                 int maxPoints)
        {
            if (points.size() <= maxPoints || maxPoints <= 0)
                return points;

            double minX = std::numeric_limits<double>::max();
            double maxX = std::numeric_limits<double>::lowest();
            for (const QPointF &point : points)
            {
                minX = std::min(minX, point.x());
                maxX = std::max(maxX, point.x());
            }

            if (!std::isfinite(minX) || !std::isfinite(maxX) || minX == maxX)
            {
                QList<QPointF> reduced;
                const int step = std::max(qsizetype(1), points.size() / maxPoints);
                for (int idx = 0; idx < points.size(); idx += step)
                    reduced.append(points[idx]);
                return reduced;
            }

            struct Bucket
            {
                bool hasPoint = false;
                QPointF minYPoint;
                QPointF maxYPoint;
            };

            const int bucketCount = std::max(1, maxPoints / 2);
            QVector<Bucket> buckets(bucketCount);
            const double bucketWidth = (maxX - minX) / bucketCount;

            for (const QPointF &point : points)
            {
                int bucketIndex = static_cast<int>((point.x() - minX) / std::max(bucketWidth, 1e-12));
                bucketIndex = std::clamp(bucketIndex, 0, bucketCount - 1);
                Bucket &bucket = buckets[bucketIndex];
                if (!bucket.hasPoint)
                {
                    bucket.hasPoint = true;
                    bucket.minYPoint = point;
                    bucket.maxYPoint = point;
                }
                else
                {
                    if (point.y() < bucket.minYPoint.y())
                        bucket.minYPoint = point;
                    if (point.y() > bucket.maxYPoint.y())
                        bucket.maxYPoint = point;
                }
            }

            QList<QPointF> reduced;
            reduced.reserve(maxPoints);
            for (const Bucket &bucket : buckets)
            {
                if (!bucket.hasPoint)
                    continue;
                reduced.append(bucket.minYPoint);
                if (bucket.maxYPoint != bucket.minYPoint)
                    reduced.append(bucket.maxYPoint);
            }

            std::sort(reduced.begin(), reduced.end(),
                      [](const QPointF &lhs, const QPointF &rhs)
                      {
                          if (lhs.x() == rhs.x())
                              return lhs.y() < rhs.y();
                          return lhs.x() < rhs.x();
                      });

            if (reduced.size() <= maxPoints)
                return reduced;

            QList<QPointF> capped;
            capped.reserve(maxPoints);
            const double step = static_cast<double>(reduced.size() - 1) / std::max(1, maxPoints - 1);
            for (int idx = 0; idx < maxPoints; ++idx)
            {
                capped.append(reduced[static_cast<int>(std::round(idx * step))]);
            }
            return capped;
        }
    } // namespace

    // AnalysisWorker implementation
    AnalysisWorker::AnalysisWorker(const CombinedData &data, const AnalysisRunConfig &config)
        : m_inputData(data), m_config(config), m_shouldStop(false)
    {
    }

    void AnalysisWorker::runAnalysis()
    {
        try
        {
            QElapsedTimer analysisTimer;
            analysisTimer.start();
            emit progressChanged(0, "Starting StatMaker analysis...");

            CombinedData filteredData = m_inputData;
            if (m_config.ppmThreshold > 0.0)
            {
                QVector<QString> filteredGeneNames;
                QVector<QVector<double>> filteredCountMatrix;
                QVector<QVector<int>> filteredRawCountMatrix;

                for (int geneIdx = 0; geneIdx < filteredData.geneNames.size(); ++geneIdx)
                {
                    double maxPpm = 0.0;
                    for (double value : filteredData.countMatrix[geneIdx])
                        maxPpm = std::max(maxPpm, value);

                    if (maxPpm >= m_config.ppmThreshold)
                    {
                        filteredGeneNames.append(filteredData.geneNames[geneIdx]);
                        filteredCountMatrix.append(filteredData.countMatrix[geneIdx]);
                        filteredRawCountMatrix.append(filteredData.rawCountMatrix[geneIdx]);
                    }
                }

                emit debugMessage(QString("PPM pre-filter kept %1 of %2 genes at threshold %3")
                                      .arg(filteredGeneNames.size())
                                      .arg(filteredData.geneNames.size())
                                      .arg(m_config.ppmThreshold));

                filteredData.geneNames = filteredGeneNames;
                filteredData.countMatrix = filteredCountMatrix;
                filteredData.rawCountMatrix = filteredRawCountMatrix;
            }

            if (filteredData.geneNames.isEmpty())
                throw std::runtime_error("No genes remain after PPM filtering.");

            // Convert Qt data to Eigen format
            auto [counts, metadata] = convertToEigenFormat(filteredData);

            emit progressChanged(10, "Data conversion completed");

            // Check if we should stop
            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 1: Create DeseqDataSet
            emit progressChanged(15, "Creating DeseqDataSet...");
            DeseqDataSet dds(counts, metadata, "~condition", true);

            // Step 2: Fit size factors
            emit progressChanged(25, "Fitting size factors...");
            dds.fitSizeFactors();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                    return;
            }

            // Step 3: Fit dispersions
            emit progressChanged(35, "Fitting genewise dispersions...");
            dds.fitGenewiseDispersions();

            emit progressChanged(45, "Fitting dispersion trend...");
            dds.fitDispersionTrend();

            emit progressChanged(55, "Fitting dispersion prior...");
            dds.fitDispersionPrior();

            emit progressChanged(65, "Fitting MAP dispersions...");
            dds.fitMAPDispersions();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 4: Fit log fold changes
            emit progressChanged(70, "Fitting log fold changes...");
            dds.fitLFC();

            // Step 5: Handle outliers
            emit progressChanged(75, "Calculating Cooks distances...");
            dds.calculateCooks();
            dds.refit();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 6: Statistical testing
            emit progressChanged(80, "Running statistical tests...");

            // Determine the number of unique conditions
            QSet<QString> uniqueGroups;
            for (const auto &group : filteredData.groupNames)
                uniqueGroups.insert(group);
            QList<QString> groupList = uniqueGroups.values();
            std::sort(groupList.begin(), groupList.end());
            int numGroups = groupList.size();
            const QString primaryContrastLabel =
                numGroups >= 2 ? QString("%1 vs %2").arg(groupList[1], groupList[0])
                               : QStringLiteral("analysis");

            std::vector<std::string> geneNamesStd;
            geneNamesStd.reserve(filteredData.geneNames.size());
            for (const auto &gene : filteredData.geneNames)
                geneNamesStd.push_back(gene.toStdString());

            for (size_t geneIdx = 0; geneIdx < geneNamesStd.size(); ++geneIdx)
            {
                if (geneNameHasInvalidBytes(geneNamesStd[geneIdx]))
                {
                    throw std::runtime_error(
                        QString("Invalid control byte detected in gene name at row %1")
                            .arg(static_cast<qulonglong>(geneIdx))
                            .toStdString());
                }
            }

            // For the primary contrast (condition 1 vs condition 0)
            Eigen::VectorXd contrast(numGroups);
            contrast.setZero();
            contrast(0) = 0.0;
            if (numGroups >= 2)
                contrast(1) = 1.0;

            DeseqStats ds(dds, contrast, m_config.pValueThreshold, true, true);

            emit progressChanged(85, "Running Wald test...");
            ds.runWaldTest();

            emit progressChanged(90, "Applying Cook's filtering...");
            ds.cooksFiltering();

            emit progressChanged(92, "Applying independent filtering...");
            ds.independentFiltering();

            // Step 7: Generate results
            emit progressChanged(94, "Generating results...");
            Eigen::MatrixXd results = ds.summary();

            // For three-way comparisons, run additional contrasts
            std::vector<Eigen::MatrixXd> allContrastResults;
            std::vector<std::string> allContrastLabels;
            std::vector<PairwiseContrast> pairwiseBaitContrasts;

            // Store first contrast result
            allContrastResults.push_back(results);
            if (numGroups >= 2)
                allContrastLabels.push_back(
                    groupList[1].toStdString() + " vs " + groupList[0].toStdString());

            if (numGroups >= 3)
            {
                // Run additional contrasts for each non-reference group vs reference
                for (int g = 2; g < numGroups; ++g)
                {
                    emit progressChanged(94 + g - 1, QString("Running contrast %1 vs %2...")
                                                         .arg(groupList[g])
                                                         .arg(groupList[0]));

                    Eigen::VectorXd extraContrast(numGroups);
                    extraContrast.setZero();
                    extraContrast(g) = 1.0;

                    DeseqStats extraDs(dds, extraContrast, m_config.pValueThreshold, true, true);
                    extraDs.runWaldTest();
                    extraDs.cooksFiltering();
                    extraDs.independentFiltering();

                    allContrastResults.push_back(extraDs.summary());
                    allContrastLabels.push_back(
                        groupList[g].toStdString() + " vs " + groupList[0].toStdString());
                }

                // Also run pairwise contrasts between non-reference groups
                for (int g = 2; g < numGroups; ++g)
                {
                    for (int h = 1; h < g; ++h)
                    {
                        emit progressChanged(95, QString("Running contrast %1 vs %2...")
                                                     .arg(groupList[g])
                                                     .arg(groupList[h]));

                        Eigen::VectorXd pairContrast(numGroups);
                        pairContrast.setZero();
                        pairContrast(g) = 1.0;
                        pairContrast(h) = -1.0;

                        DeseqStats pairDs(dds, pairContrast, m_config.pValueThreshold, true, true);
                        pairDs.runWaldTest();
                        pairDs.cooksFiltering();
                        pairDs.independentFiltering();

                        Eigen::MatrixXd pairSummary = pairDs.summary();
                        allContrastResults.push_back(pairSummary);
                        allContrastLabels.push_back(
                            groupList[g].toStdString() + " vs " + groupList[h].toStdString());

                        PairwiseContrast pairwiseContrast;
                        pairwiseContrast.bait_numerator = groupList[g].toStdString();
                        pairwiseContrast.bait_denominator = groupList[h].toStdString();
                        pairwiseContrast.results = pairSummary;
                        pairwiseContrast.gene_names = geneNamesStd;
                        pairwiseBaitContrasts.push_back(pairwiseContrast);
                    }
                }
            }

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Prepare results structure
            AnalysisResults analysisResults;
            analysisResults.results = results;

            // Store multi-contrast results for three-way comparisons
            analysisResults.contrastResults = allContrastResults;
            analysisResults.contrastLabels = allContrastLabels;
            analysisResults.activeContrast = 0;
            analysisResults.activeContrastLabel = primaryContrastLabel;
            analysisResults.createdAt = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);

            // Store dispersion data for visualization
            analysisResults.genewiseDispersions = dds.getGenewiseDispersions();
            analysisResults.fittedDispersions = dds.getFittedDispersions();
            analysisResults.baseMeans = dds.getNormedMeans();

            // Convert gene names and sample names to std::vector<std::string>
            analysisResults.geneNames = geneNamesStd;

            analysisResults.sampleNames.reserve(filteredData.sampleNames.size());
            for (const auto &sample : filteredData.sampleNames)
            {
                analysisResults.sampleNames.push_back(sample.toStdString());
            }

            // Y2H-SCORES: Enrichment scoring
            emit progressChanged(96, "Computing Y2H enrichment scores...");
            QElapsedTimer enrichmentTimer;
            enrichmentTimer.start();
            deseq2::EnrichmentScorer enrichmentScorer;
            QStringList groupNameList;
            for (const QString &groupName : filteredData.groupNames)
                groupNameList << groupName;
            const std::string baitName = databaseBaitName(groupNameList).toStdString();
            auto enrichmentResults = enrichmentScorer.compute(
                results,
                geneNamesStd,
                baitName,
                m_config.enrichmentPValueThreshold,
                m_config.enrichmentFoldChangeThreshold);
            emit debugMessage(QString("Y2H enrichment completed in %1 ms (%2 genes scored)")
                                  .arg(enrichmentTimer.elapsed())
                                  .arg(enrichmentResults.size()));

            // Y2H-SCORES: Specificity scoring where multi-bait contrasts exist
            std::vector<deseq2::SpecificityResult> specificityResults;
            if (!pairwiseBaitContrasts.empty())
            {
                emit progressChanged(97, "Computing Y2H specificity scores...");
                QElapsedTimer specificityTimer;
                specificityTimer.start();
                deseq2::SpecificityScorer specificityScorer;
                specificityResults = specificityScorer.compute(
                    pairwiseBaitContrasts,
                    std::max(1, numGroups - 1),
                    m_config.specificityPValueThreshold,
                    m_config.enrichmentFoldChangeThreshold);
                emit debugMessage(QString("Y2H specificity completed in %1 ms (%2 gene-bait scores)")
                                      .arg(specificityTimer.elapsed())
                                      .arg(specificityResults.size()));
            }
            else
            {
                emit debugMessage("Y2H specificity skipped: no pairwise bait contrasts available.");
            }

            // Y2H-SCORES: In-frame scoring from junction SQLite inputs when selected/non-selected files are available
            std::vector<deseq2::InFrameResult> inFrameResults;
            if (!m_config.junctionFiles.isEmpty())
            {
                struct JunctionAccumulator
                {
                    std::string transcriptId;
                    std::string gene;
                    int inFrameSelected = 0;
                    int totalSelected = 0;
                    int inFrameNonSelected = 0;
                    int totalNonSelected = 0;
                };

                std::map<std::string, JunctionAccumulator> accumulators;
                int matchedJunctionFiles = 0;

                for (const QString &junctionFile : m_config.junctionFiles)
                {
                    const QString junctionSample = QFileInfo(junctionFile).baseName();
                    QString sampleGroup;

                    for (int sampleIdx = 0; sampleIdx < filteredData.sampleNames.size(); ++sampleIdx)
                    {
                        const QString sampleName = filteredData.sampleNames[sampleIdx];
                        if (junctionSample.compare(sampleName, Qt::CaseInsensitive) == 0 ||
                            junctionSample.startsWith(sampleName + ".", Qt::CaseInsensitive) ||
                            sampleName.startsWith(junctionSample + ".", Qt::CaseInsensitive))
                        {
                            sampleGroup = filteredData.groupNames[sampleIdx];
                            break;
                        }
                    }

                    if (sampleGroup.isEmpty())
                    {
                        sampleGroup = inferGroupNameForSample(junctionSample);
                    }

                    const QString lowerGroup = sampleGroup.toLower();
                    const bool isSelected =
                        lowerGroup.contains("selected") && !lowerGroup.contains("non");
                    const bool isNonSelected =
                        lowerGroup.contains("non-selected") || lowerGroup.contains("nonselected") ||
                        lowerGroup.contains("non_selection") || lowerGroup.contains("vector") ||
                        lowerGroup == "control";

                    if (!isSelected && !isNonSelected)
                        continue;

                    QString connName = QString("statmaker_if_%1").arg(QUuid::createUuid().toString(QUuid::Id128));
                    {
                        QSqlDatabase junctionDb = QSqlDatabase::addDatabase("QSQLITE", connName);
                        junctionDb.setDatabaseName(junctionFile);
                        junctionDb.setConnectOptions("QSQLITE_OPEN_READONLY");
                        if (!junctionDb.open())
                        {
                            emit debugMessage(QString("Failed to open junction database %1: %2")
                                                  .arg(junctionFile, junctionDb.lastError().text()));
                        }
                        else
                        {
                            QSqlQuery junctionQuery(junctionDb);
                            const QString sql =
                                "SELECT gene, COALESCE(NULLIF(refseq, ''), gene) AS transcript_id, "
                                "COUNT(DISTINCT CASE WHEN frame = '+0_frame' THEN read END) AS in_frame_reads, "
                                "COUNT(DISTINCT read) AS total_reads "
                                "FROM maps GROUP BY gene, transcript_id";
                            if (junctionQuery.exec(sql))
                            {
                                matchedJunctionFiles++;
                                while (junctionQuery.next())
                                {
                                    const QString gene = junctionQuery.value(0).toString().trimmed();
                                    const QString transcriptId = junctionQuery.value(1).toString().trimmed();
                                    const int inFrameReads = junctionQuery.value(2).toInt();
                                    const int totalReads = junctionQuery.value(3).toInt();
                                    if (gene.isEmpty() || totalReads <= 0)
                                        continue;

                                    const std::string key =
                                        gene.toStdString() + "\t" + transcriptId.toStdString();
                                    auto &acc = accumulators[key];
                                    acc.gene = gene.toStdString();
                                    acc.transcriptId = transcriptId.toStdString();
                                    if (isSelected)
                                    {
                                        acc.inFrameSelected += inFrameReads;
                                        acc.totalSelected += totalReads;
                                    }
                                    else if (isNonSelected)
                                    {
                                        acc.inFrameNonSelected += inFrameReads;
                                        acc.totalNonSelected += totalReads;
                                    }
                                }
                            }
                            else
                            {
                                emit debugMessage(QString("Failed to query junction database %1: %2")
                                                      .arg(junctionFile, junctionQuery.lastError().text()));
                            }
                            junctionDb.close();
                        }
                    }
                    QSqlDatabase::removeDatabase(connName);
                }

                if (!accumulators.empty())
                {
                    emit progressChanged(98, "Computing Y2H in-frame scores...");
                    QElapsedTimer inFrameTimer;
                    inFrameTimer.start();
                    std::vector<deseq2::JunctionData> junctionData;
                    junctionData.reserve(accumulators.size());
                    for (const auto &[key, acc] : accumulators)
                    {
                        Q_UNUSED(key);
                        deseq2::JunctionData row;
                        row.transcript_id = acc.transcriptId;
                        row.gene = acc.gene;
                        row.in_frame_reads_selected = acc.inFrameSelected;
                        row.total_reads_selected = acc.totalSelected;
                        row.in_frame_reads_nonselected = acc.inFrameNonSelected;
                        row.total_reads_nonselected = acc.totalNonSelected;
                        junctionData.push_back(row);
                    }

                    deseq2::InFrameScorer inFrameScorer;
                    inFrameResults = inFrameScorer.compute(junctionData, baitName);
                    emit debugMessage(QString("Y2H in-frame completed in %1 ms (%2 genes scored from %3 junction files)")
                                          .arg(inFrameTimer.elapsed())
                                          .arg(inFrameResults.size())
                                          .arg(matchedJunctionFiles));
                }
                else
                {
                    emit debugMessage("Y2H in-frame skipped: no selected/non-selected junction datasets matched the loaded samples.");
                }
            }
            else
            {
                emit debugMessage("Y2H in-frame skipped: no junction databases configured.");
            }

            // Y2H-SCORES: Borda aggregation across all available metrics
            emit progressChanged(99, "Computing Y2H Borda aggregation...");
            QElapsedTimer bordaTimer;
            bordaTimer.start();
            deseq2::BordaAggregator borda;
            auto y2hResults = borda.aggregate(enrichmentResults, specificityResults, inFrameResults);
            emit debugMessage(QString("Y2H Borda aggregation completed in %1 ms (%2 results)")
                                  .arg(bordaTimer.elapsed())
                                  .arg(y2hResults.size()));

            analysisResults.enrichmentScores = enrichmentResults;
            analysisResults.specificityScores = specificityResults;
            analysisResults.inFrameScores = inFrameResults;
            analysisResults.y2hScores = y2hResults;

            // Calculate summary statistics
            analysisResults.totalGenes = results.rows();
            analysisResults.significantGenes = 0;
            analysisResults.upregulatedGenes = 0;
            analysisResults.downregulatedGenes = 0;
            analysisResults.pValueThreshold = m_config.pValueThreshold;
            analysisResults.log2FCThreshold = 0.0; // Default, can be customized later

            for (int i = 0; i < results.rows(); ++i)
            {
                double padj = results(i, 5); // Adjusted p-value column (FIXED: was 4, should be 5)
                double lfc = results(i, 1);  // Log2 fold change column

                if (!std::isnan(padj) && padj < m_config.pValueThreshold)
                {
                    analysisResults.significantGenes++;
                    if (lfc > 0)
                    {
                        analysisResults.upregulatedGenes++;
                    }
                    else if (lfc < 0)
                    {
                        analysisResults.downregulatedGenes++;
                    }
                }
            }

            analysisResults.isValid = true;
            analysisResults.errorMessage = "";

            emit progressChanged(100, "Analysis completed successfully!");
            emit debugMessage(QString("Analysis completed in %1 ms").arg(analysisTimer.elapsed()));
            emit analysisFinished(analysisResults);
            emit finished(); // Signal that worker is completely done
        }
        catch (const std::exception &e)
        {
            emit analysisError(QString("Analysis failed: %1").arg(e.what()));
            emit finished(); // Signal that worker is done even on error
        }
        catch (...)
        {
            emit analysisError("Analysis failed: Unknown error occurred");
            emit finished(); // Signal that worker is done even on error
        }
    }

    void AnalysisWorker::stopAnalysis()
    {
        QMutexLocker locker(&m_stopMutex);
        m_shouldStop = true;
    }

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> AnalysisWorker::convertToEigenFormat(const CombinedData &data)
    {
        // Convert count matrix (Qt format: genes x samples) to Eigen format (samples x genes)
        int numSamples = data.sampleNames.size();
        int numGenes = data.geneNames.size();

        Eigen::MatrixXd counts(numSamples, numGenes);

        for (int gene = 0; gene < numGenes; ++gene)
        {
            for (int sample = 0; sample < numSamples; ++sample)
            {
                // Use raw integer counts directly (from xlsx or back-computed from PPM)
                int integerCount = data.rawCountMatrix[gene][sample];
                counts(sample, gene) = static_cast<double>(integerCount);
            }
        }

        // Create metadata matrix based on group names
        Eigen::MatrixXd metadata(numSamples, 1);

        // Get unique group names to assign numeric codes
        QSet<QString> uniqueGroups;
        for (const auto &group : data.groupNames)
        {
            uniqueGroups.insert(group);
        }

        QList<QString> groupList = uniqueGroups.values();
        std::sort(groupList.begin(), groupList.end());

        // Assign numeric codes to groups
        for (int sample = 0; sample < numSamples; ++sample)
        {
            QString groupName = data.groupNames[sample];
            int groupIndex = groupList.indexOf(groupName);
            metadata(sample, 0) = static_cast<double>(groupIndex);
        }

        return std::make_tuple(counts, metadata);
    }

    MainWindow::MainWindow(QWidget *parent)
        : QMainWindow(parent),
          m_geneCountHandler(std::make_unique<GeneCountHandler>()),
          m_combinedData(),
          // Initialize all pointer members to nullptr
          m_tabWidget(nullptr),
          m_inputTab(nullptr),
          m_resultsTab(nullptr),
          m_visualizationTab(nullptr),
          m_analysisTab(nullptr),
          m_addPpmFilesButton(nullptr),
          m_clearFilesButton(nullptr),
          m_generateDataButton(nullptr),
          m_exportGeneratedDataButton(nullptr),
          m_geneCountFilesTable(nullptr),
          m_groupNameEdit(nullptr),
          m_assignGroupButton(nullptr),
          m_autoAssignButton(nullptr),
          m_countsPreviewTable(nullptr),
          m_metadataPreviewTable(nullptr),
          m_runAnalysisButton(nullptr),
          m_stopAnalysisButton(nullptr),
          m_resetAnalysisButton(nullptr),
          m_analysisProgressBar(nullptr),
          m_totalGenesLabel(nullptr),
          m_significantGenesLabel(nullptr),
          m_upregulatedLabel(nullptr),
          m_downregulatedLabel(nullptr),
          m_resultsTable(nullptr),
          m_pValueThresholdSpinBox(nullptr),
          m_log2FCThresholdSpinBox(nullptr),
          m_applyFilterButton(nullptr),
          m_clearFilterButton(nullptr),
          m_exportResultsButton(nullptr),
          m_saveGeneListsButton(nullptr),
          m_plotTypeCombo(nullptr),
          m_colorSchemeCombo(nullptr),
          m_pointSizeSpinBox(nullptr),
          m_plotWidget(nullptr),
          m_savePlotButton(nullptr),
          m_exportDataButton(nullptr),
          m_printPlotButton(nullptr),
          m_progressOutputText(nullptr),
          m_clearOutputButton(nullptr),
          m_saveOutputButton(nullptr),
          m_actionOpenCounts(nullptr),
          m_actionOpenMetadata(nullptr),
          m_actionSaveResults(nullptr),
          m_actionExportVisualizations(nullptr),
          m_actionExit(nullptr),
          m_actionRunAnalysis(nullptr),
          m_actionStopAnalysis(nullptr),
          m_actionResetAnalysis(nullptr),
          m_actionInputTab(nullptr),
          m_actionAnalysisTab(nullptr),
          m_actionResultsTab(nullptr),
          m_actionVisualizationTab(nullptr),
          m_actionAbout(nullptr),
          m_actionUserGuide(nullptr),
          m_baitGroupingCombo(nullptr),
          m_enrichPvalSpinBox(nullptr),
          m_enrichFcSpinBox(nullptr),
          m_specPvalSpinBox(nullptr),
          m_junctionFilesEdit(nullptr),
          m_browseJunctionBtn(nullptr),
          m_contrastSelectorCombo(nullptr)
    {
        // Initialize last used directory to Documents folder
        m_lastUsedDirectory = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
        // Load saved last used directory from settings
        loadLastUsedDirectory();
        qRegisterMetaType<deseq2::AnalysisResults>("deseq2::AnalysisResults");
        // Setup UI from the UI file
        m_ui = std::make_unique<Ui::MainWindow>();
        m_ui->setupUi(this);
        // Connect UI elements to member variables
        connectUiElements();
        setupConnections();
        setupMenusAndToolbar();
        initializeTables();
        updateUiState();
        addProgressMessage("StatMaker++ started successfully.");

        // Initialize analysis results
        m_analysisResults.isValid = false;
    }

    MainWindow::~MainWindow()
    {
        // Stop analysis if running
        if (m_analysisThread && m_analysisThread->isRunning())
        {
            if (m_analysisWorker)
            {
                m_analysisWorker->stopAnalysis();
            }
            m_analysisThread->quit();
            m_analysisThread->wait(5000); // Wait up to 5 seconds
        }
    }

    void MainWindow::loadWorkingDirectory(const QString &workdir)
    {
        if (workdir.trimmed().isEmpty())
            return;

        QDir dir(QFileInfo(workdir).absoluteFilePath());
        if (!dir.exists())
            dir.mkpath(".");

        m_workingDirectory = dir.absolutePath();
        m_lastUsedDirectory = m_workingDirectory;
        saveLastUsedDirectory();

        QDir analyzedDir(dir.filePath("analyzed_files"));
        if (!analyzedDir.exists())
            dir.mkpath("analyzed_files");

        addProgressMessage(QString("Loading working directory: %1").arg(m_workingDirectory));

        populateAutoDiscoveredFiles();
        autoDetectJunctionFiles();

        const QString sqlitePath = resultsDatabasePathForWorkdir(m_workingDirectory);
        if (QFile::exists(sqlitePath))
        {
            loadResultsFromSqlite(sqlitePath);
        }
    }

    void MainWindow::loadInputFiles(const QStringList &filePaths, bool clearExisting)
    {
        if (clearExisting)
        {
            m_geneCountHandler->clearAllFiles();
            if (m_geneCountFilesTable)
                m_geneCountFilesTable->setRowCount(0);
            if (m_countsPreviewTable)
            {
                m_countsPreviewTable->clear();
                m_countsPreviewTable->setRowCount(0);
                m_countsPreviewTable->setColumnCount(0);
            }
            if (m_metadataPreviewTable)
            {
                m_metadataPreviewTable->clear();
                m_metadataPreviewTable->setRowCount(0);
                m_metadataPreviewTable->setColumnCount(0);
            }
            m_combinedData = CombinedData();
        }

        int addedCount = 0;
        QStringList normalizedFiles;
        for (const QString &path : filePaths)
        {
            const QString absolutePath = QFileInfo(path.trimmed()).absoluteFilePath();
            if (absolutePath.isEmpty() || !QFile::exists(absolutePath))
                continue;
            if (QFileInfo(absolutePath).fileName().compare("statmaker_results.sqlite", Qt::CaseInsensitive) == 0)
                continue;
            normalizedFiles << absolutePath;
        }
        normalizedFiles.removeDuplicates();
        std::sort(normalizedFiles.begin(), normalizedFiles.end());

        if (!normalizedFiles.isEmpty())
            updateLastUsedDirectory(normalizedFiles.first());

        for (const QString &absolutePath : normalizedFiles)
        {
            if (m_geneCountHandler->addGeneCountFile(absolutePath))
                addedCount++;
        }

        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);
        m_resultsColumnsSized = false;
        invalidatePlotCache();
        updateUiState();

        if (addedCount == 0)
            return;

        addProgressMessage(QString("Loaded %1 StatMaker input file(s).").arg(addedCount));

        QMap<QString, QString> patterns;
        patterns["NON"] = "Non-Selected";
        patterns["non-selected"] = "Non-Selected";
        patterns["nonselected"] = "Non-Selected";
        patterns["_NON_"] = "Non-Selected";
        patterns["SEL"] = "Selected";
        patterns["selected"] = "Selected";
        patterns["_SEL_"] = "Selected";
        patterns["vector"] = "Vector Control";
        patterns["control"] = "Vector Control";
        patterns["bait"] = "Bait";
        const int autoAssigned = m_geneCountHandler->autoAssignGroups(patterns);
        if (autoAssigned > 0)
        {
            addProgressMessage(QString("Auto-assigned groups for %1 input file(s).").arg(autoAssigned));
            m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);
        }

        bool allAssigned = true;
        for (const auto &data : m_geneCountHandler->getGeneCountFiles())
        {
            if (data.groupName.isEmpty())
            {
                allAssigned = false;
                break;
            }
        }

        if (allAssigned && !m_geneCountHandler->getGeneCountFiles().isEmpty())
        {
            m_combinedData = m_geneCountHandler->generateCountMatrixAndMetadata();
            if (m_combinedData.isValid)
            {
                updateCountsPreviewWithConvertedValues(m_countsPreviewTable, m_combinedData);
                m_geneCountHandler->updateMetadataPreviewTable(m_metadataPreviewTable, m_combinedData);
                updateUiState();
                addProgressMessage(QString("Prepared count matrix automatically: %1 genes x %2 samples.")
                                       .arg(m_combinedData.geneNames.size())
                                       .arg(m_combinedData.sampleNames.size()));
            }
            else
            {
                addProgressMessage(QString("Failed to prepare count matrix automatically: %1")
                                       .arg(m_combinedData.errorMessage));
            }
        }

        if (m_workingDirectory.isEmpty())
            m_workingDirectory = resolveWorkingDirectory();
    }

    void MainWindow::populateAutoDiscoveredFiles()
    {
        if (m_workingDirectory.isEmpty())
            return;

        QStringList sqliteInputs = discoverFiles(
            QDir(m_workingDirectory).filePath("analyzed_files"),
            {"*.sqlite", "*.db"});
        sqliteInputs.erase(
            std::remove_if(sqliteInputs.begin(), sqliteInputs.end(),
                           [](const QString &path)
                           {
                               return QFileInfo(path).fileName().compare(
                                          "statmaker_results.sqlite", Qt::CaseInsensitive) == 0;
                           }),
            sqliteInputs.end());

        QStringList csvInputs;
        if (sqliteInputs.isEmpty())
        {
            csvInputs = discoverFiles(
                QDir(m_workingDirectory).filePath("gene_count_summary"),
                {"*.csv"});
        }

        const QStringList inputFiles = !sqliteInputs.isEmpty() ? sqliteInputs : csvInputs;
        if (inputFiles.isEmpty())
        {
            addProgressMessage("No GeneCount inputs were discovered in the working directory.");
            return;
        }

        loadInputFiles(inputFiles, true);
    }

    void MainWindow::autoDetectJunctionFiles()
    {
        if (m_workingDirectory.isEmpty() || !m_junctionFilesEdit)
            return;

        QStringList junctionFiles = discoverFiles(
            QDir(m_workingDirectory).filePath("analyzed_files"),
            {"*.sqlite", "*.db"});
        junctionFiles.erase(
            std::remove_if(junctionFiles.begin(), junctionFiles.end(),
                           [](const QString &path)
                           {
                               return QFileInfo(path).fileName().compare(
                                          "statmaker_results.sqlite", Qt::CaseInsensitive) == 0;
                           }),
            junctionFiles.end());
        junctionFiles.removeDuplicates();
        std::sort(junctionFiles.begin(), junctionFiles.end());

        if (!junctionFiles.isEmpty())
        {
            m_junctionFilesEdit->setText(joinPaths(junctionFiles));
            addProgressMessage(QString("Auto-detected %1 junction SQLite file(s).").arg(junctionFiles.size()));
        }
    }

    QString MainWindow::resolveWorkingDirectory() const
    {
        if (!m_workingDirectory.isEmpty())
            return m_workingDirectory;

        const auto &files = m_geneCountHandler->getGeneCountFiles();
        if (files.isEmpty())
            return m_lastUsedDirectory;

        QDir dir(QFileInfo(files.first().filePath).absolutePath());
        const QString dirName = dir.dirName();
        if (dirName == "gene_count_summary" || dirName == "analyzed_files")
            dir.cdUp();
        return dir.absolutePath();
    }

    QString MainWindow::resultsDatabasePathForWorkdir(const QString &workdir) const
    {
        return QDir(workdir).filePath("analyzed_files/statmaker_results.sqlite");
    }

    void MainWindow::invalidatePlotCache()
    {
        m_lastPlotCacheKey.clear();
        m_resultsRevision++;
    }

    void MainWindow::connectUiElements()
    {
        // Connect UI elements from the UI file to member variables
        m_tabWidget = m_ui->tabWidget;
        m_inputTab = m_ui->inputTab;
        m_resultsTab = m_ui->resultsTab;
        m_visualizationTab = m_ui->visualizationTab;
        m_analysisTab = nullptr; // Not present in UI file
        // Analysis settings
        m_analysisModeCombo = m_ui->analysisModeCombo;
        m_ppmThresholdSpinBox_y2h = m_ui->ppmThresholdSpinBox;
        m_groupTypeCombo = m_ui->groupTypeCombo;

        // Input tab elements
        m_addPpmFilesButton = m_ui->addPpmFilesButton;
        m_clearFilesButton = m_ui->clearFilesButton;
        m_generateDataButton = m_ui->generateDataButton;
        m_exportGeneratedDataButton = m_ui->exportGeneratedDataButton;
        m_geneCountFilesTable = m_ui->geneCountFilesTable;
        m_groupNameEdit = m_ui->groupNameEdit;
        m_assignGroupButton = m_ui->assignGroupButton;
        m_autoAssignButton = m_ui->autoAssignButton;
        m_countsPreviewTable = m_ui->countsPreviewTable;
        m_metadataPreviewTable = m_ui->metadataPreviewTable;

        // Analysis tab elements - connect to UI elements if they exist
        m_runAnalysisButton = m_ui->runAnalysisButton;
        m_stopAnalysisButton = m_ui->stopAnalysisButton;
        m_resetAnalysisButton = m_ui->resetAnalysisButton;

        // Add progress bar if not already in UI
        if (!m_analysisProgressBar)
        {
            m_analysisProgressBar = new QProgressBar(this);
            m_analysisProgressBar->setVisible(false);
            statusBar()->addPermanentWidget(m_analysisProgressBar);
        }

        // Results tab elements
        m_totalGenesLabel = m_ui->totalGenesLabel;
        m_significantGenesLabel = m_ui->significantGenesLabel;
        m_upregulatedLabel = m_ui->upregulatedLabel;
        m_downregulatedLabel = m_ui->downregulatedLabel;
        m_resultsTable = m_ui->resultsTable;
        m_pValueThresholdSpinBox = m_ui->pValueThresholdSpinBox;
        m_log2FCThresholdSpinBox = m_ui->log2FCThresholdSpinBox;
        m_applyFilterButton = m_ui->applyFilterButton;
        m_clearFilterButton = m_ui->clearFilterButton;
        m_exportResultsButton = m_ui->exportResultsButton;
        m_saveGeneListsButton = m_ui->saveGeneListsButton;
        // Visualization tab elements
        m_plotTypeCombo = m_ui->plotTypeCombo;
        m_colorSchemeCombo = m_ui->colorSchemeCombo;
        m_pointSizeSpinBox = m_ui->pointSizeSpinBox;
        // Replace the plain QWidget from UI with a QChartView (Qt Charts)
        {
            QWidget *placeholder = m_ui->plotWidget;
            m_plotWidget = new QChartView(placeholder->parentWidget());
            m_plotWidget->setMinimumSize(placeholder->minimumSize());
            m_plotWidget->setSizePolicy(placeholder->sizePolicy());
            m_plotWidget->setRenderHint(QPainter::Antialiasing);

            // Create an initial empty chart
            QChart *chart = new QChart();
            chart->setTitle("No data - run analysis first");
            m_plotWidget->setChart(chart);

            // Replace in the layout
            QLayout *parentLayout = placeholder->parentWidget()->layout();
            if (parentLayout)
            {
                QLayoutItem *item = parentLayout->replaceWidget(placeholder, m_plotWidget);
                delete item;
            }
            placeholder->hide();
            placeholder->deleteLater();
        }
        m_savePlotButton = m_ui->savePlotButton;
        m_exportDataButton = m_ui->exportDataButton;
        m_printPlotButton = m_ui->printPlotButton;
        // Progress output elements
        m_progressOutputText = m_ui->progressOutputText;
        m_clearOutputButton = m_ui->clearOutputButton;
        m_saveOutputButton = m_ui->saveOutputButton;

        // Y2H-SCORES Settings elements
        m_baitGroupingCombo = m_ui->baitGroupingCombo;
        m_enrichPvalSpinBox = m_ui->enrichPvalSpinBox;
        m_enrichFcSpinBox = m_ui->enrichFcSpinBox;
        m_specPvalSpinBox = m_ui->specPvalSpinBox;

        // In-frame scoring (junction files) elements
        m_junctionFilesEdit = m_ui->junctionFilesEdit;
        m_browseJunctionBtn = m_ui->browseJunctionBtn;

        // Contrast selector
        m_contrastSelectorCombo = m_ui->contrastSelectorCombo;
    }

    void MainWindow::setupConnections()
    {
        // File operations
        if (m_addPpmFilesButton)
        {
            connect(m_addPpmFilesButton, &QPushButton::clicked, this, &MainWindow::onAddGeneCountFiles);
        }
        if (m_clearFilesButton)
        {
            connect(m_clearFilesButton, &QPushButton::clicked, this, &MainWindow::onClearFiles);
        }
        if (m_generateDataButton)
        {
            connect(m_generateDataButton, &QPushButton::clicked, this, &MainWindow::onGenerateData);
        }
        if (m_exportGeneratedDataButton)
        {
            connect(m_exportGeneratedDataButton, &QPushButton::clicked, this, &MainWindow::onExportGeneratedData);
        }

        // Group assignment
        if (m_assignGroupButton)
        {
            connect(m_assignGroupButton, &QPushButton::clicked, this, &MainWindow::onAssignGroup);
        }
        if (m_autoAssignButton)
        {
            connect(m_autoAssignButton, &QPushButton::clicked, this, &MainWindow::onAutoAssignGroups);
        }

        // Analysis control
        if (m_runAnalysisButton)
        {
            connect(m_runAnalysisButton, &QPushButton::clicked, this, &MainWindow::onRunAnalysis);
        }
        if (m_stopAnalysisButton)
        {
            connect(m_stopAnalysisButton, &QPushButton::clicked, this, &MainWindow::onStopAnalysis);
        }
        if (m_resetAnalysisButton)
        {
            connect(m_resetAnalysisButton, &QPushButton::clicked, this, &MainWindow::onResetAnalysis);
        }

        // Results filtering
        if (m_applyFilterButton)
        {
            connect(m_applyFilterButton, &QPushButton::clicked, this, &MainWindow::onApplyFilter);
        }
        if (m_clearFilterButton)
        {
            connect(m_clearFilterButton, &QPushButton::clicked, this, &MainWindow::onClearFilter);
        }

        // Export operations
        if (m_exportResultsButton)
        {
            connect(m_exportResultsButton, &QPushButton::clicked, this, &MainWindow::onExportResults);
        }
        if (m_saveGeneListsButton)
        {
            connect(m_saveGeneListsButton, &QPushButton::clicked, this, &MainWindow::onSaveGeneLists);
        }

        // Visualization
        if (m_plotTypeCombo)
        {
            connect(m_plotTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                    this, &MainWindow::onPlotTypeChanged);
        }
        if (m_savePlotButton)
        {
            connect(m_savePlotButton, &QPushButton::clicked, this, &MainWindow::onSavePlot);
        }
        if (m_exportDataButton)
        {
            connect(m_exportDataButton, &QPushButton::clicked, this, &MainWindow::onExportPlotData);
        }
        if (m_printPlotButton)
        {
            connect(m_printPlotButton, &QPushButton::clicked, this, &MainWindow::onPrintPlot);
        }

        // Progress output
        if (m_clearOutputButton)
        {
            connect(m_clearOutputButton, &QPushButton::clicked, this, &MainWindow::onClearOutput);
        }
        if (m_saveOutputButton)
        {
            connect(m_saveOutputButton, &QPushButton::clicked, this, &MainWindow::onSaveOutput);
        }

        // Junction files browse button
        if (m_browseJunctionBtn)
        {
            connect(m_browseJunctionBtn, &QPushButton::clicked, this, &MainWindow::onBrowseJunctionFiles);
        }

        // Contrast selector
        if (m_contrastSelectorCombo)
        {
            connect(m_contrastSelectorCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                    this, &MainWindow::onContrastChanged);
        }
    }

    void MainWindow::setupMenusAndToolbar()
    {
        // Create menu bar
        QMenuBar *menuBar = this->menuBar();

        // File menu
        QMenu *fileMenu = menuBar->addMenu("&File");
        m_actionOpenCounts = fileMenu->addAction("Add Gene Count Files");
        m_actionOpenCounts->setShortcut(QKeySequence::Open);
        connect(m_actionOpenCounts, &QAction::triggered, this, &MainWindow::onOpenCounts);

        m_actionOpenMetadata = fileMenu->addAction("Generate Data");
        m_actionOpenMetadata->setShortcut(QKeySequence("Ctrl+G"));
        connect(m_actionOpenMetadata, &QAction::triggered, this, &MainWindow::onOpenMetadata);

        fileMenu->addSeparator();

        m_actionSaveResults = fileMenu->addAction("Save Results");
        m_actionSaveResults->setShortcut(QKeySequence::Save);
        connect(m_actionSaveResults, &QAction::triggered, this, &MainWindow::onSaveResults);

        m_actionExportVisualizations = fileMenu->addAction("Export Visualizations");
        m_actionExportVisualizations->setShortcut(QKeySequence("Ctrl+E"));
        connect(m_actionExportVisualizations, &QAction::triggered, this, &MainWindow::onExportVisualizations);

        fileMenu->addSeparator();

        m_actionExit = fileMenu->addAction("Exit");
        m_actionExit->setShortcut(QKeySequence::Quit);
        connect(m_actionExit, &QAction::triggered, this, &MainWindow::onExit);

        // Analysis menu
        QMenu *analysisMenu = menuBar->addMenu("&Analysis");
        m_actionRunAnalysis = analysisMenu->addAction("Run Analysis");
        m_actionRunAnalysis->setShortcut(QKeySequence("F5"));
        connect(m_actionRunAnalysis, &QAction::triggered, this, &MainWindow::onRunAnalysisAction);

        m_actionStopAnalysis = analysisMenu->addAction("Stop Analysis");
        m_actionStopAnalysis->setShortcut(QKeySequence("F6"));
        connect(m_actionStopAnalysis, &QAction::triggered, this, &MainWindow::onStopAnalysisAction);

        m_actionResetAnalysis = analysisMenu->addAction("Reset Analysis");
        m_actionResetAnalysis->setShortcut(QKeySequence("F7"));
        connect(m_actionResetAnalysis, &QAction::triggered, this, &MainWindow::onResetAnalysisAction);

        // View menu
        QMenu *viewMenu = menuBar->addMenu("&View");
        m_actionInputTab = viewMenu->addAction("Input Tab");
        m_actionInputTab->setShortcut(QKeySequence("Ctrl+1"));
        connect(m_actionInputTab, &QAction::triggered, this, &MainWindow::onInputTab);

        m_actionAnalysisTab = viewMenu->addAction("Analysis Tab");
        m_actionAnalysisTab->setShortcut(QKeySequence("Ctrl+2"));
        connect(m_actionAnalysisTab, &QAction::triggered, this, &MainWindow::onAnalysisTab);

        m_actionResultsTab = viewMenu->addAction("Results Tab");
        m_actionResultsTab->setShortcut(QKeySequence("Ctrl+3"));
        connect(m_actionResultsTab, &QAction::triggered, this, &MainWindow::onResultsTab);

        m_actionVisualizationTab = viewMenu->addAction("Visualization Tab");
        m_actionVisualizationTab->setShortcut(QKeySequence("Ctrl+4"));
        connect(m_actionVisualizationTab, &QAction::triggered, this, &MainWindow::onVisualizationTab);

        // Help menu
        QMenu *helpMenu = menuBar->addMenu("&Help");
        m_actionAbout = helpMenu->addAction("About");
        connect(m_actionAbout, &QAction::triggered, this, &MainWindow::onAbout);

        m_actionUserGuide = helpMenu->addAction("User Guide");
        m_actionUserGuide->setShortcut(QKeySequence::HelpContents);
        connect(m_actionUserGuide, &QAction::triggered, this, &MainWindow::onUserGuide);

        // Create toolbar
        QToolBar *toolBar = addToolBar("Main Toolbar");
        toolBar->addAction(m_actionOpenCounts);
        toolBar->addSeparator();
        toolBar->addAction(m_actionRunAnalysis);
        toolBar->addAction(m_actionStopAnalysis);
        toolBar->addSeparator();
        toolBar->addAction(m_actionSaveResults);
        toolBar->addAction(m_actionExportVisualizations);
    }

    void MainWindow::initializeTables()
    {
        // Initialize gene count files table
        if (m_geneCountFilesTable)
        {
            m_geneCountFilesTable->setColumnCount(5);
            QStringList headers = {"File Path", "Sample Name", "Group", "Total Reads", "Status"};
            m_geneCountFilesTable->setHorizontalHeaderLabels(headers);
            m_geneCountFilesTable->horizontalHeader()->setStretchLastSection(true);
            m_geneCountFilesTable->setAlternatingRowColors(true);
            m_geneCountFilesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_geneCountFilesTable->setSortingEnabled(true);
            m_geneCountFilesTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
        }

        // Initialize counts preview table
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->setAlternatingRowColors(true);
            m_countsPreviewTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_countsPreviewTable->setSortingEnabled(true);
        }

        // Initialize metadata preview table
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->setAlternatingRowColors(true);
            m_metadataPreviewTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_metadataPreviewTable->setSortingEnabled(true);
        }

        // Initialize results table
        if (m_resultsTable)
        {
            m_resultsTable->setAlternatingRowColors(true);
            m_resultsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_resultsTable->setSortingEnabled(true);
            m_resultsTable->setContextMenuPolicy(Qt::CustomContextMenu);
            connect(m_resultsTable, &QTableWidget::customContextMenuRequested,
                    this, &MainWindow::onResultsContextMenu);
        }
    }

    void MainWindow::updateUiState()
    {
        // Update button states based on current data
        bool hasFiles = !m_geneCountHandler->getGeneCountFiles().isEmpty();
        bool hasValidData = m_combinedData.isValid;

        if (m_generateDataButton)
        {
            m_generateDataButton->setEnabled(hasFiles);
        }
        if (m_exportGeneratedDataButton)
        {
            m_exportGeneratedDataButton->setEnabled(hasValidData);
        }

        // Update analysis buttons
        bool canRunAnalysis = hasValidData;
        if (m_runAnalysisButton)
        {
            m_runAnalysisButton->setEnabled(canRunAnalysis);
        }
        if (m_actionRunAnalysis)
        {
            m_actionRunAnalysis->setEnabled(canRunAnalysis);
        }

        // Update group assignment buttons
        if (m_assignGroupButton)
        {
            m_assignGroupButton->setEnabled(hasFiles);
        }
        if (m_autoAssignButton)
        {
            m_autoAssignButton->setEnabled(hasFiles);
        }
    }

    void MainWindow::addProgressMessage(const QString &message)
    {
        if (m_progressOutputText)
        {
            QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
            QString formattedMessage = QString("[%1] %2").arg(timestamp).arg(message);
            m_progressOutputText->append(formattedMessage);

            // Auto-scroll to bottom
            QTextCursor cursor = m_progressOutputText->textCursor();
            cursor.movePosition(QTextCursor::End);
            m_progressOutputText->setTextCursor(cursor);
        }
    }

    void MainWindow::updateLastUsedDirectory(const QString &filePath)
    {
        if (!filePath.isEmpty())
        {
            QFileInfo fileInfo(filePath);
            QString directory = fileInfo.absolutePath();
            if (!directory.isEmpty())
            {
                m_lastUsedDirectory = directory;
            }
        }
    }

    QString MainWindow::getLastUsedDirectory() const
    {
        return m_lastUsedDirectory;
    }

    void MainWindow::saveLastUsedDirectory()
    {
        QSettings settings;
        settings.setValue("LastUsedDirectory", m_lastUsedDirectory);
    }

    void MainWindow::loadLastUsedDirectory()
    {
        QSettings settings;
        QString savedDirectory = settings.value("LastUsedDirectory").toString();
        if (!savedDirectory.isEmpty() && QDir(savedDirectory).exists())
        {
            m_lastUsedDirectory = savedDirectory;
        }
    }

    void MainWindow::clearAllData()
    {
        // Clear tables
        if (m_geneCountFilesTable)
        {
            m_geneCountFilesTable->setRowCount(0);
        }
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->setRowCount(0);
            m_countsPreviewTable->setColumnCount(0);
        }
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->setRowCount(0);
            m_metadataPreviewTable->setColumnCount(0);
        }
        if (m_resultsTable)
        {
            m_resultsTable->setRowCount(0);
            m_resultsTable->setColumnCount(0);
        }

        // Reset labels
        if (m_totalGenesLabel)
        {
            m_totalGenesLabel->setText("Total Genes: 0");
        }
        if (m_significantGenesLabel)
        {
            m_significantGenesLabel->setText("Significant: 0");
        }
        if (m_upregulatedLabel)
        {
            m_upregulatedLabel->setText("Upregulated: 0");
        }
        if (m_downregulatedLabel)
        {
            m_downregulatedLabel->setText("Downregulated: 0");
        }

        // Clear group name
        if (m_groupNameEdit)
        {
            m_groupNameEdit->clear();
        }

        m_resultsColumnsSized = false;
        invalidatePlotCache();
        updateUiState();
        addProgressMessage("All data cleared.");
    }

    // Slot implementations - File operations
    void MainWindow::onAddGeneCountFiles()
    {
        QStringList fileNames = QFileDialog::getOpenFileNames(
            static_cast<QWidget *>(this),
            "Select Gene Count Files",
            m_lastUsedDirectory,
            "Gene Count Files (*.db *.sqlite *.csv);;SQLite Databases (*.db *.sqlite);;CSV Files (*.csv);;All Files (*)");

        if (!fileNames.isEmpty())
        {
            // Update last used directory from the first selected file
            updateLastUsedDirectory(fileNames.first());

            int addedCount = 0;
            for (const QString &fileName : fileNames)
            {
                if (m_geneCountHandler->addGeneCountFile(fileName))
                {
                    addedCount++;
                }
            }

            addProgressMessage(QString("Successfully added %1 out of %2 gene count files.").arg(addedCount).arg(fileNames.size()));

            // Update the gene count files table
            m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

            // Update UI state
            updateUiState();
        }
    }

    void MainWindow::onClearFiles()
    {
        m_geneCountHandler->clearAllFiles();
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        // Clear preview tables
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->clear();
            m_countsPreviewTable->setRowCount(0);
            m_countsPreviewTable->setColumnCount(0);
        }
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->clear();
            m_metadataPreviewTable->setRowCount(0);
            m_metadataPreviewTable->setColumnCount(0);
        }

        m_combinedData = CombinedData();
        updateUiState();
        addProgressMessage("All gene count files cleared.");
    }

    void MainWindow::onGenerateData()
    {
        addProgressMessage("Generating count matrix and metadata...");

        // Check if we have files loaded
        if (m_geneCountHandler->getGeneCountFiles().isEmpty())
        {
            addProgressMessage("No gene count files loaded. Please add files first.");
            return;
        }

        // Check if all files have group assignments
        bool allAssigned = true;
        for (const auto &data : m_geneCountHandler->getGeneCountFiles())
        {
            if (data.groupName.isEmpty())
            {
                allAssigned = false;
                break;
            }
        }

        if (!allAssigned)
        {
            addProgressMessage("Some files do not have group assignments. Please assign groups to all files.");
            return;
        }

        // Generate count matrix and metadata
        m_combinedData = m_geneCountHandler->generateCountMatrixAndMetadata();

        if (m_combinedData.isValid)
        {
            if (m_workingDirectory.isEmpty())
                m_workingDirectory = resolveWorkingDirectory();

            // Update preview tables with converted integer counts for count matrix
            updateCountsPreviewWithConvertedValues(m_countsPreviewTable, m_combinedData);
            m_geneCountHandler->updateMetadataPreviewTable(m_metadataPreviewTable, m_combinedData);

            addProgressMessage(QString("Successfully generated count matrix with %1 genes and %2 samples.").arg(m_combinedData.geneNames.size()).arg(m_combinedData.sampleNames.size()));
            updateUiState();
        }
        else
        {
            addProgressMessage(QString("Failed to generate data: %1").arg(m_combinedData.errorMessage));
        }
    }

    void MainWindow::onExportGeneratedData()
    {
        if (!m_combinedData.isValid)
        {
            addProgressMessage("No data to export. Please generate data first.");
            return;
        }

        QString baseFileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Export Generated Data",
            m_lastUsedDirectory,
            "CSV Files (*.csv);;All Files (*)");

        if (!baseFileName.isEmpty())
        {
            // Update last used directory
            updateLastUsedDirectory(baseFileName);

            // Remove extension if present
            QFileInfo fileInfo(baseFileName);
            QString basePath = fileInfo.absolutePath() + "/" + fileInfo.baseName();

            // Export count matrix (converted from PPM to counts) and metadata
            bool countMatrixSuccess = exportCountMatrixToCSV(basePath + "_count_matrix.csv", m_combinedData);
            bool metadataSuccess = m_geneCountHandler->exportMetadataToCSV(basePath + "_metadata.csv", m_combinedData);

            // Provide detailed feedback
            if (countMatrixSuccess)
            {
                addProgressMessage(QString("✓ Count matrix exported to: %1_count_matrix.csv").arg(basePath));
            }
            else
            {
                addProgressMessage("✗ Failed to export count matrix.");
            }

            if (metadataSuccess)
            {
                addProgressMessage(QString("✓ Metadata exported to: %1_metadata.csv").arg(basePath));
            }
            else
            {
                addProgressMessage("✗ Failed to export metadata.");
            }

            // Summary message
            if (countMatrixSuccess && metadataSuccess)
            {
                addProgressMessage(QString("✓ Export completed successfully! Both files saved in: %1").arg(fileInfo.absolutePath()));

                // Show success dialog with file locations
                QString message = QString("Export completed successfully!\n\n"
                                          "Files saved:\n"
                                          "• Count Matrix: %1_count_matrix.csv\n"
                                          "• Metadata: %1_metadata.csv\n\n"
                                          "Both files are required for StatMaker analysis.")
                                      .arg(fileInfo.baseName());

                QMessageBox::information(static_cast<QWidget *>(this), "Export Successful", message);
            }
            else if (countMatrixSuccess || metadataSuccess)
            {
                addProgressMessage("⚠ Partial export completed. Some files failed to export.");

                QString message = QString("Partial export completed.\n\n"
                                          "Successfully exported:\n");
                if (countMatrixSuccess)
                    message += "• Count Matrix\n";
                if (metadataSuccess)
                    message += "• Metadata\n";
                message += "\nFailed to export:\n";
                if (!countMatrixSuccess)
                    message += "• Count Matrix\n";
                if (!metadataSuccess)
                    message += "• Metadata\n";
                message += "\nBoth files are required for StatMaker analysis.";

                QMessageBox::warning(static_cast<QWidget *>(this), "Partial Export", message);
            }
            else
            {
                addProgressMessage("✗ Export failed. No files were saved.");
                QMessageBox::critical(static_cast<QWidget *>(this), "Export Failed",
                                      "Failed to export both count matrix and metadata.\n"
                                      "Please check file permissions and try again.");
            }
        }
    }

    // Slot implementations - Group assignment
    void MainWindow::onAssignGroup()
    {
        if (!m_geneCountFilesTable)
            return;

        // Get group name from combo or custom text
        QString groupName;
        if (m_groupTypeCombo) {
            QString selected = m_groupTypeCombo->currentText();
            if (selected == "Custom...") {
                if (!m_groupNameEdit || m_groupNameEdit->text().isEmpty()) {
                    QMessageBox::warning(this, "Warning", "Please enter a custom group name.");
                    return;
                }
                groupName = m_groupNameEdit->text();
            } else {
                groupName = selected;
            }
        } else if (m_groupNameEdit && !m_groupNameEdit->text().isEmpty()) {
            groupName = m_groupNameEdit->text();
        } else {
            QMessageBox::warning(this, "Warning", "Please select or enter a group name.");
            return;
        }
        QList<QTableWidgetItem *> selectedItems = m_geneCountFilesTable->selectedItems();

        if (selectedItems.isEmpty())
        {
            QMessageBox::warning(static_cast<QWidget *>(this), "Warning", "Please select files to assign group.");
            return;
        }

        int assignedCount = 0;
        QSet<int> processedRows;

        for (QTableWidgetItem *item : selectedItems)
        {
            int row = item->row();
            if (processedRows.contains(row))
            {
                continue;
            }
            processedRows.insert(row);

            QTableWidgetItem *fileNameItem = m_geneCountFilesTable->item(row, 0);
            if (fileNameItem)
            {
                QString fileName = fileNameItem->text();
                if (m_geneCountHandler->assignGroupToFile(fileName, groupName))
                {
                    assignedCount++;
                }
            }
        }

        // Update the table to show the new group assignments
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        addProgressMessage(QString("Assigned group '%1' to %2 files.").arg(groupName).arg(assignedCount));
    }

    void MainWindow::onAutoAssignGroups()
    {
        addProgressMessage("Auto-assigning groups by filename...");

        // Define patterns for auto-assignment based on the example files
        QMap<QString, QString> patterns;
        patterns["NON"] = "Non-Selected";
        patterns["non-selected"] = "Non-Selected";
        patterns["nonselected"] = "Non-Selected";
        patterns["_NON_"] = "Non-Selected";
        patterns["SEL"] = "Selected";
        patterns["selected"] = "Selected";
        patterns["_SEL_"] = "Selected";
        patterns["vector"] = "Vector Control";
        patterns["control"] = "Vector Control";
        patterns["bait"] = "Bait";

        int assignedCount = m_geneCountHandler->autoAssignGroups(patterns);

        // Update the table to show the new group assignments
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        addProgressMessage(QString("Auto-assigned groups to %1 files.").arg(assignedCount));
    }

    // Slot implementations - Analysis control
    void MainWindow::onRunAnalysis()
    {
        // Validate that we have data
        if (!m_combinedData.isValid)
        {
            QMessageBox::warning(this, "No Data",
                                 "Please load gene count files and generate data before running analysis.");
            addProgressMessage("Cannot run analysis: No valid data available.");
            return;
        }

        // Check that we have at least 2 different groups
        QSet<QString> uniqueGroups;
        for (const auto &group : m_combinedData.groupNames)
        {
            uniqueGroups.insert(group);
        }

        if (uniqueGroups.size() < 2)
        {
            QMessageBox::warning(this, "Insufficient Groups",
                                 "Analysis requires at least 2 different groups. Please assign groups to your samples.");
            addProgressMessage("Cannot run analysis: Need at least 2 different groups.");
            return;
        }

        // Prevent multiple simultaneous analyses
        if (m_analysisThread && m_analysisThread->isRunning())
        {
            QMessageBox::information(this, "Analysis Running",
                                     "Analysis is already running. Please wait for it to complete or stop it first.");
            return;
        }

        if (m_workingDirectory.isEmpty())
            m_workingDirectory = resolveWorkingDirectory();

        addProgressMessage(QString("Starting StatMaker analysis in %1").arg(resolveWorkingDirectory()));

        if (m_junctionFilesEdit && m_junctionFilesEdit->text().isEmpty())
            autoDetectJunctionFiles();

        // Populate the contrast selector for three-way comparisons (Task 5)
        if (m_contrastSelectorCombo && uniqueGroups.size() >= 3)
        {
            m_contrastSelectorCombo->blockSignals(true);
            m_contrastSelectorCombo->clear();

            QList<QString> groupList = uniqueGroups.values();
            std::sort(groupList.begin(), groupList.end());

            // Generate all pairwise contrasts: each group vs group 0 (reference)
            for (int i = 1; i < groupList.size(); ++i)
            {
                m_contrastSelectorCombo->addItem(
                    QString("%1 vs %2").arg(groupList[i]).arg(groupList[0]));
            }
            // Also add non-reference contrasts for three-way
            if (groupList.size() >= 3)
            {
                for (int i = 2; i < groupList.size(); ++i)
                {
                    for (int j = 1; j < i; ++j)
                    {
                        m_contrastSelectorCombo->addItem(
                            QString("%1 vs %2").arg(groupList[i]).arg(groupList[j]));
                    }
                }
            }

            m_contrastSelectorCombo->blockSignals(false);
            addProgressMessage(QString("Three-way comparison: %1 contrasts available.").arg(m_contrastSelectorCombo->count()));
        }
        else if (m_contrastSelectorCombo)
        {
            // Two-way: single contrast
            m_contrastSelectorCombo->blockSignals(true);
            m_contrastSelectorCombo->clear();
            QList<QString> groupList = uniqueGroups.values();
            std::sort(groupList.begin(), groupList.end());
            if (groupList.size() == 2)
            {
                m_contrastSelectorCombo->addItem(
                    QString("%1 vs %2").arg(groupList[1]).arg(groupList[0]));
            }
            m_contrastSelectorCombo->blockSignals(false);
        }

        // Get analysis parameters
        AnalysisRunConfig config;
        config.pValueThreshold = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        config.ppmThreshold = m_ppmThresholdSpinBox_y2h ? m_ppmThresholdSpinBox_y2h->value() : 0.0;
        config.enrichmentPValueThreshold = m_enrichPvalSpinBox ? m_enrichPvalSpinBox->value() : 1.0;
        config.enrichmentFoldChangeThreshold = m_enrichFcSpinBox ? m_enrichFcSpinBox->value() : 0.0;
        config.specificityPValueThreshold = m_specPvalSpinBox ? m_specPvalSpinBox->value() : 1.0;
        config.baitGroupSize = 10;
        config.junctionFiles = m_junctionFilesEdit ? splitPaths(m_junctionFilesEdit->text()) : QStringList{};
        config.workingDirectory = resolveWorkingDirectory();

        // Update UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(false);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(true);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(false);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(true);

        // Show progress bar
        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(true);
            m_analysisProgressBar->setValue(0);
        }

        // Create worker and thread
        m_analysisWorker = std::make_unique<AnalysisWorker>(m_combinedData, config);
        m_analysisThread = std::make_unique<QThread>();

        // Move worker to thread
        m_analysisWorker->moveToThread(m_analysisThread.get());

        // Connect signals with explicit connection types for cross-thread communication
        connect(m_analysisThread.get(), &QThread::started,
                m_analysisWorker.get(), &AnalysisWorker::runAnalysis);
        connect(m_analysisWorker.get(), &AnalysisWorker::progressChanged,
                this, &MainWindow::onAnalysisProgress, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::analysisFinished,
                this, &MainWindow::onAnalysisFinished, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::analysisError,
                this, &MainWindow::onAnalysisError, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::debugMessage,
                this, &MainWindow::onDebugMessage, Qt::QueuedConnection);

        // Proper thread lifecycle management
        connect(m_analysisWorker.get(), &AnalysisWorker::finished,
                m_analysisThread.get(), &QThread::quit);
        connect(m_analysisThread.get(), &QThread::finished,
                this, &MainWindow::onThreadFinished);

        // Force signal processing by connecting finished signal to a dummy slot
        connect(m_analysisWorker.get(), &AnalysisWorker::finished,
                this, [this]()
                { addProgressMessage("DEBUG: Worker finished signal received"); });

        // Start analysis (large stack for Eigen matrix operations on many genes)
        m_analysisThread->setStackSize(64 * 1024 * 1024); // 64MB
        m_analysisThread->start();
    }

    void MainWindow::onStopAnalysis()
    {
        if (m_analysisWorker)
        {
            addProgressMessage("Stopping analysis...");
            m_analysisWorker->stopAnalysis();
        }

        if (m_analysisThread && m_analysisThread->isRunning())
        {
            m_analysisThread->quit();
            if (!m_analysisThread->wait(5000))
            {
                m_analysisThread->terminate();
                m_analysisThread->wait();
            }
        }

        // Reset UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(true);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(false);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(true);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(false);

        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(false);
        }

        addProgressMessage("Analysis stopped by user.");
    }

    void MainWindow::onResetAnalysis()
    {
        // Stop any running analysis
        onStopAnalysis();

        // Clear results
        m_analysisResults.isValid = false;
        m_analysisResults.results = Eigen::MatrixXd();
        m_analysisResults.geneNames.clear();
        m_resultsSqlitePath.clear();
        m_resultsColumnsSized = false;
        invalidatePlotCache();

        // Clear results table
        if (m_resultsTable)
        {
            m_resultsTable->setRowCount(0);
            m_resultsTable->setColumnCount(0);
        }

        // Reset summary labels
        updateResultsSummary(m_analysisResults);

        addProgressMessage("Analysis results cleared.");
    }

    void MainWindow::onAnalysisProgress(int progress, const QString &message)
    {
        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setValue(progress);
        }
        addProgressMessage(message);
    }

    void MainWindow::onAnalysisFinished(const AnalysisResults &results)
    {
        QString validationError;
        if (!validateAnalysisResults(results, &validationError))
        {
            onAnalysisError(QString("Analysis completed with invalid results: %1").arg(validationError));
            return;
        }

        // Store results
        m_analysisResults = results;
        invalidatePlotCache();
        m_resultsColumnsSized = false;

        // Update UI
        updateResultsTable(results);
        updateResultsSummary(results);

        // Switch to results tab
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab (FIXED: was 2, should be 1)
        }

        addProgressMessage(QString("Analysis completed! Found %1 significant genes out of %2 total genes.")
                               .arg(results.significantGenes)
                               .arg(results.totalGenes));

        // Write results to SQLite
        writeResultsToSqlite(results);

        // Refresh visualization plot with new results
        refreshPlot();
    }

    void MainWindow::writeResultsToSqlite(const AnalysisResults &results)
    {
        QElapsedTimer timer;
        timer.start();

        const QString workdir = resolveWorkingDirectory();
        QDir baseDir(workdir.isEmpty() ? m_lastUsedDirectory : workdir);
        if (!baseDir.exists())
            baseDir.mkpath(".");
        if (!baseDir.exists("analyzed_files"))
            baseDir.mkpath("analyzed_files");

        const QString outputPath = resultsDatabasePathForWorkdir(baseDir.absolutePath());
        m_resultsSqlitePath = outputPath;

        const QString activeBait = results.activeContrastLabel.section(" vs ", 0, 0).trimmed();
        auto baitMatches = [&](const std::string &bait)
        {
            return activeBait.isEmpty() || QString::fromStdString(bait).compare(activeBait, Qt::CaseInsensitive) == 0;
        };

        QString connName = QString("statmaker_results_write_%1").arg(QUuid::createUuid().toString(QUuid::Id128));
        bool opened = false;
        {
            QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
            db.setDatabaseName(outputPath);
            opened = db.open();
            if (!opened)
            {
                addProgressMessage("Error: cannot write results to " + outputPath);
            }
            if (opened)
            {
                QSqlQuery q(db);
                q.exec("PRAGMA journal_mode = WAL");
                q.exec("DROP TABLE IF EXISTS analysis_results");
                q.exec("DROP TABLE IF EXISTS analysis_summary");
                q.exec("DROP TABLE IF EXISTS contrasts");
                q.exec("DROP TABLE IF EXISTS run_parameters");
                q.exec("DROP TABLE IF EXISTS input_files");
                q.exec("CREATE TABLE analysis_results ("
                       "gene TEXT PRIMARY KEY, "
                       "bait TEXT, "
                       "base_mean REAL, "
                       "log2_fold_change REAL, "
                       "lfc_se REAL, "
                       "stat REAL, "
                       "pvalue REAL, "
                       "padj REAL, "
                       "enrichment_call TEXT, "
                       "enrichment_score REAL, "
                       "specificity_score REAL, "
                       "in_frame_score REAL, "
                       "borda_score REAL, "
                       "in_frame_transcripts TEXT, "
                       "active_contrast_label TEXT, "
                       "created_at TEXT)");
                q.exec("CREATE TABLE analysis_summary (key TEXT PRIMARY KEY, value TEXT)");
                q.exec("CREATE TABLE contrasts (contrast_index INTEGER PRIMARY KEY, label TEXT, is_active INTEGER)");
                q.exec("CREATE TABLE run_parameters (key TEXT PRIMARY KEY, value TEXT)");
                q.exec("CREATE TABLE input_files ("
                       "path TEXT PRIMARY KEY, "
                       "sample_name TEXT, "
                       "group_name TEXT)");

                // Build lookup maps from Y2H-SCORES results by gene name
                std::map<std::string, double> enrichmentScoreMap;
                for (const auto &er : results.enrichmentScores)
                {
                    if (baitMatches(er.bait))
                        enrichmentScoreMap[er.gene] = er.total_score;
                }
                std::map<std::string, double> specificityScoreMap;
                for (const auto &spec : results.specificityScores)
                {
                    if (baitMatches(spec.bait))
                        specificityScoreMap[spec.gene] = spec.total_score;
                }
                std::map<std::string, double> inFrameScoreMap;
                std::map<std::string, QString> inFrameTranscriptMap;
                for (const auto &inFrame : results.inFrameScores)
                {
                    if (baitMatches(inFrame.bait))
                    {
                        inFrameScoreMap[inFrame.gene] = inFrame.freq_score;
                        inFrameTranscriptMap[inFrame.gene] = QString::fromStdString(inFrame.transcripts);
                    }
                }
                std::map<std::string, double> bordaScoreMap;
                for (const auto &y2h : results.y2hScores)
                {
                    if (baitMatches(y2h.bait))
                    {
                        bordaScoreMap[y2h.gene] = y2h.borda_score;
                        if (!QString::fromStdString(y2h.in_frame_transcripts).isEmpty())
                            inFrameTranscriptMap[y2h.gene] = QString::fromStdString(y2h.in_frame_transcripts);
                    }
                }

                db.transaction();
                q.prepare("INSERT INTO analysis_results VALUES "
                          "(:gene, :bait, :basemean, :lfc, :se, :stat, :pval, :padj, :enrich_call, "
                          ":enrichment_score, :specificity_score, :in_frame_score, :borda_score, "
                          ":in_frame_transcripts, :active_contrast_label, :created_at)");

                for (int i = 0; i < results.results.rows(); i++)
                {
                    double padj = results.results(i, 5);
                    double lfc = results.results(i, 1);
                    QString enrichment = "NS";
                    if (!std::isnan(padj) && padj < results.pValueThreshold)
                    {
                        enrichment = (lfc > 0) ? "Enriched" : "Depleted";
                    }

                    const std::string &geneName = results.geneNames[i];
                    double escore = 0.0;
                    auto esIt = enrichmentScoreMap.find(geneName);
                    if (esIt != enrichmentScoreMap.end())
                    {
                        escore = esIt->second;
                    }
                    double specScore = 0.0;
                    auto specIt = specificityScoreMap.find(geneName);
                    if (specIt != specificityScoreMap.end())
                        specScore = specIt->second;
                    double ifScore = 0.0;
                    auto ifIt = inFrameScoreMap.find(geneName);
                    if (ifIt != inFrameScoreMap.end())
                        ifScore = ifIt->second;
                    double bscore = 0.0;
                    auto bsIt = bordaScoreMap.find(geneName);
                    if (bsIt != bordaScoreMap.end())
                    {
                        bscore = bsIt->second;
                    }

                    q.bindValue(":gene", sanitizeGeneName(geneName));
                    q.bindValue(":bait", activeBait);
                    q.bindValue(":basemean", results.results(i, 0));
                    q.bindValue(":lfc", lfc);
                    q.bindValue(":se", results.results(i, 2));
                    q.bindValue(":stat", results.results(i, 3));
                    q.bindValue(":pval", results.results(i, 4));
                    q.bindValue(":padj", padj);
                    q.bindValue(":enrich_call", enrichment);
                    q.bindValue(":enrichment_score", escore);
                    q.bindValue(":specificity_score", specScore);
                    q.bindValue(":in_frame_score", ifScore);
                    q.bindValue(":borda_score", bscore);
                    q.bindValue(":in_frame_transcripts",
                                inFrameTranscriptMap.count(geneName)
                                    ? inFrameTranscriptMap[geneName]
                                    : QString());
                    q.bindValue(":active_contrast_label", results.activeContrastLabel);
                    q.bindValue(":created_at", results.createdAt);
                    q.exec();
                }

                // Summary table
                q.prepare("INSERT INTO analysis_summary VALUES (:key, :value)");
                q.bindValue(":key", "total_genes"); q.bindValue(":value", QString::number(results.totalGenes)); q.exec();
                q.bindValue(":key", "significant_genes"); q.bindValue(":value", QString::number(results.significantGenes)); q.exec();
                q.bindValue(":key", "upregulated"); q.bindValue(":value", QString::number(results.upregulatedGenes)); q.exec();
                q.bindValue(":key", "downregulated"); q.bindValue(":value", QString::number(results.downregulatedGenes)); q.exec();
                q.bindValue(":key", "pvalue_threshold"); q.bindValue(":value", QString::number(results.pValueThreshold)); q.exec();
                q.bindValue(":key", "active_contrast_label"); q.bindValue(":value", results.activeContrastLabel); q.exec();
                q.bindValue(":key", "working_directory"); q.bindValue(":value", baseDir.absolutePath()); q.exec();

                q.prepare("INSERT INTO contrasts VALUES (:idx, :label, :active)");
                for (int idx = 0; idx < static_cast<int>(results.contrastLabels.size()); ++idx)
                {
                    q.bindValue(":idx", idx);
                    q.bindValue(":label", QString::fromStdString(results.contrastLabels[idx]));
                    q.bindValue(":active", idx == results.activeContrast ? 1 : 0);
                    q.exec();
                }

                q.prepare("INSERT INTO run_parameters VALUES (:key, :value)");
                q.bindValue(":key", "output_path"); q.bindValue(":value", outputPath); q.exec();
                q.bindValue(":key", "created_at"); q.bindValue(":value", results.createdAt); q.exec();
                q.bindValue(":key", "active_bait"); q.bindValue(":value", activeBait); q.exec();
                q.bindValue(":key", "junction_files"); q.bindValue(":value", m_junctionFilesEdit ? m_junctionFilesEdit->text() : QString()); q.exec();

                q.prepare("INSERT INTO input_files VALUES (:path, :sample_name, :group_name)");
                for (const auto &inputFile : m_geneCountHandler->getGeneCountFiles())
                {
                    q.bindValue(":path", inputFile.filePath);
                    q.bindValue(":sample_name", inputFile.sampleName);
                    q.bindValue(":group_name", inputFile.groupName);
                    q.exec();
                }

                db.commit();
                q.clear();
                db.close();
            }
        }
        QSqlDatabase::removeDatabase(connName);

        if (!opened)
            return;

        addProgressMessage(QString("Results saved to %1 (%2 ms)").arg(outputPath).arg(timer.elapsed()));
    }

    bool MainWindow::loadResultsFromSqlite(const QString &sqlitePath)
    {
        if (!QFile::exists(sqlitePath))
            return false;

        QElapsedTimer timer;
        timer.start();

        AnalysisResults loadedResults;
        loadedResults.pValueThreshold = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        loadedResults.log2FCThreshold = m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 0.0;
        std::vector<std::array<double, kDeseqColumnCount>> rows;

        const QString connName = QString("statmaker_results_read_%1").arg(QUuid::createUuid().toString(QUuid::Id128));
        bool opened = false;
        bool queryOk = false;
        {
            QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
            db.setDatabaseName(sqlitePath);
            db.setConnectOptions("QSQLITE_OPEN_READONLY");
            opened = db.open();
            if (opened)
            {
                QSqlQuery q(db);
                queryOk = q.exec("SELECT gene, bait, base_mean, log2_fold_change, lfc_se, stat, pvalue, padj, "
                                 "enrichment_score, specificity_score, in_frame_score, borda_score, "
                                 "in_frame_transcripts, active_contrast_label, created_at "
                                 "FROM analysis_results ORDER BY padj ASC, gene ASC");
                if (queryOk)
                {
                    while (q.next())
                    {
                        const QString gene = q.value(0).toString();
                        const QString bait = q.value(1).toString();

                        rows.push_back({{
                            q.value(2).toDouble(),
                            q.value(3).toDouble(),
                            q.value(4).toDouble(),
                            q.value(5).toDouble(),
                            q.value(6).toDouble(),
                            q.value(7).toDouble(),
                        }});
                        loadedResults.geneNames.push_back(gene.toStdString());

                        EnrichmentResult enrichment;
                        enrichment.gene = gene.toStdString();
                        enrichment.bait = bait.toStdString();
                        enrichment.total_score = q.value(8).toDouble();
                        if (enrichment.total_score > 0.0)
                            loadedResults.enrichmentScores.push_back(enrichment);

                        SpecificityResult specificity;
                        specificity.gene = gene.toStdString();
                        specificity.bait = bait.toStdString();
                        specificity.total_score = q.value(9).toDouble();
                        if (specificity.total_score > 0.0)
                            loadedResults.specificityScores.push_back(specificity);

                        InFrameResult inFrame;
                        inFrame.gene = gene.toStdString();
                        inFrame.bait = bait.toStdString();
                        inFrame.freq_score = q.value(10).toDouble();
                        inFrame.transcripts = q.value(12).toString().toStdString();
                        if (inFrame.freq_score > 0.0 || !inFrame.transcripts.empty())
                            loadedResults.inFrameScores.push_back(inFrame);

                        Y2HScore y2h;
                        y2h.gene = gene.toStdString();
                        y2h.bait = bait.toStdString();
                        y2h.enrichment_score = enrichment.total_score;
                        y2h.specificity_score = specificity.total_score;
                        y2h.in_frame_score = inFrame.freq_score;
                        y2h.in_frame_transcripts = inFrame.transcripts;
                        y2h.borda_score = q.value(11).toDouble();
                        y2h.sum_scores =
                            y2h.enrichment_score + y2h.specificity_score + y2h.in_frame_score;
                        if (y2h.borda_score > 0.0)
                            loadedResults.y2hScores.push_back(y2h);

                        if (loadedResults.activeContrastLabel.isEmpty())
                            loadedResults.activeContrastLabel = q.value(13).toString();
                        if (loadedResults.createdAt.isEmpty())
                            loadedResults.createdAt = q.value(14).toString();
                    }
                    q.clear();

                    loadedResults.results = Eigen::MatrixXd(static_cast<int>(rows.size()), kDeseqColumnCount);
                    for (int row = 0; row < static_cast<int>(rows.size()); ++row)
                    {
                        for (int col = 0; col < kDeseqColumnCount; ++col)
                            loadedResults.results(row, col) = rows[row][col];
                    }

                    if (q.exec("SELECT key, value FROM analysis_summary"))
                    {
                        while (q.next())
                        {
                            const QString key = q.value(0).toString();
                            const QString value = q.value(1).toString();
                            if (key == "total_genes")
                                loadedResults.totalGenes = value.toInt();
                            else if (key == "significant_genes")
                                loadedResults.significantGenes = value.toInt();
                            else if (key == "upregulated")
                                loadedResults.upregulatedGenes = value.toInt();
                            else if (key == "downregulated")
                                loadedResults.downregulatedGenes = value.toInt();
                            else if (key == "pvalue_threshold")
                                loadedResults.pValueThreshold = value.toDouble();
                            else if (key == "active_contrast_label" && loadedResults.activeContrastLabel.isEmpty())
                                loadedResults.activeContrastLabel = value;
                        }
                    }
                    q.clear();

                    if (q.exec("SELECT contrast_index, label, is_active FROM contrasts ORDER BY contrast_index"))
                    {
                        while (q.next())
                        {
                            loadedResults.contrastLabels.push_back(q.value(1).toString().toStdString());
                            if (q.value(2).toInt() == 1)
                            {
                                loadedResults.activeContrast = q.value(0).toInt();
                                loadedResults.activeContrastLabel = q.value(1).toString();
                            }
                        }
                    }
                    q.clear();
                }
                db.close();
            }
        }
        QSqlDatabase::removeDatabase(connName);

        if (!opened || !queryOk)
            return false;

        loadedResults.totalGenes = loadedResults.totalGenes > 0 ? loadedResults.totalGenes : loadedResults.results.rows();
        loadedResults.isValid = loadedResults.results.rows() > 0 && loadedResults.geneNames.size() == static_cast<size_t>(loadedResults.results.rows());
        loadedResults.errorMessage.clear();

        QString validationError;
        if (!validateAnalysisResults(loadedResults, &validationError))
        {
            addProgressMessage(QString("Existing StatMaker results were ignored: %1").arg(validationError));
            return false;
        }

        m_resultsSqlitePath = sqlitePath;
        m_analysisResults = loadedResults;
        invalidatePlotCache();
        m_resultsColumnsSized = false;
        updateResultsTable(m_analysisResults);
        updateResultsSummary(m_analysisResults);
        addProgressMessage(QString("Loaded previous StatMaker results from %1 (%2 ms)")
                               .arg(sqlitePath)
                               .arg(timer.elapsed()));
        return true;
    }

    bool MainWindow::validateAnalysisResults(const AnalysisResults &results, QString *errorMessage) const
    {
        if (!results.isValid)
        {
            if (errorMessage)
                *errorMessage = results.errorMessage.isEmpty()
                                    ? QStringLiteral("Results are marked invalid.")
                                    : results.errorMessage;
            return false;
        }

        if (results.results.rows() <= 0 || results.results.cols() < kDeseqColumnCount)
        {
            if (errorMessage)
                *errorMessage = QStringLiteral("Result matrix is empty or malformed.");
            return false;
        }

        if (results.geneNames.size() != static_cast<size_t>(results.results.rows()))
        {
            if (errorMessage)
            {
                *errorMessage = QStringLiteral("Gene name count (%1) does not match result rows (%2).")
                                    .arg(results.geneNames.size())
                                    .arg(results.results.rows());
            }
            return false;
        }

        for (size_t idx = 0; idx < results.geneNames.size(); ++idx)
        {
            if (geneNameHasInvalidBytes(results.geneNames[idx]))
            {
                if (errorMessage)
                {
                    *errorMessage = QStringLiteral("Invalid control bytes detected in gene name at row %1.")
                                        .arg(static_cast<qulonglong>(idx));
                }
                return false;
            }
        }

        return true;
    }

    void MainWindow::onAnalysisError(const QString &errorMessage)
    {
        // Show error message
        QMessageBox::critical(this, "Analysis Error", errorMessage);
        addProgressMessage("Analysis failed: " + errorMessage);
    }

    void MainWindow::onDebugMessage(const QString &message)
    {
        addProgressMessage("DEBUG: " + message);
    }

    void MainWindow::onThreadFinished()
    {
        // Reset UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(true);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(false);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(true);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(false);

        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(false);
        }
    }

    void MainWindow::updateResultsTable(const AnalysisResults &results)
    {
        if (!m_resultsTable || !results.isValid)
        {
            return;
        }

        QElapsedTimer timer;
        timer.start();

        const auto &matrix = results.results;
        const auto &geneNames = results.geneNames;
        const QString activeBait = results.activeContrastLabel.section(" vs ", 0, 0).trimmed();
        auto baitMatches = [&](const std::string &bait)
        {
            return activeBait.isEmpty() || QString::fromStdString(bait).compare(activeBait, Qt::CaseInsensitive) == 0;
        };

        std::map<std::string, double> enrichmentScoreMap;
        for (const auto &entry : results.enrichmentScores)
        {
            if (baitMatches(entry.bait))
                enrichmentScoreMap[entry.gene] = entry.total_score;
        }
        std::map<std::string, double> specificityScoreMap;
        for (const auto &entry : results.specificityScores)
        {
            if (baitMatches(entry.bait))
                specificityScoreMap[entry.gene] = entry.total_score;
        }
        std::map<std::string, double> inFrameScoreMap;
        std::map<std::string, QString> inFrameTranscriptMap;
        for (const auto &entry : results.inFrameScores)
        {
            if (baitMatches(entry.bait))
            {
                inFrameScoreMap[entry.gene] = entry.freq_score;
                inFrameTranscriptMap[entry.gene] = QString::fromStdString(entry.transcripts);
            }
        }
        std::map<std::string, double> bordaScoreMap;
        for (const auto &entry : results.y2hScores)
        {
            if (baitMatches(entry.bait))
            {
                bordaScoreMap[entry.gene] = entry.borda_score;
                if (!QString::fromStdString(entry.in_frame_transcripts).isEmpty())
                    inFrameTranscriptMap[entry.gene] = QString::fromStdString(entry.in_frame_transcripts);
            }
        }

        // Setup table
        m_resultsTable->setSortingEnabled(false);
        m_resultsTable->clearContents();
        m_resultsTable->setRowCount(matrix.rows());
        m_resultsTable->setColumnCount(13);

        // Set column headers
        QStringList headers;
        headers << "Gene" << "baseMean" << "log2FoldChange" << "lfcSE" << "stat" << "pvalue" << "padj"
                << "Enrichment Call" << "Enrichment Score" << "Specificity Score"
                << "In-Frame Score" << "Borda Score" << "In-Frame Transcripts";
        m_resultsTable->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < matrix.rows(); ++row)
        {
            const std::string &geneKey = geneNames[row];
            const QString geneDisplay = sanitizeGeneName(geneKey);

            // Gene name
            QTableWidgetItem *geneItem = new QTableWidgetItem(geneDisplay);
            geneItem->setFlags(geneItem->flags() & ~Qt::ItemIsEditable);
            m_resultsTable->setItem(row, 0, geneItem);

            // Results columns
            for (int col = 0; col < matrix.cols(); ++col)
            {
                QTableWidgetItem *item = new QTableWidgetItem(QString::number(matrix(row, col), 'g', 6));
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);

                // Color significant genes - padj is column 5 in matrix (0-indexed)
                if (col == 5 && matrix(row, 5) < results.pValueThreshold) // padj column (5th column in matrix, 0-indexed)
                {
                    item->setBackground(QColor(255, 255, 200)); // Light yellow background
                }

                m_resultsTable->setItem(row, col + 1, item);
            }

            const double padj = matrix(row, 5);
            const double lfc = matrix(row, 1);
            QString enrichmentCall = "NS";
            if (!std::isnan(padj) && padj < results.pValueThreshold)
                enrichmentCall = lfc > 0 ? "Enriched" : "Depleted";

            auto makeReadOnlyItem = [](const QString &value)
            {
                QTableWidgetItem *item = new QTableWidgetItem(value);
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);
                return item;
            };

            const auto enrichmentIt = enrichmentScoreMap.find(geneKey);
            const auto specificityIt = specificityScoreMap.find(geneKey);
            const auto inFrameIt = inFrameScoreMap.find(geneKey);
            const auto bordaIt = bordaScoreMap.find(geneKey);
            const auto transcriptsIt = inFrameTranscriptMap.find(geneKey);

            m_resultsTable->setItem(row, kResultsColumnEnrichmentCall, makeReadOnlyItem(enrichmentCall));
            m_resultsTable->setItem(row, kResultsColumnEnrichmentScore,
                                    makeReadOnlyItem(enrichmentIt != enrichmentScoreMap.end()
                                                         ? QString::number(enrichmentIt->second, 'g', 6)
                                                         : QString()));
            m_resultsTable->setItem(row, kResultsColumnSpecificityScore,
                                    makeReadOnlyItem(specificityIt != specificityScoreMap.end()
                                                         ? QString::number(specificityIt->second, 'g', 6)
                                                         : QString()));
            m_resultsTable->setItem(row, kResultsColumnInFrameScore,
                                    makeReadOnlyItem(inFrameIt != inFrameScoreMap.end()
                                                         ? QString::number(inFrameIt->second, 'g', 6)
                                                         : QString()));
            m_resultsTable->setItem(row, kResultsColumnBordaScore,
                                    makeReadOnlyItem(bordaIt != bordaScoreMap.end()
                                                         ? QString::number(bordaIt->second, 'g', 6)
                                                         : QString()));
            m_resultsTable->setItem(row, kResultsColumnInFrameTranscripts,
                                    makeReadOnlyItem(transcriptsIt != inFrameTranscriptMap.end()
                                                         ? transcriptsIt->second
                                                         : QString()));

            // Progress update every 1000 rows
            if (row % 1000 == 0 && row > 0)
            {
                addProgressMessage(QString("Populated %1/%2 result rows...")
                                       .arg(row)
                                       .arg(matrix.rows()));
            }
        }

        if (!m_resultsColumnsSized || matrix.rows() <= 2000)
        {
            m_resultsTable->resizeColumnsToContents();
            m_resultsColumnsSized = true;
        }

        // Enable sorting
        m_resultsTable->setSortingEnabled(true);

        // Sort by adjusted p-value by default
        m_resultsTable->sortItems(kResultsColumnPadj, Qt::AscendingOrder);
        addProgressMessage(QString("Results table updated in %1 ms (%2 rows)")
                               .arg(timer.elapsed())
                               .arg(matrix.rows()));
    }

    void MainWindow::updateResultsSummary(const AnalysisResults &results)
    {
        if (results.isValid)
        {
            if (m_totalGenesLabel)
            {
                m_totalGenesLabel->setText(QString("Total Genes: %1").arg(results.totalGenes));
            }
            if (m_significantGenesLabel)
            {
                m_significantGenesLabel->setText(QString("Significant: %1").arg(results.significantGenes));
            }
            if (m_upregulatedLabel)
            {
                m_upregulatedLabel->setText(QString("Upregulated: %1").arg(results.upregulatedGenes));
            }
            if (m_downregulatedLabel)
            {
                m_downregulatedLabel->setText(QString("Downregulated: %1").arg(results.downregulatedGenes));
            }
        }
        else
        {
            // Reset to defaults
            if (m_totalGenesLabel)
                m_totalGenesLabel->setText("Total Genes: 0");
            if (m_significantGenesLabel)
                m_significantGenesLabel->setText("Significant: 0");
            if (m_upregulatedLabel)
                m_upregulatedLabel->setText("Upregulated: 0");
            if (m_downregulatedLabel)
                m_downregulatedLabel->setText("Downregulated: 0");
        }
    }

    // Slot implementations - Results filtering
    void MainWindow::onApplyFilter()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        applyCurrentFilter();
    }

    void MainWindow::onClearFilter()
    {
        if (m_pValueThresholdSpinBox)
        {
            m_pValueThresholdSpinBox->setValue(0.05);
        }
        if (m_log2FCThresholdSpinBox)
        {
            m_log2FCThresholdSpinBox->setValue(0.0);
        }

        if (m_analysisResults.isValid && m_resultsTable)
        {
            // Show all rows instead of repopulating the table
            for (int row = 0; row < m_resultsTable->rowCount(); ++row)
            {
                m_resultsTable->setRowHidden(row, false);
            }
            addProgressMessage("Filter cleared - showing all results.");
        }
    }

    void MainWindow::applyCurrentFilter()
    {
        if (!m_resultsTable || !m_analysisResults.isValid)
        {
            return;
        }

        double pValueThreshold = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        double log2FCThreshold = m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 0.0;

        int visibleRows = 0;
        for (int row = 0; row < m_resultsTable->rowCount(); ++row)
        {
            // Get p-value and log2FC
            double padj = m_resultsTable->item(row, 6)->text().toDouble();             // padj column
            double log2FC = std::abs(m_resultsTable->item(row, 2)->text().toDouble()); // log2FoldChange column

            bool passesFilter = (padj < pValueThreshold) && (log2FC > log2FCThreshold);

            m_resultsTable->setRowHidden(row, !passesFilter);
            if (passesFilter)
                visibleRows++;
        }

        addProgressMessage(QString("Applied filter: p-value < %1, |log2FC| > %2. Showing %3 genes.")
                               .arg(pValueThreshold)
                               .arg(log2FCThreshold)
                               .arg(visibleRows));
    }

    // Slot implementations - Export operations
    void MainWindow::onExportResults()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Export Results",
            m_lastUsedDirectory + "/statmaker_results.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportResultsToCSV(fileName, m_analysisResults, false))
            {
                addProgressMessage("Results exported successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to export results.");
                addProgressMessage("Failed to export results.");
            }
        }
    }

    void MainWindow::onSaveGeneLists()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Save Significant Genes",
            m_lastUsedDirectory + "/statmaker_significant_genes.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportResultsToCSV(fileName, m_analysisResults, true))
            {
                addProgressMessage("Significant genes saved successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to save gene lists.");
                addProgressMessage("Failed to save gene lists.");
            }
        }
    }

    void MainWindow::onExportCountMatrix()
    {
        if (!m_combinedData.isValid)
        {
            QMessageBox::information(this, "No Data", "Please load and generate data first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Export Count Matrix",
            m_lastUsedDirectory + "/count_matrix.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportCountMatrixToCSV(fileName, m_combinedData))
            {
                addProgressMessage("Count matrix exported successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to export count matrix.");
                addProgressMessage("Failed to export count matrix.");
            }
        }
    }

    bool MainWindow::exportResultsToCSV(const QString &filePath, const AnalysisResults &results, bool significantOnly)
    {
        // Export from SQLite if available, otherwise from memory
        if (!results.isValid && !m_resultsSqlitePath.isEmpty() && QFile::exists(m_resultsSqlitePath)) {
            return exportSqliteToCSV(filePath, significantOnly);
        }

        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            return false;

        QTextStream stream(&file);
        stream << "Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,enrichment_call,"
                  "enrichment_score,specificity_score,in_frame_score,borda_score,in_frame_transcripts,"
                  "active_contrast_label,created_at\n";

        const QString activeBait = results.activeContrastLabel.section(" vs ", 0, 0).trimmed();
        auto baitMatches = [&](const std::string &bait)
        {
            return activeBait.isEmpty() || QString::fromStdString(bait).compare(activeBait, Qt::CaseInsensitive) == 0;
        };
        std::map<std::string, double> enrichmentScoreMap;
        for (const auto &entry : results.enrichmentScores)
        {
            if (baitMatches(entry.bait))
                enrichmentScoreMap[entry.gene] = entry.total_score;
        }
        std::map<std::string, double> specificityScoreMap;
        for (const auto &entry : results.specificityScores)
        {
            if (baitMatches(entry.bait))
                specificityScoreMap[entry.gene] = entry.total_score;
        }
        std::map<std::string, double> inFrameScoreMap;
        std::map<std::string, QString> inFrameTranscriptMap;
        for (const auto &entry : results.inFrameScores)
        {
            if (baitMatches(entry.bait))
            {
                inFrameScoreMap[entry.gene] = entry.freq_score;
                inFrameTranscriptMap[entry.gene] = QString::fromStdString(entry.transcripts);
            }
        }
        std::map<std::string, double> bordaScoreMap;
        for (const auto &entry : results.y2hScores)
        {
            if (baitMatches(entry.bait))
            {
                bordaScoreMap[entry.gene] = entry.borda_score;
                if (!QString::fromStdString(entry.in_frame_transcripts).isEmpty())
                    inFrameTranscriptMap[entry.gene] = QString::fromStdString(entry.in_frame_transcripts);
            }
        }

        for (int row = 0; row < results.results.rows(); ++row) {
            double padj = results.results(row, 5);
            if (significantOnly && padj >= results.pValueThreshold)
                continue;

            double lfc = results.results(row, 1);
            QString enrichment = "ns";
            if (!std::isnan(padj) && padj < results.pValueThreshold)
                enrichment = (lfc > 0) ? "enriched" : "depleted";

            const std::string &geneKey = results.geneNames[row];
            stream << csvEscape(sanitizeGeneName(geneKey));
            for (int col = 0; col < results.results.cols(); ++col)
                stream << "," << results.results(row, col);
            stream << "," << csvEscape(enrichment)
                   << "," << enrichmentScoreMap[geneKey]
                   << "," << specificityScoreMap[geneKey]
                   << "," << inFrameScoreMap[geneKey]
                   << "," << bordaScoreMap[geneKey]
                   << "," << csvEscape(inFrameTranscriptMap[geneKey])
                   << "," << csvEscape(results.activeContrastLabel)
                   << "," << csvEscape(results.createdAt)
                   << "\n";
        }
        return true;
    }

    bool MainWindow::exportSqliteToCSV(const QString &csvPath, bool significantOnly)
    {
        QString connName = "deseq2_csv_export";
        bool success = false;
        {
            QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
            db.setDatabaseName(m_resultsSqlitePath);
            if (!db.open()) return false;

            QFile file(csvPath);
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                db.close();
                return false;
            }

            QTextStream stream(&file);
            stream << "Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,enrichment_call,"
                      "enrichment_score,specificity_score,in_frame_score,borda_score,"
                      "in_frame_transcripts,active_contrast_label,created_at\n";

            QString query = "SELECT gene, base_mean, log2_fold_change, lfc_se, "
                           "stat, pvalue, padj, enrichment_call, enrichment_score, "
                           "specificity_score, in_frame_score, borda_score, in_frame_transcripts, "
                           "active_contrast_label, created_at "
                           "FROM analysis_results";
            if (significantOnly) {
                query += " WHERE enrichment_call != 'NS'";
            }
            query += " ORDER BY padj ASC";

            QSqlQuery q(db);
            q.exec(query);
            while (q.next()) {
                stream << csvEscape(q.value(0).toString());
                for (int i = 1; i <= 14; i++)
                {
                    if (i == 7 || i == 12 || i == 13 || i == 14)
                        stream << "," << csvEscape(q.value(i).toString());
                    else
                        stream << "," << q.value(i).toString();
                }
                stream << "\n";
            }

            q.clear();
            db.close();
            success = true;
        }
        QSqlDatabase::removeDatabase(connName);
        return success;
    }

    bool MainWindow::exportCountMatrixToCSV(const QString &filePath, const CombinedData &data)
    {
        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return false;
        }

        QTextStream stream(&file);

        // Write header with sample names
        stream << "Gene";
        for (const auto &sampleName : data.sampleNames)
        {
            stream << "," << sampleName;
        }
        stream << "\n";

        // Write data - use raw integer counts
        for (int gene = 0; gene < data.geneNames.size(); ++gene)
        {
            // Gene name
            stream << data.geneNames[gene];

            // Count values for each sample
            for (int sample = 0; sample < data.sampleNames.size(); ++sample)
            {
                stream << "," << data.rawCountMatrix[gene][sample];
            }
            stream << "\n";
        }

        return true;
    }

    // Slot implementations - Visualization
    void MainWindow::onPlotTypeChanged()
    {
        if (!m_plotTypeCombo)
            return;
        refreshPlot();
    }

    void MainWindow::refreshPlot()
    {
        if (!m_plotTypeCombo || !m_analysisResults.isValid)
            return;

        QElapsedTimer timer;
        timer.start();

        QString plotType = m_plotTypeCombo->currentText();
        const QString cacheKey = QString("%1|%2|%3|%4|%5|%6")
                                     .arg(plotType)
                                     .arg(m_analysisResults.activeContrast)
                                     .arg(m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05, 0, 'g', 6)
                                     .arg(m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 0.0, 0, 'g', 6)
                                     .arg(m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5)
                                     .arg(m_resultsRevision);

        if (cacheKey == m_lastPlotCacheKey && m_plotWidget && m_plotWidget->chart())
        {
            addProgressMessage(QString("Plot refresh skipped for %1: cached view is still current.").arg(plotType));
            return;
        }

        if (plotType == "Volcano Plot")
            plotVolcanoPlot();
        else if (plotType == "MA Plot")
            plotMAPlot();
        else if (plotType == "Dispersion Plot")
            plotDispersionPlot();
        else if (plotType == "Three-way Scatter")
            plotThreeWayScatter();

        m_lastPlotCacheKey = cacheKey;
        addProgressMessage(QString("Rendered %1 in %2 ms").arg(plotType).arg(timer.elapsed()));
    }

    void MainWindow::plotMAPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        const auto &res = m_analysisResults.results;
        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;
        qreal markerSize = static_cast<qreal>(ptSize);

        auto *sigSeries = new QScatterSeries();
        sigSeries->setName("Significant");
        sigSeries->setMarkerSize(markerSize);
        sigSeries->setColor(QColor(220, 50, 50));
        sigSeries->setBorderColor(QColor(220, 50, 50));

        auto *nsSeries = new QScatterSeries();
        nsSeries->setName("Non-significant");
        nsSeries->setMarkerSize(markerSize);
        nsSeries->setColor(QColor(150, 150, 150));
        nsSeries->setBorderColor(QColor(150, 150, 150));

        QList<QPointF> sigPoints;
        QList<QPointF> nsPoints;

        for (int i = 0; i < res.rows(); ++i)
        {
            double baseMean = res(i, 0);
            double lfc = res(i, 1);
            double padj = res(i, 5);

            if (baseMean <= 0 || !std::isfinite(lfc))
                continue;
            double x = std::log10(baseMean);

            if (padj < pThresh)
                sigPoints.append(QPointF(x, lfc));
            else
                nsPoints.append(QPointF(x, lfc));
        }

        if (nsPoints.size() > kInteractivePointLimit)
            nsPoints = downsamplePointsByXBucket(nsPoints, kInteractivePointLimit);

        nsSeries->replace(nsPoints);
        sigSeries->replace(sigPoints);

        // Build the chart
        QChart *chart = new QChart();
        chart->setTitle("MA Plot");
        chart->addSeries(nsSeries);
        chart->addSeries(sigSeries);

        auto *axisX = new QValueAxis();
        axisX->setTitleText("log10(baseMean)");
        chart->addAxis(axisX, Qt::AlignBottom);
        nsSeries->attachAxis(axisX);
        sigSeries->attachAxis(axisX);

        auto *axisY = new QValueAxis();
        axisY->setTitleText("log2 Fold Change");
        chart->addAxis(axisY, Qt::AlignLeft);
        nsSeries->attachAxis(axisY);
        sigSeries->attachAxis(axisY);

        chart->legend()->setVisible(true);
        chart->legend()->setAlignment(Qt::AlignBottom);

        // Replace chart on the view (QChartView takes ownership of previous chart)
        m_plotWidget->setChart(chart);
    }

    void MainWindow::plotVolcanoPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        const auto &res = m_analysisResults.results;
        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        double fcThresh = m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 1.0;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;
        qreal markerSize = static_cast<qreal>(ptSize);

        auto *sigSeries = new QScatterSeries();
        sigSeries->setName("Significant");
        sigSeries->setMarkerSize(markerSize);
        sigSeries->setColor(QColor(220, 50, 50));
        sigSeries->setBorderColor(QColor(220, 50, 50));

        auto *nsSeries = new QScatterSeries();
        nsSeries->setName("Non-significant");
        nsSeries->setMarkerSize(markerSize);
        nsSeries->setColor(QColor(150, 150, 150));
        nsSeries->setBorderColor(QColor(150, 150, 150));

        QList<QPointF> sigPoints;
        QList<QPointF> nsPoints;

        double maxNegLog10P = 0.0;
        double maxAbsLfc = 0.0;

        for (int i = 0; i < res.rows(); ++i)
        {
            double lfc = res(i, 1);
            double pval = res(i, 4);
            double padj = res(i, 5);

            if (pval <= 0 || !std::isfinite(pval) || !std::isfinite(lfc))
                continue;
            double negLog10P = -std::log10(pval);

            if (padj < pThresh && std::abs(lfc) > fcThresh)
                sigPoints.append(QPointF(lfc, negLog10P));
            else
                nsPoints.append(QPointF(lfc, negLog10P));

            maxNegLog10P = std::max(maxNegLog10P, negLog10P);
            maxAbsLfc = std::max(maxAbsLfc, std::abs(lfc));
        }

        if (nsPoints.size() > kInteractivePointLimit)
            nsPoints = downsamplePointsByXBucket(nsPoints, kInteractivePointLimit);

        nsSeries->replace(nsPoints);
        sigSeries->replace(sigPoints);

        // Threshold lines: horizontal at -log10(pThresh), vertical at +/- fcThresh
        double hLine = -std::log10(pThresh);

        auto *pThreshLine = new QLineSeries();
        pThreshLine->setName(QString("p = %1").arg(pThresh));
        pThreshLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));
        pThreshLine->append(-maxAbsLfc * 1.2, hLine);
        pThreshLine->append(maxAbsLfc * 1.2, hLine);

        auto *fcPosLine = new QLineSeries();
        fcPosLine->setName(QString("log2FC = %1").arg(fcThresh));
        fcPosLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));
        fcPosLine->append(fcThresh, 0);
        fcPosLine->append(fcThresh, maxNegLog10P * 1.1);

        auto *fcNegLine = new QLineSeries();
        fcNegLine->setName(QString("log2FC = -%1").arg(fcThresh));
        fcNegLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));
        fcNegLine->append(-fcThresh, 0);
        fcNegLine->append(-fcThresh, maxNegLog10P * 1.1);

        // Build the chart
        QChart *chart = new QChart();
        chart->setTitle("Volcano Plot");
        chart->addSeries(nsSeries);
        chart->addSeries(sigSeries);
        chart->addSeries(pThreshLine);
        chart->addSeries(fcPosLine);
        chart->addSeries(fcNegLine);

        auto *axisX = new QValueAxis();
        axisX->setTitleText("log2 Fold Change");
        chart->addAxis(axisX, Qt::AlignBottom);
        nsSeries->attachAxis(axisX);
        sigSeries->attachAxis(axisX);
        pThreshLine->attachAxis(axisX);
        fcPosLine->attachAxis(axisX);
        fcNegLine->attachAxis(axisX);

        auto *axisY = new QValueAxis();
        axisY->setTitleText("-log10(p-value)");
        chart->addAxis(axisY, Qt::AlignLeft);
        nsSeries->attachAxis(axisY);
        sigSeries->attachAxis(axisY);
        pThreshLine->attachAxis(axisY);
        fcPosLine->attachAxis(axisY);
        fcNegLine->attachAxis(axisY);

        chart->legend()->setVisible(true);
        chart->legend()->setAlignment(Qt::AlignBottom);

        // Hide threshold line legend markers (keep only Significant/Non-significant visible)
        for (auto *marker : chart->legend()->markers())
        {
            if (marker->series() == pThreshLine ||
                marker->series() == fcPosLine ||
                marker->series() == fcNegLine)
                marker->setVisible(false);
        }

        m_plotWidget->setChart(chart);
    }

    void MainWindow::plotDispersionPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        const auto &baseMeans = m_analysisResults.baseMeans;
        const auto &genewise = m_analysisResults.genewiseDispersions;
        const auto &fitted = m_analysisResults.fittedDispersions;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;
        qreal markerSize = static_cast<qreal>(ptSize);

        bool hasData = baseMeans.size() > 0 && genewise.size() > 0 && fitted.size() > 0;
        if (!hasData)
        {
            addProgressMessage("No dispersion data available for plotting.");
            return;
        }

        // Genewise dispersions scatter
        auto *geneSeries = new QScatterSeries();
        geneSeries->setName("Genewise");
        geneSeries->setMarkerSize(markerSize);
        geneSeries->setColor(QColor(70, 130, 220));
        geneSeries->setBorderColor(QColor(70, 130, 220));

        QList<QPointF> genePoints;

        for (int i = 0; i < baseMeans.size(); ++i)
        {
            if (baseMeans(i) > 0 && genewise(i) > 0)
            {
                genePoints.append(QPointF(std::log10(baseMeans(i)), std::log10(genewise(i))));
            }
        }

        if (genePoints.size() > kInteractivePointLimit)
            genePoints = downsamplePointsByXBucket(genePoints, kInteractivePointLimit);
        geneSeries->replace(genePoints);

        // Fitted trend line (sorted by x for connected line)
        std::vector<std::pair<double, double>> fitPoints;
        for (int i = 0; i < baseMeans.size(); ++i)
        {
            if (baseMeans(i) > 0 && fitted(i) > 0)
            {
                fitPoints.push_back({std::log10(baseMeans(i)), std::log10(fitted(i))});
            }
        }
        std::sort(fitPoints.begin(), fitPoints.end());

        auto *fitSeries = new QLineSeries();
        fitSeries->setName("Fitted trend");
        fitSeries->setPen(QPen(QColor(220, 50, 50), 2));

        QList<QPointF> fitLinePoints;
        fitLinePoints.reserve(static_cast<int>(fitPoints.size()));
        for (const auto &pt : fitPoints)
            fitLinePoints.append(QPointF(pt.first, pt.second));
        fitSeries->replace(fitLinePoints);

        // Build the chart
        QChart *chart = new QChart();
        chart->setTitle("Dispersion Plot");
        chart->addSeries(geneSeries);
        chart->addSeries(fitSeries);

        auto *axisX = new QValueAxis();
        axisX->setTitleText("log10(baseMean)");
        chart->addAxis(axisX, Qt::AlignBottom);
        geneSeries->attachAxis(axisX);
        fitSeries->attachAxis(axisX);

        auto *axisY = new QValueAxis();
        axisY->setTitleText("log10(dispersion)");
        chart->addAxis(axisY, Qt::AlignLeft);
        geneSeries->attachAxis(axisY);
        fitSeries->attachAxis(axisY);

        chart->legend()->setVisible(true);
        chart->legend()->setAlignment(Qt::AlignBottom);

        m_plotWidget->setChart(chart);
    }

    void MainWindow::plotThreeWayScatter()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        // Requires at least 2 contrasts from a three-way comparison
        if (m_analysisResults.contrastResults.size() < 2)
        {
            addProgressMessage("Three-way scatter requires a three-way comparison (at least 2 contrasts). "
                               "Run a three-way analysis first.");

            // Show an informative empty chart
            QChart *chart = new QChart();
            chart->setTitle("Three-way Scatter: No multi-contrast data available");
            m_plotWidget->setChart(chart);
            return;
        }

        const auto &res1 = m_analysisResults.contrastResults[0]; // contrast 1: log2FC on x-axis
        const auto &res2 = m_analysisResults.contrastResults[1]; // contrast 2: log2FC on y-axis

        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;
        qreal markerSize = static_cast<qreal>(ptSize);

        // Three categories: significant in both, significant in one only, neither
        auto *bothSeries = new QScatterSeries();
        bothSeries->setName("Significant in both");
        bothSeries->setMarkerSize(markerSize);
        bothSeries->setColor(QColor(220, 50, 50));   // Red
        bothSeries->setBorderColor(QColor(220, 50, 50));

        auto *oneSeries = new QScatterSeries();
        oneSeries->setName("Significant in one");
        oneSeries->setMarkerSize(markerSize);
        oneSeries->setColor(QColor(50, 100, 220));    // Blue
        oneSeries->setBorderColor(QColor(50, 100, 220));

        auto *nsSeries = new QScatterSeries();
        nsSeries->setName("Not significant");
        nsSeries->setMarkerSize(markerSize);
        nsSeries->setColor(QColor(180, 180, 180));    // Gray
        nsSeries->setBorderColor(QColor(180, 180, 180));

        QList<QPointF> bothPoints;
        QList<QPointF> onePoints;
        QList<QPointF> nsPoints;

        int nGenes = std::min(static_cast<int>(res1.rows()), static_cast<int>(res2.rows()));
        double maxAbsX = 0.0;
        double maxAbsY = 0.0;

        for (int i = 0; i < nGenes; ++i)
        {
            double lfc1 = res1(i, 1); // log2FC from contrast 1
            double lfc2 = res2(i, 1); // log2FC from contrast 2
            double padj1 = res1(i, 5); // adjusted p-value from contrast 1
            double padj2 = res2(i, 5); // adjusted p-value from contrast 2

            if (!std::isfinite(lfc1) || !std::isfinite(lfc2))
                continue;

            bool sig1 = !std::isnan(padj1) && padj1 < pThresh;
            bool sig2 = !std::isnan(padj2) && padj2 < pThresh;

            if (sig1 && sig2)
                bothPoints.append(QPointF(lfc1, lfc2));
            else if (sig1 || sig2)
                onePoints.append(QPointF(lfc1, lfc2));
            else
                nsPoints.append(QPointF(lfc1, lfc2));

            maxAbsX = std::max(maxAbsX, std::abs(lfc1));
            maxAbsY = std::max(maxAbsY, std::abs(lfc2));
        }

        if (nsPoints.size() > kInteractivePointLimit)
            nsPoints = downsamplePointsByXBucket(nsPoints, kInteractivePointLimit);
        nsSeries->replace(nsPoints);
        oneSeries->replace(onePoints);
        bothSeries->replace(bothPoints);

        // Reference lines at x=0 and y=0
        double xExtent = maxAbsX * 1.2;
        double yExtent = maxAbsY * 1.2;

        auto *hRefLine = new QLineSeries();
        hRefLine->setName("y = 0");
        hRefLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));
        hRefLine->append(-xExtent, 0.0);
        hRefLine->append(xExtent, 0.0);

        auto *vRefLine = new QLineSeries();
        vRefLine->setName("x = 0");
        vRefLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));
        vRefLine->append(0.0, -yExtent);
        vRefLine->append(0.0, yExtent);

        // Build chart
        QChart *chart = new QChart();
        chart->setTitle("Three-way Scatter");
        chart->addSeries(nsSeries);
        chart->addSeries(oneSeries);
        chart->addSeries(bothSeries);
        chart->addSeries(hRefLine);
        chart->addSeries(vRefLine);

        // Determine axis labels from contrast labels if available
        QString xLabel = "log2FC (Contrast 1)";
        QString yLabel = "log2FC (Contrast 2)";
        if (m_analysisResults.contrastLabels.size() >= 2)
        {
            xLabel = QString("log2FC (%1)").arg(QString::fromStdString(m_analysisResults.contrastLabels[0]));
            yLabel = QString("log2FC (%1)").arg(QString::fromStdString(m_analysisResults.contrastLabels[1]));
        }

        auto *axisX = new QValueAxis();
        axisX->setTitleText(xLabel);
        chart->addAxis(axisX, Qt::AlignBottom);
        nsSeries->attachAxis(axisX);
        oneSeries->attachAxis(axisX);
        bothSeries->attachAxis(axisX);
        hRefLine->attachAxis(axisX);
        vRefLine->attachAxis(axisX);

        auto *axisY = new QValueAxis();
        axisY->setTitleText(yLabel);
        chart->addAxis(axisY, Qt::AlignLeft);
        nsSeries->attachAxis(axisY);
        oneSeries->attachAxis(axisY);
        bothSeries->attachAxis(axisY);
        hRefLine->attachAxis(axisY);
        vRefLine->attachAxis(axisY);

        chart->legend()->setVisible(true);
        chart->legend()->setAlignment(Qt::AlignBottom);

        // Hide reference line legend markers
        for (auto *marker : chart->legend()->markers())
        {
            if (marker->series() == hRefLine || marker->series() == vRefLine)
                marker->setVisible(false);
        }

        m_plotWidget->setChart(chart);
    }

    void MainWindow::onSavePlot()
    {
        if (!m_plotWidget || !m_plotWidget->chart())
            return;

        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Save Plot",
            m_lastUsedDirectory,
            "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            bool ok = false;
            if (fileName.endsWith(".svg", Qt::CaseInsensitive))
            {
                // SVG export using QSvgGenerator
                QSvgGenerator generator;
                generator.setFileName(fileName);
                generator.setSize(QSize(800, 600));
                generator.setViewBox(QRect(0, 0, 800, 600));
                generator.setTitle("StatMaker++ Plot");
                generator.setDescription("Exported from StatMaker++ analysis");
                QPainter painter(&generator);
                m_plotWidget->render(&painter);
                painter.end();
                ok = true;
            }
            else if (fileName.endsWith(".pdf", Qt::CaseInsensitive))
            {
                QPrinter printer(QPrinter::HighResolution);
                printer.setOutputFormat(QPrinter::PdfFormat);
                printer.setOutputFileName(fileName);
                QPainter painter(&printer);
                m_plotWidget->render(&painter);
                painter.end();
                ok = true;
            }
            else
            {
                // PNG export with configurable DPI
                bool dpiOk = false;
                int dpi = QInputDialog::getInt(
                    this,
                    "PNG Export DPI",
                    "Enter export DPI (72 = screen, 150 = print, 300 = publication):",
                    150,   // default
                    72,    // min
                    600,   // max
                    1,     // step
                    &dpiOk);

                if (!dpiOk)
                    return; // User cancelled

                // Calculate scale factor relative to screen DPI (72)
                double scaleFactor = static_cast<double>(dpi) / 72.0;
                int w = static_cast<int>(800 * scaleFactor);
                int h = static_cast<int>(600 * scaleFactor);

                QPixmap pixmap(w, h);
                pixmap.fill(Qt::white);
                QPainter painter(&pixmap);
                painter.setRenderHint(QPainter::Antialiasing);
                m_plotWidget->render(&painter, QRect(0, 0, w, h));
                painter.end();
                ok = pixmap.save(fileName, "PNG");

                if (ok)
                    addProgressMessage(QString("Plot saved at %1 DPI (%2x%3 px) to: %4")
                                           .arg(dpi).arg(w).arg(h).arg(fileName));
            }

            if (ok && !fileName.endsWith(".png", Qt::CaseInsensitive))
                addProgressMessage("Plot saved to: " + fileName);
            else if (!ok)
                addProgressMessage("Failed to save plot.");
        }
    }

    void MainWindow::onExportPlotData()
    {
        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Export Plot Data",
            m_lastUsedDirectory,
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (m_analysisResults.isValid)
            {
                exportResultsToCSV(fileName, m_analysisResults, false);
                addProgressMessage("Plot data exported to: " + fileName);
            }
            else
            {
                addProgressMessage("No analysis results to export.");
            }
        }
    }

    void MainWindow::onPrintPlot()
    {
        if (!m_plotWidget || !m_plotWidget->chart())
            return;

        // Save to a temp PDF and inform the user
        QString tempPath = QDir::tempPath() + "/statmaker_plot.pdf";
        QPrinter printer(QPrinter::HighResolution);
        printer.setOutputFormat(QPrinter::PdfFormat);
        printer.setOutputFileName(tempPath);
        QPainter painter(&printer);
        m_plotWidget->render(&painter);
        painter.end();

        addProgressMessage("Plot saved to temporary PDF: " + tempPath);
    }

    // Slot implementations - Progress output
    void MainWindow::onClearOutput()
    {
        if (m_progressOutputText)
        {
            m_progressOutputText->clear();
        }
    }

    void MainWindow::onSaveOutput()
    {
        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Save Output",
            m_lastUsedDirectory,
            "Text Files (*.txt);;All Files (*)");

        if (!fileName.isEmpty() && m_progressOutputText)
        {
            // Update last used directory
            updateLastUsedDirectory(fileName);

            QFile file(fileName);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QTextStream stream(&file);
                stream << m_progressOutputText->toPlainText();
                addProgressMessage("Output saved successfully.");
            }
            else
            {
                QMessageBox::warning(static_cast<QWidget *>(this), "Error", "Could not save output file.");
            }
        }
    }

    // Slot implementations - Menu actions
    void MainWindow::onOpenCounts()
    {
        onAddGeneCountFiles();
    }

    void MainWindow::onOpenMetadata()
    {
        onGenerateData();
    }

    void MainWindow::onSaveResults()
    {
        onExportResults();
    }

    void MainWindow::onExportVisualizations()
    {
        onSavePlot();
    }

    void MainWindow::onExit()
    {
        close();
    }

    void MainWindow::onRunAnalysisAction()
    {
        onRunAnalysis();
    }

    void MainWindow::onStopAnalysisAction()
    {
        onStopAnalysis();
    }

    void MainWindow::onResetAnalysisAction()
    {
        onResetAnalysis();
    }

    void MainWindow::onInputTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(0);
        }
    }

    void MainWindow::onAnalysisTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab (no analysis tab exists)
        }
    }

    void MainWindow::onResultsTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab
        }
    }

    void MainWindow::onVisualizationTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(2); // Visualization tab
        }
        refreshPlot();
    }

    void MainWindow::onAbout()
    {
        QMessageBox::about(static_cast<QWidget *>(this), "About StatMaker++",
                           "StatMaker++ Y2H Statistical Analysis\n\n"
                           "Version 1.0.0\n"
                           "A Qt6-based interface for DESeq2 and Y2H-SCORES analysis.");
    }

    void MainWindow::onUserGuide()
    {
        QMessageBox::information(static_cast<QWidget *>(this), "User Guide",
                                 "User guide functionality will be implemented in future versions.");
    }

    void MainWindow::updateCountsPreviewWithConvertedValues(QTableWidget *table, const CombinedData &data)
    {
        if (!table || !data.isValid)
        {
            return;
        }

        const int maxRows = std::min(20, static_cast<int>(data.geneNames.size()));
        const int maxCols = std::min(10, static_cast<int>(data.sampleNames.size()));

        table->clear();
        table->setRowCount(maxRows);
        table->setColumnCount(maxCols + 1); // +1 for gene names

        QStringList headers = {"Gene"};
        for (int i = 0; i < maxCols; ++i)
        {
            headers.append(data.sampleNames[i]);
        }
        table->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < maxRows; ++row)
        {
            // Set gene name in first column
            table->setItem(row, 0, new QTableWidgetItem(data.geneNames[row]));

            for (int col = 0; col < maxCols; ++col)
            {
                // Use raw integer counts directly
                int integerCount = data.rawCountMatrix[row][col];

                QTableWidgetItem *item = new QTableWidgetItem(QString::number(integerCount));
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);
                table->setItem(row, col + 1, item);
            }
        }

        table->resizeColumnsToContents();
    }

    // =========================================================================
    // Results table context menu slots
    // =========================================================================

    void MainWindow::onResultsContextMenu(const QPoint &pos)
    {
        QTableWidgetItem *item = m_resultsTable->itemAt(pos);
        if (!item)
            return;

        QMenu menu(this);
        menu.addAction("Open in MultiQuery++", this, &MainWindow::onOpenInMultiQuery);
        menu.addAction("Open in ReadDepth++", this, &MainWindow::onOpenInReadDepth);
        menu.addSeparator();
        menu.addAction("Copy Gene Name", this, &MainWindow::onCopyGeneName);
        menu.exec(m_resultsTable->viewport()->mapToGlobal(pos));
    }

    void MainWindow::onOpenInMultiQuery()
    {
        int row = m_resultsTable->currentRow();
        if (row < 0 || !m_resultsTable->item(row, 0))
            return;

        QString geneName = m_resultsTable->item(row, 0)->text();

        // Resolve path relative to the bundled app location:
        // StatMaker++.app/Contents/MacOS (applicationDirPath)
        //   -> cdUp -> Contents
        //   -> cdUp -> StatMaker++.app (which lives in DEEPN++.app/Contents/Resources/)
        // Sibling app: MultiQuery++.app/Contents/MacOS/MultiQuery++
        QDir appDir(QCoreApplication::applicationDirPath());
        appDir.cdUp(); // Contents
        appDir.cdUp(); // StatMaker++.app -> now in DEEPN++.app/Contents/Resources/
        QString multiQueryPath = appDir.absoluteFilePath(
            "MultiQuery++.app/Contents/MacOS/MultiQuery++");

        if (!QFile::exists(multiQueryPath))
        {
            QMessageBox::warning(this, "Application Not Found",
                                 "Could not find MultiQuery++ at:\n" + multiQueryPath);
            return;
        }

        const QString workdir = resolveWorkingDirectory();
        QStringList datasets = discoverFiles(QDir(workdir).filePath("analyzed_files"), {"*.sqlite", "*.db"});
        datasets.erase(
            std::remove_if(datasets.begin(), datasets.end(),
                           [](const QString &path)
                           {
                               return QFileInfo(path).fileName().compare(
                                          "statmaker_results.sqlite", Qt::CaseInsensitive) == 0;
                           }),
            datasets.end());

        QStringList arguments;
        if (!workdir.isEmpty())
            arguments << "--workdir" << workdir;
        if (!datasets.isEmpty())
            arguments << "--datasets" << datasets.join(",");
        arguments << "--gene" << geneName;

        QProcess::startDetached(multiQueryPath, arguments);
    }

    void MainWindow::onOpenInReadDepth()
    {
        int row = m_resultsTable->currentRow();
        if (row < 0 || !m_resultsTable->item(row, 0))
            return;

        QString geneName = m_resultsTable->item(row, 0)->text();

        // Same path resolution as MultiQuery++
        QDir appDir(QCoreApplication::applicationDirPath());
        appDir.cdUp(); // Contents
        appDir.cdUp(); // StatMaker++.app -> now in DEEPN++.app/Contents/Resources/
        QString readDepthPath = appDir.absoluteFilePath(
            "ReadDepth++.app/Contents/MacOS/ReadDepth++");

        if (!QFile::exists(readDepthPath))
        {
            QMessageBox::warning(this, "Application Not Found",
                                 "Could not find ReadDepth++ at:\n" + readDepthPath);
            return;
        }

        const QString workdir = resolveWorkingDirectory();
        QStringList datasets = discoverFiles(QDir(workdir).filePath("analyzed_files"), {"*.sqlite", "*.db"});
        datasets.erase(
            std::remove_if(datasets.begin(), datasets.end(),
                           [](const QString &path)
                           {
                               return QFileInfo(path).fileName().compare(
                                          "statmaker_results.sqlite", Qt::CaseInsensitive) == 0;
                           }),
            datasets.end());

        QStringList arguments;
        if (!workdir.isEmpty())
            arguments << "--workdir" << workdir;
        if (!datasets.isEmpty())
            arguments << "--datasets" << datasets.join(",");
        arguments << "--gene" << geneName;

        QProcess::startDetached(readDepthPath, arguments);
    }

    void MainWindow::onCopyGeneName()
    {
        int row = m_resultsTable->currentRow();
        if (row >= 0 && m_resultsTable->item(row, 0))
        {
            QString geneName = m_resultsTable->item(row, 0)->text();
            QApplication::clipboard()->setText(geneName);
        }
    }

    // =========================================================================
    // Junction files browsing (Task 4)
    // =========================================================================

    void MainWindow::onBrowseJunctionFiles()
    {
        QStringList fileNames = QFileDialog::getOpenFileNames(
            this,
            "Select Junction SQLite Files",
            m_lastUsedDirectory,
            "SQLite Files (*.sqlite *.db);;All Files (*)");

        if (!fileNames.isEmpty())
        {
            updateLastUsedDirectory(fileNames.first());

            // Display selected files as semicolon-separated list
            if (m_junctionFilesEdit)
            {
                m_junctionFilesEdit->setText(fileNames.join("; "));
            }

            addProgressMessage(QString("Selected %1 junction SQLite file(s).").arg(fileNames.size()));
        }
    }

    // =========================================================================
    // Contrast selector for three-way comparisons (Task 5)
    // =========================================================================

    void MainWindow::onContrastChanged()
    {
        if (!m_contrastSelectorCombo || !m_analysisResults.isValid)
            return;

        int contrastIdx = m_contrastSelectorCombo->currentIndex();

        // Check if we have multi-contrast results
        if (m_analysisResults.contrastResults.empty())
        {
            addProgressMessage("No multi-contrast results available. "
                               "Run a three-way comparison to use the contrast selector.");
            return;
        }

        if (contrastIdx < 0 || contrastIdx >= static_cast<int>(m_analysisResults.contrastResults.size()))
        {
            addProgressMessage(QString("Invalid contrast index: %1").arg(contrastIdx));
            return;
        }

        addProgressMessage(QString("Switching to contrast: %1").arg(m_contrastSelectorCombo->currentText()));

        // Update the active contrast
        m_analysisResults.activeContrast = contrastIdx;
        if (contrastIdx < static_cast<int>(m_analysisResults.contrastLabels.size()))
            m_analysisResults.activeContrastLabel = QString::fromStdString(m_analysisResults.contrastLabels[contrastIdx]);

        // Replace the main results matrix with the selected contrast's results
        m_analysisResults.results = m_analysisResults.contrastResults[contrastIdx];

        // Recalculate summary statistics
        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        m_analysisResults.significantGenes = 0;
        m_analysisResults.upregulatedGenes = 0;
        m_analysisResults.downregulatedGenes = 0;

        const auto &res = m_analysisResults.results;
        for (int i = 0; i < res.rows(); ++i)
        {
            double padj = res(i, 5);
            double lfc = res(i, 1);
            if (!std::isnan(padj) && padj < pThresh)
            {
                m_analysisResults.significantGenes++;
                if (lfc > 0)
                    m_analysisResults.upregulatedGenes++;
                else if (lfc < 0)
                    m_analysisResults.downregulatedGenes++;
            }
        }

        // Refresh UI
        m_resultsColumnsSized = false;
        invalidatePlotCache();
        updateResultsTable(m_analysisResults);
        updateResultsSummary(m_analysisResults);
        refreshPlot();

        addProgressMessage(QString("Contrast '%1': %2 significant genes (%3 up, %4 down)")
                               .arg(m_contrastSelectorCombo->currentText())
                               .arg(m_analysisResults.significantGenes)
                               .arg(m_analysisResults.upregulatedGenes)
                               .arg(m_analysisResults.downregulatedGenes));
    }

} // namespace deseq2

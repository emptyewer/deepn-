#include "batch_runner.h"
#include "depth_calculator.h"
#include "export_engine.h"
#include "sqlite_junction_loader.h"

#include <QDebug>
#include <QDir>
#include <QFileInfo>

namespace deepn {

BatchRunner::BatchRunner(QObject* parent)
    : QObject(parent)
{
}

BatchRunner::~BatchRunner()
{
    stop();
}

void BatchRunner::start(const Config& config)
{
    if (m_thread && m_thread->isRunning()) {
        qDebug() << "BatchRunner already running";
        return;
    }

    m_config = config;
    m_shouldStop = false;

    m_thread = new QThread(this);

    // Use moveToThread pattern: create a worker lambda on the thread
    connect(m_thread, &QThread::started, this, &BatchRunner::runBatch);
    connect(m_thread, &QThread::finished, m_thread, &QThread::deleteLater);

    // We run the batch on the current thread but in a worker thread context
    // Actually, use a simpler approach: run in a separate thread via QThread::create
    if (m_thread) {
        delete m_thread;
        m_thread = nullptr;
    }

    m_thread = QThread::create([this]() {
        runBatch();
    });
    connect(m_thread, &QThread::finished, m_thread, &QThread::deleteLater);
    m_thread->start();
}

void BatchRunner::stop()
{
    m_shouldStop = true;
    if (m_thread && m_thread->isRunning()) {
        m_thread->quit();
        m_thread->wait(5000);
    }
    m_thread = nullptr;
}

bool BatchRunner::isRunning() const
{
    return m_thread && m_thread->isRunning();
}

void BatchRunner::runBatch()
{
    // Validate config
    if (m_config.dbPaths.isEmpty()) {
        emit error("No database paths specified");
        emit finished(false);
        return;
    }

    if (m_config.geneList.isEmpty()) {
        emit error("No genes specified");
        emit finished(false);
        return;
    }

    // Create output directory
    QDir outDir(m_config.outputDir);
    if (!outDir.exists()) {
        if (!outDir.mkpath(".")) {
            emit error(QStringLiteral("Cannot create output directory: %1").arg(m_config.outputDir));
            emit finished(false);
            return;
        }
    }

    // Limit to topN genes
    QStringList genesToProcess = m_config.geneList;
    if (m_config.topN > 0 && genesToProcess.size() > m_config.topN) {
        genesToProcess = genesToProcess.mid(0, m_config.topN);
    }

    int totalTasks = genesToProcess.size() * m_config.dbPaths.size();
    int currentTask = 0;

    QVector<GeneJunctionProfile> allProfiles;

    for (const QString& dbPath : m_config.dbPaths) {
        if (m_shouldStop)
            break;

        SqliteJunctionLoader loader;
        if (!loader.open(dbPath)) {
            emit error(QStringLiteral("Cannot open database: %1 - %2").arg(dbPath, loader.lastError()));
            continue;
        }

        QString dbBaseName = QFileInfo(dbPath).completeBaseName();

        for (const QString& gene : genesToProcess) {
            if (m_shouldStop)
                break;

            currentTask++;
            emit progressChanged(currentTask, totalTasks, gene);

            // Load junction data
            GeneJunctionProfile profile = loader.loadGeneJunctions(gene);
            profile.datasetLabel = dbBaseName;

            if (profile.sites.isEmpty()) {
                qDebug() << "No junctions found for gene:" << gene << "in" << dbBaseName;
                continue;
            }

            allProfiles.append(profile);

            if (m_config.exportCSV) {
                // Export junction CSV
                QString junctionFile = outDir.filePath(
                    QStringLiteral("%1_%2_junctions.csv").arg(gene, dbBaseName));
                if (!ExportEngine::exportJunctionCSV(junctionFile, profile)) {
                    qDebug() << "Failed to export junction CSV:" << ExportEngine::lastError();
                }

                // Export collapsed CSV
                auto collapsed = SqliteJunctionLoader::collapseByPosition(profile.sites);
                QString collapsedFile = outDir.filePath(
                    QStringLiteral("%1_%2_collapsed.csv").arg(gene, dbBaseName));
                if (!ExportEngine::exportCollapsedCSV(collapsedFile, collapsed)) {
                    qDebug() << "Failed to export collapsed CSV:" << ExportEngine::lastError();
                }

                // Export depth CSV if annotation available
                if (profile.annotation.isValid()) {
                    DepthCalculator calc;
                    DepthProfile depth = calc.calculate(dbPath, gene, profile.annotation);
                    if (!depth.points.isEmpty()) {
                        QString depthFile = outDir.filePath(
                            QStringLiteral("%1_%2_depth.csv").arg(gene, dbBaseName));
                        if (!ExportEngine::exportDepthCSV(depthFile, depth)) {
                            qDebug() << "Failed to export depth CSV:" << ExportEngine::lastError();
                        }
                    }
                }
            }
        }

        loader.close();
    }

    // Export batch summary
    if (m_config.exportCSV && !allProfiles.isEmpty()) {
        QString summaryFile = outDir.filePath("batch_summary.csv");
        if (!ExportEngine::exportBatchSummaryCSV(summaryFile, allProfiles)) {
            qDebug() << "Failed to export batch summary:" << ExportEngine::lastError();
        }
    }

    if (m_shouldStop) {
        qDebug() << "Batch processing stopped by user";
        emit finished(false);
    } else {
        qDebug() << "Batch processing complete:" << currentTask << "tasks";
        emit finished(true);
    }
}

}  // namespace deepn

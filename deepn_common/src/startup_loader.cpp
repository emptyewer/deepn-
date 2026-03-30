#include "startup_loader.h"
#include "csv_utils.h"

#include <QDir>
#include <QFileInfo>

namespace deepn {

StartupLoader::StartupLoader(QObject* parent)
    : QObject(parent)
{
    moveToThread(&m_thread);
    connect(&m_thread, &QThread::started, this, &StartupLoader::doWork);
}

void StartupLoader::start(const Config& config)
{
    m_config = config;
    m_thread.start();
}

void StartupLoader::doWork()
{
    // 1. Load gene reference
    if (!m_config.generefPath.isEmpty()) {
        emit message("Loading gene reference: " + QFileInfo(m_config.generefPath).fileName());
        if (m_config.generefPath.endsWith(".sqlite", Qt::CaseInsensitive)) {
            m_annotationDB.loadFromSqlite(m_config.generefPath);
        } else {
            m_annotationDB.loadFromFasta(m_config.generefPath);
        }
    }

    // 2. Load DESeq2 results
    if (!m_config.resultsCSV.isEmpty()) {
        emit message("Loading DESeq2 results...");
        m_deseq2Results = parseDESeq2CSV(m_config.resultsCSV);
    }

    // 3. Discover datasets in working directory
    if (!m_config.workdir.isEmpty()) {
        emit message("Discovering datasets...");
        QDir analyzedDir(QDir(m_config.workdir).filePath("analyzed_files"));
        if (analyzedDir.exists()) {
            QFileInfoList dbFiles = analyzedDir.entryInfoList(
                {"*.sqlite", "*.db"}, QDir::Files, QDir::Name);
            for (const QFileInfo& fi : dbFiles) {
                if (fi.fileName().compare("statmaker_results.sqlite", Qt::CaseInsensitive) == 0)
                    continue;
                m_discoveredDatasets.append(fi.absoluteFilePath());
            }
        }

        // Discover DESeq2 results if not already provided
        if (m_deseq2Results.isEmpty()) {
            QDir dir(m_config.workdir);
            QStringList csvDirs = {m_config.workdir, dir.filePath("results")};
            for (const QString& csvDir : csvDirs) {
                QDir d(csvDir);
                if (!d.exists()) continue;
                QFileInfoList csvFiles = d.entryInfoList({"*.csv"}, QDir::Files);
                for (const QFileInfo& fi : csvFiles) {
                    if (fi.fileName().contains("deseq2", Qt::CaseInsensitive)) {
                        emit message("Loading DESeq2 results: " + fi.fileName());
                        m_deseq2Results = parseDESeq2CSV(fi.absoluteFilePath());
                        break;
                    }
                }
                if (!m_deseq2Results.isEmpty()) break;
            }
        }
    }

    // 4. Add explicit dataset paths
    for (const QString& path : m_config.datasetPaths) {
        QString trimmed = path.trimmed();
        if (!trimmed.isEmpty() && !m_discoveredDatasets.contains(trimmed)) {
            m_discoveredDatasets.append(trimmed);
        }
    }

    m_initialGene = m_config.geneName;

    emit message(QString("Ready (%1 genes, %2 datasets).")
                     .arg(m_annotationDB.count())
                     .arg(m_discoveredDatasets.size()));
    emit finished();
    m_thread.quit();
}

}  // namespace deepn

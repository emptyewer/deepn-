#ifndef DEEPN_STARTUP_LOADER_H
#define DEEPN_STARTUP_LOADER_H

#include "data_structures.h"
#include "gene_annotation_db.h"

#include <QObject>
#include <QStringList>
#include <QThread>

namespace deepn {

class StartupLoader : public QObject {
    Q_OBJECT
public:
    struct Config {
        QString generefPath;      // .sqlite or .fasta
        QString workdir;
        QString resultsCSV;
        QStringList datasetPaths; // comma-separated from CLI
        QString geneName;
        QString junctionsCSV;     // ReadDepth++ only
    };

    explicit StartupLoader(QObject* parent = nullptr);

    void start(const Config& config);

    // Results accessible after finished()
    GeneAnnotationDB& annotationDB() { return m_annotationDB; }
    QVector<DESeq2Result>& deseq2Results() { return m_deseq2Results; }
    QStringList& discoveredDatasets() { return m_discoveredDatasets; }
    QString& initialGene() { return m_initialGene; }

signals:
    void message(const QString& msg);
    void finished();

private slots:
    void doWork();

private:
    Config m_config;
    GeneAnnotationDB m_annotationDB;
    QVector<DESeq2Result> m_deseq2Results;
    QStringList m_discoveredDatasets;
    QString m_initialGene;
    QThread m_thread;
};

}  // namespace deepn

#endif  // DEEPN_STARTUP_LOADER_H

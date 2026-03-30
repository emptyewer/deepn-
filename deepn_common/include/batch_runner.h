#ifndef DEEPN_BATCH_RUNNER_H
#define DEEPN_BATCH_RUNNER_H

#include "data_structures.h"

#include <QObject>
#include <QStringList>
#include <QThread>

namespace deepn {

class BatchRunner : public QObject {
    Q_OBJECT
public:
    struct Config {
        QStringList dbPaths;         // SQLite databases to process
        QStringList geneList;        // genes to process (from DESeq2 results)
        QString outputDir;
        int topN = 50;
        double pValueCutoff = 0.05;
        double log2FCCutoff = 1.0;
        bool exportCSV = true;
        bool exportFigures = false;
    };

    explicit BatchRunner(QObject* parent = nullptr);
    ~BatchRunner() override;

    void start(const Config& config);
    void stop();
    bool isRunning() const;

signals:
    void progressChanged(int current, int total, const QString& currentGene);
    void finished(bool success);
    void error(const QString& message);

private:
    void runBatch();
    QThread* m_thread = nullptr;
    Config m_config;
    bool m_shouldStop = false;
};

}  // namespace deepn

#endif  // DEEPN_BATCH_RUNNER_H

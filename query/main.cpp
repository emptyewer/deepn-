#include "mainwindow.h"

#include <loading_overlay.h>
#include <startup_loader.h>

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("MultiQuery++");
    app.setOrganizationName("DEEPN++");
    app.setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Junction Analysis & Visualization Dashboard");
    parser.addHelpOption();
    parser.addVersionOption();

    QCommandLineOption workdirOpt("workdir", "Working directory.", "path");
    QCommandLineOption datasetsOpt("datasets", "Comma-separated SQLite paths.", "paths");
    QCommandLineOption generefOpt("generef", "Gene reference file (.sqlite or .fasta).", "path");
    QCommandLineOption resultsOpt("results", "DESeq2 results CSV.", "path");
    QCommandLineOption geneOpt("gene", "Initial gene to display.", "name");
    parser.addOption(workdirOpt);
    parser.addOption(datasetsOpt);
    parser.addOption(generefOpt);
    parser.addOption(resultsOpt);
    parser.addOption(geneOpt);
    parser.process(app);

    MainWindow w;
    w.show();

    auto* overlay = w.loadingOverlay();
    overlay->show("Loading MultiQuery++...");

    // Configure background loader
    deepn::StartupLoader loader;
    deepn::StartupLoader::Config config;
    config.generefPath = parser.value(generefOpt);
    config.workdir = parser.value(workdirOpt);
    config.resultsCSV = parser.value(resultsOpt);
    config.geneName = parser.value(geneOpt);
    if (parser.isSet(datasetsOpt)) {
        config.datasetPaths = parser.value(datasetsOpt).split(',', Qt::SkipEmptyParts);
    }

    // Stream messages to overlay
    QObject::connect(&loader, &deepn::StartupLoader::message,
                     overlay, &deepn::LoadingOverlay::addMessage);

    // When loading finishes, apply results to the window on the main thread
    QObject::connect(&loader, &deepn::StartupLoader::finished, &w, [&]() {
        // Load annotation DB on main thread (fast — metadata only, no sequences)
        if (!config.generefPath.isEmpty()) {
            if (config.generefPath.endsWith(".sqlite", Qt::CaseInsensitive)) {
                w.loadGeneReferenceSqlite(config.generefPath);
            } else {
                w.loadGeneReference(config.generefPath);
            }
        }

        w.setDESeq2Results(loader.deseq2Results());

        for (const QString& path : loader.discoveredDatasets()) {
            w.loadDataset(path);
        }

        w.refreshGeneSelector();

        if (!loader.initialGene().isEmpty()) {
            w.displayGene(loader.initialGene());
        }

        overlay->hide();
    });

    loader.start(config);

    return app.exec();
}

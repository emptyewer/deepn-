#include "mainwindow.h"

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

    QCommandLineOption workdirOpt(
        "workdir",
        "Working directory containing experiment data.",
        "path");
    parser.addOption(workdirOpt);

    QCommandLineOption datasetsOpt(
        "datasets",
        "Comma-separated SQLite database paths.",
        "paths");
    parser.addOption(datasetsOpt);

    QCommandLineOption generefOpt(
        "generef",
        "Path to FASTA gene reference file.",
        "path");
    parser.addOption(generefOpt);

    QCommandLineOption resultsOpt(
        "results",
        "DESeq2 results CSV path.",
        "path");
    parser.addOption(resultsOpt);

    QCommandLineOption geneOpt(
        "gene",
        "Initial gene to display.",
        "name");
    parser.addOption(geneOpt);

    parser.process(app);

    MainWindow w;
    w.show();

    // Apply CLI arguments
    if (parser.isSet(generefOpt)) {
        w.loadGeneReference(parser.value(generefOpt));
    }

    if (parser.isSet(resultsOpt)) {
        w.loadDESeq2Results(parser.value(resultsOpt));
    }

    if (parser.isSet(workdirOpt)) {
        w.loadWorkingDirectory(parser.value(workdirOpt));
    }

    if (parser.isSet(datasetsOpt)) {
        const QStringList paths = parser.value(datasetsOpt).split(',', Qt::SkipEmptyParts);
        for (const QString& path : paths) {
            w.loadDataset(path.trimmed());
        }
    }

    if (parser.isSet(geneOpt)) {
        w.displayGene(parser.value(geneOpt));
    }

    return app.exec();
}

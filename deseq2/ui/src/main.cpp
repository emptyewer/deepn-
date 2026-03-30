#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QStandardPaths>
#include <QDebug>
#include "main_window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Set application properties
    app.setApplicationName("StatMaker++");
    app.setApplicationVersion("1.0.0");
    app.setOrganizationName("DEEPN++");
    app.setOrganizationDomain("deepn-plus.com");

    // Set application style
    app.setStyle("Fusion");

    // Create and show main window
    deseq2::MainWindow mainWindow;
    mainWindow.show();

    QCommandLineParser parser;
    parser.setApplicationDescription("StatMaker++ Y2H statistical analysis");
    parser.addHelpOption();
    parser.addVersionOption();

    QCommandLineOption workdirOpt(
        "workdir",
        "Working directory containing experiment data.",
        "path");
    parser.addOption(workdirOpt);

    QCommandLineOption datasetsOpt(
        "datasets",
        "Comma-separated GeneCount or SQLite dataset paths.",
        "paths");
    parser.addOption(datasetsOpt);

    parser.addPositionalArgument("files", "Optional GeneCount or SQLite input files.");
    parser.process(app);

    if (parser.isSet(workdirOpt))
    {
        mainWindow.loadWorkingDirectory(parser.value(workdirOpt));
    }

    if (parser.isSet(datasetsOpt))
    {
        const QStringList files =
            parser.value(datasetsOpt).split(',', Qt::SkipEmptyParts);
        mainWindow.loadInputFiles(files);
    }

    const QStringList positionalFiles = parser.positionalArguments();
    if (!positionalFiles.isEmpty())
    {
        mainWindow.loadInputFiles(positionalFiles);
    }

    // Start event loop
    return app.exec();
}

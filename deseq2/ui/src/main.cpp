#include <QApplication>
#include <QDir>
#include <QStandardPaths>
#include <QDebug>
#include "main_window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Set application properties
    app.setApplicationName("DESeq2 GUI");
    app.setApplicationVersion("1.0.0");
    app.setOrganizationName("DeepN Plus");
    app.setOrganizationDomain("deepn-plus.com");

    // Set application style
    app.setStyle("Fusion");

    // Create and show main window
    deseq2::MainWindow mainWindow;
    mainWindow.show();

    // Start event loop
    return app.exec();
}
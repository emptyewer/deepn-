#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include "signals.h"
#include "datastructs.h"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
Q_OBJECT

public:
    explicit MainWindow(int argc, char *argv[], QWidget *parent = nullptr);

    ~MainWindow() override;

private slots:

    void on_install_btn_clicked();

    void on_run_btn_clicked();

private:
    Ui::MainWindow *ui;
    QList<QString> files;
    Stat *stat = new Stat();
    Signals *sig = Signals::getCommonInstance();

    static void customMessageHandler(QtMsgType type, const QMessageLogContext &context, const QString &msg);

    static void handleMessage(QtMsgType type, const QMessageLogContext &context, const QString &msg);

    void loadFiles();

    void setupSlots();

    void displayStatus(const QString &status);

    bool copyFolderToFrameworks(const QString &sourcePath, const QString &destinationPath);

    static bool pathExists(const QString &path);

    bool checkRFramework();

    void organizeFiles();

    void writeRInputFile(const QString &outputFile);
};

#endif // MAINWINDOW_H
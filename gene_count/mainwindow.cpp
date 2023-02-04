#include "mainwindow.h"

#include <seqan/sequence.h>

#include <QDebug>
#include <QFile>

#include "qdir.h"
#include "ui_mainwindow.h"

QString appendPath(const QString& path1, const QString& path2) {
  return QDir::cleanPath(path1 + QDir::separator() + path2);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  setupPython();
  setupScriptPath();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::setupPython() {
  PythonQt::init();
  // get a smart pointer to the __main__ module of the Python interpreter
  python = PythonQt::self()->getMainModule();
}

void MainWindow::setupScriptPath() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    script_path = appendPath(
        appendPath(application_directory.path(), "Contents"), "Scripts");
  }
}

bool MainWindow::readPythonScript(QString file_name) {
  QString s = appendPath(script_path, file_name);
  QFile file(s);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    return false;
  }
  QTextStream in(&file);
  python.evalScript(in.readAll());
  return true;
}

void MainWindow::on_debug_btn_clicked() {
  readPythonScript("gene_count.py");
  QVariantList args;
  args << 42 << 47;
  QVariant result = python.call("multiply", args);
  qDebug() << args;
  qDebug() << result;
}

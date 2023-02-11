#include "mainwindow.h"

#include <QCoreApplication>
#include <QDebug>
#include <QSysInfo>

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  connect(QCoreApplication::instance(), &QCoreApplication::aboutToQuit,
          &process, &QProcess::kill);
}

MainWindow::~MainWindow() { delete ui; }

QString MainWindow::appendPath(const QString& path1, const QString& path2) {
  return QDir::cleanPath(path1 + QDir::separator() + path2);
}

void MainWindow::on_gene_count_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString gene_count_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    gene_count_path = appendPath(
        application_directory.path(),
        "Contents/Resources/GeneCount++.app/Contents/MacOS/GeneCount++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    gene_count_path =
        appendPath(application_directory.path(), "gene_count/GeneCount++.exe");
  } else {
    gene_count_path =
        appendPath(application_directory.path(), "gene_count/GeneCount++");
  }
  ui->status_text->appendPlainText(gene_count_path);
  QStringList arguments;
  arguments << "/c C:/Users/firstname secondname/desktop/mybatchfile.bat 2";
  process.start(QDir::toNativeSeparators(gene_count_path), arguments);
  process.waitForFinished(-1);
}

void MainWindow::on_junction_make_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString junction_make_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    junction_make_path = appendPath(
        application_directory.path(),
        "/Contents/Resources/JunctionMake++.app/Contents/MacOS/JunctionMake++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    junction_make_path = appendPath(application_directory.path(),
                                    "junction_make/JunctionMake++.exe");
  } else {
    junction_make_path = appendPath(application_directory.path(),
                                    "junction_make/JunctionMake++");
  }
  QStringList arguments;
  arguments << "/c C:/Users/firstname secondname/desktop/mybatchfile.bat 2";
  process.startDetached(QDir::toNativeSeparators(junction_make_path),
                        arguments);
}

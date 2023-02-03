#include "mainwindow.h"

#include <seqan/sequence.h>

#include <QDebug>
#include <QFile>

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  initializePython();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::initializePython() {
}

void MainWindow::on_pushButton_clicked() {
  // add a QObject as variable of name "example" to the namespace of the
  // __main__ module
  qDebug() << "App path : " << qApp->applicationDirPath();
  QFile f(qApp->applicationDirPath() + "/gene_count.py");
  if (!f.open(QFile::ReadOnly | QFile::Text))
    qDebug() << "Error File";
  QTextStream in(&f);
  QString input = in.readAll();
}

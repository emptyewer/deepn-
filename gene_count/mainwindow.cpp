#include "mainwindow.h"

// #include <seqan/sequence.h>

#include <QDebug>
#include <QFile>

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
}

MainWindow::~MainWindow() { delete ui; }

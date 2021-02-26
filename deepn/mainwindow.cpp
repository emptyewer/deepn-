#include<QDebug>

#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);
}

MainWindow::~MainWindow()
{
  delete ui;
}

QString appendPath(const QString& path1, const QString& path2)
{
    return QDir::cleanPath(path1 + QDir::separator() + path2);
}
void MainWindow::on_gene_count_btn_clicked()
{
    QDir application_directory = QDir(QCoreApplication::applicationDirPath());
    application_directory.cdUp();
    application_directory.cdUp();
    QString gene_count_path = appendPath(application_directory.path(), "Contents/Apps/GeneCount++.app/Contents/MacOS/GeneCount++");
    QProcess::execute(QDir::toNativeSeparators(gene_count_path));
}


void MainWindow::on_junction_make_btn_clicked()
{
    QDir application_directory = QDir(QCoreApplication::applicationDirPath());
    application_directory.cdUp();
    application_directory.cdUp();
    QString junction_make_path = appendPath(application_directory.path(), "/Contents/Apps/JunctionMake++.app/Contents/MacOS/JunctionMake++");
    QProcess::execute(QDir::toNativeSeparators(junction_make_path));
}

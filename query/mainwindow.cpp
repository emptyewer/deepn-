#include "mainwindow.h"
#include "dockwidget.h"
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

void MainWindow::on_add_clicked()
{
    DockWidget *wdgt = new DockWidget(this);
    ui->widgets_layout->addWidget(wdgt);
}

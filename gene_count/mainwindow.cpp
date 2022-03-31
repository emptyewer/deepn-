#include "mainwindow.h"
#include "PythonQt.h"
#include "PythonQt_QtAll.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QFile>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  initializePython();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::initializePython() {
  // init PythonQt and Python
  PythonQt::init();
}

void MainWindow::on_pushButton_clicked() {
  // get a smart pointer to the __main__ module of the Python interpreter
  PythonQtObjectPtr context = PythonQt::self()->getMainModule();
  // add a QObject as variable of name "example" to the namespace of the
  // __main__ module
  qDebug() << "App path : " << qApp->applicationDirPath();
  QFile f(qApp->applicationDirPath() + "/gene_count.py");
  if (!f.open(QFile::ReadOnly | QFile::Text))
    qDebug() << "Error File";
  QTextStream in(&f);
  QString input = in.readAll();

  // do something
  context.evalScript(input);
  QVariantList args;
  args << 42 << 47;
  QVariant result = context.call("multiply", args);
  qDebug() << result;
}

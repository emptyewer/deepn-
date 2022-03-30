#include "mainwindow.h"
#include "pybind11/embed.h"
#include "ui_mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  initializePython();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::initializePython() {
  py::initialize_interpreter();
  py::object scope = py::module_::import("__main__").attr("__dict__");
  py::eval_file("gene_count.py", scope);
  py::finalize_interpreter();
}

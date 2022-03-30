#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "pybind11/pybind11.h"
#include <Python.h>
#include <QMainWindow>
#include <pybind11/eval.h>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

namespace py = pybind11;

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

private:
  Ui::MainWindow *ui;
  void initializePython();
};
#endif // MAINWINDOW_H

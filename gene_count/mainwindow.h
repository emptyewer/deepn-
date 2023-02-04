#ifndef MAINWINDOW_H_GENECOUNT
#define MAINWINDOW_H_GENECOUNT
#include <QMainWindow>

#include "PythonQt.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

 private slots:
  void on_debug_btn_clicked();

 private:
  Ui::MainWindow *ui;
  void setupPython();
  PythonQtObjectPtr python;
  QString script_path;
  bool readPythonScript(QString file_name);
  void setupScriptPath();
};
#endif  // MAINWINDOW_H_GENECOUNT

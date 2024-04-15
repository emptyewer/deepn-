#ifndef MAINWINDOW_H_GENECOUNT
#define MAINWINDOW_H_GENECOUNT
#include <QMainWindow>
#include <QThread>
#include <QTimer>
#include <iostream>

#include "PythonQt.h"
#include "datastructs.h"
#include "signals.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow(int argc, char *argv[], QWidget *parent = nullptr);
  ~MainWindow();

 signals:

 private slots:
  void on_debug_btn_clicked();
  void geneCountFinished();
  void updateGeneCountProgress();

 private:
  Ui::MainWindow *ui;
  PythonQtObjectPtr python;
  QString script_path;
  Signals *sig = Signals::getCommonInstance();
  QMap<QString, GCStat> statistics = {};
  QStringList files = {};
  bool readPythonScript(QString file_name);
  void setupPython();
  void setupScriptPath();
  QString appendPath(const QString &path1, const QString &path2);
  void setupSlots();
  void launchGeneCount(QString file);
};

#endif  // MAINWINDOW_H_GENECOUNT

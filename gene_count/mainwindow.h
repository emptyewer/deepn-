#ifndef MAINWINDOW_H_GENECOUNT
#define MAINWINDOW_H_GENECOUNT
#include <QMainWindow>
#include <QThread>
#include <QTimer>

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
  void on_start_btn_clicked();
  void geneCountFinished();
  void updateGeneCountProgress();

 private:
  Ui::MainWindow *ui;
  QString script_path;
  Signals *sig = Signals::getCommonInstance();
  QMap<QString, GCStat> statistics = {};
  QStringList files = {};
  QStringList pendingFiles = {};
  int activeWorkers = 0;
  int maxWorkers = 1;
  QString appendPath(const QString &path1, const QString &path2);
  void setupSlots();
  void launchGeneCount(QString file);
  void launchNext();
};

#endif  // MAINWINDOW_H_GENECOUNT

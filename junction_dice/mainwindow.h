#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "datastructs.h"
#include "signals.h"

QT_BEGIN_NAMESPACE
    namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow(int argc, char *argv[], QWidget *parent = nullptr);
  ~MainWindow();

 private slots:
  void on_dice_btn_clicked();
  void junctionDiceFinished();
  void updateJunctionDiceProgress();

 private:
  Ui::MainWindow *ui;
  QStringList files = {};
  QStringList pendingFiles = {};
  QString junction_sequence = "";
  QString db_path = "";
  Signals *sig = Signals::getCommonInstance();
  QMap<QString, JDStat> statistics = {};
  int activeWorkers = 0;
  int maxWorkers = 1;
  void launchJunctionDice(QString file);
  void launchNext();
  void setupSlots();
  int showDialog(QString message, QString accept, QString cancel);
};
#endif // MAINWINDOW_H

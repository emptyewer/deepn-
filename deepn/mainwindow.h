#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDir>
#include <QFileSystemWatcher>
#include <QListWidgetItem>
#include <QMainWindow>
#include <QProcess>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

 signals:

 private slots:
  void gatherFilesToggleButtons();
  void on_gene_count_btn_clicked();
  void on_select_folder_btn_clicked();
  void on_db_list_wgt_currentItemChanged(QListWidgetItem *current,
                                         QListWidgetItem *previous);
  void on_junction_dice_btn_clicked();
  void on_deseq2_btn_clicked();
  void on_query_blast_btn_clicked();
  void on_read_depth_btn_clicked();
  void on_actionDB_Path_triggered();

 private:
  Ui::MainWindow *ui;
  QList<QMap<QString, QVariant>> data;
  QFileSystemWatcher watcher;
  QProcess process;
  QDir parentDir = QDir();
  QMap<QString, QStringList> files = {};
  QString appendPath(const QString &path1, const QString &path2);
  void loadDatabase();
  void initializeWorkDir();
  void monitorSubDirs(QDir parentDir);
  void createSubDirs();
  void gatherFiles(QDir dir, QString key);
  void readJsonFile(QString filename);
};
#endif // MAINWINDOW_H

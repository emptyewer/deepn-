#include "mainwindow.h"

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QtSql>

#include "gcworker.h"
#include "ui_mainwindow.h"

QString MainWindow::appendPath(const QString& path1, const QString& path2) {
  return QDir::cleanPath(path1 + QDir::separator() + path2);
}

MainWindow::MainWindow(int argc, char* argv[], QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  this->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint |
                       Qt::WindowStaysOnTopHint);
  setupSlots();
  qRegisterMetaType<GCStat>("GCStat");
  for (int i = 1; i < argc; i++) {
      files << argv[i];
  }
  //  files << "/Users/vkrishnamani/Downloads/DEEPN_Example_Data/mapped_files/"
  //           "Piper_48_1_lane2_20200221000_S97_L002_R1_001.blat.txt";
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::launchGeneCount(QString file) {
  GCStat stat = GCStat();
  QFileInfo fi(file);
  QThread* thread = new QThread;
  stat.input = file;
  statistics[fi.baseName()] = stat;
  GCWorker* worker = new GCWorker(&statistics[fi.baseName()]);
  worker->moveToThread(thread);
  connect(worker, &GCWorker::finished, thread, &QThread::quit);
  connect(thread, &QThread::finished, worker, &QObject::deleteLater);
  connect(thread, &QThread::finished, thread, &QObject::deleteLater);
  thread->start();
  QMetaObject::invokeMethod(worker, "run");
}

void MainWindow::updateGeneCountProgress() {
  ui->gc_output->clear();
  foreach (GCStat stat, statistics.values()) {
    if (stat.running == false) {
      ui->gc_output->appendPlainText(QString(">>> %1 <<<").arg(stat.input));
      ui->gc_output->appendPlainText(QString("Finished Counting %2 in %1 secs ")
                                         .arg(stat.elapsedTime)
                                         .arg(stat.readCount));
    } else {
      ui->gc_output->appendPlainText(QString("*** %1 ***").arg(stat.input));
      ui->gc_output->appendPlainText(
          QString("Elapsed Time : %1 secs for %2 reads")
              .arg(stat.elapsedTime)
              .arg(stat.readCount));
      ui->gc_output->appendPlainText(
          QString("Current Read : %1").arg(stat.readName));
    }
    //    ui->gc_output->appendPlainText(stat.counter.getStats());
  }
}

void MainWindow::geneCountFinished() {
  activeWorkers--;
  if (activeWorkers <= 0) {
    ui->start_btn->setEnabled(true);
    ui->exit_btn->setEnabled(true);
  }
}

void MainWindow::setupSlots() {
  connect(sig, &Signals::gc_update_progress_sig, this,
          &MainWindow::updateGeneCountProgress);
  connect(sig, &Signals::gc_finished_sig, this, &MainWindow::geneCountFinished);
}

void MainWindow::on_start_btn_clicked() {
  ui->gc_output->appendPlainText("Starting Gene Count...");
  ui->start_btn->setEnabled(false);
  ui->exit_btn->setEnabled(false);
  activeWorkers = files.length();
  for (int i = 0; i < files.length(); i++) {
    launchGeneCount(files.at(i));
  }
}

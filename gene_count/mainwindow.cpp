#include "mainwindow.h"

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QScrollBar>
#include <QStatusBar>
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
  QScrollBar* sb = ui->gc_output->verticalScrollBar();
  int scrollPos = sb->value();
  bool wasAtBottom = (scrollPos >= sb->maximum() - 4);
  ui->gc_output->clear();
  bool first = true;
  foreach (GCStat stat, statistics.values()) {
    if (!first) {
      ui->gc_output->appendPlainText(QString(""));
      ui->gc_output->appendPlainText(QString("────────────────────────────────"));
      ui->gc_output->appendPlainText(QString(""));
    }
    first = false;
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
  }
  if (wasAtBottom) {
    sb->setValue(sb->maximum());
  } else {
    sb->setValue(scrollPos);
  }
}

void MainWindow::geneCountFinished() {
  activeWorkers--;
  launchNext();
  int finished = files.length() - activeWorkers - pendingFiles.length();
  if (activeWorkers <= 0 && pendingFiles.isEmpty()) {
    ui->start_btn->setEnabled(true);
    ui->exit_btn->setEnabled(true);
    statusBar()->showMessage(QString("All %1 file(s) completed.").arg(files.length()));
  } else {
    statusBar()->showMessage(QString("%1 of %2 completed, %3 running, %4 queued")
                                 .arg(finished).arg(files.length())
                                 .arg(activeWorkers).arg(pendingFiles.length()));
  }
}

void MainWindow::launchNext() {
  while (activeWorkers < maxWorkers && !pendingFiles.isEmpty()) {
    launchGeneCount(pendingFiles.takeFirst());
    activeWorkers++;
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
  maxWorkers = qMax(1, QThread::idealThreadCount() - 2);
  pendingFiles = files;
  activeWorkers = 0;
  while (activeWorkers < maxWorkers && !pendingFiles.isEmpty()) {
    launchGeneCount(pendingFiles.takeFirst());
    activeWorkers++;
  }
  statusBar()->showMessage(QString("Processing %1 file(s), %2 concurrent...")
                               .arg(files.length()).arg(activeWorkers));
}

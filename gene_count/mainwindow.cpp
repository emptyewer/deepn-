#include "mainwindow.h"

#include <QDebug>
#include <QFile>
#include <QProcess>

#include "gcworker.h"
#include "qdir.h"
#include "ui_mainwindow.h"

QString MainWindow::appendPath(const QString& path1, const QString& path2) {
  return QDir::cleanPath(path1 + QDir::separator() + path2);
}

MainWindow::MainWindow(int argc, char* argv[], QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  qDebug() << "#### " << argc << argv[1];
  ui->setupUi(this);
  setupPython();
  setupScriptPath();
  setupSlots();
  qRegisterMetaType<GCStat>("GCStat");
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::launchGeneCount() {
  QThread* workerThread = new QThread();
  GCWorker* worker = new GCWorker(
      "/Users/vkrishnamani/Downloads/DEEPN_Example_Data/Mm10_mappedreads.sam");

  worker->moveToThread(workerThread);
  connect(sig, &Signals::gc_update_progress_sig, this,
          &MainWindow::updateGeneCountProgress);
  connect(sig, &Signals::gc_finished_sig, this, &MainWindow::geneCountFinished);
  connect(workerThread, &QThread::started, worker, &GCWorker::doWork);
  workerThread->start();
}

void MainWindow::updateGeneCountProgress(GCStat stat) {
  ui->gc_output->clear();
  if (stat.running == false) {
    ui->gc_output->appendPlainText(QString(">>> %1 <<<").arg(stat.filename));
    ui->gc_output->appendPlainText(
        QString("Finished Counting a %2 reads in %1 secs")
            .arg(stat.elapsedTime)
            .arg(stat.readCount));
  } else {
    ui->gc_output->appendPlainText(QString("*** %1 ***").arg(stat.filename));
    ui->gc_output->appendPlainText(
        QString("Elapsed Time : %1 secs").arg(stat.elapsedTime));
    ui->gc_output->appendPlainText(QString("Current Read: %1 (%2)")
                                       .arg(stat.readCount)
                                       .arg(stat.currentRefSeq));
  }
}

void MainWindow::geneCountFinished() {}

void MainWindow::setupSlots() {}

void MainWindow::setupPython() {
  PythonQt::init();
  // get a smart pointer to the __main__ module of the Python interpreter
  python = PythonQt::self()->getMainModule();
}

void MainWindow::setupScriptPath() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    script_path = appendPath(
        appendPath(application_directory.path(), "Contents"), "Scripts");
    PythonQt::self()->addSysPath(script_path);
  }
}

bool MainWindow::readPythonScript(QString file_name) {
  QString s = appendPath(script_path, file_name);
  QFile file(s);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    return false;
  }
  QTextStream in(&file);
  python.evalScript(in.readAll());
  return true;
}

void MainWindow::on_debug_btn_clicked() {
  readPythonScript("gene_count.py");
  QVariantList args;
  args << 8 << 4;
  QVariant result = python.call("multiply", args);
  QVariant pwd = python.call("get_current_dir");
  qDebug() << result << pwd;
  ui->gc_output->appendPlainText(result.toString());
  ui->gc_output->appendPlainText(pwd.toString());
  ui->gc_output->appendPlainText(python.call("test_utils").toString());
  ui->gc_output->appendPlainText("Emitting...");
  launchGeneCount();
}


#include "mainwindow.h"

#include <QDebug>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QFile>
#include <QLabel>
#include <QScrollBar>
#include <QStatusBar>
#include <QVBoxLayout>

#include "jdworker.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(int argc, char* argv[], QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  this->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint |
                       Qt::WindowStaysOnTopHint);
  setupSlots();
  qRegisterMetaType<JDStat>("JDStat");

  for (int i = 1; i < argc - 2; i++) {
    files << argv[i];
  }
  junction_sequence = argv[argc - 2];
  db_path = argv[argc - 1];

  //  files << "/Users/vkrishnamani/Downloads/DEEPN_Example_Data/fastq/"
  //           "Piper_48_1_lane2_20200221000_S97_L002_R1_001.fastq.gz";
  //  junction_sequence = "CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC";
  //  db_path = "hg38GeneList2023.unique.fasta";

  QFileInfo db(db_path);
  if (db.path() == ".") {
    QDir application_directory = QDir(QCoreApplication::applicationDirPath());
    application_directory.cdUp();
    db_path = QDir::cleanPath(application_directory.path() + QDir::separator() +
                              "Data" + QDir::separator() + db_path);
  }
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::setupSlots() {
  connect(sig, &Signals::jd_update_progress_sig, this,
          &MainWindow::updateJunctionDiceProgress);
  connect(sig, &Signals::jd_finished_sig, this,
          &MainWindow::junctionDiceFinished);
}

int MainWindow::showDialog(QString message, QString accept, QString cancel) {
  QDialog dialog;
  QVBoxLayout layout(&dialog);
  QLabel label(message);
  layout.addWidget(&label);

  QDialogButtonBox buttonBox(&dialog);
  buttonBox.addButton(accept, QDialogButtonBox::AcceptRole);
  buttonBox.addButton(cancel, QDialogButtonBox::RejectRole);
  layout.addWidget(&buttonBox);

  QObject::connect(&buttonBox, &QDialogButtonBox::accepted, &dialog,
                   &QDialog::accept);
  QObject::connect(&buttonBox, &QDialogButtonBox::rejected, &dialog,
                   &QDialog::reject);

  int result = dialog.exec();
  return result;
}

void MainWindow::launchJunctionDice(QString file) {
  QFileInfo fi(file);
  JDStat stat = JDStat();
  QThread* thread = new QThread;
  thread->start();
  stat.dstat.input = file;
  if (ui->blat->isChecked()) {
    stat.mstat.algo = blat;
  } else if (ui->blastn->isChecked()) {
    stat.mstat.algo = blastn;
  } else if (ui->megablast->isChecked()) {
    stat.mstat.algo = megablast;
  }
  stat.mstat.db = db_path;
  stat.dstat.matchLength = ui->jseq_match_len->value();
  stat.dstat.junction = junction_sequence;
  statistics[fi.baseName()] = stat;
  JDWorker* worker = new JDWorker(&statistics[fi.baseName()], files.count());
  if (QFileInfo(statistics[fi.baseName()].dstat.output).exists()) {
    QString msg =
        QString(
            "Jiced Output file %1 already exists...\nDo you want to "
            "re-dice and overwrite or use the current file for mapping?")
            .arg(QFileInfo(statistics[fi.baseName()].dstat.output).fileName());
    int result = showDialog(msg, "Overwrite", "Use Current");
    if (result == QDialog::Accepted) {
      statistics[fi.baseName()].readdice = false;
      statistics[fi.baseName()].dicing = true;
    } else {
      statistics[fi.baseName()].readdice = true;
      statistics[fi.baseName()].dicing = false;
    }
  }
  worker->moveToThread(thread);
  connect(worker, &JDWorker::finished, thread, &QThread::quit);
  connect(thread, &QThread::finished, worker, &QObject::deleteLater);
  connect(thread, &QThread::finished, thread, &QObject::deleteLater);
  QMetaObject::invokeMethod(worker, "run");
}

void MainWindow::on_dice_btn_clicked() {
  //  readPythonScript("gene_count.py");
  //  QVariantList args;
  //  args << 8 << 4;
  //  QVariant result = python.call("multiply", args);
  //  QVariant pwd = python.call("get_current_dir");
  //  qDebug() << result << pwd;
  //  ui->gc_output->appendPlainText(result.toString());
  //  ui->gc_output->appendPlainText(pwd.toString());
  //  ui->gc_output->appendPlainText(python.call("test_utils").toString());
  //  ui->gc_output->appendPlainText("Emitting...");
  ui->jd_output->appendPlainText("Starting Junction Dice...");
  maxWorkers = qMax(1, QThread::idealThreadCount() - 2);
  pendingFiles = files;
  activeWorkers = 0;
  ui->dice_btn->setEnabled(false);
  ui->exit_btn->setEnabled(false);
  ui->jseq_match_len->setEnabled(false);
  ui->blastChoices->setEnabled(false);
  // Launch up to maxWorkers initially
  while (activeWorkers < maxWorkers && !pendingFiles.isEmpty()) {
    launchJunctionDice(pendingFiles.takeFirst());
    activeWorkers++;
  }
  statusBar()->showMessage(QString("Processing %1 file(s), %2 concurrent...")
                               .arg(files.length()).arg(activeWorkers));
}

void MainWindow::junctionDiceFinished() {
  activeWorkers--;
  launchNext();
  int finished = files.length() - activeWorkers - pendingFiles.length();
  if (activeWorkers <= 0 && pendingFiles.isEmpty()) {
    ui->exit_btn->setEnabled(true);
    ui->dice_btn->setEnabled(true);
    ui->jseq_match_len->setEnabled(true);
    ui->blastChoices->setEnabled(true);
    statusBar()->showMessage(QString("All %1 file(s) completed.").arg(files.length()));
  } else {
    statusBar()->showMessage(QString("%1 of %2 completed, %3 running, %4 queued")
                                 .arg(finished).arg(files.length())
                                 .arg(activeWorkers).arg(pendingFiles.length()));
  }
}

void MainWindow::launchNext() {
  while (activeWorkers < maxWorkers && !pendingFiles.isEmpty()) {
    launchJunctionDice(pendingFiles.takeFirst());
    activeWorkers++;
  }
}

void MainWindow::updateJunctionDiceProgress() {
  QScrollBar* sb = ui->jd_output->verticalScrollBar();
  int scrollPos = sb->value();
  bool wasAtBottom = (scrollPos >= sb->maximum() - 4);
  ui->jd_output->clear();
  foreach (JDStat stat, statistics.values()) {
    ui->jd_output->appendPlainText(QString("--------------------------------"));
    if (stat.dicing == false) {
      ui->jd_output->appendPlainText(
          QString("Junction Diced File %1\n")
              .arg(QFileInfo(stat.dstat.output).fileName()));
      ui->jd_output->appendPlainText(QString("Read %1 MB in %2 secs")
                                         .arg(stat.dstat.fileSize / 1000000)
                                         .arg(stat.dstat.elapsedTime));
      ui->jd_output->appendPlainText(
          QString("Diced %1 forward and %2 reverse junction sequences in %3 "
                  "total reads")
              .arg(stat.dstat.forwardMatches)
              .arg(stat.dstat.reverseMatches)
              .arg(stat.dstat.totalReads));
      ui->jd_output->appendPlainText(
          QString("Wrote %1 reads to diced output file")
              .arg(stat.dstat.writtenReads));
      ui->jd_output->appendPlainText(
          QString("Skipped %1 reads\n").arg(stat.dstat.removedReads));
    } else {
      ui->jd_output->appendPlainText(QString("<<< %1 (%2 MB)\n")
                                         .arg(stat.dstat.input)
                                         .arg(stat.dstat.fileSize / 1000000));
      ui->jd_output->appendPlainText(
          QString("Junction Sequence: %1").arg(stat.dstat.junction));
      ui->jd_output->appendPlainText(QString("%1bp 3'-end Match Sequence: %2")
                                         .arg(stat.dstat.matchLength)
                                         .arg(stat.dstat.match));
      ui->jd_output->appendPlainText(
          QString("\nElapsed Time : %1 secs").arg(stat.dstat.elapsedTime));
      ui->jd_output->appendPlainText(
          QString("Diced %1 forward and %2 reverse junction sequences in %3 "
                  "total reads\n")
              .arg(stat.dstat.forwardMatches)
              .arg(stat.dstat.reverseMatches)
              .arg(stat.dstat.totalReads));
    }

    if (stat.blasting == true) {
      ui->jd_output->appendPlainText(
          QString("\nReference DB : %1")
              .arg(QFileInfo(stat.mstat.db).fileName()));
      if (stat.mstat.algo == blat) {
        ui->jd_output->appendPlainText(
            QString("Algorithm : BLAT (Blast Like Alignment Tool)"));
      } else if (stat.mstat.algo == blastn) {
        ui->jd_output->appendPlainText(
            QString("Algorithm : BLASTN (Nucleotide Blast)"));
      } else if (stat.mstat.algo == megablast) {
        ui->jd_output->appendPlainText(QString("Algorithm : MEGABLAST "));
      }
      ui->jd_output->appendPlainText(
          QString("<<< %1").arg(QFileInfo(stat.dstat.output).fileName()));
      ui->jd_output->appendPlainText(
          QString(">>> %1").arg(QFileInfo(stat.mstat.output).fileName()));
      ui->jd_output->appendPlainText(
          QString("\nCurrent Read : %1").arg(stat.mstat.currentRead));
      ui->jd_output->appendPlainText(QString("%1% mapped in %2 secs...")
                                         .arg(stat.mstat.percentComplete)
                                         .arg(stat.mstat.elapsedTime));
    }

    if (stat.blasting == false && stat.dicing == false &&
        stat.mstat.percentComplete > 50) {
      ui->jd_output->appendPlainText(
          QString("\n>>> %1").arg(QFileInfo(stat.mstat.output).fileName()));
      ui->jd_output->appendPlainText(
          QString("Reference DB : %1")
              .arg(QFileInfo(stat.mstat.db).fileName()));
      ui->jd_output->appendPlainText(QString("%1% mapped in %2 secs...")
                                         .arg(stat.mstat.percentComplete)
                                         .arg(stat.mstat.elapsedTime));
      ui->jd_output->appendPlainText(
          QString("\nDICING and MAPPING Completed!"));
      ui->jd_output->appendPlainText(
          QString("--------------------------------\n"));
    }
  }
  if (wasAtBottom) {
    sb->setValue(sb->maximum());
  } else {
    sb->setValue(scrollPos);
  }
}

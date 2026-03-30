#include "mainwindow.h"

#include <QCoreApplication>
#include <QDebug>
#include <QFileDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QStandardPaths>
#include <QSysInfo>
#include <QUuid>

#include "darkpaint.h"
#include "ui_mainwindow.h"

QMap<QString, QStringList> SUBDIR_FILES = {
    {"fastq", {"*.fastq.gz", "*.fastq"}},
    {"junction_diced_fasta", {"*.jdice.fasta"}},
    {"gene_count_summary", {"*.csv"}},
    {"query_files", {"*.bq"}},
    {"mapped_files", {"*.megablast.txt", "*.blat.txt", "*.blastn.txt"}},
    {"analyzed_files", {"*.sqlite"}}};

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  connect(QCoreApplication::instance(), &QCoreApplication::aboutToQuit,
          &process, &QProcess::kill);
  QObject::connect(&watcher, &QFileSystemWatcher::directoryChanged, this,
                   &MainWindow::gatherFilesToggleButtons);
  loadDatabase();
}

MainWindow::~MainWindow() { delete ui; }

QString MainWindow::appendPath(const QString& path1, const QString& path2) {
  return QDir::cleanPath(path1 + QDir::separator() + path2);
}

void MainWindow::gatherFiles(QDir dir, QString key) {
  files[key].clear();
  dir.setFilter(QDir::Files | QDir::NoDotAndDotDot);
  QFileInfoList fileInfoList = dir.entryInfoList(SUBDIR_FILES[key]);
  foreach (QFileInfo fileInfo, fileInfoList) {
    files[key] << fileInfo.filePath();
  }
}

void MainWindow::gatherFilesToggleButtons() {
  foreach (QString subDirName, SUBDIR_FILES.keys()) {
    QDir sub(parentDir.absoluteFilePath(subDirName));
    gatherFiles(sub, subDirName);
  }
  if (files["fastq"].length() > 0) {
    ui->junction_dice_btn->setEnabled(true);
  } else {
    ui->junction_dice_btn->setEnabled(false);
  }

  if (files["mapped_files"].length() > 0) {
    ui->gene_count_btn->setEnabled(true);
  } else {
    ui->gene_count_btn->setEnabled(false);
  }

  if (files["analyzed_files"].length() > 0) {
    ui->query_blast_btn->setEnabled(true);
  } else {
    ui->query_blast_btn->setEnabled(false);
  }

  if (files["analyzed_files"].length() > 0) {
    ui->read_depth_btn->setEnabled(true);
  } else {
    ui->read_depth_btn->setEnabled(false);
  }

  if (files["gene_count_summary"].length() >= 2 ||
      files["analyzed_files"].length() >= 2) {
    ui->deseq2_btn->setEnabled(true);
  } else {
    ui->deseq2_btn->setEnabled(false);
  }
}

void MainWindow::monitorSubDirs(QDir parentDir) {
  foreach (QString subDirName, SUBDIR_FILES.keys()) {
    QDir sub(parentDir.absoluteFilePath(subDirName));
    watcher.addPath(sub.absolutePath());
  }
}

void MainWindow::createSubDirs() {
  foreach (QString subdir, SUBDIR_FILES.keys()) {
    QDir sub(parentDir.absoluteFilePath(subdir));
    if (!sub.exists()) {
      ui->status_text->appendPlainText(
          QString("Creating subfolder : %1").arg(sub.path()));
      sub.mkpath(".");
    }
  }
}

void MainWindow::initializeWorkDir() {
  if (parentDir.absolutePath() != "") {
    ui->status_text->appendPlainText(
        QString("Setting Work Directory : %1").arg(parentDir.absolutePath()));
    if (!parentDir.exists()) {
      parentDir.mkdir(".");
    }
    createSubDirs();
    monitorSubDirs(parentDir);
    gatherFilesToggleButtons();
  }
}

void MainWindow::readJsonFile(QString filename) {
  // Open the file
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qWarning() << "Could not open file for reading:" << filename;
    return;
  }

  // Read the JSON data
  QByteArray jsonData = file.readAll();
  QJsonDocument doc = QJsonDocument::fromJson(jsonData);

  // Check if the JSON data is valid
  if (doc.isNull()) {
    qWarning() << "Invalid JSON data in file:" << filename;
    return;
  }

  // Convert the JSON data to a list of dictionaries
  QJsonArray jsonArray = doc.array();
  for (auto it = jsonArray.begin(); it != jsonArray.end(); ++it) {
    QJsonObject jsonObj = it->toObject();
    QMap<QString, QVariant> map;
    for (auto iter = jsonObj.begin(); iter != jsonObj.end(); ++iter) {
      map.insert(iter.key(), iter.value().toVariant());
    }
    data.append(map);
  }
}

void MainWindow::loadDatabase() {
  QString data_path;
  bool sel = false;
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    data_path = QDir::cleanPath(application_directory.path() +
                                QDir::separator() + "Contents/Data/deepn.json");
  }
  readJsonFile(data_path);
  for (const auto& map : qAsConst(data)) {
    QListWidgetItem* item = new QListWidgetItem();
    item->setData(Qt::UserRole, map);
    item->setText(map["display_name"].toString());
    ui->db_list_wgt->addItem(item);
    // Iterate over the key-value pairs in the map
    if (sel) {
      ui->db_list_wgt->setCurrentItem(item);
      sel = false;
    }
  }
  QApplication::processEvents();
}

void MainWindow::on_gene_count_btn_clicked() {
  //  Overlay* overlay_ = new Overlay(this);
  //  overlay_->resize(size());
  //  overlay_->setVisible(true);
  //  QApplication::processEvents();

  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString gene_count_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    gene_count_path = appendPath(
        application_directory.path(),
        "Contents/Resources/GeneCount++.app/Contents/MacOS/GeneCount++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    gene_count_path =
        appendPath(application_directory.path(), "gene_count/GeneCount++.exe");
  } else {
    gene_count_path =
        appendPath(application_directory.path(), "gene_count/GeneCount++");
  }
  ui->status_text->appendPlainText(gene_count_path);
  QStringList arguments;
  arguments << files["mapped_files"];
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(gene_count_path), arguments);
  process.waitForFinished(-1);
}

void MainWindow::on_db_list_wgt_currentItemChanged(QListWidgetItem* current,
                                                   QListWidgetItem* previous) {
  QMap<QString, QVariant> data = current->data(Qt::UserRole).toMap();
  ui->junction_sequence_txt->setText(data["junction_sequence"].toString());
  ui->database_path->setText(data["map_db"].toString());
}

void MainWindow::on_select_folder_btn_clicked() {
  QString workDir = QFileDialog::getExistingDirectory(
      this, tr("Select DEEPN++ Work Directory"), QDir::homePath(),
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  ui->folder_path_lbl->setText(workDir);
  parentDir = QDir(workDir);
  initializeWorkDir();
}

void MainWindow::on_junction_dice_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString junction_dice_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    junction_dice_path = appendPath(application_directory.path(),
                                    "Contents/Resources/JunctionDice++.app/"
                                    "Contents/MacOS/JunctionDice++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    junction_dice_path = appendPath(application_directory.path(),
                                    "junction_dice/JunctionDice++.exe");
  } else {
    junction_dice_path = appendPath(application_directory.path(),
                                    "junction_dice/JunctionDice++");
  }
  ui->status_text->appendPlainText(junction_dice_path);
  QStringList arguments;
  arguments << files["fastq"] << ui->junction_sequence_txt->text()
            << ui->database_path->text();
  qDebug() << arguments;
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(junction_dice_path), arguments);
  process.waitForFinished(-1);
}

// void MainWindow::on_junction_make_btn_clicked() {
//   QDir application_directory =
//   QDir(QCoreApplication::applicationDirPath()); QString junction_make_path;
//   if (QSysInfo::productType() == "osx" || QSysInfo::productType() ==
//   "macos")
//   {
//     application_directory.cdUp();
//     application_directory.cdUp();
//     junction_make_path = appendPath(
//         application_directory.path(),
//         "/Contents/Resources/JunctionMake++.app/Contents/MacOS/JunctionMake++");
//   } else if (QSysInfo::productType() == "windows" ||
//              QSysInfo::productType() == "winrt") {
//     junction_make_path = appendPath(application_directory.path(),
//                                     "junction_make/JunctionMake++.exe");
//   } else {
//     junction_make_path = appendPath(application_directory.path(),
//                                     "junction_make/JunctionMake++");
//   }
//   QStringList arguments;
//   arguments << "/c C:/Users/firstname secondname/desktop/mybatchfile.bat
//   2"; process.start(QDir::toNativeSeparators(junction_make_path),
//   arguments); process.waitForFinished(-1);
// }

void MainWindow::on_deseq2_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString deseq2_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    deseq2_path = appendPath(
        application_directory.path(),
        "Contents/Resources/StatMaker++.app/Contents/MacOS/StatMaker++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    deseq2_path =
        appendPath(application_directory.path(), "deseq2/StatMaker++.exe");
  } else {
    deseq2_path =
        appendPath(application_directory.path(), "deseq2/StatMaker++");
  }
  ui->status_text->appendPlainText(deseq2_path);
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath();
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(deseq2_path), arguments);
  process.waitForFinished(-1);
}

void MainWindow::on_query_blast_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString query_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    query_path = appendPath(
        application_directory.path(),
        "Contents/Resources/MultiQuery++.app/Contents/MacOS/MultiQuery++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    query_path =
        appendPath(application_directory.path(), "query/MultiQuery++.exe");
  } else {
    query_path =
        appendPath(application_directory.path(), "query/MultiQuery++");
  }
  ui->status_text->appendPlainText(query_path);
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath()
            << "--datasets" << files["analyzed_files"].join(",");
  // Pass gene reference database path if available
  QListWidgetItem* currentDb = ui->db_list_wgt->currentItem();
  if (currentDb) {
    QMap<QString, QVariant> dbData = currentDb->data(Qt::UserRole).toMap();
    arguments << "--generef" << dbData["map_db"].toString();
  }
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(query_path), arguments);
  process.waitForFinished(-1);
}

void MainWindow::on_read_depth_btn_clicked() {
  QDir application_directory = QDir(QCoreApplication::applicationDirPath());
  QString depth_path;
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    application_directory.cdUp();
    application_directory.cdUp();
    depth_path = appendPath(
        application_directory.path(),
        "Contents/Resources/ReadDepth++.app/Contents/MacOS/ReadDepth++");
  } else if (QSysInfo::productType() == "windows" ||
             QSysInfo::productType() == "winrt") {
    depth_path =
        appendPath(application_directory.path(), "read_depth/ReadDepth++.exe");
  } else {
    depth_path =
        appendPath(application_directory.path(), "read_depth/ReadDepth++");
  }
  ui->status_text->appendPlainText(depth_path);
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath()
            << "--datasets" << files["analyzed_files"].join(",");
  QListWidgetItem* currentDb = ui->db_list_wgt->currentItem();
  if (currentDb) {
    QMap<QString, QVariant> dbData = currentDb->data(Qt::UserRole).toMap();
    arguments << "--generef" << dbData["map_db"].toString();
  }
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(depth_path), arguments);
  process.waitForFinished(-1);
}

void MainWindow::on_actionDB_Path_triggered() {}

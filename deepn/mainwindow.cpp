#include "mainwindow.h"

#include <QCoreApplication>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QStandardPaths>
#include <QSysInfo>
#include <QUuid>

#include "darkpaint.h"
#include "ui_mainwindow.h"

namespace {
// Check if a SQLite database contains specific tables
bool sqliteHasTables(const QString& dbPath, const QStringList& requiredTables) {
    QString connName = "probe_" + QUuid::createUuid().toString(QUuid::Id128);
    bool ok = false;
    {
        QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connName);
        db.setDatabaseName(dbPath);
        db.setConnectOptions("QSQLITE_OPEN_READONLY");
        if (db.open()) {
            int found = 0;
            for (const QString& table : requiredTables) {
                QSqlQuery q(db);
                if (q.exec(QString("SELECT 1 FROM sqlite_master WHERE type='table' AND name='%1'").arg(table))) {
                    if (q.next()) found++;
                }
            }
            ok = (found == requiredTables.size());
            db.close();
        }
    }
    QSqlDatabase::removeDatabase(connName);
    return ok;
}

// Count how many SQLite files in a list contain the required tables
int countValidSqliteFiles(const QStringList& files, const QStringList& requiredTables) {
    int count = 0;
    for (const QString& f : files) {
        if (sqliteHasTables(f, requiredTables)) count++;
    }
    return count;
}
}  // namespace

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
  connect(&process, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
          this, &MainWindow::onSubprocessFinished);
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

  // All analysis buttons require a reference genome to be selected
  bool hasGenomeRef = (ui->db_list_wgt->currentItem() != nullptr);

  ui->junction_dice_btn->setEnabled(hasGenomeRef && files["fastq"].length() > 0);
  ui->gene_count_btn->setEnabled(hasGenomeRef && files["mapped_files"].length() > 0);

  // MultiQuery++ and ReadDepth++ need SQLite files with a 'maps' table
  int withMaps = countValidSqliteFiles(files["analyzed_files"], {"maps"});
  ui->query_blast_btn->setEnabled(hasGenomeRef && withMaps > 0);
  ui->read_depth_btn->setEnabled(hasGenomeRef && withMaps > 0);

  // StatMaker++ needs >=2 SQLite files with 'maps' table
  ui->deseq2_btn->setEnabled(hasGenomeRef && withMaps >= 2);
}

void MainWindow::monitorSubDirs(QDir parentDir) {
  foreach (QString subDirName, SUBDIR_FILES.keys()) {
    QDir sub(parentDir.absoluteFilePath(subDirName));
    watcher.addPath(sub.absolutePath());
  }
}

void MainWindow::createSubDirs() {
  // Only create input directories; output directories (gene_count_summary,
  // query_files, junction_diced_fasta, mapped_files, analyzed_files)
  // are created by the pipeline modules when they produce output.
  static const QStringList inputDirs = {"fastq"};
  for (const QString& subdir : inputDirs) {
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
  }
  // No genome selected by default — user must choose before analysis
  ui->db_list_wgt->setCurrentRow(-1);
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
  QStringList arguments;
  arguments << files["mapped_files"];
  launchSubprocess(gene_count_path, arguments, "GeneCount++");
}

void MainWindow::on_db_list_wgt_currentItemChanged(QListWidgetItem* current,
                                                   QListWidgetItem* previous) {
  Q_UNUSED(previous)
  if (!current) {
    ui->junction_sequence_txt->clear();
    ui->database_path->clear();
    gatherFilesToggleButtons();
    return;
  }
  QMap<QString, QVariant> data = current->data(Qt::UserRole).toMap();
  ui->junction_sequence_txt->setText(data["junction_sequence"].toString());
  ui->database_path->setText(data["map_db"].toString());
  gatherFilesToggleButtons();
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
  QStringList arguments;
  arguments << files["fastq"] << ui->junction_sequence_txt->text()
            << ui->database_path->text();
  qDebug() << arguments;
  launchSubprocess(junction_dice_path, arguments, "JunctionDice++");
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
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath();
  launchSubprocess(deseq2_path, arguments, "StatMaker++");
}

QString MainWindow::resolveGeneRefPath() const {
  // Resolve the gene reference database filename to a full path.
  // Prefer the pre-built .sqlite annotation DB; fall back to FASTA.
  QListWidgetItem* currentDb = ui->db_list_wgt->currentItem();
  if (!currentDb) return {};

  QString dbName = currentDb->data(Qt::UserRole).toMap()["map_db"].toString();
  if (dbName.isEmpty()) return {};

  // Build list of directories to search for gene reference data.
  // On macOS, data is bundled inside the nested JunctionDice++ app.
  QStringList searchDirs;
  QDir appDir(QCoreApplication::applicationDirPath());
  if (QSysInfo::productType() == "osx" || QSysInfo::productType() == "macos") {
    appDir.cdUp();  // Contents/MacOS -> Contents
    // Nested child app data (primary location for bundled builds)
    searchDirs << appDir.filePath("Resources/JunctionDice++.app/Contents/Data");
    // Orchestrator's own Data dir
    searchDirs << appDir.filePath("Data");
  } else {
    searchDirs << appDir.filePath("Data");
  }
  // Also try working directory
  searchDirs << parentDir.absolutePath();

  for (const QString& dir : searchDirs) {
    // Try pre-built .sqlite first
    QString sqlitePath = QDir(dir).filePath(dbName + ".sqlite");
    if (QFileInfo::exists(sqlitePath)) return sqlitePath;
    // Try FASTA
    QString fastaPath = QDir(dir).filePath(dbName);
    if (QFileInfo::exists(fastaPath)) return fastaPath;
  }

  return dbName;  // Return as-is, let the child app handle the error
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
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath()
            << "--datasets" << files["analyzed_files"].join(",");
  QString geneRef = resolveGeneRefPath();
  if (!geneRef.isEmpty()) {
    arguments << "--generef" << geneRef;
  }
  launchSubprocess(query_path, arguments, "MultiQuery++");
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
  QStringList arguments;
  arguments << "--workdir" << parentDir.absolutePath()
            << "--datasets" << files["analyzed_files"].join(",");
  QString geneRef = resolveGeneRefPath();
  if (!geneRef.isEmpty()) {
    arguments << "--generef" << geneRef;
  }
  launchSubprocess(depth_path, arguments, "ReadDepth++");
}

void MainWindow::on_actionDB_Path_triggered() {}

void MainWindow::launchSubprocess(const QString& execPath, const QStringList& arguments,
                                   const QString& displayName) {
  if (process.state() != QProcess::NotRunning) {
    QMessageBox::warning(this, "Process Running",
                         "Another subprocess is already running. Please wait for it to finish.");
    return;
  }
  ui->status_text->appendPlainText(execPath + " " + arguments.join(" "));
  showOverlay(QString("Running %1...").arg(displayName));
  process.setProcessChannelMode(QProcess::ForwardedChannels);
  process.start(QDir::toNativeSeparators(execPath), arguments);
}

void MainWindow::onSubprocessFinished(int exitCode, QProcess::ExitStatus exitStatus) {
  Q_UNUSED(exitCode)
  Q_UNUSED(exitStatus)
  hideOverlay();
  gatherFilesToggleButtons();
}

void MainWindow::showOverlay(const QString& message) {
  if (!m_overlay) {
    m_overlay = new QWidget(this);
    m_overlay->setStyleSheet("background-color: rgba(0, 0, 0, 230);");
    m_overlay->setAttribute(Qt::WA_StyledBackground, true);
    m_overlay->setAutoFillBackground(true);

    m_overlayLabel = new QLabel(m_overlay);
    m_overlayLabel->setAlignment(Qt::AlignCenter);
    m_overlayLabel->setStyleSheet(
        "color: white; font-size: 18px; font-weight: bold; background: transparent;");
  }
  m_overlayLabel->setText(message);
  m_overlay->setGeometry(0, 0, width(), height());
  m_overlayLabel->setGeometry(0, 0, width(), height());
  m_overlay->raise();
  m_overlay->show();
}

void MainWindow::hideOverlay() {
  if (m_overlay) {
    m_overlay->hide();
  }
}

void MainWindow::resizeEvent(QResizeEvent* event) {
  QMainWindow::resizeEvent(event);
  if (m_overlay && m_overlay->isVisible()) {
    m_overlay->setGeometry(0, 0, width(), height());
    m_overlayLabel->setGeometry(0, 0, width(), height());
  }
}

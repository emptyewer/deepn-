#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "junction_plot_widget.h"
#include "junction_table_model.h"
#include "gene_detail_panel.h"
#include "comparison_manager.h"

#include <gene_selector_widget.h>
#include <mrna_track_widget.h>
#include <export_engine.h>
#include <batch_runner.h>

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDirIterator>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QKeyEvent>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QShortcut>
#include <QSortFilterProxyModel>
#include <QSplitter>
#include <QStatusBar>
#include <QTableView>
#include <QTextStream>
#include <QToolBar>
#include <QToolButton>
#include <QVBoxLayout>

// ────────────────────────────────────────────────────────────────────
// Construction / Destruction
// ────────────────────────────────────────────────────────────────────

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle("MultiQuery++");
    resize(1200, 800);

    setupUI();
    setupMenuBar();
    setupToolBar();
    connectSignals();

    statusBar()->showMessage("Ready -- load a working directory or dataset to begin.");
}

MainWindow::~MainWindow()
{
    m_primaryLoader.close();
    m_secondaryLoader.close();
    delete ui;
}

// ────────────────────────────────────────────────────────────────────
// UI Construction (all programmatic)
// ────────────────────────────────────────────────────────────────────

void MainWindow::setupUI()
{
    auto* central = new QWidget(this);
    auto* rootLayout = new QVBoxLayout(central);
    rootLayout->setContentsMargins(4, 4, 4, 4);
    rootLayout->setSpacing(4);

    // ── Top bar ──────────────────────────────────────────────────
    auto* topBar = new QHBoxLayout;
    topBar->setSpacing(8);

    m_geneSelector = new deepn::GeneSelectorWidget(this);
    topBar->addWidget(m_geneSelector, 1);

    topBar->addWidget(new QLabel("Dataset:", this));
    m_datasetCombo = new QComboBox(this);
    m_datasetCombo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    m_datasetCombo->setMinimumWidth(180);
    topBar->addWidget(m_datasetCombo);

    m_compareCheck = new QCheckBox("Compare", this);
    m_compareCheck->setToolTip("Enable side-by-side dataset comparison");
    topBar->addWidget(m_compareCheck);

    m_secondaryLabel = new QLabel("vs:", this);
    m_secondaryLabel->setVisible(false);
    topBar->addWidget(m_secondaryLabel);

    m_secondaryDatasetCombo = new QComboBox(this);
    m_secondaryDatasetCombo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    m_secondaryDatasetCombo->setMinimumWidth(180);
    m_secondaryDatasetCombo->setVisible(false);
    topBar->addWidget(m_secondaryDatasetCombo);

    rootLayout->addLayout(topBar);

    // ── Main content area ────────────────────────────────────────
    m_mainSplitter = new QSplitter(Qt::Horizontal, this);

    // Left side: charts + track + table
    m_leftSplitter = new QSplitter(Qt::Vertical, this);

    // Chart area -- contains either single plot or split comparison
    m_chartSplitter = new QSplitter(Qt::Horizontal, this);

    // Primary chart container (plot + mRNA track, stacked)
    m_primaryChartContainer = new QWidget(this);
    auto* primaryLayout = new QVBoxLayout(m_primaryChartContainer);
    primaryLayout->setContentsMargins(0, 0, 0, 0);
    primaryLayout->setSpacing(0);

    m_primaryPlot = new JunctionPlotWidget(this);
    m_primaryPlot->setMinimumHeight(200);
    primaryLayout->addWidget(m_primaryPlot, 3);

    m_primaryTrack = new deepn::MRNATrackWidget(this);
    m_primaryTrack->setFixedHeight(50);
    primaryLayout->addWidget(m_primaryTrack);

    m_chartSplitter->addWidget(m_primaryChartContainer);

    // Secondary chart container (hidden by default)
    m_secondaryChartContainer = new QWidget(this);
    auto* secondaryLayout = new QVBoxLayout(m_secondaryChartContainer);
    secondaryLayout->setContentsMargins(0, 0, 0, 0);
    secondaryLayout->setSpacing(0);

    m_secondaryPlot = new JunctionPlotWidget(this);
    m_secondaryPlot->setMinimumHeight(200);
    secondaryLayout->addWidget(m_secondaryPlot, 3);

    m_secondaryTrack = new deepn::MRNATrackWidget(this);
    m_secondaryTrack->setFixedHeight(50);
    secondaryLayout->addWidget(m_secondaryTrack);

    m_chartSplitter->addWidget(m_secondaryChartContainer);
    m_secondaryChartContainer->setVisible(false);

    m_leftSplitter->addWidget(m_chartSplitter);

    // Junction table
    m_tableModel = new JunctionTableModel(this);
    m_tableView = new QTableView(this);
    m_tableView->setModel(m_tableModel);
    m_tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_tableView->setSelectionMode(QAbstractItemView::SingleSelection);
    m_tableView->setSortingEnabled(true);
    m_tableView->setAlternatingRowColors(true);
    m_tableView->horizontalHeader()->setStretchLastSection(true);
    m_tableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    m_tableView->verticalHeader()->setDefaultSectionSize(22);
    m_tableView->setMinimumHeight(120);

    m_leftSplitter->addWidget(m_tableView);
    m_leftSplitter->setStretchFactor(0, 3);
    m_leftSplitter->setStretchFactor(1, 2);

    m_mainSplitter->addWidget(m_leftSplitter);

    // Right side: gene detail panel
    m_detailPanel = new GeneDetailPanel(this);
    m_detailPanel->setMinimumWidth(220);
    m_detailPanel->setMaximumWidth(350);
    m_mainSplitter->addWidget(m_detailPanel);

    m_mainSplitter->setStretchFactor(0, 4);
    m_mainSplitter->setStretchFactor(1, 1);

    rootLayout->addWidget(m_mainSplitter, 1);

    // ── Bottom bar ───────────────────────────────────────────────
    auto* bottomBar = new QHBoxLayout;
    bottomBar->setSpacing(8);

    m_collapseCheck = new QCheckBox("Collapse by Position", this);
    m_collapseCheck->setToolTip("Merge junctions at the same position (sum PPM across QueryStart variants)");
    bottomBar->addWidget(m_collapseCheck);

    m_inFrameOnlyCheck = new QCheckBox("In-Frame Only", this);
    m_inFrameOnlyCheck->setToolTip("Show only in-frame (+0_frame) junctions");
    bottomBar->addWidget(m_inFrameOnlyCheck);

    bottomBar->addStretch(1);

    auto* exportCSVBtn = new QPushButton("Export CSV", this);
    exportCSVBtn->setToolTip("Export junction results table to CSV");
    connect(exportCSVBtn, &QPushButton::clicked, this, &MainWindow::onExportCSV);
    bottomBar->addWidget(exportCSVBtn);

    auto* exportFigBtn = new QPushButton("Export Figure", this);
    exportFigBtn->setToolTip("Export junction plot as SVG/PDF/PNG");
    connect(exportFigBtn, &QPushButton::clicked, this, &MainWindow::onExportFigure);
    bottomBar->addWidget(exportFigBtn);

    auto* batchBtn = new QPushButton("Batch", this);
    batchBtn->setToolTip("Run batch analysis on top N genes from DESeq2 results");
    connect(batchBtn, &QPushButton::clicked, this, &MainWindow::onBatchRun);
    bottomBar->addWidget(batchBtn);

    m_batchProgress = new QProgressBar(this);
    m_batchProgress->setVisible(false);
    m_batchProgress->setMaximumWidth(200);
    bottomBar->addWidget(m_batchProgress);

    rootLayout->addLayout(bottomBar);

    setCentralWidget(central);

    // ── Comparison manager ───────────────────────────────────────
    m_comparisonMgr = new ComparisonManager(this);
    m_comparisonMgr->setPrimaryPlot(m_primaryPlot);
    m_comparisonMgr->setSecondaryPlot(m_secondaryPlot);
    m_comparisonMgr->setPrimaryTrack(m_primaryTrack);
    m_comparisonMgr->setSecondaryTrack(m_secondaryTrack);
}

void MainWindow::setupMenuBar()
{
    // File menu
    QMenu* fileMenu = menuBar()->addMenu("&File");

    QAction* openDirAction = fileMenu->addAction("Open Working Directory...");
    openDirAction->setShortcut(QKeySequence("Ctrl+O"));
    connect(openDirAction, &QAction::triggered, this, [this]() {
        QString dir = QFileDialog::getExistingDirectory(this, "Open Working Directory", m_workdir);
        if (!dir.isEmpty()) {
            loadWorkingDirectory(dir);
        }
    });

    QAction* openDbAction = fileMenu->addAction("Open Database...");
    openDbAction->setShortcut(QKeySequence("Ctrl+Shift+O"));
    connect(openDbAction, &QAction::triggered, this, [this]() {
        QString path = QFileDialog::getOpenFileName(this, "Open SQLite Database",
                                                     m_workdir, "SQLite databases (*.sqlite)");
        if (!path.isEmpty()) {
            loadDataset(path);
        }
    });

    QAction* openRefAction = fileMenu->addAction("Open Gene Reference...");
    connect(openRefAction, &QAction::triggered, this, [this]() {
        QString path = QFileDialog::getOpenFileName(this, "Open Gene Reference FASTA",
                                                     m_workdir, "FASTA files (*.fa *.fasta *.fna);;All files (*)");
        if (!path.isEmpty()) {
            loadGeneReference(path);
        }
    });

    QAction* openResultsAction = fileMenu->addAction("Load DESeq2 Results...");
    connect(openResultsAction, &QAction::triggered, this, [this]() {
        QString path = QFileDialog::getOpenFileName(this, "Open DESeq2 Results CSV",
                                                     m_workdir, "CSV files (*.csv);;All files (*)");
        if (!path.isEmpty()) {
            loadDESeq2Results(path);
        }
    });

    fileMenu->addSeparator();

    QAction* exportCSVAction = fileMenu->addAction("Export CSV...");
    exportCSVAction->setShortcut(QKeySequence("Ctrl+E"));
    connect(exportCSVAction, &QAction::triggered, this, &MainWindow::onExportCSV);

    QAction* exportFigAction = fileMenu->addAction("Export Figure...");
    exportFigAction->setShortcut(QKeySequence("Ctrl+Shift+E"));
    connect(exportFigAction, &QAction::triggered, this, &MainWindow::onExportFigure);

    fileMenu->addSeparator();

    QAction* quitAction = fileMenu->addAction("Quit");
    quitAction->setShortcut(QKeySequence::Quit);
    connect(quitAction, &QAction::triggered, this, &QWidget::close);

    // View menu
    QMenu* viewMenu = menuBar()->addMenu("&View");

    QAction* collapseAction = viewMenu->addAction("Collapse by Position");
    collapseAction->setCheckable(true);
    connect(collapseAction, &QAction::toggled, m_collapseCheck, &QCheckBox::setChecked);
    connect(m_collapseCheck, &QCheckBox::toggled, collapseAction, &QAction::setChecked);

    QAction* inFrameAction = viewMenu->addAction("In-Frame Only");
    inFrameAction->setCheckable(true);
    connect(inFrameAction, &QAction::toggled, m_inFrameOnlyCheck, &QCheckBox::setChecked);
    connect(m_inFrameOnlyCheck, &QCheckBox::toggled, inFrameAction, &QAction::setChecked);

    viewMenu->addSeparator();

    QAction* compareAction = viewMenu->addAction("Compare Datasets");
    compareAction->setCheckable(true);
    connect(compareAction, &QAction::toggled, m_compareCheck, &QCheckBox::setChecked);
    connect(m_compareCheck, &QCheckBox::toggled, compareAction, &QAction::setChecked);

    // Analysis menu
    QMenu* analysisMenu = menuBar()->addMenu("&Analysis");

    QAction* batchAction = analysisMenu->addAction("Batch Analysis...");
    batchAction->setShortcut(QKeySequence("Ctrl+B"));
    connect(batchAction, &QAction::triggered, this, &MainWindow::onBatchRun);
}

void MainWindow::setupToolBar()
{
    // Keyboard shortcuts for navigation
    auto* prevShortcut = new QShortcut(QKeySequence(Qt::Key_Left), this);
    connect(prevShortcut, &QShortcut::activated, this, [this]() {
        m_geneSelector->goToPrev();
    });

    auto* nextShortcut = new QShortcut(QKeySequence(Qt::Key_Right), this);
    connect(nextShortcut, &QShortcut::activated, this, [this]() {
        m_geneSelector->goToNext();
    });

    auto* escShortcut = new QShortcut(QKeySequence(Qt::Key_Escape), this);
    connect(escShortcut, &QShortcut::activated, this, [this]() {
        m_primaryPlot->highlightPosition(-1);
        if (m_secondaryPlot->isVisible()) {
            m_secondaryPlot->highlightPosition(-1);
        }
        m_tableView->clearSelection();
        m_detailPanel->setSelectedJunction({});
    });
}

void MainWindow::connectSignals()
{
    // Gene selector
    connect(m_geneSelector, &deepn::GeneSelectorWidget::geneSelected,
            this, &MainWindow::onGeneSelected);

    // Primary plot clicked
    connect(m_primaryPlot, &JunctionPlotWidget::junctionClicked,
            this, &MainWindow::onJunctionClicked);

    // Secondary plot clicked
    connect(m_secondaryPlot, &JunctionPlotWidget::junctionClicked,
            this, &MainWindow::onJunctionClicked);

    // Table selection
    connect(m_tableView->selectionModel(), &QItemSelectionModel::currentRowChanged,
            this, &MainWindow::onTableRowSelected);

    // Filters
    connect(m_collapseCheck, &QCheckBox::toggled,
            this, &MainWindow::onCollapseToggled);
    connect(m_inFrameOnlyCheck, &QCheckBox::toggled,
            this, &MainWindow::onInFrameOnlyToggled);

    // Compare mode
    connect(m_compareCheck, &QCheckBox::toggled,
            this, &MainWindow::onCompareToggled);

    // Dataset combos
    connect(m_datasetCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onDatasetChanged);
    connect(m_secondaryDatasetCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onSecondaryDatasetChanged);

    // Sync chart zoom to mRNA track
    connect(m_primaryPlot, &JunctionPlotWidget::rangeChanged,
            this, [this](qreal xMin, qreal xMax) {
        m_primaryTrack->setVisibleRange(static_cast<int>(xMin), static_cast<int>(xMax));
    });

    connect(m_secondaryPlot, &JunctionPlotWidget::rangeChanged,
            this, [this](qreal xMin, qreal xMax) {
        m_secondaryTrack->setVisibleRange(static_cast<int>(xMin), static_cast<int>(xMax));
    });
}

// ────────────────────────────────────────────────────────────────────
// Data Loading
// ────────────────────────────────────────────────────────────────────

void MainWindow::loadWorkingDirectory(const QString& workdir)
{
    m_workdir = workdir;
    autoDiscoverFiles();
    statusBar()->showMessage(QString("Loaded working directory: %1").arg(workdir));
}

void MainWindow::loadDataset(const QString& dbPath)
{
    if (m_datasetPaths.contains(dbPath)) {
        // Already loaded, just select it
        int idx = m_datasetPaths.indexOf(dbPath);
        m_datasetCombo->setCurrentIndex(idx);
        return;
    }

    m_datasetPaths.append(dbPath);
    QFileInfo fi(dbPath);
    m_datasetCombo->addItem(fi.fileName(), dbPath);
    m_secondaryDatasetCombo->addItem(fi.fileName(), dbPath);

    // If this is the first dataset, select it and load
    if (m_datasetPaths.size() == 1) {
        m_datasetCombo->setCurrentIndex(0);
    }

    statusBar()->showMessage(QString("Added dataset: %1").arg(fi.fileName()));
}

void MainWindow::loadGeneReference(const QString& fastaPath)
{
    if (!m_annotationDB.loadFromFasta(fastaPath)) {
        QMessageBox::warning(this, "Error Loading Gene Reference",
                             QString("Failed to load gene reference:\n%1")
                                 .arg(m_annotationDB.lastError()));
        return;
    }

    // Populate gene selector with all gene names
    QStringList genes = m_annotationDB.allGeneNames();
    m_geneSelector->setGeneList(genes);

    statusBar()->showMessage(QString("Loaded %1 gene annotations from %2")
                                 .arg(m_annotationDB.count())
                                 .arg(QFileInfo(fastaPath).fileName()));
}

void MainWindow::loadDESeq2Results(const QString& csvPath)
{
    m_deseq2Results = parseDESeq2CSV(csvPath);
    if (m_deseq2Results.isEmpty()) {
        QMessageBox::warning(this, "Error Loading Results",
                             "No valid results found in CSV file.");
        return;
    }

    m_geneSelector->setDESeq2Results(m_deseq2Results);

    statusBar()->showMessage(QString("Loaded %1 DESeq2 results from %2")
                                 .arg(m_deseq2Results.size())
                                 .arg(QFileInfo(csvPath).fileName()));
}

void MainWindow::displayGene(const QString& gene)
{
    m_geneSelector->selectGene(gene);
}

// ────────────────────────────────────────────────────────────────────
// Auto-discovery
// ────────────────────────────────────────────────────────────────────

void MainWindow::autoDiscoverFiles()
{
    if (m_workdir.isEmpty()) return;

    QDir dir(m_workdir);

    // Discover SQLite databases in analyzed_files/
    QDir analyzedDir(dir.filePath("analyzed_files"));
    if (analyzedDir.exists()) {
        QStringList filters;
        filters << "*.sqlite";
        QFileInfoList dbFiles = analyzedDir.entryInfoList(filters, QDir::Files, QDir::Name);
        for (const QFileInfo& fi : dbFiles) {
            loadDataset(fi.absoluteFilePath());
        }
    }

    // Discover gene reference FASTA in data/ or at workdir level
    QStringList fastaExtensions = {"*.fa", "*.fasta", "*.fna"};
    auto findFasta = [&](const QString& searchDir) -> QString {
        QDir d(searchDir);
        if (!d.exists()) return {};
        QFileInfoList fastaFiles = d.entryInfoList(fastaExtensions, QDir::Files, QDir::Name);
        if (!fastaFiles.isEmpty()) return fastaFiles.first().absoluteFilePath();
        return {};
    };

    QString fastaPath = findFasta(dir.filePath("data"));
    if (fastaPath.isEmpty()) {
        fastaPath = findFasta(m_workdir);
    }
    if (!fastaPath.isEmpty() && m_annotationDB.count() == 0) {
        loadGeneReference(fastaPath);
    }

    // Discover DESeq2 results CSV
    QStringList csvExtensions = {"*.csv"};
    auto findCSV = [&](const QString& searchDir, const QString& pattern) -> QString {
        QDir d(searchDir);
        if (!d.exists()) return {};
        QFileInfoList csvFiles = d.entryInfoList(csvExtensions, QDir::Files, QDir::Name);
        for (const QFileInfo& fi : csvFiles) {
            if (fi.fileName().contains(pattern, Qt::CaseInsensitive)) {
                return fi.absoluteFilePath();
            }
        }
        if (!csvFiles.isEmpty()) return csvFiles.first().absoluteFilePath();
        return {};
    };

    if (m_deseq2Results.isEmpty()) {
        QString resultsPath = findCSV(m_workdir, "deseq2");
        if (resultsPath.isEmpty()) {
            resultsPath = findCSV(dir.filePath("results"), "deseq2");
        }
        if (!resultsPath.isEmpty()) {
            loadDESeq2Results(resultsPath);
        }
    }

    // If gene selector has items and primary loader is open, display first gene
    if (m_geneSelector->currentGene().isEmpty() && m_datasetPaths.size() > 0) {
        // Load available genes from the primary dataset
        if (m_primaryLoader.open(m_datasetPaths.first())) {
            QStringList availGenes = m_primaryLoader.availableGenes();
            if (!availGenes.isEmpty()) {
                if (m_geneSelector->currentGene().isEmpty()) {
                    m_geneSelector->selectGene(availGenes.first());
                }
            }
        }
    }
}

void MainWindow::populateDatasetCombo()
{
    m_datasetCombo->clear();
    m_secondaryDatasetCombo->clear();
    for (const QString& path : m_datasetPaths) {
        QFileInfo fi(path);
        m_datasetCombo->addItem(fi.fileName(), path);
        m_secondaryDatasetCombo->addItem(fi.fileName(), path);
    }
}

// ────────────────────────────────────────────────────────────────────
// Slot implementations
// ────────────────────────────────────────────────────────────────────

void MainWindow::onGeneSelected(const QString& gene)
{
    if (gene.isEmpty()) return;

    // Check if primary loader is ready
    int idx = m_datasetCombo->currentIndex();
    if (idx < 0 || idx >= m_datasetPaths.size()) return;

    // Ensure the loader is connected to the right database
    QString dbPath = m_datasetPaths.at(idx);
    if (m_primaryLoader.databasePath() != dbPath) {
        m_primaryLoader.close();
        if (!m_primaryLoader.open(dbPath)) {
            statusBar()->showMessage(QString("Failed to open database: %1").arg(m_primaryLoader.lastError()));
            return;
        }
    }

    // Get annotation for the gene
    deepn::GeneAnnotation annotation = m_annotationDB.findByGeneName(gene);
    if (!annotation.isValid()) {
        // Try by refseq
        annotation = m_annotationDB.findByRefseq(gene);
    }

    // Load junction profile
    m_currentProfile = m_primaryLoader.loadGeneJunctions(gene, annotation);
    m_currentProfile.datasetLabel = m_datasetCombo->currentText();

    refreshDisplay();

    // If comparison mode is active, also load secondary
    if (m_compareCheck->isChecked()) {
        refreshSecondaryDisplay();
    }

    statusBar()->showMessage(QString("Displaying %1: %2 junctions")
                                 .arg(gene)
                                 .arg(m_currentProfile.sites.size()));
}

void MainWindow::onJunctionClicked(int position)
{
    // Find the junction in the table and select it
    for (int row = 0; row < m_tableModel->rowCount(); ++row) {
        deepn::JunctionSite site = m_tableModel->junctionAt(row);
        if (site.position == position) {
            QModelIndex idx = m_tableModel->index(row, 0);
            m_tableView->selectionModel()->setCurrentIndex(idx, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
            m_tableView->scrollTo(idx);
            break;
        }
    }

    // Highlight on both plots
    m_primaryPlot->highlightPosition(position);
    if (m_compareCheck->isChecked()) {
        m_secondaryPlot->highlightPosition(position);
    }
}

void MainWindow::onTableRowSelected(const QModelIndex& current, const QModelIndex& /*previous*/)
{
    if (!current.isValid()) return;

    deepn::JunctionSite site = m_tableModel->junctionAt(current.row());
    m_detailPanel->setSelectedJunction(site);

    // Highlight on chart
    m_primaryPlot->highlightPosition(site.position);
    if (m_compareCheck->isChecked()) {
        m_secondaryPlot->highlightPosition(site.position);
    }
}

void MainWindow::onCollapseToggled(bool checked)
{
    m_collapseByPosition = checked;
    refreshDisplay();
    if (m_compareCheck->isChecked()) {
        refreshSecondaryDisplay();
    }
}

void MainWindow::onInFrameOnlyToggled(bool checked)
{
    m_inFrameOnly = checked;
    refreshDisplay();
    if (m_compareCheck->isChecked()) {
        refreshSecondaryDisplay();
    }
}

void MainWindow::onCompareToggled(bool checked)
{
    m_secondaryChartContainer->setVisible(checked);
    m_secondaryLabel->setVisible(checked);
    m_secondaryDatasetCombo->setVisible(checked);
    m_comparisonMgr->setEnabled(checked);

    if (checked && !m_geneSelector->currentGene().isEmpty()) {
        refreshSecondaryDisplay();
    }
}

void MainWindow::onDatasetChanged(int index)
{
    if (index < 0 || index >= m_datasetPaths.size()) return;

    m_primaryLoader.close();
    QString dbPath = m_datasetPaths.at(index);
    if (!m_primaryLoader.open(dbPath)) {
        statusBar()->showMessage(QString("Failed to open: %1").arg(m_primaryLoader.lastError()));
        return;
    }

    // Reload current gene if one is selected
    QString gene = m_geneSelector->currentGene();
    if (!gene.isEmpty()) {
        onGeneSelected(gene);
    }
}

void MainWindow::onSecondaryDatasetChanged(int index)
{
    if (index < 0 || index >= m_datasetPaths.size()) return;

    m_secondaryLoader.close();
    QString dbPath = m_datasetPaths.at(index);
    if (!m_secondaryLoader.open(dbPath)) {
        statusBar()->showMessage(QString("Failed to open secondary: %1").arg(m_secondaryLoader.lastError()));
        return;
    }

    if (m_compareCheck->isChecked() && !m_geneSelector->currentGene().isEmpty()) {
        refreshSecondaryDisplay();
    }
}

void MainWindow::onExportCSV()
{
    if (m_currentProfile.sites.isEmpty()) {
        QMessageBox::information(this, "No Data", "No junction data to export. Select a gene first.");
        return;
    }

    QString defaultName = QString("%1_junctions.csv").arg(m_currentProfile.annotation.geneName);
    QString filePath = QFileDialog::getSaveFileName(this, "Export Junction CSV",
                                                     QDir(m_workdir).filePath(defaultName),
                                                     "CSV files (*.csv)");
    if (filePath.isEmpty()) return;

    bool ok = false;
    if (m_collapseByPosition) {
        auto collapsed = deepn::SqliteJunctionLoader::collapseByPosition(m_currentProfile.sites);
        ok = deepn::ExportEngine::exportCollapsedCSV(filePath, collapsed);
    } else {
        ok = deepn::ExportEngine::exportJunctionCSV(filePath, m_currentProfile);
    }

    if (ok) {
        statusBar()->showMessage(QString("Exported CSV to %1").arg(filePath));
    } else {
        QMessageBox::warning(this, "Export Failed",
                             QString("Failed to export CSV:\n%1").arg(deepn::ExportEngine::lastError()));
    }
}

void MainWindow::onExportFigure()
{
    if (m_currentProfile.sites.isEmpty()) {
        QMessageBox::information(this, "No Data", "No junction data to export. Select a gene first.");
        return;
    }

    QString defaultName = QString("%1_junction_plot").arg(m_currentProfile.annotation.geneName);
    QString filePath = QFileDialog::getSaveFileName(this, "Export Figure",
                                                     QDir(m_workdir).filePath(defaultName),
                                                     "SVG (*.svg);;PDF (*.pdf);;PNG (*.png)");
    if (filePath.isEmpty()) return;

    bool ok = false;
    QSize figSize(800, 400);

    if (filePath.endsWith(".svg", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigureSVG(filePath, m_primaryPlot, figSize);
    } else if (filePath.endsWith(".pdf", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigurePDF(filePath, m_primaryPlot, figSize);
    } else if (filePath.endsWith(".png", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigurePNG(filePath, m_primaryPlot, figSize, 300);
    } else {
        // Default to SVG
        filePath += ".svg";
        ok = deepn::ExportEngine::exportFigureSVG(filePath, m_primaryPlot, figSize);
    }

    if (ok) {
        statusBar()->showMessage(QString("Exported figure to %1").arg(filePath));
    } else {
        QMessageBox::warning(this, "Export Failed",
                             QString("Failed to export figure:\n%1").arg(deepn::ExportEngine::lastError()));
    }
}

void MainWindow::onBatchRun()
{
    if (m_datasetPaths.isEmpty()) {
        QMessageBox::information(this, "No Data", "No datasets loaded. Open a working directory first.");
        return;
    }

    bool ok = false;
    int topN = QInputDialog::getInt(this, "Batch Analysis",
                                     "Number of top genes to process:",
                                     50, 1, 1000, 1, &ok);
    if (!ok) return;

    QString outDir = QFileDialog::getExistingDirectory(this, "Select Output Directory",
                                                        m_workdir);
    if (outDir.isEmpty()) return;

    // Build gene list from DESeq2 results or available genes
    QStringList geneList;
    if (!m_deseq2Results.isEmpty()) {
        for (int i = 0; i < qMin(topN, m_deseq2Results.size()); ++i) {
            geneList.append(m_deseq2Results.at(i).gene);
        }
    } else {
        // Fall back to available genes from primary loader
        geneList = m_primaryLoader.availableGenes();
        if (geneList.size() > topN) {
            geneList = geneList.mid(0, topN);
        }
    }

    deepn::BatchRunner::Config config;
    config.dbPaths = m_datasetPaths;
    config.geneList = geneList;
    config.outputDir = outDir;
    config.topN = topN;
    config.exportCSV = true;
    config.exportFigures = true;

    auto* runner = new deepn::BatchRunner(this);

    m_batchProgress->setVisible(true);
    m_batchProgress->setRange(0, geneList.size());
    m_batchProgress->setValue(0);

    connect(runner, &deepn::BatchRunner::progressChanged,
            this, [this](int current, int total, const QString& gene) {
        m_batchProgress->setRange(0, total);
        m_batchProgress->setValue(current);
        statusBar()->showMessage(QString("Batch: processing %1 (%2/%3)").arg(gene).arg(current).arg(total));
    });

    connect(runner, &deepn::BatchRunner::finished,
            this, [this, runner](bool success) {
        m_batchProgress->setVisible(false);
        if (success) {
            statusBar()->showMessage("Batch analysis complete.");
            QMessageBox::information(this, "Batch Complete",
                                     "Batch analysis finished successfully.");
        } else {
            statusBar()->showMessage("Batch analysis failed.");
        }
        runner->deleteLater();
    });

    connect(runner, &deepn::BatchRunner::error,
            this, [this](const QString& message) {
        QMessageBox::warning(this, "Batch Error", message);
    });

    runner->start(config);
}

// ────────────────────────────────────────────────────────────────────
// Display Refresh
// ────────────────────────────────────────────────────────────────────

void MainWindow::refreshDisplay()
{
    const deepn::GeneAnnotation& annot = m_currentProfile.annotation;

    // Update annotation-based views
    if (annot.isValid()) {
        m_detailPanel->setAnnotation(annot);
        m_primaryTrack->setAnnotation(annot);
        m_primaryTrack->setVisibleRange(0, annot.mRNALength);
    }

    // Apply filters to get display data
    QVector<deepn::JunctionSite> filtered;
    applyFilters(m_currentProfile, filtered);

    if (m_collapseByPosition) {
        auto collapsed = deepn::SqliteJunctionLoader::collapseByPosition(filtered);
        m_primaryPlot->setCollapsedView(collapsed);
        m_tableModel->setCollapsedJunctions(collapsed);
    } else {
        deepn::GeneJunctionProfile displayProfile = m_currentProfile;
        displayProfile.sites = filtered;
        m_primaryPlot->setProfile(displayProfile);
        m_tableModel->setJunctions(filtered);
    }

    // Profile summary
    int inFrameCount = 0;
    double maxPpm = 0.0;
    for (const auto& site : m_currentProfile.sites) {
        if (site.isInFrame()) ++inFrameCount;
        if (site.ppm > maxPpm) maxPpm = site.ppm;
    }
    m_detailPanel->setProfileSummary(m_currentProfile.sites.size(), inFrameCount, maxPpm);

    // Resize table columns to content
    m_tableView->resizeColumnsToContents();
}

void MainWindow::refreshSecondaryDisplay()
{
    if (!m_compareCheck->isChecked()) return;

    int idx = m_secondaryDatasetCombo->currentIndex();
    if (idx < 0 || idx >= m_datasetPaths.size()) return;

    QString dbPath = m_datasetPaths.at(idx);
    if (m_secondaryLoader.databasePath() != dbPath) {
        m_secondaryLoader.close();
        if (!m_secondaryLoader.open(dbPath)) return;
    }

    QString gene = m_geneSelector->currentGene();
    if (gene.isEmpty()) return;

    deepn::GeneAnnotation annotation = m_annotationDB.findByGeneName(gene);
    if (!annotation.isValid()) {
        annotation = m_annotationDB.findByRefseq(gene);
    }

    m_secondaryProfile = m_secondaryLoader.loadGeneJunctions(gene, annotation);
    m_secondaryProfile.datasetLabel = m_secondaryDatasetCombo->currentText();

    if (annotation.isValid()) {
        m_secondaryTrack->setAnnotation(annotation);
        m_secondaryTrack->setVisibleRange(0, annotation.mRNALength);
    }

    QVector<deepn::JunctionSite> filtered;
    applyFilters(m_secondaryProfile, filtered);

    if (m_collapseByPosition) {
        auto collapsed = deepn::SqliteJunctionLoader::collapseByPosition(filtered);
        m_secondaryPlot->setCollapsedView(collapsed);
    } else {
        deepn::GeneJunctionProfile displayProfile = m_secondaryProfile;
        displayProfile.sites = filtered;
        m_secondaryPlot->setProfile(displayProfile);
    }

    m_comparisonMgr->syncRanges();
}

void MainWindow::applyFilters(const deepn::GeneJunctionProfile& profile,
                               QVector<deepn::JunctionSite>& filtered) const
{
    filtered.clear();
    filtered.reserve(profile.sites.size());

    for (const auto& site : profile.sites) {
        if (m_inFrameOnly && !site.isInFrame()) {
            continue;
        }
        filtered.append(site);
    }
}

QVector<deepn::DESeq2Result> MainWindow::parseDESeq2CSV(const QString& path) const
{
    QVector<deepn::DESeq2Result> results;

    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return results;
    }

    QTextStream in(&file);
    QString header = in.readLine();  // skip header

    // Detect separator
    QChar sep = ',';
    if (header.contains('\t')) sep = '\t';

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        QStringList fields = line.split(sep);
        if (fields.size() < 7) continue;

        deepn::DESeq2Result r;
        r.gene = fields[0].trimmed().remove('"');
        r.baseMean = fields[1].toDouble();
        r.log2FoldChange = fields[2].toDouble();
        r.lfcSE = fields[3].toDouble();
        r.stat = fields[4].toDouble();
        r.pvalue = fields[5].toDouble();
        r.padj = fields[6].toDouble();
        if (fields.size() > 7) {
            r.enrichment = fields[7].trimmed().remove('"');
        }

        results.append(r);
    }

    return results;
}

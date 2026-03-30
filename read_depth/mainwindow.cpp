#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "batch_runner.h"

#include <QApplication>
#include <QButtonGroup>
#include <QDebug>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QProgressDialog>
#include <QStatusBar>
#include <QTextStream>
#include <QToolBar>
#include <QVBoxLayout>

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle("ReadDepth++");
    resize(1100, 700);

    setupUI();
    setupMenuBar();
    setupStatusBar();
    connectSignals();
}

MainWindow::~MainWindow()
{
    delete ui;
}

// ---------------------------------------------------------------------------
// UI Setup
// ---------------------------------------------------------------------------

void MainWindow::setupUI()
{
    auto* central = new QWidget(this);
    setCentralWidget(central);

    auto* mainLayout = new QVBoxLayout(central);
    mainLayout->setContentsMargins(6, 6, 6, 6);
    mainLayout->setSpacing(4);

    // === Top bar row 1: Gene selector + dataset ===
    auto* topBar1 = new QHBoxLayout();
    topBar1->setSpacing(8);

    m_geneSelector = new deepn::GeneSelectorWidget(this);
    topBar1->addWidget(m_geneSelector, 1);

    auto* datasetLabel = new QLabel("Dataset:", this);
    topBar1->addWidget(datasetLabel);

    m_datasetCombo = new QComboBox(this);
    m_datasetCombo->setMinimumWidth(180);
    topBar1->addWidget(m_datasetCombo);

    mainLayout->addLayout(topBar1);

    // === Top bar row 2: Interval controls + comparison controls ===
    auto* topBar2 = new QHBoxLayout();
    topBar2->setSpacing(8);

    auto* widthLabel = new QLabel("Interval Width:", this);
    topBar2->addWidget(widthLabel);

    m_intervalWidthSpin = new QSpinBox(this);
    m_intervalWidthSpin->setRange(10, 100);
    m_intervalWidthSpin->setValue(m_intervalWidth);
    m_intervalWidthSpin->setSuffix(" bp");
    m_intervalWidthSpin->setSingleStep(5);
    m_intervalWidthSpin->setFixedWidth(90);
    topBar2->addWidget(m_intervalWidthSpin);

    auto* spacingLabel = new QLabel("Interval Spacing:", this);
    topBar2->addWidget(spacingLabel);

    m_intervalSpacingSpin = new QSpinBox(this);
    m_intervalSpacingSpin->setRange(25, 500);
    m_intervalSpacingSpin->setValue(m_intervalSpacing);
    m_intervalSpacingSpin->setSuffix(" nt");
    m_intervalSpacingSpin->setSingleStep(25);
    m_intervalSpacingSpin->setFixedWidth(90);
    topBar2->addWidget(m_intervalSpacingSpin);

    topBar2->addSpacing(20);

    m_compareCheck = new QCheckBox("Compare", this);
    topBar2->addWidget(m_compareCheck);

    // Comparison mode panel (hidden by default)
    m_comparisonPanel = new QWidget(this);
    auto* compLayout = new QHBoxLayout(m_comparisonPanel);
    compLayout->setContentsMargins(0, 0, 0, 0);
    compLayout->setSpacing(4);

    m_splitRadio = new QRadioButton("Split", m_comparisonPanel);
    m_overlayRadio = new QRadioButton("Overlay", m_comparisonPanel);
    m_splitRadio->setChecked(true);

    auto* modeGroup = new QButtonGroup(this);
    modeGroup->addButton(m_splitRadio);
    modeGroup->addButton(m_overlayRadio);

    compLayout->addWidget(m_splitRadio);
    compLayout->addWidget(m_overlayRadio);

    auto* secLabel = new QLabel("vs:", m_comparisonPanel);
    compLayout->addWidget(secLabel);

    m_secondaryDatasetCombo = new QComboBox(m_comparisonPanel);
    m_secondaryDatasetCombo->setMinimumWidth(180);
    compLayout->addWidget(m_secondaryDatasetCombo);

    m_comparisonPanel->setVisible(false);
    topBar2->addWidget(m_comparisonPanel);

    topBar2->addStretch();

    mainLayout->addLayout(topBar2);

    // === Main area ===
    m_mainSplitter = new QSplitter(Qt::Vertical, this);

    // Chart area (horizontal splitter for split mode)
    m_chartSplitter = new QSplitter(Qt::Horizontal, this);

    // Primary panel
    m_primaryPanel = new QWidget(this);
    auto* primaryLayout = new QVBoxLayout(m_primaryPanel);
    primaryLayout->setContentsMargins(0, 0, 0, 0);
    primaryLayout->setSpacing(0);

    m_primaryPlot = new DepthPlotWidget(m_primaryPanel);
    m_primaryPlot->setMinimumHeight(200);
    primaryLayout->addWidget(m_primaryPlot, 1);

    m_primaryTrack = new deepn::MRNATrackWidget(m_primaryPanel);
    primaryLayout->addWidget(m_primaryTrack);

    m_chartSplitter->addWidget(m_primaryPanel);

    // Secondary panel (hidden by default)
    m_secondaryPanel = new QWidget(this);
    auto* secondaryLayout = new QVBoxLayout(m_secondaryPanel);
    secondaryLayout->setContentsMargins(0, 0, 0, 0);
    secondaryLayout->setSpacing(0);

    m_secondaryPlot = new DepthPlotWidget(m_secondaryPanel);
    m_secondaryPlot->setMinimumHeight(200);
    secondaryLayout->addWidget(m_secondaryPlot, 1);

    m_secondaryTrack = new deepn::MRNATrackWidget(m_secondaryPanel);
    secondaryLayout->addWidget(m_secondaryTrack);

    m_chartSplitter->addWidget(m_secondaryPanel);
    m_secondaryPanel->setVisible(false);

    m_mainSplitter->addWidget(m_chartSplitter);

    // Boundary info panel
    auto* infoPanel = new QWidget(this);
    auto* infoLayout = new QVBoxLayout(infoPanel);
    infoLayout->setContentsMargins(8, 4, 8, 4);

    m_boundaryLabel = new QLabel(this);
    m_boundaryLabel->setTextFormat(Qt::RichText);
    m_boundaryLabel->setWordWrap(true);
    m_boundaryLabel->setText("Load a gene to see 3' boundary analysis.");
    infoLayout->addWidget(m_boundaryLabel);

    m_insertExtentLabel = new QLabel(this);
    m_insertExtentLabel->setTextFormat(Qt::RichText);
    m_insertExtentLabel->setWordWrap(true);
    m_insertExtentLabel->setVisible(false);
    infoLayout->addWidget(m_insertExtentLabel);

    m_mainSplitter->addWidget(infoPanel);

    // Set splitter proportions: chart area gets 80%, info panel 20%
    m_mainSplitter->setStretchFactor(0, 4);
    m_mainSplitter->setStretchFactor(1, 1);

    mainLayout->addWidget(m_mainSplitter, 1);

    // === Bottom button bar ===
    auto* bottomBar = new QHBoxLayout();
    bottomBar->setSpacing(8);

    auto* resetZoomBtn = new QPushButton("Reset Zoom", this);
    bottomBar->addWidget(resetZoomBtn);
    connect(resetZoomBtn, &QPushButton::clicked, this, &MainWindow::onResetZoom);

    bottomBar->addStretch();

    auto* exportCSVBtn = new QPushButton("Export CSV", this);
    bottomBar->addWidget(exportCSVBtn);
    connect(exportCSVBtn, &QPushButton::clicked, this, &MainWindow::onExportCSV);

    auto* exportFigBtn = new QPushButton("Export Figure", this);
    bottomBar->addWidget(exportFigBtn);
    connect(exportFigBtn, &QPushButton::clicked, this, &MainWindow::onExportFigure);

    auto* batchBtn = new QPushButton("Batch Report", this);
    bottomBar->addWidget(batchBtn);
    connect(batchBtn, &QPushButton::clicked, this, &MainWindow::onBatchRun);

    mainLayout->addLayout(bottomBar);

    // Sync controller for split mode
    m_syncController = new deepn::SyncController(this);
}

void MainWindow::setupMenuBar()
{
    QMenuBar* mb = menuBar();

    // File menu
    QMenu* fileMenu = mb->addMenu("&File");

    fileMenu->addAction("Open &Working Directory...", this, [this]() {
        QString dir = QFileDialog::getExistingDirectory(
            this, "Open Working Directory", m_workdir);
        if (!dir.isEmpty()) {
            loadWorkingDirectory(dir);
        }
    });

    fileMenu->addAction("Open &Dataset...", this, [this]() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open SQLite Dataset", m_workdir,
            "SQLite Database (*.sqlite *.db);;All Files (*)");
        if (!path.isEmpty()) {
            loadDataset(path);
        }
    });

    fileMenu->addAction("Open Gene &Reference...", this, [this]() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open Gene Reference FASTA", m_workdir,
            "FASTA Files (*.fa *.fasta *.fna);;All Files (*)");
        if (!path.isEmpty()) {
            loadGeneReference(path);
        }
    });

    fileMenu->addAction("Load DESeq2 &Results...", this, [this]() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open DESeq2 Results CSV", m_workdir,
            "CSV Files (*.csv *.tsv);;All Files (*)");
        if (!path.isEmpty()) {
            loadDESeq2Results(path);
        }
    });

    fileMenu->addAction("Load &Junction Data...", this, [this]() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open Junction CSV (from MultiQuery++)", m_workdir,
            "CSV Files (*.csv);;All Files (*)");
        if (!path.isEmpty()) {
            loadJunctionData(path);
        }
    });

    fileMenu->addSeparator();

    fileMenu->addAction("Export &CSV...", this, &MainWindow::onExportCSV);
    fileMenu->addAction("Export &Figure...", this, &MainWindow::onExportFigure);

    fileMenu->addSeparator();

    fileMenu->addAction("&Quit", QKeySequence::Quit, qApp, &QApplication::quit);

    // View menu
    QMenu* viewMenu = mb->addMenu("&View");
    viewMenu->addAction("Reset &Zoom", QKeySequence("Ctrl+0"), this, &MainWindow::onResetZoom);
    viewMenu->addAction("Zoom &In", QKeySequence::ZoomIn, this, [this]() {
        // Simulate zoom in by adjusting axis
        if (m_primaryPlot && m_primaryPlot->depthChart()) {
            auto axes = m_primaryPlot->depthChart()->axes(Qt::Horizontal);
            if (!axes.isEmpty()) {
                auto* xAxis = qobject_cast<QValueAxis*>(axes.first());
                if (xAxis) {
                    qreal range = xAxis->max() - xAxis->min();
                    qreal center = (xAxis->max() + xAxis->min()) / 2.0;
                    qreal newRange = range * 0.8;
                    if (newRange < 200.0) newRange = 200.0;
                    xAxis->setRange(center - newRange / 2, center + newRange / 2);
                }
            }
        }
    });
    viewMenu->addAction("Zoom &Out", QKeySequence::ZoomOut, this, [this]() {
        if (m_primaryPlot && m_primaryPlot->depthChart()) {
            auto axes = m_primaryPlot->depthChart()->axes(Qt::Horizontal);
            if (!axes.isEmpty()) {
                auto* xAxis = qobject_cast<QValueAxis*>(axes.first());
                if (xAxis) {
                    qreal range = xAxis->max() - xAxis->min();
                    qreal center = (xAxis->max() + xAxis->min()) / 2.0;
                    qreal newRange = range * 1.25;
                    xAxis->setRange(center - newRange / 2, center + newRange / 2);
                }
            }
        }
    });

    // Tools menu
    QMenu* toolsMenu = mb->addMenu("&Tools");
    toolsMenu->addAction("&Batch Report...", this, &MainWindow::onBatchRun);
}

void MainWindow::setupStatusBar()
{
    m_cursorLabel = new QLabel("Position: -- | Depth: -- reads", this);
    statusBar()->addWidget(m_cursorLabel, 1);

    m_zoomLabel = new QLabel("Zoom: 100%", this);
    statusBar()->addPermanentWidget(m_zoomLabel);
}

void MainWindow::connectSignals()
{
    connect(m_geneSelector, &deepn::GeneSelectorWidget::geneSelected,
            this, &MainWindow::onGeneSelected);

    connect(m_intervalWidthSpin, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &MainWindow::onIntervalWidthChanged);
    connect(m_intervalSpacingSpin, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &MainWindow::onIntervalSpacingChanged);

    connect(m_compareCheck, &QCheckBox::toggled,
            this, &MainWindow::onCompareToggled);
    connect(m_overlayRadio, &QRadioButton::toggled,
            this, &MainWindow::onOverlayToggled);
    connect(m_splitRadio, &QRadioButton::toggled,
            this, &MainWindow::onSplitToggled);

    connect(m_datasetCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onDatasetChanged);
    connect(m_secondaryDatasetCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onSecondaryDatasetChanged);

    connect(m_primaryPlot, &DepthPlotWidget::cursorMoved,
            this, &MainWindow::onCursorMoved);
    connect(m_primaryPlot, &DepthPlotWidget::rangeChanged,
            this, &MainWindow::onRangeChanged);

    connect(m_secondaryPlot, &DepthPlotWidget::cursorMoved,
            this, &MainWindow::onCursorMoved);
}

// ---------------------------------------------------------------------------
// Data loading
// ---------------------------------------------------------------------------

void MainWindow::loadWorkingDirectory(const QString& workdir)
{
    m_workdir = workdir;
    setWindowTitle(QStringLiteral("ReadDepth++ - %1").arg(QDir(workdir).dirName()));
    autoDiscoverFiles();
}

void MainWindow::loadDataset(const QString& dbPath)
{
    if (!m_datasetPaths.contains(dbPath)) {
        m_datasetPaths.append(dbPath);
        populateDatasetCombos();
    }

    // If this is the first dataset, select it
    if (m_datasetPaths.size() == 1) {
        m_datasetCombo->setCurrentIndex(0);
    }

    statusBar()->showMessage(QStringLiteral("Loaded dataset: %1")
                                 .arg(QFileInfo(dbPath).fileName()), 3000);
}

void MainWindow::loadGeneReference(const QString& fastaPath)
{
    statusBar()->showMessage("Loading gene reference...");
    QApplication::processEvents();

    if (!m_annotationDB.loadFromFasta(fastaPath)) {
        QMessageBox::warning(this, "Error",
                             QStringLiteral("Failed to load gene reference:\n%1")
                                 .arg(m_annotationDB.lastError()));
        statusBar()->showMessage("Gene reference load failed", 3000);
        return;
    }

    QStringList genes = m_annotationDB.allGeneNames();
    m_geneSelector->setGeneList(genes);

    statusBar()->showMessage(QStringLiteral("Loaded %1 genes from reference")
                                 .arg(genes.size()), 3000);
}

void MainWindow::loadDESeq2Results(const QString& csvPath)
{
    m_deseq2Results = parseDESeq2CSV(csvPath);
    if (m_deseq2Results.isEmpty()) {
        QMessageBox::warning(this, "Error",
                             "Failed to parse DESeq2 results or file is empty.");
        return;
    }

    m_geneSelector->setDESeq2Results(m_deseq2Results);

    statusBar()->showMessage(QStringLiteral("Loaded %1 DESeq2 results")
                                 .arg(m_deseq2Results.size()), 3000);
}

void MainWindow::loadJunctionData(const QString& csvPath)
{
    m_junctionData = parseJunctionCSV(csvPath);
    if (m_junctionData.isEmpty()) {
        QMessageBox::warning(this, "Warning",
                             "No junction data found in CSV file.");
        return;
    }

    statusBar()->showMessage(QStringLiteral("Loaded junction data for %1 genes")
                                 .arg(m_junctionData.size()), 3000);

    // If currently displaying a gene, refresh to show junction overlay
    if (!m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

void MainWindow::displayGene(const QString& gene)
{
    m_geneSelector->selectGene(gene);
}

// ---------------------------------------------------------------------------
// Slots
// ---------------------------------------------------------------------------

void MainWindow::onGeneSelected(const QString& gene)
{
    m_currentGene = gene;
    refreshDisplay();
}

void MainWindow::onIntervalWidthChanged(int value)
{
    m_intervalWidth = value;
    if (!m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

void MainWindow::onIntervalSpacingChanged(int value)
{
    m_intervalSpacing = value;
    if (!m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

void MainWindow::onCompareToggled(bool checked)
{
    m_comparisonPanel->setVisible(checked);
    updateComparisonMode();
}

void MainWindow::onOverlayToggled(bool checked)
{
    if (checked) {
        updateComparisonMode();
    }
}

void MainWindow::onSplitToggled(bool checked)
{
    if (checked) {
        updateComparisonMode();
    }
}

void MainWindow::onDatasetChanged(int index)
{
    Q_UNUSED(index);
    if (!m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

void MainWindow::onSecondaryDatasetChanged(int index)
{
    Q_UNUSED(index);
    if (m_compareCheck->isChecked() && !m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

void MainWindow::onCursorMoved(int position, int depth)
{
    m_cursorLabel->setText(QStringLiteral("Position: %1 nt | Depth: %2 reads")
                              .arg(QLocale().toString(position))
                              .arg(QLocale().toString(depth)));
}

void MainWindow::onRangeChanged(qreal xMin, qreal xMax)
{
    // Update mRNA track widgets to match chart zoom
    m_primaryTrack->setVisibleRange(static_cast<int>(xMin), static_cast<int>(xMax));

    if (m_secondaryPanel->isVisible()) {
        m_secondaryTrack->setVisibleRange(static_cast<int>(xMin), static_cast<int>(xMax));
    }

    // Update zoom level display
    if (!m_primaryProfile.points.isEmpty()) {
        qreal fullRange = m_primaryProfile.annotation.mRNALength;
        qreal currentRange = xMax - xMin;
        if (currentRange > 0 && fullRange > 0) {
            int zoomPct = static_cast<int>(fullRange / currentRange * 100.0);
            m_zoomLabel->setText(QStringLiteral("Zoom: %1%").arg(zoomPct));
        }
    }

    // Sync if in split mode
    if (m_compareCheck->isChecked() && m_splitRadio->isChecked()) {
        m_syncController->setXRange(xMin, xMax);
    }
}

void MainWindow::onExportCSV()
{
    if (m_primaryProfile.points.isEmpty()) {
        QMessageBox::information(this, "Export", "No depth data to export. Load a gene first.");
        return;
    }

    QString defaultName = QStringLiteral("%1_depth.csv").arg(m_currentGene);
    QString filePath = QFileDialog::getSaveFileName(
        this, "Export Depth CSV",
        m_workdir.isEmpty() ? defaultName : m_workdir + "/" + defaultName,
        "CSV Files (*.csv);;All Files (*)");

    if (filePath.isEmpty()) return;

    if (deepn::ExportEngine::exportDepthCSV(filePath, m_primaryProfile)) {
        statusBar()->showMessage(QStringLiteral("Exported to %1").arg(filePath), 3000);
    } else {
        QMessageBox::warning(this, "Export Error",
                             QStringLiteral("Failed to export CSV:\n%1")
                                 .arg(deepn::ExportEngine::lastError()));
    }
}

void MainWindow::onExportFigure()
{
    if (m_primaryProfile.points.isEmpty()) {
        QMessageBox::information(this, "Export", "No depth data to export. Load a gene first.");
        return;
    }

    QString defaultName = QStringLiteral("%1_depth").arg(m_currentGene);
    QString filePath = QFileDialog::getSaveFileName(
        this, "Export Figure",
        m_workdir.isEmpty() ? defaultName : m_workdir + "/" + defaultName,
        "SVG (*.svg);;PDF (*.pdf);;PNG (*.png);;All Files (*)");

    if (filePath.isEmpty()) return;

    QSize exportSize(800, 400);
    bool ok = false;

    if (filePath.endsWith(".svg", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigureSVG(filePath, m_primaryPlot, exportSize);
    } else if (filePath.endsWith(".pdf", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigurePDF(filePath, m_primaryPlot, exportSize);
    } else if (filePath.endsWith(".png", Qt::CaseInsensitive)) {
        ok = deepn::ExportEngine::exportFigurePNG(filePath, m_primaryPlot, exportSize, 300);
    } else {
        // Default to PNG
        if (!filePath.endsWith(".png", Qt::CaseInsensitive)) {
            filePath += ".png";
        }
        ok = deepn::ExportEngine::exportFigurePNG(filePath, m_primaryPlot, exportSize, 300);
    }

    if (ok) {
        statusBar()->showMessage(QStringLiteral("Exported figure to %1").arg(filePath), 3000);
    } else {
        QMessageBox::warning(this, "Export Error",
                             QStringLiteral("Failed to export figure:\n%1")
                                 .arg(deepn::ExportEngine::lastError()));
    }
}

void MainWindow::onBatchRun()
{
    if (m_datasetPaths.isEmpty()) {
        QMessageBox::information(this, "Batch", "No datasets loaded. Load at least one dataset first.");
        return;
    }

    QStringList geneList;
    if (!m_deseq2Results.isEmpty()) {
        for (const auto& r : m_deseq2Results) {
            geneList.append(r.gene);
        }
    } else {
        geneList = m_annotationDB.allGeneNames();
    }

    if (geneList.isEmpty()) {
        QMessageBox::information(this, "Batch", "No gene list available for batch processing.");
        return;
    }

    // Ask for number of top genes
    bool ok = false;
    int topN = QInputDialog::getInt(this, "Batch Report",
                                     "Number of top genes to process:",
                                     50, 1, geneList.size(), 1, &ok);
    if (!ok) return;

    // Ask for output directory
    QString outputDir = QFileDialog::getExistingDirectory(
        this, "Select Output Directory for Batch Report",
        m_workdir.isEmpty() ? QDir::homePath() : m_workdir);
    if (outputDir.isEmpty()) return;

    // Configure and run batch
    deepn::BatchRunner::Config config;
    config.dbPaths = m_datasetPaths;
    config.geneList = geneList;
    config.outputDir = outputDir;
    config.topN = topN;
    config.exportCSV = true;
    config.exportFigures = false;

    auto* runner = new deepn::BatchRunner(this);
    auto* progressDlg = new QProgressDialog("Processing genes...", "Cancel", 0, topN, this);
    progressDlg->setWindowModality(Qt::WindowModal);
    progressDlg->setMinimumDuration(0);

    connect(runner, &deepn::BatchRunner::progressChanged,
            this, [progressDlg](int current, int total, const QString& gene) {
                progressDlg->setMaximum(total);
                progressDlg->setValue(current);
                progressDlg->setLabelText(QStringLiteral("Processing: %1 (%2/%3)")
                                              .arg(gene).arg(current).arg(total));
            });

    connect(runner, &deepn::BatchRunner::finished,
            this, [this, progressDlg, runner, outputDir](bool success) {
                progressDlg->close();
                if (success) {
                    QMessageBox::information(this, "Batch Complete",
                                             QStringLiteral("Batch report saved to:\n%1").arg(outputDir));
                } else {
                    QMessageBox::warning(this, "Batch Error", "Batch processing failed or was cancelled.");
                }
                runner->deleteLater();
                progressDlg->deleteLater();
            });

    connect(progressDlg, &QProgressDialog::canceled, runner, &deepn::BatchRunner::stop);

    runner->start(config);
}

void MainWindow::onResetZoom()
{
    m_primaryPlot->resetZoom();
    if (m_secondaryPanel->isVisible()) {
        m_secondaryPlot->resetZoom();
    }
}

// ---------------------------------------------------------------------------
// Core display logic
// ---------------------------------------------------------------------------

void MainWindow::refreshDisplay()
{
    if (m_currentGene.isEmpty()) return;

    int primaryIdx = m_datasetCombo->currentIndex();
    if (primaryIdx < 0 || primaryIdx >= m_datasetPaths.size()) {
        m_boundaryLabel->setText("No dataset selected.");
        return;
    }

    QString primaryDbPath = m_datasetPaths[primaryIdx];

    // Get gene annotation
    deepn::GeneAnnotation annotation = m_annotationDB.findByGeneName(m_currentGene);
    if (!annotation.isValid()) {
        // Try by refseq
        annotation = m_annotationDB.findByRefseq(m_currentGene);
    }
    if (!annotation.isValid()) {
        m_boundaryLabel->setText(QStringLiteral("Gene '%1' not found in reference database.")
                                    .arg(m_currentGene));
        return;
    }

    statusBar()->showMessage(QStringLiteral("Calculating depth profile for %1...")
                                 .arg(m_currentGene));
    QApplication::processEvents();

    // Calculate primary depth profile
    m_primaryProfile = m_depthCalc.calculate(
        primaryDbPath, m_currentGene, annotation, m_intervalWidth, m_intervalSpacing);

    if (m_primaryProfile.points.isEmpty()) {
        m_boundaryLabel->setText(QStringLiteral("No depth data for gene '%1' in dataset '%2'.")
                                    .arg(m_currentGene,
                                         QFileInfo(primaryDbPath).fileName()));
        m_primaryPlot->setProfile(m_primaryProfile);
        statusBar()->showMessage("No data found.", 3000);
        return;
    }

    m_primaryProfile.datasetLabel = QFileInfo(primaryDbPath).baseName();

    // Detect boundary
    m_primaryBoundary = deepn::DepthCalculator::detect3PrimeBoundary(m_primaryProfile);

    // Set primary plot
    m_primaryPlot->setProfile(m_primaryProfile);
    m_primaryPlot->setBoundary(m_primaryBoundary);
    m_primaryTrack->setAnnotation(annotation);

    // Junction overlay
    if (m_junctionData.contains(m_currentGene)) {
        m_primaryPlot->setJunctionOverlay(m_junctionData[m_currentGene]);
    } else {
        m_primaryPlot->clearJunctionOverlay();
    }

    // Handle comparison mode
    bool comparing = m_compareCheck->isChecked();
    if (comparing) {
        int secIdx = m_secondaryDatasetCombo->currentIndex();
        if (secIdx >= 0 && secIdx < m_datasetPaths.size()) {
            QString secDbPath = m_datasetPaths[secIdx];
            m_secondaryProfile = m_depthCalc.calculate(
                secDbPath, m_currentGene, annotation, m_intervalWidth, m_intervalSpacing);
            m_secondaryProfile.datasetLabel = QFileInfo(secDbPath).baseName();
            m_secondaryBoundary = deepn::DepthCalculator::detect3PrimeBoundary(m_secondaryProfile);

            if (m_overlayRadio->isChecked()) {
                // Overlay mode: show both on primary plot
                m_primaryPlot->setOverlayProfile(m_secondaryProfile);
                m_secondaryPanel->setVisible(false);
            } else {
                // Split mode: show in secondary panel
                m_primaryPlot->clearOverlay();
                m_secondaryPlot->setProfile(m_secondaryProfile);
                m_secondaryPlot->setBoundary(m_secondaryBoundary);
                m_secondaryTrack->setAnnotation(annotation);

                if (m_junctionData.contains(m_currentGene)) {
                    m_secondaryPlot->setJunctionOverlay(m_junctionData[m_currentGene]);
                } else {
                    m_secondaryPlot->clearJunctionOverlay();
                }

                m_secondaryPanel->setVisible(true);

                // Set up sync controller
                m_syncController->clear();
                m_syncController->addChart(m_primaryPlot->depthChart());
                m_syncController->addChart(m_secondaryPlot->depthChart());
            }
        }
    } else {
        m_primaryPlot->clearOverlay();
        m_secondaryPanel->setVisible(false);
        m_syncController->clear();
    }

    // Update boundary info
    updateBoundaryInfo();

    statusBar()->showMessage(QStringLiteral("Displaying %1 (%2 intervals, %3 total reads)")
                                 .arg(m_currentGene)
                                 .arg(m_primaryProfile.points.size())
                                 .arg(m_primaryProfile.totalReads), 5000);
}

void MainWindow::updateBoundaryInfo()
{
    // Find 5' junction from junction data (if available)
    int fivePrimeJunction = 0;
    if (m_junctionData.contains(m_currentGene)) {
        const auto& junctions = m_junctionData[m_currentGene];
        // Find the most prominent in-frame in-ORF junction
        double bestPpm = 0.0;
        for (const auto& jnc : junctions) {
            if (jnc.isInFrame() && jnc.cdsClass == "in_orf" && jnc.ppm > bestPpm) {
                bestPpm = jnc.ppm;
                fivePrimeJunction = jnc.position;
            }
        }
    }

    // Analyze primary boundary
    auto detailedResult = BoundaryDetector::analyze(m_primaryProfile, fivePrimeJunction);
    m_boundaryLabel->setText(detailedResult.summary);

    if (detailedResult.hasInsertExtent) {
        m_insertExtentLabel->setText(
            BoundaryDetector::formatInsertExtent(
                detailedResult.insertExtent, m_primaryProfile.annotation));
        m_insertExtentLabel->setVisible(true);
    } else {
        m_insertExtentLabel->setVisible(false);
    }
}

void MainWindow::updateComparisonMode()
{
    bool comparing = m_compareCheck->isChecked();

    if (!comparing) {
        m_secondaryPanel->setVisible(false);
        m_primaryPlot->clearOverlay();
        m_syncController->clear();
    }

    if (!m_currentGene.isEmpty()) {
        refreshDisplay();
    }
}

// ---------------------------------------------------------------------------
// File discovery
// ---------------------------------------------------------------------------

void MainWindow::autoDiscoverFiles()
{
    if (m_workdir.isEmpty()) return;

    QDir dir(m_workdir);

    // Discover SQLite databases in analyzed_files/
    QDir analyzedDir(dir.filePath("analyzed_files"));
    if (analyzedDir.exists()) {
        QStringList sqliteFiles = analyzedDir.entryList({"*.sqlite", "*.db"}, QDir::Files);
        for (const QString& f : sqliteFiles) {
            QString fullPath = analyzedDir.filePath(f);
            if (!m_datasetPaths.contains(fullPath)) {
                m_datasetPaths.append(fullPath);
            }
        }
    }

    // Also check for sqlite files directly in workdir
    QStringList rootSqlite = dir.entryList({"*.sqlite", "*.db"}, QDir::Files);
    for (const QString& f : rootSqlite) {
        QString fullPath = dir.filePath(f);
        if (!m_datasetPaths.contains(fullPath)) {
            m_datasetPaths.append(fullPath);
        }
    }

    populateDatasetCombos();

    // Discover gene reference FASTA files in data/ or gene_dictionary/
    QStringList refDirs = {"data", "gene_dictionary", "data/gene_dictionary"};
    for (const QString& refDirName : refDirs) {
        QDir refDir(dir.filePath(refDirName));
        if (refDir.exists()) {
            QStringList fastaFiles = refDir.entryList({"*.fa", "*.fasta", "*.fna"}, QDir::Files);
            if (!fastaFiles.isEmpty()) {
                loadGeneReference(refDir.filePath(fastaFiles.first()));
                break;
            }
        }
    }

    // Discover DESeq2 results
    QStringList resultsDirs = {".", "deseq2_results", "gene_count_summary"};
    for (const QString& resDirName : resultsDirs) {
        QDir resDir(dir.filePath(resDirName));
        if (resDir.exists()) {
            QStringList csvFiles = resDir.entryList({"*deseq2*.csv", "*results*.csv"}, QDir::Files);
            if (!csvFiles.isEmpty()) {
                loadDESeq2Results(resDir.filePath(csvFiles.first()));
                break;
            }
        }
    }

    // Discover junction data from MultiQuery++
    QStringList junctionDirs = {".", "multiquery_results"};
    for (const QString& jncDirName : junctionDirs) {
        QDir jncDir(dir.filePath(jncDirName));
        if (jncDir.exists()) {
            QStringList jncFiles = jncDir.entryList({"*junction*.csv", "*collapsed*.csv"}, QDir::Files);
            if (!jncFiles.isEmpty()) {
                loadJunctionData(jncDir.filePath(jncFiles.first()));
                break;
            }
        }
    }
}

void MainWindow::populateDatasetCombos()
{
    // Block signals while repopulating
    m_datasetCombo->blockSignals(true);
    m_secondaryDatasetCombo->blockSignals(true);

    int prevPrimary = m_datasetCombo->currentIndex();
    int prevSecondary = m_secondaryDatasetCombo->currentIndex();

    m_datasetCombo->clear();
    m_secondaryDatasetCombo->clear();

    for (const QString& path : m_datasetPaths) {
        QString label = QFileInfo(path).fileName();
        m_datasetCombo->addItem(label, path);
        m_secondaryDatasetCombo->addItem(label, path);
    }

    // Restore selection
    if (prevPrimary >= 0 && prevPrimary < m_datasetCombo->count()) {
        m_datasetCombo->setCurrentIndex(prevPrimary);
    } else if (m_datasetCombo->count() > 0) {
        m_datasetCombo->setCurrentIndex(0);
    }

    if (prevSecondary >= 0 && prevSecondary < m_secondaryDatasetCombo->count()) {
        m_secondaryDatasetCombo->setCurrentIndex(prevSecondary);
    } else if (m_secondaryDatasetCombo->count() > 1) {
        m_secondaryDatasetCombo->setCurrentIndex(1);
    } else if (m_secondaryDatasetCombo->count() > 0) {
        m_secondaryDatasetCombo->setCurrentIndex(0);
    }

    m_datasetCombo->blockSignals(false);
    m_secondaryDatasetCombo->blockSignals(false);
}

// ---------------------------------------------------------------------------
// CSV Parsing
// ---------------------------------------------------------------------------

QVector<deepn::DESeq2Result> MainWindow::parseDESeq2CSV(const QString& csvPath)
{
    QVector<deepn::DESeq2Result> results;

    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "Cannot open DESeq2 CSV:" << csvPath;
        return results;
    }

    QTextStream in(&file);
    QString headerLine = in.readLine();
    if (headerLine.isEmpty()) {
        return results;
    }

    // Parse header to find column indices
    QStringList headers = headerLine.split(',');
    auto findCol = [&headers](const QStringList& names) -> int {
        for (const QString& name : names) {
            for (int i = 0; i < headers.size(); ++i) {
                if (headers[i].trimmed().compare(name, Qt::CaseInsensitive) == 0 ||
                    headers[i].trimmed().replace('"', "").compare(name, Qt::CaseInsensitive) == 0) {
                    return i;
                }
            }
        }
        return -1;
    };

    int geneCol = findCol({"gene", "Gene", "gene_name", "geneName", ""});
    int baseMeanCol = findCol({"baseMean", "basemean", "base_mean"});
    int l2fcCol = findCol({"log2FoldChange", "log2foldchange", "l2fc", "log2FC"});
    int lfcSECol = findCol({"lfcSE", "lfcse", "lfc_se"});
    int statCol = findCol({"stat", "Stat"});
    int pvalCol = findCol({"pvalue", "pval", "p_value"});
    int padjCol = findCol({"padj", "p_adj", "adjusted_pvalue"});
    int enrichCol = findCol({"enrichment", "Enrichment", "direction"});

    // If gene column is -1, try the first column (common for R output where
    // the gene names are row names)
    if (geneCol < 0) geneCol = 0;

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        // Handle quoted CSV fields
        QStringList fields;
        bool inQuotes = false;
        QString current;
        for (QChar ch : line) {
            if (ch == '"') {
                inQuotes = !inQuotes;
            } else if (ch == ',' && !inQuotes) {
                fields.append(current.trimmed());
                current.clear();
            } else {
                current += ch;
            }
        }
        fields.append(current.trimmed());

        deepn::DESeq2Result r;
        auto safeField = [&fields](int idx) -> QString {
            if (idx >= 0 && idx < fields.size()) {
                return fields[idx].replace('"', "").trimmed();
            }
            return {};
        };

        r.gene = safeField(geneCol);
        if (r.gene.isEmpty()) continue;

        if (baseMeanCol >= 0) r.baseMean = safeField(baseMeanCol).toDouble();
        if (l2fcCol >= 0) r.log2FoldChange = safeField(l2fcCol).toDouble();
        if (lfcSECol >= 0) r.lfcSE = safeField(lfcSECol).toDouble();
        if (statCol >= 0) r.stat = safeField(statCol).toDouble();
        if (pvalCol >= 0) r.pvalue = safeField(pvalCol).toDouble();
        if (padjCol >= 0) r.padj = safeField(padjCol).toDouble();
        if (enrichCol >= 0) r.enrichment = safeField(enrichCol);

        results.append(r);
    }

    return results;
}

QMap<QString, QVector<deepn::JunctionSite>> MainWindow::parseJunctionCSV(const QString& csvPath)
{
    QMap<QString, QVector<deepn::JunctionSite>> result;

    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "Cannot open junction CSV:" << csvPath;
        return result;
    }

    QTextStream in(&file);
    QString headerLine = in.readLine();
    if (headerLine.isEmpty()) {
        return result;
    }

    QStringList headers = headerLine.split(',');
    auto findCol = [&headers](const QStringList& names) -> int {
        for (const QString& name : names) {
            for (int i = 0; i < headers.size(); ++i) {
                if (headers[i].trimmed().compare(name, Qt::CaseInsensitive) == 0 ||
                    headers[i].trimmed().replace('"', "").compare(name, Qt::CaseInsensitive) == 0) {
                    return i;
                }
            }
        }
        return -1;
    };

    int geneCol = findCol({"Gene", "gene", "gene_name"});
    int posCol = findCol({"Position", "position", "rstart"});
    int posEndCol = findCol({"PositionEnd", "positionEnd", "rend"});
    int ppmCol = findCol({"PPM", "ppm", "TotalPPM"});
    int frameCol = findCol({"Frame", "frame", "DominantFrame"});
    int cdsCol = findCol({"CDS_Class", "cds_class", "location"});
    int refseqCol = findCol({"RefSeq", "refseq"});

    if (geneCol < 0 || posCol < 0) {
        qDebug() << "Junction CSV missing Gene or Position columns";
        return result;
    }

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        QStringList fields = line.split(',');
        auto safeField = [&fields](int idx) -> QString {
            if (idx >= 0 && idx < fields.size()) {
                return fields[idx].trimmed().replace('"', "");
            }
            return {};
        };

        QString gene = safeField(geneCol);
        if (gene.isEmpty()) continue;

        deepn::JunctionSite site;
        site.geneName = gene;
        site.position = safeField(posCol).toInt();
        if (posEndCol >= 0) site.positionEnd = safeField(posEndCol).toInt();
        if (ppmCol >= 0) site.ppm = safeField(ppmCol).toDouble();
        if (frameCol >= 0) site.frame = safeField(frameCol);
        if (cdsCol >= 0) site.cdsClass = safeField(cdsCol);
        if (refseqCol >= 0) site.refseq = safeField(refseqCol);

        result[gene].append(site);
    }

    return result;
}

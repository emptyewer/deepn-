#include "main_window.h"
#include "ui_mainwindow.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QFileDialog>
#include <QMessageBox>
#include <QHeaderView>
#include <QStandardPaths>
#include <QDateTime>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QSettings>
#include <QDir>
#include <QProgressBar>
#include <QThread>
#include <QMutex>
#include <QMutexLocker>
#include <stdexcept>

namespace deseq2
{

    // AnalysisWorker implementation
    AnalysisWorker::AnalysisWorker(const CombinedData &data, double pValueThreshold)
        : m_inputData(data), m_pValueThreshold(pValueThreshold), m_shouldStop(false)
    {
    }

    void AnalysisWorker::runAnalysis()
    {
        try
        {
            emit progressChanged(0, "Starting DESeq2 analysis...");

            // Convert Qt data to Eigen format
            auto [counts, metadata] = convertToEigenFormat(m_inputData);

            emit progressChanged(10, "Data conversion completed");

            // Check if we should stop
            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 1: Create DeseqDataSet
            emit progressChanged(15, "Creating DeseqDataSet...");
            DeseqDataSet dds(counts, metadata, "~condition", true);

            // Step 2: Fit size factors
            emit progressChanged(25, "Fitting size factors...");
            dds.fitSizeFactors();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                    return;
            }

            // Step 3: Fit dispersions
            emit progressChanged(35, "Fitting genewise dispersions...");
            dds.fitGenewiseDispersions();

            emit progressChanged(45, "Fitting dispersion trend...");
            dds.fitDispersionTrend();

            emit progressChanged(55, "Fitting dispersion prior...");
            dds.fitDispersionPrior();

            emit progressChanged(65, "Fitting MAP dispersions...");
            dds.fitMAPDispersions();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 4: Fit log fold changes
            emit progressChanged(70, "Fitting log fold changes...");
            dds.fitLFC();

            // Step 5: Handle outliers
            emit progressChanged(75, "Calculating Cooks distances...");
            dds.calculateCooks();
            dds.refit();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Step 6: Statistical testing
            emit progressChanged(80, "Running statistical tests...");

            // Create contrast for the analysis (assuming first condition is reference)
            Eigen::VectorXd contrast(2);
            contrast << 0.0, 1.0; // Test second condition vs first

            DeseqStats ds(dds, contrast, m_pValueThreshold, true, true);

            emit progressChanged(85, "Running Wald test...");
            ds.runWaldTest();

            emit progressChanged(90, "Applying Cook's filtering...");
            ds.cooksFiltering();

            emit progressChanged(95, "Applying independent filtering...");
            ds.independentFiltering();

            // Step 7: Generate results
            emit progressChanged(98, "Generating results...");
            Eigen::MatrixXd results = ds.summary();

            {
                QMutexLocker locker(&m_stopMutex);
                if (m_shouldStop)
                {
                    emit debugMessage("Analysis stopped by user");
                    emit finished();
                    return;
                }
            }

            // Prepare results structure
            AnalysisResults analysisResults;
            analysisResults.results = results;

            // Store dispersion data for visualization
            analysisResults.genewiseDispersions = dds.getGenewiseDispersions();
            analysisResults.fittedDispersions = dds.getFittedDispersions();
            analysisResults.baseMeans = dds.getNormedMeans();

            // Convert gene names and sample names
            analysisResults.geneNames.reserve(m_inputData.geneNames.size());
            for (const auto &gene : m_inputData.geneNames)
            {
                analysisResults.geneNames.push_back(gene.toStdString());
            }

            analysisResults.sampleNames.reserve(m_inputData.sampleNames.size());
            for (const auto &sample : m_inputData.sampleNames)
            {
                analysisResults.sampleNames.push_back(sample.toStdString());
            }

            // Calculate summary statistics
            analysisResults.totalGenes = results.rows();
            analysisResults.significantGenes = 0;
            analysisResults.upregulatedGenes = 0;
            analysisResults.downregulatedGenes = 0;
            analysisResults.pValueThreshold = m_pValueThreshold;
            analysisResults.log2FCThreshold = 0.0; // Default, can be customized later

            for (int i = 0; i < results.rows(); ++i)
            {
                double padj = results(i, 5); // Adjusted p-value column (FIXED: was 4, should be 5)
                double lfc = results(i, 1);  // Log2 fold change column

                if (padj < m_pValueThreshold)
                {
                    analysisResults.significantGenes++;
                    if (lfc > 0)
                    {
                        analysisResults.upregulatedGenes++;
                    }
                    else if (lfc < 0)
                    {
                        analysisResults.downregulatedGenes++;
                    }
                }
            }

            analysisResults.isValid = true;
            analysisResults.errorMessage = "";

            emit progressChanged(100, "Analysis completed successfully!");
            emit debugMessage("About to emit analysisFinished signal");
            emit debugMessage("analysisResults.isValid = " + QString(analysisResults.isValid ? "true" : "false"));
            emit debugMessage("analysisResults.totalGenes = " + QString::number(analysisResults.totalGenes));
            emit debugMessage("analysisResults.significantGenes = " + QString::number(analysisResults.significantGenes));
            emit analysisFinished(&analysisResults);
            emit debugMessage("analysisFinished signal emitted");
            emit finished(); // Signal that worker is completely done
        }
        catch (const std::exception &e)
        {
            emit analysisError(QString("Analysis failed: %1").arg(e.what()));
            emit finished(); // Signal that worker is done even on error
        }
        catch (...)
        {
            emit analysisError("Analysis failed: Unknown error occurred");
            emit finished(); // Signal that worker is done even on error
        }
    }

    void AnalysisWorker::stopAnalysis()
    {
        QMutexLocker locker(&m_stopMutex);
        m_shouldStop = true;
    }

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> AnalysisWorker::convertToEigenFormat(const CombinedData &data)
    {
        // Convert count matrix (Qt format: genes x samples) to Eigen format (samples x genes)
        int numSamples = data.sampleNames.size();
        int numGenes = data.geneNames.size();

        Eigen::MatrixXd counts(numSamples, numGenes);

        for (int gene = 0; gene < numGenes; ++gene)
        {
            for (int sample = 0; sample < numSamples; ++sample)
            {
                // Use raw integer counts directly (from xlsx or back-computed from PPM)
                int integerCount = data.rawCountMatrix[gene][sample];
                counts(sample, gene) = static_cast<double>(integerCount);
            }
        }

        // Create metadata matrix based on group names
        Eigen::MatrixXd metadata(numSamples, 1);

        // Get unique group names to assign numeric codes
        QSet<QString> uniqueGroups;
        for (const auto &group : data.groupNames)
        {
            uniqueGroups.insert(group);
        }

        QList<QString> groupList = uniqueGroups.values();
        std::sort(groupList.begin(), groupList.end());

        // Assign numeric codes to groups
        for (int sample = 0; sample < numSamples; ++sample)
        {
            QString groupName = data.groupNames[sample];
            int groupIndex = groupList.indexOf(groupName);
            metadata(sample, 0) = static_cast<double>(groupIndex);
        }

        return std::make_tuple(counts, metadata);
    }

    MainWindow::MainWindow(QWidget *parent)
        : QMainWindow(parent),
          m_geneCountHandler(std::make_unique<GeneCountHandler>()),
          m_combinedData(),
          // Initialize all pointer members to nullptr
          m_tabWidget(nullptr),
          m_inputTab(nullptr),
          m_resultsTab(nullptr),
          m_visualizationTab(nullptr),
          m_analysisTab(nullptr),
          m_addPpmFilesButton(nullptr),
          m_clearFilesButton(nullptr),
          m_generateDataButton(nullptr),
          m_exportGeneratedDataButton(nullptr),
          m_geneCountFilesTable(nullptr),
          m_groupNameEdit(nullptr),
          m_assignGroupButton(nullptr),
          m_autoAssignButton(nullptr),
          m_countsPreviewTable(nullptr),
          m_metadataPreviewTable(nullptr),
          m_runAnalysisButton(nullptr),
          m_stopAnalysisButton(nullptr),
          m_resetAnalysisButton(nullptr),
          m_analysisProgressBar(nullptr),
          m_totalGenesLabel(nullptr),
          m_significantGenesLabel(nullptr),
          m_upregulatedLabel(nullptr),
          m_downregulatedLabel(nullptr),
          m_resultsTable(nullptr),
          m_pValueThresholdSpinBox(nullptr),
          m_log2FCThresholdSpinBox(nullptr),
          m_applyFilterButton(nullptr),
          m_clearFilterButton(nullptr),
          m_exportResultsButton(nullptr),
          m_saveGeneListsButton(nullptr),
          m_plotTypeCombo(nullptr),
          m_colorSchemeCombo(nullptr),
          m_pointSizeSpinBox(nullptr),
          m_plotWidget(nullptr),
          m_savePlotButton(nullptr),
          m_exportDataButton(nullptr),
          m_printPlotButton(nullptr),
          m_progressOutputText(nullptr),
          m_clearOutputButton(nullptr),
          m_saveOutputButton(nullptr),
          m_actionOpenCounts(nullptr),
          m_actionOpenMetadata(nullptr),
          m_actionSaveResults(nullptr),
          m_actionExportVisualizations(nullptr),
          m_actionExit(nullptr),
          m_actionRunAnalysis(nullptr),
          m_actionStopAnalysis(nullptr),
          m_actionResetAnalysis(nullptr),
          m_actionInputTab(nullptr),
          m_actionAnalysisTab(nullptr),
          m_actionResultsTab(nullptr),
          m_actionVisualizationTab(nullptr),
          m_actionAbout(nullptr),
          m_actionUserGuide(nullptr)
    {
        // Initialize last used directory to Documents folder
        m_lastUsedDirectory = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
        // Load saved last used directory from settings
        loadLastUsedDirectory();
        // Setup UI from the UI file
        m_ui = std::make_unique<Ui::MainWindow>();
        m_ui->setupUi(this);
        // Connect UI elements to member variables
        connectUiElements();
        setupConnections();
        setupMenusAndToolbar();
        initializeTables();
        updateUiState();
        addProgressMessage("DESeq2 GUI Application started successfully.");

        // Initialize analysis results
        m_analysisResults.isValid = false;
    }

    MainWindow::~MainWindow()
    {
        // Stop analysis if running
        if (m_analysisThread && m_analysisThread->isRunning())
        {
            if (m_analysisWorker)
            {
                m_analysisWorker->stopAnalysis();
            }
            m_analysisThread->quit();
            m_analysisThread->wait(5000); // Wait up to 5 seconds
        }
    }

    void MainWindow::connectUiElements()
    {
        // Connect UI elements from the UI file to member variables
        m_tabWidget = m_ui->tabWidget;
        m_inputTab = m_ui->inputTab;
        m_resultsTab = m_ui->resultsTab;
        m_visualizationTab = m_ui->visualizationTab;
        m_analysisTab = nullptr; // Not present in UI file
        // Input tab elements
        m_addPpmFilesButton = m_ui->addPpmFilesButton;
        m_clearFilesButton = m_ui->clearFilesButton;
        m_generateDataButton = m_ui->generateDataButton;
        m_exportGeneratedDataButton = m_ui->exportGeneratedDataButton;
        m_geneCountFilesTable = m_ui->geneCountFilesTable;
        m_groupNameEdit = m_ui->groupNameEdit;
        m_assignGroupButton = m_ui->assignGroupButton;
        m_autoAssignButton = m_ui->autoAssignButton;
        m_countsPreviewTable = m_ui->countsPreviewTable;
        m_metadataPreviewTable = m_ui->metadataPreviewTable;

        // Analysis tab elements - connect to UI elements if they exist
        m_runAnalysisButton = m_ui->runAnalysisButton;
        m_stopAnalysisButton = m_ui->stopAnalysisButton;
        m_resetAnalysisButton = m_ui->resetAnalysisButton;

        // Add progress bar if not already in UI
        if (!m_analysisProgressBar)
        {
            m_analysisProgressBar = new QProgressBar(this);
            m_analysisProgressBar->setVisible(false);
            statusBar()->addPermanentWidget(m_analysisProgressBar);
        }

        // Results tab elements
        m_totalGenesLabel = m_ui->totalGenesLabel;
        m_significantGenesLabel = m_ui->significantGenesLabel;
        m_upregulatedLabel = m_ui->upregulatedLabel;
        m_downregulatedLabel = m_ui->downregulatedLabel;
        m_resultsTable = m_ui->resultsTable;
        m_pValueThresholdSpinBox = m_ui->pValueThresholdSpinBox;
        m_log2FCThresholdSpinBox = m_ui->log2FCThresholdSpinBox;
        m_applyFilterButton = m_ui->applyFilterButton;
        m_clearFilterButton = m_ui->clearFilterButton;
        m_exportResultsButton = m_ui->exportResultsButton;
        m_saveGeneListsButton = m_ui->saveGeneListsButton;
        // Visualization tab elements
        m_plotTypeCombo = m_ui->plotTypeCombo;
        m_colorSchemeCombo = m_ui->colorSchemeCombo;
        m_pointSizeSpinBox = m_ui->pointSizeSpinBox;
        // Replace the plain QWidget from UI with a QCustomPlot
        {
            QWidget *placeholder = m_ui->plotWidget;
            m_plotWidget = new QCustomPlot(placeholder->parentWidget());
            m_plotWidget->setMinimumSize(placeholder->minimumSize());
            m_plotWidget->setSizePolicy(placeholder->sizePolicy());

            // Replace in the layout
            QLayout *parentLayout = placeholder->parentWidget()->layout();
            if (parentLayout)
            {
                QLayoutItem *item = parentLayout->replaceWidget(placeholder, m_plotWidget);
                delete item;
            }
            placeholder->hide();
            placeholder->deleteLater();
        }
        m_savePlotButton = m_ui->savePlotButton;
        m_exportDataButton = m_ui->exportDataButton;
        m_printPlotButton = m_ui->printPlotButton;
        // Progress output elements
        m_progressOutputText = m_ui->progressOutputText;
        m_clearOutputButton = m_ui->clearOutputButton;
        m_saveOutputButton = m_ui->saveOutputButton;
    }

    void MainWindow::setupConnections()
    {
        // File operations
        if (m_addPpmFilesButton)
        {
            connect(m_addPpmFilesButton, &QPushButton::clicked, this, &MainWindow::onAddGeneCountFiles);
        }
        if (m_clearFilesButton)
        {
            connect(m_clearFilesButton, &QPushButton::clicked, this, &MainWindow::onClearFiles);
        }
        if (m_generateDataButton)
        {
            connect(m_generateDataButton, &QPushButton::clicked, this, &MainWindow::onGenerateData);
        }
        if (m_exportGeneratedDataButton)
        {
            connect(m_exportGeneratedDataButton, &QPushButton::clicked, this, &MainWindow::onExportGeneratedData);
        }

        // Group assignment
        if (m_assignGroupButton)
        {
            connect(m_assignGroupButton, &QPushButton::clicked, this, &MainWindow::onAssignGroup);
        }
        if (m_autoAssignButton)
        {
            connect(m_autoAssignButton, &QPushButton::clicked, this, &MainWindow::onAutoAssignGroups);
        }

        // Analysis control
        if (m_runAnalysisButton)
        {
            connect(m_runAnalysisButton, &QPushButton::clicked, this, &MainWindow::onRunAnalysis);
        }
        if (m_stopAnalysisButton)
        {
            connect(m_stopAnalysisButton, &QPushButton::clicked, this, &MainWindow::onStopAnalysis);
        }
        if (m_resetAnalysisButton)
        {
            connect(m_resetAnalysisButton, &QPushButton::clicked, this, &MainWindow::onResetAnalysis);
        }

        // Results filtering
        if (m_applyFilterButton)
        {
            connect(m_applyFilterButton, &QPushButton::clicked, this, &MainWindow::onApplyFilter);
        }
        if (m_clearFilterButton)
        {
            connect(m_clearFilterButton, &QPushButton::clicked, this, &MainWindow::onClearFilter);
        }

        // Export operations
        if (m_exportResultsButton)
        {
            connect(m_exportResultsButton, &QPushButton::clicked, this, &MainWindow::onExportResults);
        }
        if (m_saveGeneListsButton)
        {
            connect(m_saveGeneListsButton, &QPushButton::clicked, this, &MainWindow::onSaveGeneLists);
        }

        // Visualization
        if (m_plotTypeCombo)
        {
            connect(m_plotTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                    this, &MainWindow::onPlotTypeChanged);
        }
        if (m_savePlotButton)
        {
            connect(m_savePlotButton, &QPushButton::clicked, this, &MainWindow::onSavePlot);
        }
        if (m_exportDataButton)
        {
            connect(m_exportDataButton, &QPushButton::clicked, this, &MainWindow::onExportPlotData);
        }
        if (m_printPlotButton)
        {
            connect(m_printPlotButton, &QPushButton::clicked, this, &MainWindow::onPrintPlot);
        }

        // Progress output
        if (m_clearOutputButton)
        {
            connect(m_clearOutputButton, &QPushButton::clicked, this, &MainWindow::onClearOutput);
        }
        if (m_saveOutputButton)
        {
            connect(m_saveOutputButton, &QPushButton::clicked, this, &MainWindow::onSaveOutput);
        }
    }

    void MainWindow::setupMenusAndToolbar()
    {
        // Create menu bar
        QMenuBar *menuBar = this->menuBar();

        // File menu
        QMenu *fileMenu = menuBar->addMenu("&File");
        m_actionOpenCounts = fileMenu->addAction("Add Gene Count Files");
        m_actionOpenCounts->setShortcut(QKeySequence::Open);
        connect(m_actionOpenCounts, &QAction::triggered, this, &MainWindow::onOpenCounts);

        m_actionOpenMetadata = fileMenu->addAction("Generate Data");
        m_actionOpenMetadata->setShortcut(QKeySequence("Ctrl+G"));
        connect(m_actionOpenMetadata, &QAction::triggered, this, &MainWindow::onOpenMetadata);

        fileMenu->addSeparator();

        m_actionSaveResults = fileMenu->addAction("Save Results");
        m_actionSaveResults->setShortcut(QKeySequence::Save);
        connect(m_actionSaveResults, &QAction::triggered, this, &MainWindow::onSaveResults);

        m_actionExportVisualizations = fileMenu->addAction("Export Visualizations");
        m_actionExportVisualizations->setShortcut(QKeySequence("Ctrl+E"));
        connect(m_actionExportVisualizations, &QAction::triggered, this, &MainWindow::onExportVisualizations);

        fileMenu->addSeparator();

        m_actionExit = fileMenu->addAction("Exit");
        m_actionExit->setShortcut(QKeySequence::Quit);
        connect(m_actionExit, &QAction::triggered, this, &MainWindow::onExit);

        // Analysis menu
        QMenu *analysisMenu = menuBar->addMenu("&Analysis");
        m_actionRunAnalysis = analysisMenu->addAction("Run Analysis");
        m_actionRunAnalysis->setShortcut(QKeySequence("F5"));
        connect(m_actionRunAnalysis, &QAction::triggered, this, &MainWindow::onRunAnalysisAction);

        m_actionStopAnalysis = analysisMenu->addAction("Stop Analysis");
        m_actionStopAnalysis->setShortcut(QKeySequence("F6"));
        connect(m_actionStopAnalysis, &QAction::triggered, this, &MainWindow::onStopAnalysisAction);

        m_actionResetAnalysis = analysisMenu->addAction("Reset Analysis");
        m_actionResetAnalysis->setShortcut(QKeySequence("F7"));
        connect(m_actionResetAnalysis, &QAction::triggered, this, &MainWindow::onResetAnalysisAction);

        // View menu
        QMenu *viewMenu = menuBar->addMenu("&View");
        m_actionInputTab = viewMenu->addAction("Input Tab");
        m_actionInputTab->setShortcut(QKeySequence("Ctrl+1"));
        connect(m_actionInputTab, &QAction::triggered, this, &MainWindow::onInputTab);

        m_actionAnalysisTab = viewMenu->addAction("Analysis Tab");
        m_actionAnalysisTab->setShortcut(QKeySequence("Ctrl+2"));
        connect(m_actionAnalysisTab, &QAction::triggered, this, &MainWindow::onAnalysisTab);

        m_actionResultsTab = viewMenu->addAction("Results Tab");
        m_actionResultsTab->setShortcut(QKeySequence("Ctrl+3"));
        connect(m_actionResultsTab, &QAction::triggered, this, &MainWindow::onResultsTab);

        m_actionVisualizationTab = viewMenu->addAction("Visualization Tab");
        m_actionVisualizationTab->setShortcut(QKeySequence("Ctrl+4"));
        connect(m_actionVisualizationTab, &QAction::triggered, this, &MainWindow::onVisualizationTab);

        // Help menu
        QMenu *helpMenu = menuBar->addMenu("&Help");
        m_actionAbout = helpMenu->addAction("About");
        connect(m_actionAbout, &QAction::triggered, this, &MainWindow::onAbout);

        m_actionUserGuide = helpMenu->addAction("User Guide");
        m_actionUserGuide->setShortcut(QKeySequence::HelpContents);
        connect(m_actionUserGuide, &QAction::triggered, this, &MainWindow::onUserGuide);

        // Create toolbar
        QToolBar *toolBar = addToolBar("Main Toolbar");
        toolBar->addAction(m_actionOpenCounts);
        toolBar->addSeparator();
        toolBar->addAction(m_actionRunAnalysis);
        toolBar->addAction(m_actionStopAnalysis);
        toolBar->addSeparator();
        toolBar->addAction(m_actionSaveResults);
        toolBar->addAction(m_actionExportVisualizations);
    }

    void MainWindow::initializeTables()
    {
        // Initialize gene count files table
        if (m_geneCountFilesTable)
        {
            m_geneCountFilesTable->setColumnCount(5);
            QStringList headers = {"File Path", "Sample Name", "Group", "Total Reads", "Status"};
            m_geneCountFilesTable->setHorizontalHeaderLabels(headers);
            m_geneCountFilesTable->horizontalHeader()->setStretchLastSection(true);
            m_geneCountFilesTable->setAlternatingRowColors(true);
            m_geneCountFilesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_geneCountFilesTable->setSortingEnabled(true);
            m_geneCountFilesTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
        }

        // Initialize counts preview table
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->setAlternatingRowColors(true);
            m_countsPreviewTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_countsPreviewTable->setSortingEnabled(true);
        }

        // Initialize metadata preview table
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->setAlternatingRowColors(true);
            m_metadataPreviewTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_metadataPreviewTable->setSortingEnabled(true);
        }

        // Initialize results table
        if (m_resultsTable)
        {
            m_resultsTable->setAlternatingRowColors(true);
            m_resultsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
            m_resultsTable->setSortingEnabled(true);
        }
    }

    void MainWindow::updateUiState()
    {
        // Update button states based on current data
        bool hasFiles = !m_geneCountHandler->getGeneCountFiles().isEmpty();
        bool hasValidData = m_combinedData.isValid;

        if (m_generateDataButton)
        {
            m_generateDataButton->setEnabled(hasFiles);
        }
        if (m_exportGeneratedDataButton)
        {
            m_exportGeneratedDataButton->setEnabled(hasValidData);
        }

        // Update analysis buttons
        bool canRunAnalysis = hasValidData;
        if (m_runAnalysisButton)
        {
            m_runAnalysisButton->setEnabled(canRunAnalysis);
        }
        if (m_actionRunAnalysis)
        {
            m_actionRunAnalysis->setEnabled(canRunAnalysis);
        }

        // Update group assignment buttons
        if (m_assignGroupButton)
        {
            m_assignGroupButton->setEnabled(hasFiles);
        }
        if (m_autoAssignButton)
        {
            m_autoAssignButton->setEnabled(hasFiles);
        }
    }

    void MainWindow::addProgressMessage(const QString &message)
    {
        if (m_progressOutputText)
        {
            QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
            QString formattedMessage = QString("[%1] %2").arg(timestamp).arg(message);
            m_progressOutputText->append(formattedMessage);

            // Auto-scroll to bottom
            QTextCursor cursor = m_progressOutputText->textCursor();
            cursor.movePosition(QTextCursor::End);
            m_progressOutputText->setTextCursor(cursor);
        }
    }

    void MainWindow::updateLastUsedDirectory(const QString &filePath)
    {
        if (!filePath.isEmpty())
        {
            QFileInfo fileInfo(filePath);
            QString directory = fileInfo.absolutePath();
            if (!directory.isEmpty())
            {
                m_lastUsedDirectory = directory;
            }
        }
    }

    QString MainWindow::getLastUsedDirectory() const
    {
        return m_lastUsedDirectory;
    }

    void MainWindow::saveLastUsedDirectory()
    {
        QSettings settings;
        settings.setValue("LastUsedDirectory", m_lastUsedDirectory);
    }

    void MainWindow::loadLastUsedDirectory()
    {
        QSettings settings;
        QString savedDirectory = settings.value("LastUsedDirectory").toString();
        if (!savedDirectory.isEmpty() && QDir(savedDirectory).exists())
        {
            m_lastUsedDirectory = savedDirectory;
        }
    }

    void MainWindow::clearAllData()
    {
        // Clear tables
        if (m_geneCountFilesTable)
        {
            m_geneCountFilesTable->setRowCount(0);
        }
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->setRowCount(0);
            m_countsPreviewTable->setColumnCount(0);
        }
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->setRowCount(0);
            m_metadataPreviewTable->setColumnCount(0);
        }
        if (m_resultsTable)
        {
            m_resultsTable->setRowCount(0);
            m_resultsTable->setColumnCount(0);
        }

        // Reset labels
        if (m_totalGenesLabel)
        {
            m_totalGenesLabel->setText("Total Genes: 0");
        }
        if (m_significantGenesLabel)
        {
            m_significantGenesLabel->setText("Significant: 0");
        }
        if (m_upregulatedLabel)
        {
            m_upregulatedLabel->setText("Upregulated: 0");
        }
        if (m_downregulatedLabel)
        {
            m_downregulatedLabel->setText("Downregulated: 0");
        }

        // Clear group name
        if (m_groupNameEdit)
        {
            m_groupNameEdit->clear();
        }

        updateUiState();
        addProgressMessage("All data cleared.");
    }

    // Slot implementations - File operations
    void MainWindow::onAddGeneCountFiles()
    {
        QStringList fileNames = QFileDialog::getOpenFileNames(
            static_cast<QWidget *>(this),
            "Select Gene Count Files",
            m_lastUsedDirectory,
            "Gene Count Files (*.xlsx *.csv);;XLSX Files (*.xlsx);;CSV Files (*.csv);;All Files (*)");

        if (!fileNames.isEmpty())
        {
            // Update last used directory from the first selected file
            updateLastUsedDirectory(fileNames.first());

            int addedCount = 0;
            for (const QString &fileName : fileNames)
            {
                if (m_geneCountHandler->addGeneCountFile(fileName))
                {
                    addedCount++;
                }
            }

            addProgressMessage(QString("Successfully added %1 out of %2 gene count files.").arg(addedCount).arg(fileNames.size()));

            // Update the gene count files table
            m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

            // Update UI state
            updateUiState();
        }
    }

    void MainWindow::onClearFiles()
    {
        m_geneCountHandler->clearAllFiles();
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        // Clear preview tables
        if (m_countsPreviewTable)
        {
            m_countsPreviewTable->clear();
            m_countsPreviewTable->setRowCount(0);
            m_countsPreviewTable->setColumnCount(0);
        }
        if (m_metadataPreviewTable)
        {
            m_metadataPreviewTable->clear();
            m_metadataPreviewTable->setRowCount(0);
            m_metadataPreviewTable->setColumnCount(0);
        }

        m_combinedData = CombinedData();
        updateUiState();
        addProgressMessage("All gene count files cleared.");
    }

    void MainWindow::onGenerateData()
    {
        addProgressMessage("Generating count matrix and metadata...");

        // Check if we have files loaded
        if (m_geneCountHandler->getGeneCountFiles().isEmpty())
        {
            addProgressMessage("No gene count files loaded. Please add files first.");
            return;
        }

        // Check if all files have group assignments
        bool allAssigned = true;
        for (const auto &data : m_geneCountHandler->getGeneCountFiles())
        {
            if (data.groupName.isEmpty())
            {
                allAssigned = false;
                break;
            }
        }

        if (!allAssigned)
        {
            addProgressMessage("Some files do not have group assignments. Please assign groups to all files.");
            return;
        }

        // Generate count matrix and metadata
        m_combinedData = m_geneCountHandler->generateCountMatrixAndMetadata();

        if (m_combinedData.isValid)
        {
            // Update preview tables with converted integer counts for count matrix
            updateCountsPreviewWithConvertedValues(m_countsPreviewTable, m_combinedData);
            m_geneCountHandler->updateMetadataPreviewTable(m_metadataPreviewTable, m_combinedData);

            addProgressMessage(QString("Successfully generated count matrix with %1 genes and %2 samples.").arg(m_combinedData.geneNames.size()).arg(m_combinedData.sampleNames.size()));
            updateUiState();
        }
        else
        {
            addProgressMessage(QString("Failed to generate data: %1").arg(m_combinedData.errorMessage));
        }
    }

    void MainWindow::onExportGeneratedData()
    {
        if (!m_combinedData.isValid)
        {
            addProgressMessage("No data to export. Please generate data first.");
            return;
        }

        QString baseFileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Export Generated Data",
            m_lastUsedDirectory,
            "CSV Files (*.csv);;All Files (*)");

        if (!baseFileName.isEmpty())
        {
            // Update last used directory
            updateLastUsedDirectory(baseFileName);

            // Remove extension if present
            QFileInfo fileInfo(baseFileName);
            QString basePath = fileInfo.absolutePath() + "/" + fileInfo.baseName();

            // Export count matrix (converted from PPM to counts) and metadata
            bool countMatrixSuccess = exportCountMatrixToCSV(basePath + "_count_matrix.csv", m_combinedData);
            bool metadataSuccess = m_geneCountHandler->exportMetadataToCSV(basePath + "_metadata.csv", m_combinedData);

            // Provide detailed feedback
            if (countMatrixSuccess)
            {
                addProgressMessage(QString("✓ Count matrix exported to: %1_count_matrix.csv").arg(basePath));
            }
            else
            {
                addProgressMessage("✗ Failed to export count matrix.");
            }

            if (metadataSuccess)
            {
                addProgressMessage(QString("✓ Metadata exported to: %1_metadata.csv").arg(basePath));
            }
            else
            {
                addProgressMessage("✗ Failed to export metadata.");
            }

            // Summary message
            if (countMatrixSuccess && metadataSuccess)
            {
                addProgressMessage(QString("✓ Export completed successfully! Both files saved in: %1").arg(fileInfo.absolutePath()));

                // Show success dialog with file locations
                QString message = QString("Export completed successfully!\n\n"
                                          "Files saved:\n"
                                          "• Count Matrix: %1_count_matrix.csv\n"
                                          "• Metadata: %1_metadata.csv\n\n"
                                          "Both files are required for DESeq2 analysis.")
                                      .arg(fileInfo.baseName());

                QMessageBox::information(static_cast<QWidget *>(this), "Export Successful", message);
            }
            else if (countMatrixSuccess || metadataSuccess)
            {
                addProgressMessage("⚠ Partial export completed. Some files failed to export.");

                QString message = QString("Partial export completed.\n\n"
                                          "Successfully exported:\n");
                if (countMatrixSuccess)
                    message += "• Count Matrix\n";
                if (metadataSuccess)
                    message += "• Metadata\n";
                message += "\nFailed to export:\n";
                if (!countMatrixSuccess)
                    message += "• Count Matrix\n";
                if (!metadataSuccess)
                    message += "• Metadata\n";
                message += "\nBoth files are required for DESeq2 analysis.";

                QMessageBox::warning(static_cast<QWidget *>(this), "Partial Export", message);
            }
            else
            {
                addProgressMessage("✗ Export failed. No files were saved.");
                QMessageBox::critical(static_cast<QWidget *>(this), "Export Failed",
                                      "Failed to export both count matrix and metadata.\n"
                                      "Please check file permissions and try again.");
            }
        }
    }

    // Slot implementations - Group assignment
    void MainWindow::onAssignGroup()
    {
        if (!m_groupNameEdit || m_groupNameEdit->text().isEmpty())
        {
            QMessageBox::warning(static_cast<QWidget *>(this), "Warning", "Please enter a group name.");
            return;
        }

        if (!m_geneCountFilesTable)
        {
            return;
        }

        QString groupName = m_groupNameEdit->text();
        QList<QTableWidgetItem *> selectedItems = m_geneCountFilesTable->selectedItems();

        if (selectedItems.isEmpty())
        {
            QMessageBox::warning(static_cast<QWidget *>(this), "Warning", "Please select files to assign group.");
            return;
        }

        int assignedCount = 0;
        QSet<int> processedRows;

        for (QTableWidgetItem *item : selectedItems)
        {
            int row = item->row();
            if (processedRows.contains(row))
            {
                continue;
            }
            processedRows.insert(row);

            QTableWidgetItem *fileNameItem = m_geneCountFilesTable->item(row, 0);
            if (fileNameItem)
            {
                QString fileName = fileNameItem->text();
                if (m_geneCountHandler->assignGroupToFile(fileName, groupName))
                {
                    assignedCount++;
                }
            }
        }

        // Update the table to show the new group assignments
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        addProgressMessage(QString("Assigned group '%1' to %2 files.").arg(groupName).arg(assignedCount));
    }

    void MainWindow::onAutoAssignGroups()
    {
        addProgressMessage("Auto-assigning groups by filename...");

        // Define patterns for auto-assignment based on the example files
        QMap<QString, QString> patterns;
        patterns["SEL"] = "Selection";
        patterns["NON"] = "Non-Selection";
        patterns["_SEL_"] = "Selection";
        patterns["_NON_"] = "Non-Selection";

        int assignedCount = m_geneCountHandler->autoAssignGroups(patterns);

        // Update the table to show the new group assignments
        m_geneCountHandler->updateGeneCountFilesTable(m_geneCountFilesTable);

        addProgressMessage(QString("Auto-assigned groups to %1 files.").arg(assignedCount));
    }

    // Slot implementations - Analysis control
    void MainWindow::onRunAnalysis()
    {
        // Validate that we have data
        if (!m_combinedData.isValid)
        {
            QMessageBox::warning(this, "No Data",
                                 "Please load gene count files and generate data before running analysis.");
            addProgressMessage("Cannot run analysis: No valid data available.");
            return;
        }

        // Check that we have at least 2 different groups
        QSet<QString> uniqueGroups;
        for (const auto &group : m_combinedData.groupNames)
        {
            uniqueGroups.insert(group);
        }

        if (uniqueGroups.size() < 2)
        {
            QMessageBox::warning(this, "Insufficient Groups",
                                 "Analysis requires at least 2 different groups. Please assign groups to your samples.");
            addProgressMessage("Cannot run analysis: Need at least 2 different groups.");
            return;
        }

        // Prevent multiple simultaneous analyses
        if (m_analysisThread && m_analysisThread->isRunning())
        {
            QMessageBox::information(this, "Analysis Running",
                                     "Analysis is already running. Please wait for it to complete or stop it first.");
            return;
        }

        addProgressMessage("Starting DESeq2 differential expression analysis...");

        // Get analysis parameters
        double pValueThreshold = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;

        // Update UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(false);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(true);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(false);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(true);

        // Show progress bar
        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(true);
            m_analysisProgressBar->setValue(0);
        }

        // Create worker and thread
        m_analysisWorker = std::make_unique<AnalysisWorker>(m_combinedData, pValueThreshold);
        m_analysisThread = std::make_unique<QThread>();

        // Move worker to thread
        m_analysisWorker->moveToThread(m_analysisThread.get());

        // Connect signals with explicit connection types for cross-thread communication
        connect(m_analysisThread.get(), &QThread::started,
                m_analysisWorker.get(), &AnalysisWorker::runAnalysis);
        connect(m_analysisWorker.get(), &AnalysisWorker::progressChanged,
                this, &MainWindow::onAnalysisProgress, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::analysisFinished,
                this, &MainWindow::onAnalysisFinished, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::analysisError,
                this, &MainWindow::onAnalysisError, Qt::QueuedConnection);
        connect(m_analysisWorker.get(), &AnalysisWorker::debugMessage,
                this, &MainWindow::onDebugMessage, Qt::QueuedConnection);

        // Proper thread lifecycle management
        connect(m_analysisWorker.get(), &AnalysisWorker::finished,
                m_analysisThread.get(), &QThread::quit);
        connect(m_analysisThread.get(), &QThread::finished,
                this, &MainWindow::onThreadFinished);

        // Force signal processing by connecting finished signal to a dummy slot
        connect(m_analysisWorker.get(), &AnalysisWorker::finished,
                this, [this]()
                { addProgressMessage("DEBUG: Worker finished signal received"); });

        // Start analysis
        m_analysisThread->start();
    }

    void MainWindow::onStopAnalysis()
    {
        if (m_analysisWorker)
        {
            addProgressMessage("Stopping analysis...");
            m_analysisWorker->stopAnalysis();
        }

        if (m_analysisThread && m_analysisThread->isRunning())
        {
            m_analysisThread->quit();
            if (!m_analysisThread->wait(5000))
            {
                m_analysisThread->terminate();
                m_analysisThread->wait();
            }
        }

        // Reset UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(true);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(false);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(true);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(false);

        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(false);
        }

        addProgressMessage("Analysis stopped by user.");
    }

    void MainWindow::onResetAnalysis()
    {
        // Stop any running analysis
        onStopAnalysis();

        // Clear results
        m_analysisResults.isValid = false;
        m_analysisResults.results = Eigen::MatrixXd();
        m_analysisResults.geneNames.clear();

        // Clear results table
        if (m_resultsTable)
        {
            m_resultsTable->setRowCount(0);
            m_resultsTable->setColumnCount(0);
        }

        // Reset summary labels
        updateResultsSummary(m_analysisResults);

        addProgressMessage("Analysis results cleared.");
    }

    void MainWindow::onAnalysisProgress(int progress, const QString &message)
    {
        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setValue(progress);
        }
        addProgressMessage(message);
    }

    void MainWindow::onAnalysisFinished(AnalysisResults *results)
    {
        addProgressMessage("=== ANALYSIS FINISHED CALLED ===");
        addProgressMessage("Results pointer: " + QString(results ? "VALID" : "NULL"));

        if (!results)
        {
            addProgressMessage("ERROR: Results pointer is null!");
            return;
        }

        addProgressMessage("Results valid: " + QString(results->isValid ? "YES" : "NO"));
        addProgressMessage("Total genes: " + QString::number(results->totalGenes));
        addProgressMessage("Significant genes: " + QString::number(results->significantGenes));
        addProgressMessage("Matrix size: " + QString::number(results->results.rows()) + "x" + QString::number(results->results.cols()));
        addProgressMessage("Gene names count: " + QString::number(results->geneNames.size()));

        // Store results
        m_analysisResults = *results;

        // Update UI
        addProgressMessage("About to update results table...");
        updateResultsTable(*results);
        addProgressMessage("About to update results summary...");
        updateResultsSummary(*results);

        // Switch to results tab
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab (FIXED: was 2, should be 1)
            addProgressMessage("Switched to results tab");
        }
        else
        {
            addProgressMessage("ERROR: m_tabWidget is null!");
        }

        addProgressMessage(QString("Analysis completed! Found %1 significant genes out of %2 total genes.")
                               .arg(results->significantGenes)
                               .arg(results->totalGenes));

        // Refresh visualization plot with new results
        refreshPlot();
    }

    void MainWindow::onAnalysisError(const QString &errorMessage)
    {
        // Show error message
        QMessageBox::critical(this, "Analysis Error", errorMessage);
        addProgressMessage("Analysis failed: " + errorMessage);
    }

    void MainWindow::onDebugMessage(const QString &message)
    {
        addProgressMessage("DEBUG: " + message);
    }

    void MainWindow::onThreadFinished()
    {
        // Reset UI state
        if (m_runAnalysisButton)
            m_runAnalysisButton->setEnabled(true);
        if (m_stopAnalysisButton)
            m_stopAnalysisButton->setEnabled(false);
        if (m_actionRunAnalysis)
            m_actionRunAnalysis->setEnabled(true);
        if (m_actionStopAnalysis)
            m_actionStopAnalysis->setEnabled(false);

        if (m_analysisProgressBar)
        {
            m_analysisProgressBar->setVisible(false);
        }
    }

    void MainWindow::updateResultsTable(const AnalysisResults &results)
    {
        addProgressMessage("updateResultsTable called");
        addProgressMessage("m_resultsTable: " + QString(m_resultsTable ? "VALID" : "NULL"));
        addProgressMessage("results.isValid: " + QString(results.isValid ? "YES" : "NO"));

        if (!m_resultsTable || !results.isValid)
        {
            addProgressMessage("updateResultsTable: Early return - table or results invalid");
            return;
        }

        addProgressMessage("updateResultsTable: Proceeding with table update");

        const auto &matrix = results.results;
        const auto &geneNames = results.geneNames;

        // Setup table
        addProgressMessage("Setting table dimensions: " + QString::number(matrix.rows()) + " rows, 7 columns");
        m_resultsTable->setRowCount(matrix.rows());
        m_resultsTable->setColumnCount(7); // Gene name + 6 columns from DESeq2

        // Set column headers
        QStringList headers;
        headers << "Gene" << "baseMean" << "log2FoldChange" << "lfcSE" << "stat" << "pvalue" << "padj";
        m_resultsTable->setHorizontalHeaderLabels(headers);
        addProgressMessage("Set table headers");

        // Populate table
        addProgressMessage("Starting to populate table with " + QString::number(matrix.rows()) + " rows");

        for (int row = 0; row < matrix.rows(); ++row)
        {
            // Gene name
            QTableWidgetItem *geneItem = new QTableWidgetItem(QString::fromStdString(geneNames[row]));
            geneItem->setFlags(geneItem->flags() & ~Qt::ItemIsEditable);
            m_resultsTable->setItem(row, 0, geneItem);

            // Results columns
            for (int col = 0; col < matrix.cols(); ++col)
            {
                QTableWidgetItem *item = new QTableWidgetItem(QString::number(matrix(row, col), 'g', 6));
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);

                // Color significant genes - padj is column 5 in matrix (0-indexed)
                if (col == 5 && matrix(row, 5) < results.pValueThreshold) // padj column (5th column in matrix, 0-indexed)
                {
                    item->setBackground(QColor(255, 255, 200)); // Light yellow background
                }

                m_resultsTable->setItem(row, col + 1, item);
            }

            // Progress update every 1000 rows
            if (row % 1000 == 0 && row > 0)
            {
                addProgressMessage("Populated " + QString::number(row) + " rows...");
            }
        }

        addProgressMessage("Finished populating table");

        // Auto-resize columns
        addProgressMessage("Resizing columns...");
        m_resultsTable->resizeColumnsToContents();

        // Enable sorting
        addProgressMessage("Enabling sorting...");
        m_resultsTable->setSortingEnabled(true);

        // Sort by adjusted p-value by default
        addProgressMessage("Sorting by adjusted p-value...");
        m_resultsTable->sortItems(6, Qt::AscendingOrder);

        addProgressMessage("updateResultsTable completed successfully!");
    }

    void MainWindow::updateResultsSummary(const AnalysisResults &results)
    {
        if (results.isValid)
        {
            if (m_totalGenesLabel)
            {
                m_totalGenesLabel->setText(QString("Total Genes: %1").arg(results.totalGenes));
            }
            if (m_significantGenesLabel)
            {
                m_significantGenesLabel->setText(QString("Significant: %1").arg(results.significantGenes));
            }
            if (m_upregulatedLabel)
            {
                m_upregulatedLabel->setText(QString("Upregulated: %1").arg(results.upregulatedGenes));
            }
            if (m_downregulatedLabel)
            {
                m_downregulatedLabel->setText(QString("Downregulated: %1").arg(results.downregulatedGenes));
            }
        }
        else
        {
            // Reset to defaults
            if (m_totalGenesLabel)
                m_totalGenesLabel->setText("Total Genes: 0");
            if (m_significantGenesLabel)
                m_significantGenesLabel->setText("Significant: 0");
            if (m_upregulatedLabel)
                m_upregulatedLabel->setText("Upregulated: 0");
            if (m_downregulatedLabel)
                m_downregulatedLabel->setText("Downregulated: 0");
        }
    }

    // Slot implementations - Results filtering
    void MainWindow::onApplyFilter()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        applyCurrentFilter();
    }

    void MainWindow::onClearFilter()
    {
        if (m_pValueThresholdSpinBox)
        {
            m_pValueThresholdSpinBox->setValue(0.05);
        }
        if (m_log2FCThresholdSpinBox)
        {
            m_log2FCThresholdSpinBox->setValue(0.0);
        }

        if (m_analysisResults.isValid && m_resultsTable)
        {
            // Show all rows instead of repopulating the table
            for (int row = 0; row < m_resultsTable->rowCount(); ++row)
            {
                m_resultsTable->setRowHidden(row, false);
            }
            addProgressMessage("Filter cleared - showing all results.");
        }
    }

    void MainWindow::applyCurrentFilter()
    {
        if (!m_resultsTable || !m_analysisResults.isValid)
        {
            return;
        }

        double pValueThreshold = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        double log2FCThreshold = m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 0.0;

        int visibleRows = 0;
        for (int row = 0; row < m_resultsTable->rowCount(); ++row)
        {
            // Get p-value and log2FC
            double padj = m_resultsTable->item(row, 6)->text().toDouble();             // padj column
            double log2FC = std::abs(m_resultsTable->item(row, 2)->text().toDouble()); // log2FoldChange column

            bool passesFilter = (padj < pValueThreshold) && (log2FC > log2FCThreshold);

            m_resultsTable->setRowHidden(row, !passesFilter);
            if (passesFilter)
                visibleRows++;
        }

        addProgressMessage(QString("Applied filter: p-value < %1, |log2FC| > %2. Showing %3 genes.")
                               .arg(pValueThreshold)
                               .arg(log2FCThreshold)
                               .arg(visibleRows));
    }

    // Slot implementations - Export operations
    void MainWindow::onExportResults()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Export Results",
            m_lastUsedDirectory + "/deseq2_results.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportResultsToCSV(fileName, m_analysisResults, false))
            {
                addProgressMessage("Results exported successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to export results.");
                addProgressMessage("Failed to export results.");
            }
        }
    }

    void MainWindow::onSaveGeneLists()
    {
        if (!m_analysisResults.isValid)
        {
            QMessageBox::information(this, "No Results", "Please run analysis first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Save Significant Genes",
            m_lastUsedDirectory + "/significant_genes.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportResultsToCSV(fileName, m_analysisResults, true))
            {
                addProgressMessage("Significant genes saved successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to save gene lists.");
                addProgressMessage("Failed to save gene lists.");
            }
        }
    }

    void MainWindow::onExportCountMatrix()
    {
        if (!m_combinedData.isValid)
        {
            QMessageBox::information(this, "No Data", "Please load and generate data first.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this,
            "Export Count Matrix",
            m_lastUsedDirectory + "/count_matrix.csv",
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (exportCountMatrixToCSV(fileName, m_combinedData))
            {
                addProgressMessage("Count matrix exported successfully to: " + fileName);
            }
            else
            {
                QMessageBox::warning(this, "Export Error", "Failed to export count matrix.");
                addProgressMessage("Failed to export count matrix.");
            }
        }
    }

    bool MainWindow::exportResultsToCSV(const QString &filePath, const AnalysisResults &results, bool significantOnly)
    {
        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return false;
        }

        QTextStream stream(&file);

        // Write header
        stream << "Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n";

        const auto &matrix = results.results;
        const auto &geneNames = results.geneNames;

        // Write data
        for (int row = 0; row < matrix.rows(); ++row)
        {
            double padj = matrix(row, 5); // Adjusted p-value column (FIXED: was 4, should be 5)

            // Skip non-significant genes if requested
            if (significantOnly && padj >= results.pValueThreshold)
            {
                continue;
            }

            // Gene name
            stream << QString::fromStdString(geneNames[row]);

            // Results columns
            for (int col = 0; col < matrix.cols(); ++col)
            {
                stream << "," << matrix(row, col);
            }
            stream << "\n";
        }

        return true;
    }

    bool MainWindow::exportCountMatrixToCSV(const QString &filePath, const CombinedData &data)
    {
        QFile file(filePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            return false;
        }

        QTextStream stream(&file);

        // Write header with sample names
        stream << "Gene";
        for (const auto &sampleName : data.sampleNames)
        {
            stream << "," << sampleName;
        }
        stream << "\n";

        // Write data - use raw integer counts
        for (int gene = 0; gene < data.geneNames.size(); ++gene)
        {
            // Gene name
            stream << data.geneNames[gene];

            // Count values for each sample
            for (int sample = 0; sample < data.sampleNames.size(); ++sample)
            {
                stream << "," << data.rawCountMatrix[gene][sample];
            }
            stream << "\n";
        }

        return true;
    }

    // Slot implementations - Visualization
    void MainWindow::onPlotTypeChanged()
    {
        if (!m_plotTypeCombo)
            return;
        refreshPlot();
    }

    void MainWindow::refreshPlot()
    {
        if (!m_plotTypeCombo || !m_analysisResults.isValid)
            return;

        QString plotType = m_plotTypeCombo->currentText();
        addProgressMessage(QString("Drawing %1...").arg(plotType));

        if (plotType == "Volcano Plot")
            plotVolcanoPlot();
        else if (plotType == "MA Plot")
            plotMAPlot();
        else if (plotType == "Dispersion Plot")
            plotDispersionPlot();
    }

    void MainWindow::plotMAPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        m_plotWidget->clearPlottables();
        m_plotWidget->clearItems();

        const auto &res = m_analysisResults.results;
        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;

        QCPGraph *sigGraph = m_plotWidget->addGraph();
        sigGraph->setLineStyle(QCPGraph::lsNone);
        sigGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor(220, 50, 50), ptSize));
        sigGraph->setName("Significant");

        QCPGraph *nsGraph = m_plotWidget->addGraph();
        nsGraph->setLineStyle(QCPGraph::lsNone);
        nsGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor(150, 150, 150), ptSize));
        nsGraph->setName("Non-significant");

        QVector<double> sigX, sigY, nsX, nsY;

        for (int i = 0; i < res.rows(); ++i)
        {
            double baseMean = res(i, 0);
            double lfc = res(i, 1);
            double padj = res(i, 5);

            if (baseMean <= 0)
                continue;
            double x = std::log10(baseMean);

            if (padj < pThresh)
            {
                sigX.append(x);
                sigY.append(lfc);
            }
            else
            {
                nsX.append(x);
                nsY.append(lfc);
            }
        }

        sigGraph->setData(sigX, sigY);
        nsGraph->setData(nsX, nsY);

        m_plotWidget->xAxis->setLabel("log10(baseMean)");
        m_plotWidget->yAxis->setLabel("log2 Fold Change");
        m_plotWidget->legend->setVisible(true);
        m_plotWidget->rescaleAxes();
        m_plotWidget->replot();
    }

    void MainWindow::plotVolcanoPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        m_plotWidget->clearPlottables();
        m_plotWidget->clearItems();

        const auto &res = m_analysisResults.results;
        double pThresh = m_pValueThresholdSpinBox ? m_pValueThresholdSpinBox->value() : 0.05;
        double fcThresh = m_log2FCThresholdSpinBox ? m_log2FCThresholdSpinBox->value() : 1.0;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;

        QCPGraph *sigGraph = m_plotWidget->addGraph();
        sigGraph->setLineStyle(QCPGraph::lsNone);
        sigGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor(220, 50, 50), ptSize));
        sigGraph->setName("Significant");

        QCPGraph *nsGraph = m_plotWidget->addGraph();
        nsGraph->setLineStyle(QCPGraph::lsNone);
        nsGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor(150, 150, 150), ptSize));
        nsGraph->setName("Non-significant");

        QVector<double> sigX, sigY, nsX, nsY;

        for (int i = 0; i < res.rows(); ++i)
        {
            double lfc = res(i, 1);
            double pval = res(i, 4);
            double padj = res(i, 5);

            if (pval <= 0 || !std::isfinite(pval))
                continue;
            double negLog10P = -std::log10(pval);

            if (padj < pThresh && std::abs(lfc) > fcThresh)
            {
                sigX.append(lfc);
                sigY.append(negLog10P);
            }
            else
            {
                nsX.append(lfc);
                nsY.append(negLog10P);
            }
        }

        sigGraph->setData(sigX, sigY);
        nsGraph->setData(nsX, nsY);

        // Add threshold lines
        // Horizontal line at -log10(alpha)
        double hLine = -std::log10(pThresh);
        QCPItemStraightLine *pLine = new QCPItemStraightLine(m_plotWidget);
        pLine->point1->setCoords(0, hLine);
        pLine->point2->setCoords(1, hLine);
        pLine->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));

        // Vertical lines at +/- FC threshold
        QCPItemStraightLine *fcLinePos = new QCPItemStraightLine(m_plotWidget);
        fcLinePos->point1->setCoords(fcThresh, 0);
        fcLinePos->point2->setCoords(fcThresh, 1);
        fcLinePos->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));

        QCPItemStraightLine *fcLineNeg = new QCPItemStraightLine(m_plotWidget);
        fcLineNeg->point1->setCoords(-fcThresh, 0);
        fcLineNeg->point2->setCoords(-fcThresh, 1);
        fcLineNeg->setPen(QPen(Qt::darkGray, 1, Qt::DashLine));

        m_plotWidget->xAxis->setLabel("log2 Fold Change");
        m_plotWidget->yAxis->setLabel("-log10(p-value)");
        m_plotWidget->legend->setVisible(true);
        m_plotWidget->rescaleAxes();
        m_plotWidget->replot();
    }

    void MainWindow::plotDispersionPlot()
    {
        if (!m_plotWidget || !m_analysisResults.isValid)
            return;

        m_plotWidget->clearPlottables();
        m_plotWidget->clearItems();

        const auto &baseMeans = m_analysisResults.baseMeans;
        const auto &genewise = m_analysisResults.genewiseDispersions;
        const auto &fitted = m_analysisResults.fittedDispersions;
        int ptSize = m_pointSizeSpinBox ? m_pointSizeSpinBox->value() : 5;

        bool hasData = baseMeans.size() > 0 && genewise.size() > 0 && fitted.size() > 0;
        if (!hasData)
        {
            addProgressMessage("No dispersion data available for plotting.");
            return;
        }

        // Genewise dispersions scatter
        QCPGraph *geneGraph = m_plotWidget->addGraph();
        geneGraph->setLineStyle(QCPGraph::lsNone);
        geneGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor(70, 130, 220), ptSize));
        geneGraph->setName("Genewise");

        QVector<double> gX, gY;
        for (int i = 0; i < baseMeans.size(); ++i)
        {
            if (baseMeans(i) > 0 && genewise(i) > 0)
            {
                gX.append(std::log10(baseMeans(i)));
                gY.append(std::log10(genewise(i)));
            }
        }
        geneGraph->setData(gX, gY);

        // Fitted trend line (sorted by x for connected line)
        QCPGraph *fitGraph = m_plotWidget->addGraph();
        fitGraph->setLineStyle(QCPGraph::lsLine);
        fitGraph->setPen(QPen(QColor(220, 50, 50), 2));
        fitGraph->setScatterStyle(QCPScatterStyle::ssNone);
        fitGraph->setName("Fitted trend");

        // Collect valid fitted points and sort by mean
        std::vector<std::pair<double, double>> fitPoints;
        for (int i = 0; i < baseMeans.size(); ++i)
        {
            if (baseMeans(i) > 0 && fitted(i) > 0)
            {
                fitPoints.push_back({std::log10(baseMeans(i)), std::log10(fitted(i))});
            }
        }
        std::sort(fitPoints.begin(), fitPoints.end());

        QVector<double> fX, fY;
        for (const auto &pt : fitPoints)
        {
            fX.append(pt.first);
            fY.append(pt.second);
        }
        fitGraph->setData(fX, fY);

        m_plotWidget->xAxis->setLabel("log10(baseMean)");
        m_plotWidget->yAxis->setLabel("log10(dispersion)");
        m_plotWidget->legend->setVisible(true);
        m_plotWidget->rescaleAxes();
        m_plotWidget->replot();
    }

    void MainWindow::onSavePlot()
    {
        if (!m_plotWidget)
            return;

        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Save Plot",
            m_lastUsedDirectory,
            "PNG Files (*.png);;PDF Files (*.pdf);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            bool ok = false;
            if (fileName.endsWith(".pdf", Qt::CaseInsensitive))
            {
                ok = m_plotWidget->savePdf(fileName);
            }
            else
            {
                ok = m_plotWidget->savePng(fileName, 0, 0, 2.0);
            }

            if (ok)
                addProgressMessage("Plot saved to: " + fileName);
            else
                addProgressMessage("Failed to save plot.");
        }
    }

    void MainWindow::onExportPlotData()
    {
        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Export Plot Data",
            m_lastUsedDirectory,
            "CSV Files (*.csv);;All Files (*)");

        if (!fileName.isEmpty())
        {
            updateLastUsedDirectory(fileName);

            if (m_analysisResults.isValid)
            {
                exportResultsToCSV(fileName, m_analysisResults, false);
                addProgressMessage("Plot data exported to: " + fileName);
            }
            else
            {
                addProgressMessage("No analysis results to export.");
            }
        }
    }

    void MainWindow::onPrintPlot()
    {
        if (!m_plotWidget)
            return;

        // Save to a temp PDF and inform the user
        QString tempPath = QDir::tempPath() + "/deseq2_plot.pdf";
        if (m_plotWidget->savePdf(tempPath))
        {
            addProgressMessage("Plot saved to temporary PDF: " + tempPath);
        }
        else
        {
            addProgressMessage("Failed to generate print PDF.");
        }
    }

    // Slot implementations - Progress output
    void MainWindow::onClearOutput()
    {
        if (m_progressOutputText)
        {
            m_progressOutputText->clear();
        }
    }

    void MainWindow::onSaveOutput()
    {
        QString fileName = QFileDialog::getSaveFileName(
            static_cast<QWidget *>(this),
            "Save Output",
            m_lastUsedDirectory,
            "Text Files (*.txt);;All Files (*)");

        if (!fileName.isEmpty() && m_progressOutputText)
        {
            // Update last used directory
            updateLastUsedDirectory(fileName);

            QFile file(fileName);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QTextStream stream(&file);
                stream << m_progressOutputText->toPlainText();
                addProgressMessage("Output saved successfully.");
            }
            else
            {
                QMessageBox::warning(static_cast<QWidget *>(this), "Error", "Could not save output file.");
            }
        }
    }

    // Slot implementations - Menu actions
    void MainWindow::onOpenCounts()
    {
        onAddGeneCountFiles();
    }

    void MainWindow::onOpenMetadata()
    {
        onGenerateData();
    }

    void MainWindow::onSaveResults()
    {
        onExportResults();
    }

    void MainWindow::onExportVisualizations()
    {
        onSavePlot();
    }

    void MainWindow::onExit()
    {
        close();
    }

    void MainWindow::onRunAnalysisAction()
    {
        onRunAnalysis();
    }

    void MainWindow::onStopAnalysisAction()
    {
        onStopAnalysis();
    }

    void MainWindow::onResetAnalysisAction()
    {
        onResetAnalysis();
    }

    void MainWindow::onInputTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(0);
        }
    }

    void MainWindow::onAnalysisTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab (no analysis tab exists)
        }
    }

    void MainWindow::onResultsTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(1); // Results tab
        }
    }

    void MainWindow::onVisualizationTab()
    {
        if (m_tabWidget)
        {
            m_tabWidget->setCurrentIndex(2); // Visualization tab
        }
        refreshPlot();
    }

    void MainWindow::onAbout()
    {
        QMessageBox::about(static_cast<QWidget *>(this), "About DESeq2 GUI",
                           "DESeq2 Differential Expression Analysis GUI\n\n"
                           "Version 1.0.0\n"
                           "A Qt5-based interface for DESeq2 analysis.");
    }

    void MainWindow::onUserGuide()
    {
        QMessageBox::information(static_cast<QWidget *>(this), "User Guide",
                                 "User guide functionality will be implemented in future versions.");
    }

    void MainWindow::updateCountsPreviewWithConvertedValues(QTableWidget *table, const CombinedData &data)
    {
        if (!table || !data.isValid)
        {
            return;
        }

        const int maxRows = std::min(20, data.geneNames.size());
        const int maxCols = std::min(10, data.sampleNames.size());

        table->clear();
        table->setRowCount(maxRows);
        table->setColumnCount(maxCols + 1); // +1 for gene names

        QStringList headers = {"Gene"};
        for (int i = 0; i < maxCols; ++i)
        {
            headers.append(data.sampleNames[i]);
        }
        table->setHorizontalHeaderLabels(headers);

        for (int row = 0; row < maxRows; ++row)
        {
            // Set gene name in first column
            table->setItem(row, 0, new QTableWidgetItem(data.geneNames[row]));

            for (int col = 0; col < maxCols; ++col)
            {
                // Use raw integer counts directly
                int integerCount = data.rawCountMatrix[row][col];

                QTableWidgetItem *item = new QTableWidgetItem(QString::number(integerCount));
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);
                table->setItem(row, col + 1, item);
            }
        }

        table->resizeColumnsToContents();
    }

} // namespace deseq2
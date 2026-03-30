#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include <QTableWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QTextEdit>
#include <QTabWidget>
#include <QGroupBox>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QProgressBar>
#include <QThread>
#include <QMutex>
#include <QMetaType>
#include <memory>

// Statistics library includes
#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"
#include "y2h_scores.h"

// Gene count handler
#include "gene_count_handler.h"

// Qt Charts for visualizations (migrated from QCustomPlot)
#include <QChartView>
#include <QChart>

namespace Ui
{
    class MainWindow;
}

QT_BEGIN_NAMESPACE
class QWidget;
QT_END_NAMESPACE

namespace deseq2
{

    /**
     * @brief Structure to hold DESeq2 analysis results
     */
    struct AnalysisResults
    {
        Eigen::MatrixXd results;              // DESeq2 results matrix
        std::vector<std::string> geneNames;   // Gene names
        std::vector<std::string> sampleNames; // Sample names
        int totalGenes = 0;                  // Total number of genes
        int significantGenes = 0;            // Number of significant genes
        int upregulatedGenes = 0;            // Number of significant genes
        int downregulatedGenes = 0;          // Number of significant genes
        bool isValid = false;                // Whether results are valid
        QString errorMessage;                 // Error message if invalid
        double pValueThreshold = 0.05;        // P-value threshold used
        double log2FCThreshold = 0.0;         // Log2FC threshold used

        // Dispersion data for visualization
        Eigen::VectorXd genewiseDispersions;
        Eigen::VectorXd fittedDispersions;
        Eigen::VectorXd baseMeans;

        // Y2H-SCORES results
        std::vector<deseq2::EnrichmentResult> enrichmentScores;
        std::vector<deseq2::SpecificityResult> specificityScores;
        std::vector<deseq2::InFrameResult> inFrameScores;
        std::vector<deseq2::Y2HScore> y2hScores;

        // Multi-contrast results for three-way comparisons
        // contrastResults[i] holds the Nx6 results matrix for contrast i
        std::vector<Eigen::MatrixXd> contrastResults;
        // contrastLabels[i] describes each contrast (e.g. "Group 1 vs Group 0")
        std::vector<std::string> contrastLabels;
        int activeContrast = 0; // Index of currently displayed contrast
        QString activeContrastLabel;
        QString createdAt;
    };

    struct AnalysisRunConfig
    {
        double pValueThreshold = 0.05;
        double ppmThreshold = 0.0;
        double enrichmentPValueThreshold = 1.0;
        double enrichmentFoldChangeThreshold = 0.0;
        double specificityPValueThreshold = 1.0;
        int baitGroupSize = 10;
        QStringList junctionFiles;
        QString workingDirectory;
    };

    /**
     * @brief Worker class for running DESeq2 analysis in background thread
     */
    class AnalysisWorker : public QObject
    {
        Q_OBJECT

    public:
        AnalysisWorker(const CombinedData &data, const AnalysisRunConfig &config);
        ~AnalysisWorker() = default;

    public slots:
        /**
         * @brief Run the DESeq2 analysis
         */
        void runAnalysis();

        /**
         * @brief Stop the analysis
         */
        void stopAnalysis();

    signals:
        /**
         * @brief Emitted when analysis progress changes
         * @param progress Progress percentage (0-100)
         * @param message Progress message
         */
        void progressChanged(int progress, const QString &message);

        /**
         * @brief Emitted when analysis is finished
         * @param results Analysis results passed by value for safe queued delivery
         */
        void analysisFinished(deseq2::AnalysisResults results);

        /**
         * @brief Emitted when analysis encounters an error
         * @param errorMessage Error message
         */
        void analysisError(const QString &errorMessage);

        /**
         * @brief Debug signal to track analysis completion
         * @param message Debug message
         */
        void debugMessage(const QString &message);

        /**
         * @brief Emitted when worker is completely finished and thread should quit
         */
        void finished();

    private:
        CombinedData m_inputData;
        AnalysisRunConfig m_config;
        bool m_shouldStop;
        QMutex m_stopMutex;

        /**
         * @brief Convert Qt CombinedData to Eigen matrices
         * @param data Input combined data
         * @return Tuple of (counts, metadata)
         */
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> convertToEigenFormat(const CombinedData &data);
    };

    /**
     * @brief Main window class for DESeq2 differential expression analysis GUI
     *
     * This class provides the main interface for the DESeq2 analysis application,
     * including file input, data processing, analysis execution, and result visualization.
     */
    class MainWindow : public QMainWindow
    {
        Q_OBJECT

    public:
        /**
         * @brief Constructor for the main window
         * @param parent Parent widget (default: nullptr)
         */
        explicit MainWindow(QWidget *parent = nullptr);

        /**
         * @brief Destructor
         */
        ~MainWindow() override;

        // Rule of Five - disable copy operations
        MainWindow(const MainWindow &) = delete;
        MainWindow &operator=(const MainWindow &) = delete;
        MainWindow(MainWindow &&) = delete;
        MainWindow &operator=(MainWindow &&) = delete;

        void loadWorkingDirectory(const QString &workdir);
        void loadInputFiles(const QStringList &filePaths, bool clearExisting = false);

    private slots:
        // File operations
        void onAddGeneCountFiles();
        void onClearFiles();
        void onGenerateData();
        void onExportGeneratedData();

        // Group assignment
        void onAssignGroup();
        void onAutoAssignGroups();

        // Analysis control
        void onRunAnalysis();
        void onStopAnalysis();
        void onResetAnalysis();

        // Analysis worker slots
        void onAnalysisProgress(int progress, const QString &message);
        void onAnalysisFinished(const AnalysisResults &results);
        void onAnalysisError(const QString &errorMessage);
        void onDebugMessage(const QString &message);
        void onThreadFinished();

        // Results filtering
        void onApplyFilter();
        void onClearFilter();

        // Export operations
        void onExportResults();
        void onSaveGeneLists();
        void onExportCountMatrix();

        // Visualization
        void onPlotTypeChanged();
        void onSavePlot();
        void onExportPlotData();
        void onPrintPlot();

        // Progress output
        void onClearOutput();
        void onSaveOutput();

        // Results table context menu
        void onResultsContextMenu(const QPoint &pos);
        void onOpenInMultiQuery();
        void onOpenInReadDepth();
        void onCopyGeneName();

        // Junction files browsing
        void onBrowseJunctionFiles();

        // Contrast selector
        void onContrastChanged();

        // Menu actions
        void onOpenCounts();
        void onOpenMetadata();
        void onSaveResults();
        void onExportVisualizations();
        void onExit();
        void onRunAnalysisAction();
        void onStopAnalysisAction();
        void onResetAnalysisAction();
        void onInputTab();
        void onAnalysisTab();
        void onResultsTab();
        void onVisualizationTab();
        void onAbout();
        void onUserGuide();

    private:
        /**
         * @brief Connect UI elements from the UI file to member variables
         */
        void connectUiElements();

        /**
         * @brief Setup signal-slot connections
         */
        void setupConnections();

        /**
         * @brief Setup menu bar and toolbar
         */
        void setupMenusAndToolbar();

        /**
         * @brief Initialize table widgets
         */
        void initializeTables();

        /**
         * @brief Update UI state based on current data
         */
        void updateUiState();

        /**
         * @brief Add a progress message to the output
         * @param message Message to add
         */
        void addProgressMessage(const QString &message);

        /**
         * @brief Update the last used directory from a file path
         * @param filePath File path to extract directory from
         */
        void updateLastUsedDirectory(const QString &filePath);

        /**
         * @brief Get the last used directory for file dialogs
         * @return Last used directory path
         */
        QString getLastUsedDirectory() const;

        /**
         * @brief Save the last used directory to application settings
         */
        void saveLastUsedDirectory();

        /**
         * @brief Load the last used directory from application settings
         */
        void loadLastUsedDirectory();

        /**
         * @brief Clear all data and reset UI
         */
        void clearAllData();

        /**
         * @brief Update results table with analysis results
         * @param results Analysis results to display
         */
        void updateResultsTable(const AnalysisResults &results);

        /**
         * @brief Update results summary labels
         * @param results Analysis results
         */
        void updateResultsSummary(const AnalysisResults &results);

        /**
         * @brief Apply current filter settings to results
         */
        void applyCurrentFilter();

        /**
         * @brief Export analysis results to CSV file
         * @param filePath Output file path
         * @param results Results to export
         * @param significantOnly Whether to export only significant genes
         * @return True if successful, false otherwise
         */
        bool exportResultsToCSV(const QString &filePath, const AnalysisResults &results, bool significantOnly = false);
        bool exportCountMatrixToCSV(const QString &filePath, const CombinedData &data);

        /**
         * @brief Update counts preview table with converted integer counts instead of PPM values
         * @param table Table widget to update
         * @param data Combined data with PPM values to convert
         */
        void updateCountsPreviewWithConvertedValues(QTableWidget *table, const CombinedData &data);
        void populateAutoDiscoveredFiles();
        void autoDetectJunctionFiles();
        QString resolveWorkingDirectory() const;
        QString resultsDatabasePathForWorkdir(const QString &workdir) const;
        bool loadResultsFromSqlite(const QString &sqlitePath);
        bool validateAnalysisResults(const AnalysisResults &results, QString *errorMessage = nullptr) const;
        void invalidatePlotCache();

        // UI components
        QTabWidget *m_tabWidget;
        QWidget *m_inputTab;
        QWidget *m_resultsTab;
        QWidget *m_visualizationTab;
        QWidget *m_analysisTab;

        // Analysis settings
        QComboBox *m_analysisModeCombo;
        QDoubleSpinBox *m_ppmThresholdSpinBox_y2h;
        QComboBox *m_groupTypeCombo;

        // Input tab components
        QPushButton *m_addPpmFilesButton;
        QPushButton *m_clearFilesButton;
        QPushButton *m_generateDataButton;
        QPushButton *m_exportGeneratedDataButton;
        QTableWidget *m_geneCountFilesTable;
        QLineEdit *m_groupNameEdit;
        QPushButton *m_assignGroupButton;
        QPushButton *m_autoAssignButton;
        QTableWidget *m_countsPreviewTable;
        QTableWidget *m_metadataPreviewTable;

        // Analysis tab components
        QPushButton *m_runAnalysisButton;
        QPushButton *m_stopAnalysisButton;
        QPushButton *m_resetAnalysisButton;
        QProgressBar *m_analysisProgressBar;

        // Results tab components
        QLabel *m_totalGenesLabel;
        QLabel *m_significantGenesLabel;
        QLabel *m_upregulatedLabel;
        QLabel *m_downregulatedLabel;
        QTableWidget *m_resultsTable;
        QDoubleSpinBox *m_pValueThresholdSpinBox;
        QDoubleSpinBox *m_log2FCThresholdSpinBox;
        QPushButton *m_applyFilterButton;
        QPushButton *m_clearFilterButton;
        QPushButton *m_exportResultsButton;
        QPushButton *m_saveGeneListsButton;

        // Visualization tab components
        QComboBox *m_plotTypeCombo;
        QComboBox *m_colorSchemeCombo;
        QSpinBox *m_pointSizeSpinBox;
        QChartView *m_plotWidget;
        QPushButton *m_savePlotButton;
        QPushButton *m_exportDataButton;
        QPushButton *m_printPlotButton;

        // Y2H-SCORES Settings components
        QComboBox *m_baitGroupingCombo;
        QDoubleSpinBox *m_enrichPvalSpinBox;
        QDoubleSpinBox *m_enrichFcSpinBox;
        QDoubleSpinBox *m_specPvalSpinBox;

        // In-frame scoring (junction files) components
        QLineEdit *m_junctionFilesEdit;
        QPushButton *m_browseJunctionBtn;

        // Contrast selector for three-way comparisons
        QComboBox *m_contrastSelectorCombo;

        // Plot methods
        void plotMAPlot();
        void plotVolcanoPlot();
        void plotDispersionPlot();
        void plotThreeWayScatter();
        void refreshPlot();

        // Progress output components
        QTextEdit *m_progressOutputText;
        QPushButton *m_clearOutputButton;
        QPushButton *m_saveOutputButton;

        // Menu actions
        QAction *m_actionOpenCounts;
        QAction *m_actionOpenMetadata;
        QAction *m_actionSaveResults;
        QAction *m_actionExportVisualizations;
        QAction *m_actionExit;
        QAction *m_actionRunAnalysis;
        QAction *m_actionStopAnalysis;
        QAction *m_actionResetAnalysis;
        QAction *m_actionInputTab;
        QAction *m_actionAnalysisTab;
        QAction *m_actionResultsTab;
        QAction *m_actionVisualizationTab;
        QAction *m_actionAbout;
        QAction *m_actionUserGuide;

        // Data
        std::unique_ptr<GeneCountHandler> m_geneCountHandler;
        CombinedData m_combinedData;
        AnalysisResults m_analysisResults;

        // Analysis components
        std::unique_ptr<QThread> m_analysisThread;
        std::unique_ptr<AnalysisWorker> m_analysisWorker;

        // File dialog state
        QString m_lastUsedDirectory;
        QString m_workingDirectory;
        QString m_resultsSqlitePath;
        QString m_lastPlotCacheKey;
        bool m_resultsColumnsSized = false;
        quint64 m_resultsRevision = 0;

        // SQLite output
        void writeResultsToSqlite(const AnalysisResults &results);
        bool exportSqliteToCSV(const QString &csvPath, bool significantOnly = false);

        std::unique_ptr<Ui::MainWindow> m_ui;
    };

} // namespace deseq2

Q_DECLARE_METATYPE(deseq2::AnalysisResults)

#endif // MAIN_WINDOW_H

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "boundary_detector.h"
#include "depth_plot_widget.h"

#include "data_structures.h"
#include "depth_calculator.h"
#include "export_engine.h"
#include "gene_annotation_db.h"
#include "gene_selector_widget.h"
#include "mrna_track_widget.h"
#include "sync_controller.h"

#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QMainWindow>
#include <QProgressBar>
#include <QRadioButton>
#include <QSpinBox>
#include <QSplitter>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

    void loadWorkingDirectory(const QString& workdir);
    void loadDataset(const QString& dbPath);
    void loadGeneReference(const QString& fastaPath);
    void loadDESeq2Results(const QString& csvPath);
    void loadJunctionData(const QString& csvPath);
    void displayGene(const QString& gene);

private slots:
    void onGeneSelected(const QString& gene);
    void onIntervalWidthChanged(int value);
    void onIntervalSpacingChanged(int value);
    void onCompareToggled(bool checked);
    void onOverlayToggled(bool checked);
    void onSplitToggled(bool checked);
    void onDatasetChanged(int index);
    void onSecondaryDatasetChanged(int index);
    void onCursorMoved(int position, int depth);
    void onRangeChanged(qreal xMin, qreal xMax);
    void onExportCSV();
    void onExportFigure();
    void onBatchRun();
    void onResetZoom();

private:
    Ui::MainWindow* ui;

    // Shared library objects
    deepn::GeneAnnotationDB m_annotationDB;
    deepn::DepthCalculator m_depthCalc;
    deepn::DepthProfile m_primaryProfile;
    deepn::DepthProfile m_secondaryProfile;
    deepn::BoundaryResult m_primaryBoundary;
    deepn::BoundaryResult m_secondaryBoundary;

    // Junction data from MultiQuery++ (optional)
    QMap<QString, QVector<deepn::JunctionSite>> m_junctionData;

    // DESeq2 results
    QVector<deepn::DESeq2Result> m_deseq2Results;

    // UI widgets -- top bar
    deepn::GeneSelectorWidget* m_geneSelector = nullptr;
    QComboBox* m_datasetCombo = nullptr;
    QComboBox* m_secondaryDatasetCombo = nullptr;
    QCheckBox* m_compareCheck = nullptr;
    QRadioButton* m_splitRadio = nullptr;
    QRadioButton* m_overlayRadio = nullptr;
    QWidget* m_comparisonPanel = nullptr;
    QSpinBox* m_intervalWidthSpin = nullptr;
    QSpinBox* m_intervalSpacingSpin = nullptr;

    // UI widgets -- main area
    QSplitter* m_mainSplitter = nullptr;
    QSplitter* m_chartSplitter = nullptr;     // horizontal: primary | secondary
    QWidget* m_primaryPanel = nullptr;
    QWidget* m_secondaryPanel = nullptr;
    DepthPlotWidget* m_primaryPlot = nullptr;
    DepthPlotWidget* m_secondaryPlot = nullptr;
    deepn::MRNATrackWidget* m_primaryTrack = nullptr;
    deepn::MRNATrackWidget* m_secondaryTrack = nullptr;

    // Boundary info
    QLabel* m_boundaryLabel = nullptr;
    QLabel* m_insertExtentLabel = nullptr;

    // Status bar
    QLabel* m_cursorLabel = nullptr;
    QLabel* m_zoomLabel = nullptr;

    // Sync controller for split mode
    deepn::SyncController* m_syncController = nullptr;

    // State
    QStringList m_datasetPaths;
    QString m_workdir;
    int m_intervalWidth = 25;
    int m_intervalSpacing = 50;
    QString m_currentGene;

    void setupUI();
    void setupMenuBar();
    void setupStatusBar();
    void connectSignals();
    void refreshDisplay();
    void updateBoundaryInfo();
    void updateComparisonMode();
    void autoDiscoverFiles();
    void populateDatasetCombos();
    QVector<deepn::DESeq2Result> parseDESeq2CSV(const QString& csvPath);
    QMap<QString, QVector<deepn::JunctionSite>> parseJunctionCSV(const QString& csvPath);
};

#endif  // MAINWINDOW_H

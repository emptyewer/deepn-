#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QModelIndex>
#include <QVector>

#include <data_structures.h>
#include <gene_annotation_db.h>
#include <sqlite_junction_loader.h>

class QCheckBox;
class QComboBox;
class QLabel;
class QSplitter;
class QTableView;
class QToolButton;
class QProgressBar;

namespace deepn {
class GeneSelectorWidget;
class MRNATrackWidget;
}

class JunctionPlotWidget;
class JunctionTableModel;
class GeneDetailPanel;
class ComparisonManager;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

    void loadWorkingDirectory(const QString& workdir);
    void loadDataset(const QString& dbPath);
    void loadGeneReference(const QString& fastaPath);
    void loadDESeq2Results(const QString& csvPath);
    void displayGene(const QString& gene);

private slots:
    void onGeneSelected(const QString& gene);
    void onJunctionClicked(int position);
    void onTableRowSelected(const QModelIndex& current, const QModelIndex& previous);
    void onCollapseToggled(bool checked);
    void onInFrameOnlyToggled(bool checked);
    void onCompareToggled(bool checked);
    void onDatasetChanged(int index);
    void onSecondaryDatasetChanged(int index);
    void onExportCSV();
    void onExportFigure();
    void onBatchRun();

private:
    // Shared library components
    deepn::GeneAnnotationDB m_annotationDB;
    deepn::SqliteJunctionLoader m_primaryLoader;
    deepn::SqliteJunctionLoader m_secondaryLoader;
    deepn::GeneJunctionProfile m_currentProfile;
    deepn::GeneJunctionProfile m_secondaryProfile;

    // UI components
    Ui::MainWindow* ui;
    deepn::GeneSelectorWidget* m_geneSelector = nullptr;
    JunctionPlotWidget* m_primaryPlot = nullptr;
    JunctionPlotWidget* m_secondaryPlot = nullptr;
    JunctionTableModel* m_tableModel = nullptr;
    QTableView* m_tableView = nullptr;
    GeneDetailPanel* m_detailPanel = nullptr;
    deepn::MRNATrackWidget* m_primaryTrack = nullptr;
    deepn::MRNATrackWidget* m_secondaryTrack = nullptr;
    ComparisonManager* m_comparisonMgr = nullptr;

    // Layout containers
    QSplitter* m_mainSplitter = nullptr;
    QSplitter* m_leftSplitter = nullptr;
    QWidget* m_primaryChartContainer = nullptr;
    QWidget* m_secondaryChartContainer = nullptr;
    QSplitter* m_chartSplitter = nullptr;

    // Controls
    QComboBox* m_datasetCombo = nullptr;
    QComboBox* m_secondaryDatasetCombo = nullptr;
    QCheckBox* m_collapseCheck = nullptr;
    QCheckBox* m_inFrameOnlyCheck = nullptr;
    QCheckBox* m_compareCheck = nullptr;
    QLabel* m_secondaryLabel = nullptr;
    QProgressBar* m_batchProgress = nullptr;

    // State
    QStringList m_datasetPaths;
    QString m_workdir;
    QVector<deepn::DESeq2Result> m_deseq2Results;
    bool m_collapseByPosition = false;
    bool m_inFrameOnly = false;

    void setupUI();
    void setupMenuBar();
    void setupToolBar();
    void connectSignals();
    void refreshDisplay();
    void refreshSecondaryDisplay();
    void autoDiscoverFiles();
    void populateDatasetCombo();
    void applyFilters(const deepn::GeneJunctionProfile& profile,
                      QVector<deepn::JunctionSite>& filtered) const;
    QVector<deepn::DESeq2Result> parseDESeq2CSV(const QString& path) const;
};

#endif // MAINWINDOW_H

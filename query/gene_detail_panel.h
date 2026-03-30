#ifndef GENE_DETAIL_PANEL_H
#define GENE_DETAIL_PANEL_H

#include <data_structures.h>

#include <QWidget>

class QGroupBox;
class QLabel;
class QTextEdit;

class GeneDetailPanel : public QWidget
{
    Q_OBJECT

public:
    explicit GeneDetailPanel(QWidget* parent = nullptr);

    void setAnnotation(const deepn::GeneAnnotation& annotation);
    void setSelectedJunction(const deepn::JunctionSite& site);
    void setProfileSummary(int totalJunctions, int inFrameCount, double maxPpm);

private:
    // Gene info labels
    QLabel* m_geneNameLabel = nullptr;
    QLabel* m_refseqLabel = nullptr;
    QLabel* m_mrnaLengthLabel = nullptr;
    QLabel* m_cdsRangeLabel = nullptr;
    QLabel* m_chromosomeLabel = nullptr;

    // Profile summary labels
    QLabel* m_totalJunctionsLabel = nullptr;
    QLabel* m_inFrameCountLabel = nullptr;
    QLabel* m_maxPpmLabel = nullptr;

    // Selected junction detail labels
    QGroupBox* m_junctionGroup = nullptr;
    QLabel* m_jPositionLabel = nullptr;
    QLabel* m_jPpmLabel = nullptr;
    QLabel* m_jFrameLabel = nullptr;
    QLabel* m_jCdsLabel = nullptr;
    QLabel* m_jQueryStartLabel = nullptr;
    QLabel* m_jRawCountLabel = nullptr;
    QLabel* m_jCategoryLabel = nullptr;
    QTextEdit* m_jSequenceEdit = nullptr;

    // Stored state for sequence context
    deepn::GeneAnnotation m_currentAnnotation;

    void setupUI();
    QString extractSequenceContext(int position, int window = 30) const;
};

#endif // GENE_DETAIL_PANEL_H

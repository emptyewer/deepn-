#ifndef DEEPN_GENE_SELECTOR_WIDGET_H
#define DEEPN_GENE_SELECTOR_WIDGET_H

#include "data_structures.h"

#include <QCheckBox>
#include <QComboBox>
#include <QCompleter>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QStringList>
#include <QVector>
#include <QWidget>

namespace deepn {

class GeneSelectorWidget : public QWidget {
    Q_OBJECT
public:
    explicit GeneSelectorWidget(QWidget* parent = nullptr);

    // Set the gene list (from GeneAnnotationDB or DESeq2 results)
    void setGeneList(const QStringList& genes);
    void setDESeq2Results(const QVector<DESeq2Result>& results);

    // Get current selection
    QString currentGene() const;
    int currentIndex() const;

    // Navigation
    bool selectGene(const QString& gene);
    void goToNext();
    void goToPrev();

signals:
    void geneSelected(const QString& gene);
    void filterChanged();

private:
    QLineEdit* m_searchEdit;
    QCompleter* m_completer;
    QPushButton* m_prevBtn;
    QPushButton* m_nextBtn;
    QLabel* m_infoLabel;
    QComboBox* m_sortCombo;
    QCheckBox* m_filterCheck;
    QDoubleSpinBox* m_padjSpin;
    QDoubleSpinBox* m_l2fcSpin;
    QStringList m_geneList;            // currently displayed (after filter+sort)
    QStringList m_unfilteredGeneList;  // full set before filtering
    QVector<DESeq2Result> m_deseq2Results;
    QVector<DESeq2Result> m_unfilteredDeseq2;
    int m_currentIndex = -1;

    void setupUI();
    void updateNavigation();
    void onSearchTextChanged(const QString& text);
    void onSortChanged(int index);
    void onFilterChanged();
    void applySorting();
    void applyFilters();
};

}  // namespace deepn

#endif  // DEEPN_GENE_SELECTOR_WIDGET_H

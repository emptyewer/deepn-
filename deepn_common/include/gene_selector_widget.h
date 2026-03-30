#ifndef DEEPN_GENE_SELECTOR_WIDGET_H
#define DEEPN_GENE_SELECTOR_WIDGET_H

#include "data_structures.h"

#include <QComboBox>
#include <QCompleter>
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
    void selectGene(const QString& gene);
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
    QStringList m_geneList;
    QVector<DESeq2Result> m_deseq2Results;
    int m_currentIndex = -1;

    void setupUI();
    void updateNavigation();
    void onSearchTextChanged(const QString& text);
    void onSortChanged(int index);
    void applySorting();
};

}  // namespace deepn

#endif  // DEEPN_GENE_SELECTOR_WIDGET_H

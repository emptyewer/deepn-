#include "gene_selector_widget.h"

#include <QDebug>
#include <QHBoxLayout>
#include <QStringListModel>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>

namespace deepn {

GeneSelectorWidget::GeneSelectorWidget(QWidget* parent)
    : QWidget(parent)
{
    setupUI();
}

void GeneSelectorWidget::setupUI()
{
    auto* outerLayout = new QVBoxLayout(this);
    outerLayout->setContentsMargins(4, 2, 4, 2);
    outerLayout->setSpacing(2);

    // Row 1: search + navigation + sort
    auto* row1 = new QHBoxLayout;
    row1->setSpacing(6);

    auto* searchLabel = new QLabel("Search:", this);
    row1->addWidget(searchLabel);

    m_searchEdit = new QLineEdit(this);
    m_searchEdit->setPlaceholderText("Type gene name...");
    m_searchEdit->setMinimumWidth(150);
    m_searchEdit->setClearButtonEnabled(true);
    row1->addWidget(m_searchEdit);

    m_completer = new QCompleter(this);
    m_completer->setCaseSensitivity(Qt::CaseInsensitive);
    m_completer->setFilterMode(Qt::MatchContains);
    m_completer->setCompletionMode(QCompleter::PopupCompletion);
    m_searchEdit->setCompleter(m_completer);

    auto* geneLabel = new QLabel("Gene:", this);
    row1->addWidget(geneLabel);

    m_infoLabel = new QLabel("--", this);
    m_infoLabel->setMinimumWidth(100);
    m_infoLabel->setStyleSheet("font-weight: bold;");
    row1->addWidget(m_infoLabel);

    m_prevBtn = new QPushButton("<< Prev", this);
    m_prevBtn->setEnabled(false);
    m_prevBtn->setFixedWidth(80);
    row1->addWidget(m_prevBtn);

    m_nextBtn = new QPushButton("Next >>", this);
    m_nextBtn->setEnabled(false);
    m_nextBtn->setFixedWidth(80);
    row1->addWidget(m_nextBtn);

    auto* sortLabel = new QLabel("Sort:", this);
    row1->addWidget(sortLabel);

    m_sortCombo = new QComboBox(this);
    m_sortCombo->addItem("Name (A-Z)", "name_asc");
    m_sortCombo->addItem("Name (Z-A)", "name_desc");
    m_sortCombo->addItem("padj (ascending)", "padj_asc");
    m_sortCombo->addItem("padj (descending)", "padj_desc");
    m_sortCombo->addItem("log2FC (ascending)", "l2fc_asc");
    m_sortCombo->addItem("log2FC (descending)", "l2fc_desc");
    m_sortCombo->setCurrentIndex(0);
    row1->addWidget(m_sortCombo);

    row1->addStretch();
    outerLayout->addLayout(row1);

    // Row 2: filter controls (visible only when DESeq2 data loaded)
    auto* row2 = new QHBoxLayout;
    row2->setSpacing(6);

    m_filterCheck = new QCheckBox("Filter:", this);
    m_filterCheck->setChecked(false);
    row2->addWidget(m_filterCheck);

    auto* padjLabel = new QLabel("padj \u2264", this);
    row2->addWidget(padjLabel);

    m_padjSpin = new QDoubleSpinBox(this);
    m_padjSpin->setRange(0.0, 1.0);
    m_padjSpin->setDecimals(3);
    m_padjSpin->setSingleStep(0.01);
    m_padjSpin->setValue(0.05);
    m_padjSpin->setEnabled(false);
    m_padjSpin->setFixedWidth(90);
    row2->addWidget(m_padjSpin);

    auto* l2fcLabel = new QLabel("|log2FC| \u2265", this);
    row2->addWidget(l2fcLabel);

    m_l2fcSpin = new QDoubleSpinBox(this);
    m_l2fcSpin->setRange(0.0, 20.0);
    m_l2fcSpin->setDecimals(2);
    m_l2fcSpin->setSingleStep(0.5);
    m_l2fcSpin->setValue(1.0);
    m_l2fcSpin->setEnabled(false);
    m_l2fcSpin->setFixedWidth(90);
    row2->addWidget(m_l2fcSpin);

    row2->addStretch();
    outerLayout->addLayout(row2);

    // Connections
    connect(m_searchEdit, &QLineEdit::textChanged, this, &GeneSelectorWidget::onSearchTextChanged);
    connect(m_searchEdit, &QLineEdit::returnPressed, this, [this]() {
        QString text = m_searchEdit->text().trimmed();
        if (!text.isEmpty()) {
            selectGene(text);
        }
    });

    connect(m_completer, QOverload<const QString&>::of(&QCompleter::activated),
            this, [this](const QString& text) {
                selectGene(text);
            });

    connect(m_prevBtn, &QPushButton::clicked, this, &GeneSelectorWidget::goToPrev);
    connect(m_nextBtn, &QPushButton::clicked, this, &GeneSelectorWidget::goToNext);
    connect(m_sortCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &GeneSelectorWidget::onSortChanged);

    connect(m_filterCheck, &QCheckBox::toggled, this, [this](bool checked) {
        m_padjSpin->setEnabled(checked);
        m_l2fcSpin->setEnabled(checked);
        onFilterChanged();
    });
    connect(m_padjSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, [this](double) { if (m_filterCheck->isChecked()) onFilterChanged(); });
    connect(m_l2fcSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, [this](double) { if (m_filterCheck->isChecked()) onFilterChanged(); });
}

void GeneSelectorWidget::setGeneList(const QStringList& genes)
{
    m_unfilteredGeneList = genes;
    m_unfilteredDeseq2.clear();
    m_deseq2Results.clear();
    m_geneList = genes;
    m_currentIndex = -1;

    applySorting();

    auto* model = new QStringListModel(m_geneList, m_completer);
    m_completer->setModel(model);

    updateNavigation();
}

void GeneSelectorWidget::setDESeq2Results(const QVector<DESeq2Result>& results)
{
    m_unfilteredDeseq2 = results;
    m_unfilteredGeneList.clear();
    for (const auto& r : results) {
        m_unfilteredGeneList.append(r.gene);
    }

    // Apply filters if active, otherwise use full set
    applyFilters();
    applySorting();

    auto* model = new QStringListModel(m_geneList, m_completer);
    m_completer->setModel(model);

    updateNavigation();
}

QString GeneSelectorWidget::currentGene() const
{
    if (m_currentIndex >= 0 && m_currentIndex < m_geneList.size())
        return m_geneList[m_currentIndex];
    return {};
}

int GeneSelectorWidget::currentIndex() const
{
    return m_currentIndex;
}

bool GeneSelectorWidget::selectGene(const QString& gene)
{
    int idx = m_geneList.indexOf(gene);
    if (idx < 0) {
        // Try case-insensitive search
        for (int i = 0; i < m_geneList.size(); i++) {
            if (m_geneList[i].compare(gene, Qt::CaseInsensitive) == 0) {
                idx = i;
                break;
            }
        }
    }

    if (idx >= 0) {
        m_currentIndex = idx;
        updateNavigation();
        emit geneSelected(m_geneList[m_currentIndex]);
        return true;
    } else {
        qDebug() << "Gene not found:" << gene;
        return false;
    }
}

void GeneSelectorWidget::goToNext()
{
    if (m_geneList.isEmpty())
        return;

    if (m_currentIndex < m_geneList.size() - 1) {
        m_currentIndex++;
    } else {
        m_currentIndex = 0;  // wrap around
    }

    updateNavigation();
    emit geneSelected(m_geneList[m_currentIndex]);
}

void GeneSelectorWidget::goToPrev()
{
    if (m_geneList.isEmpty())
        return;

    if (m_currentIndex > 0) {
        m_currentIndex--;
    } else {
        m_currentIndex = m_geneList.size() - 1;  // wrap around
    }

    updateNavigation();
    emit geneSelected(m_geneList[m_currentIndex]);
}

void GeneSelectorWidget::updateNavigation()
{
    bool hasGenes = !m_geneList.isEmpty();
    m_prevBtn->setEnabled(hasGenes);
    m_nextBtn->setEnabled(hasGenes);

    if (m_currentIndex >= 0 && m_currentIndex < m_geneList.size()) {
        QString gene = m_geneList[m_currentIndex];
        m_infoLabel->setText(QStringLiteral("%1  (%2 of %3)")
                                 .arg(gene)
                                 .arg(m_currentIndex + 1)
                                 .arg(m_geneList.size()));

        // Don't update search text while user is typing
        if (!m_searchEdit->hasFocus()) {
            m_searchEdit->setText(gene);
        }
    } else if (hasGenes) {
        m_infoLabel->setText(QStringLiteral("-- (%1 genes)").arg(m_geneList.size()));
    } else {
        m_infoLabel->setText("--");
    }
}

void GeneSelectorWidget::onSearchTextChanged(const QString& text)
{
    Q_UNUSED(text)
    // Let the completer handle popup display; actual selection happens on
    // returnPressed or completer activated
}

void GeneSelectorWidget::onSortChanged(int index)
{
    Q_UNUSED(index)
    QString currentGene;
    if (m_currentIndex >= 0 && m_currentIndex < m_geneList.size())
        currentGene = m_geneList[m_currentIndex];

    applySorting();

    // Re-select current gene after sort
    if (!currentGene.isEmpty()) {
        m_currentIndex = m_geneList.indexOf(currentGene);
    }

    auto* model = new QStringListModel(m_geneList, m_completer);
    m_completer->setModel(model);

    updateNavigation();
    emit filterChanged();
}

void GeneSelectorWidget::applySorting()
{
    QString sortKey = m_sortCombo->currentData().toString();

    if (m_deseq2Results.isEmpty()) {
        // Simple name sorting
        if (sortKey == "name_asc") {
            std::sort(m_geneList.begin(), m_geneList.end());
        } else if (sortKey == "name_desc") {
            std::sort(m_geneList.begin(), m_geneList.end(),
                      [](const QString& a, const QString& b) { return a > b; });
        }
        // Other sort modes not applicable without DESeq2 data
        return;
    }

    // Sort DESeq2 results and rebuild gene list
    auto sorted = m_deseq2Results;

    if (sortKey == "name_asc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.gene < b.gene; });
    } else if (sortKey == "name_desc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.gene > b.gene; });
    } else if (sortKey == "padj_asc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.padj < b.padj; });
    } else if (sortKey == "padj_desc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.padj > b.padj; });
    } else if (sortKey == "l2fc_asc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.log2FoldChange < b.log2FoldChange; });
    } else if (sortKey == "l2fc_desc") {
        std::sort(sorted.begin(), sorted.end(),
                  [](const DESeq2Result& a, const DESeq2Result& b) { return a.log2FoldChange > b.log2FoldChange; });
    }

    m_deseq2Results = sorted;
    m_geneList.clear();
    for (const auto& r : sorted) {
        m_geneList.append(r.gene);
    }
}

void GeneSelectorWidget::onFilterChanged()
{
    QString currentGene;
    if (m_currentIndex >= 0 && m_currentIndex < m_geneList.size())
        currentGene = m_geneList[m_currentIndex];

    applyFilters();
    applySorting();

    auto* model = new QStringListModel(m_geneList, m_completer);
    m_completer->setModel(model);

    // Try to re-select current gene
    if (!currentGene.isEmpty()) {
        m_currentIndex = m_geneList.indexOf(currentGene);
        if (m_currentIndex < 0 && !m_geneList.isEmpty()) {
            m_currentIndex = 0;
        }
    } else if (!m_geneList.isEmpty()) {
        m_currentIndex = 0;
    }

    updateNavigation();
    if (m_currentIndex >= 0 && m_currentIndex < m_geneList.size()) {
        emit geneSelected(m_geneList[m_currentIndex]);
    }
    emit filterChanged();
}

void GeneSelectorWidget::applyFilters()
{
    if (m_unfilteredDeseq2.isEmpty()) {
        // No DESeq2 data — just copy the unfiltered gene list
        m_deseq2Results.clear();
        m_geneList = m_unfilteredGeneList;
        m_currentIndex = -1;
        return;
    }

    if (!m_filterCheck->isChecked()) {
        // Filters disabled — use full set
        m_deseq2Results = m_unfilteredDeseq2;
        m_geneList.clear();
        for (const auto& r : m_deseq2Results) {
            m_geneList.append(r.gene);
        }
        m_currentIndex = -1;
        return;
    }

    // Apply p-value and log2FC filters
    const double padjCutoff = m_padjSpin->value();
    const double l2fcCutoff = m_l2fcSpin->value();

    m_deseq2Results.clear();
    m_geneList.clear();
    m_currentIndex = -1;

    for (const auto& r : m_unfilteredDeseq2) {
        if (r.padj <= padjCutoff && std::abs(r.log2FoldChange) >= l2fcCutoff) {
            m_deseq2Results.append(r);
            m_geneList.append(r.gene);
        }
    }
}

}  // namespace deepn

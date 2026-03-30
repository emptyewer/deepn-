#include "gene_selector_widget.h"

#include <QDebug>
#include <QHBoxLayout>
#include <QStringListModel>

#include <algorithm>

namespace deepn {

GeneSelectorWidget::GeneSelectorWidget(QWidget* parent)
    : QWidget(parent)
{
    setupUI();
}

void GeneSelectorWidget::setupUI()
{
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(4, 2, 4, 2);
    layout->setSpacing(6);

    // Search label and edit
    auto* searchLabel = new QLabel("Search:", this);
    layout->addWidget(searchLabel);

    m_searchEdit = new QLineEdit(this);
    m_searchEdit->setPlaceholderText("Type gene name...");
    m_searchEdit->setMinimumWidth(150);
    m_searchEdit->setClearButtonEnabled(true);
    layout->addWidget(m_searchEdit);

    // Completer
    m_completer = new QCompleter(this);
    m_completer->setCaseSensitivity(Qt::CaseInsensitive);
    m_completer->setFilterMode(Qt::MatchContains);
    m_completer->setCompletionMode(QCompleter::PopupCompletion);
    m_searchEdit->setCompleter(m_completer);

    // Gene info label
    auto* geneLabel = new QLabel("Gene:", this);
    layout->addWidget(geneLabel);

    m_infoLabel = new QLabel("--", this);
    m_infoLabel->setMinimumWidth(100);
    m_infoLabel->setStyleSheet("font-weight: bold;");
    layout->addWidget(m_infoLabel);

    // Navigation buttons
    m_prevBtn = new QPushButton("<< Prev", this);
    m_prevBtn->setEnabled(false);
    m_prevBtn->setFixedWidth(80);
    layout->addWidget(m_prevBtn);

    m_nextBtn = new QPushButton("Next >>", this);
    m_nextBtn->setEnabled(false);
    m_nextBtn->setFixedWidth(80);
    layout->addWidget(m_nextBtn);

    // Sort combo
    auto* sortLabel = new QLabel("Sort:", this);
    layout->addWidget(sortLabel);

    m_sortCombo = new QComboBox(this);
    m_sortCombo->addItem("Name (A-Z)", "name_asc");
    m_sortCombo->addItem("Name (Z-A)", "name_desc");
    m_sortCombo->addItem("padj (ascending)", "padj_asc");
    m_sortCombo->addItem("padj (descending)", "padj_desc");
    m_sortCombo->addItem("log2FC (ascending)", "l2fc_asc");
    m_sortCombo->addItem("log2FC (descending)", "l2fc_desc");
    m_sortCombo->setCurrentIndex(0);
    layout->addWidget(m_sortCombo);

    layout->addStretch();

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
}

void GeneSelectorWidget::setGeneList(const QStringList& genes)
{
    m_geneList = genes;
    m_deseq2Results.clear();
    m_currentIndex = -1;

    applySorting();

    auto* model = new QStringListModel(m_geneList, m_completer);
    m_completer->setModel(model);

    updateNavigation();
}

void GeneSelectorWidget::setDESeq2Results(const QVector<DESeq2Result>& results)
{
    m_deseq2Results = results;
    m_geneList.clear();
    m_currentIndex = -1;

    for (const auto& r : results) {
        m_geneList.append(r.gene);
    }

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

void GeneSelectorWidget::selectGene(const QString& gene)
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
    } else {
        qDebug() << "Gene not found:" << gene;
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

}  // namespace deepn

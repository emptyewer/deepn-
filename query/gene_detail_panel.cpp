#include "gene_detail_panel.h"

#include <QFont>
#include <QGroupBox>
#include <QLabel>
#include <QScrollArea>
#include <QTextEdit>
#include <QVBoxLayout>

using namespace deepn;

// ────────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────────

GeneDetailPanel::GeneDetailPanel(QWidget* parent)
    : QWidget(parent)
{
    setupUI();
}

void GeneDetailPanel::setupUI()
{
    auto* scrollArea = new QScrollArea(this);
    scrollArea->setWidgetResizable(true);
    scrollArea->setFrameShape(QFrame::NoFrame);

    auto* scrollContent = new QWidget;
    auto* layout = new QVBoxLayout(scrollContent);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(8);

    // ── Gene Information Group ───────────────────────────────────
    auto* geneGroup = new QGroupBox("Gene Information", scrollContent);
    auto* geneLayout = new QVBoxLayout(geneGroup);
    geneLayout->setSpacing(4);

    m_geneNameLabel = new QLabel("Gene: --");
    QFont nameFont = m_geneNameLabel->font();
    nameFont.setBold(true);
    nameFont.setPointSize(nameFont.pointSize() + 2);
    m_geneNameLabel->setFont(nameFont);
    geneLayout->addWidget(m_geneNameLabel);

    m_refseqLabel = new QLabel("RefSeq: --");
    geneLayout->addWidget(m_refseqLabel);

    m_chromosomeLabel = new QLabel("Chromosome: --");
    geneLayout->addWidget(m_chromosomeLabel);

    m_mrnaLengthLabel = new QLabel("mRNA Length: --");
    geneLayout->addWidget(m_mrnaLengthLabel);

    m_cdsRangeLabel = new QLabel("CDS: --");
    geneLayout->addWidget(m_cdsRangeLabel);

    layout->addWidget(geneGroup);

    // ── Profile Summary Group ────────────────────────────────────
    auto* summaryGroup = new QGroupBox("Profile Summary", scrollContent);
    auto* summaryLayout = new QVBoxLayout(summaryGroup);
    summaryLayout->setSpacing(4);

    m_totalJunctionsLabel = new QLabel("Total Junctions: --");
    summaryLayout->addWidget(m_totalJunctionsLabel);

    m_inFrameCountLabel = new QLabel("In-Frame: --");
    summaryLayout->addWidget(m_inFrameCountLabel);

    m_maxPpmLabel = new QLabel("Max PPM: --");
    summaryLayout->addWidget(m_maxPpmLabel);

    layout->addWidget(summaryGroup);

    // ── Selected Junction Group ──────────────────────────────────
    m_junctionGroup = new QGroupBox("Selected Junction", scrollContent);
    auto* junctionLayout = new QVBoxLayout(m_junctionGroup);
    junctionLayout->setSpacing(4);

    m_jPositionLabel = new QLabel("Position: --");
    QFont posFont = m_jPositionLabel->font();
    posFont.setBold(true);
    m_jPositionLabel->setFont(posFont);
    junctionLayout->addWidget(m_jPositionLabel);

    m_jPpmLabel = new QLabel("PPM: --");
    junctionLayout->addWidget(m_jPpmLabel);

    m_jFrameLabel = new QLabel("Frame: --");
    junctionLayout->addWidget(m_jFrameLabel);

    m_jCdsLabel = new QLabel("CDS: --");
    junctionLayout->addWidget(m_jCdsLabel);

    m_jQueryStartLabel = new QLabel("QueryStart: --");
    junctionLayout->addWidget(m_jQueryStartLabel);

    m_jRawCountLabel = new QLabel("Raw Count: --");
    junctionLayout->addWidget(m_jRawCountLabel);

    m_jCategoryLabel = new QLabel("Category: --");
    junctionLayout->addWidget(m_jCategoryLabel);

    // Sequence context
    auto* seqLabel = new QLabel("Sequence Context:");
    junctionLayout->addWidget(seqLabel);

    m_jSequenceEdit = new QTextEdit;
    m_jSequenceEdit->setReadOnly(true);
    m_jSequenceEdit->setMaximumHeight(100);
    QFont monoFont("Menlo", 10);
    monoFont.setStyleHint(QFont::Monospace);
    m_jSequenceEdit->setFont(monoFont);
    m_jSequenceEdit->setLineWrapMode(QTextEdit::WidgetWidth);
    junctionLayout->addWidget(m_jSequenceEdit);

    m_junctionGroup->setEnabled(false);
    layout->addWidget(m_junctionGroup);

    layout->addStretch(1);

    scrollArea->setWidget(scrollContent);

    auto* outerLayout = new QVBoxLayout(this);
    outerLayout->setContentsMargins(0, 0, 0, 0);
    outerLayout->addWidget(scrollArea);
}

// ────────────────────────────────────────────────────────────────────
// Public API
// ────────────────────────────────────────────────────────────────────

void GeneDetailPanel::setAnnotation(const GeneAnnotation& annotation)
{
    m_currentAnnotation = annotation;

    m_geneNameLabel->setText(QString("Gene: %1").arg(
        annotation.geneName.isEmpty() ? "--" : annotation.geneName));

    m_refseqLabel->setText(QString("RefSeq: %1").arg(
        annotation.refseq.isEmpty() ? "--" : annotation.refseq));

    m_chromosomeLabel->setText(QString("Chromosome: %1").arg(
        annotation.chromosome.isEmpty() ? "--" : annotation.chromosome));

    if (annotation.mRNALength > 0) {
        m_mrnaLengthLabel->setText(QString("mRNA Length: %L1 nt").arg(annotation.mRNALength));
    } else {
        m_mrnaLengthLabel->setText("mRNA Length: --");
    }

    if (annotation.orfStart > 0 || annotation.orfEnd > 0) {
        m_cdsRangeLabel->setText(QString("CDS: %1 - %2")
                                     .arg(annotation.orfStart)
                                     .arg(annotation.orfEnd));
    } else {
        m_cdsRangeLabel->setText("CDS: --");
    }
}

void GeneDetailPanel::setSelectedJunction(const JunctionSite& site)
{
    if (site.position == 0 && site.ppm == 0.0) {
        // Deselect
        m_junctionGroup->setEnabled(false);
        m_jPositionLabel->setText("Position: --");
        m_jPpmLabel->setText("PPM: --");
        m_jFrameLabel->setText("Frame: --");
        m_jCdsLabel->setText("CDS: --");
        m_jQueryStartLabel->setText("QueryStart: --");
        m_jRawCountLabel->setText("Raw Count: --");
        m_jCategoryLabel->setText("Category: --");
        m_jSequenceEdit->clear();
        return;
    }

    m_junctionGroup->setEnabled(true);

    m_jPositionLabel->setText(QString("Position: %1").arg(site.position));
    m_jPpmLabel->setText(QString("PPM: %1").arg(site.ppm, 0, 'f', 2));
    m_jFrameLabel->setText(QString("Frame: %1").arg(site.frameLabel()));
    m_jCdsLabel->setText(QString("CDS: %1").arg(site.cdsLabel()));
    m_jQueryStartLabel->setText(QString("QueryStart: %1").arg(site.queryStart));
    m_jRawCountLabel->setText(QString("Raw Count: %1").arg(site.rawCount));

    // Category with color
    JunctionCategory cat = classifyJunction(site);
    QString catName;
    QString catColor;
    switch (cat) {
    case JunctionCategory::InOrfInFrame:
        catName = "In ORF + In-frame";
        catColor = "#1E40AF";
        break;
    case JunctionCategory::UpstreamInFrame:
        catName = "Upstream + In-frame";
        catColor = "#0D9488";
        break;
    case JunctionCategory::InOrfOutOfFrame:
        catName = "In ORF + Out-of-frame";
        catColor = "#D97706";
        break;
    case JunctionCategory::DownstreamOrBack:
        catName = "Downstream / Backwards";
        catColor = "#6B7280";
        break;
    }
    m_jCategoryLabel->setText(QString("<span style='color:%1; font-weight:bold;'>%2</span>")
                                  .arg(catColor, catName));
    m_jCategoryLabel->setTextFormat(Qt::RichText);

    // Sequence context
    QString seqContext = extractSequenceContext(site.position);
    if (!seqContext.isEmpty()) {
        m_jSequenceEdit->setPlainText(seqContext);
    } else {
        m_jSequenceEdit->setPlainText("(sequence not available)");
    }
}

void GeneDetailPanel::setProfileSummary(int totalJunctions, int inFrameCount, double maxPpm)
{
    m_totalJunctionsLabel->setText(QString("Total Junctions: %1").arg(totalJunctions));
    m_inFrameCountLabel->setText(QString("In-Frame: %1 (%2%)")
                                     .arg(inFrameCount)
                                     .arg(totalJunctions > 0
                                              ? QString::number(100.0 * inFrameCount / totalJunctions, 'f', 1)
                                              : "0"));
    m_maxPpmLabel->setText(QString("Max PPM: %1").arg(maxPpm, 0, 'f', 2));
}

// ────────────────────────────────────────────────────────────────────
// Helpers
// ────────────────────────────────────────────────────────────────────

QString GeneDetailPanel::extractSequenceContext(int position, int window) const
{
    if (m_currentAnnotation.sequence.isEmpty()) return {};
    if (position < 0 || position >= m_currentAnnotation.sequence.length()) return {};

    int start = qMax(0, position - window);
    int end = qMin(m_currentAnnotation.sequence.length(), position + window + 1);

    QString before = m_currentAnnotation.sequence.mid(start, position - start);
    QString junction = m_currentAnnotation.sequence.mid(position, 1);
    QString after = m_currentAnnotation.sequence.mid(position + 1, end - position - 1);

    // Format: lowercase flanking, uppercase junction site
    return before.toLower() + "[" + junction.toUpper() + "]" + after.toLower();
}

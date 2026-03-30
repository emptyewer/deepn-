#include "junction_plot_widget.h"

#include <QAreaSeries>
#include <QChart>
#include <QGraphicsScene>
#include <QGraphicsSimpleTextItem>
#include <QLegendMarker>
#include <QLineSeries>
#include <QMouseEvent>
#include <QToolTip>
#include <QValueAxis>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

using namespace deepn;

// ────────────────────────────────────────────────────────────────────
// Color definitions from the plan
// ────────────────────────────────────────────────────────────────────

QColor JunctionPlotWidget::colorForCategory(JunctionCategory cat)
{
    switch (cat) {
    case JunctionCategory::InOrfInFrame:     return QColor("#1E40AF");  // deep blue
    case JunctionCategory::UpstreamInFrame:  return QColor("#0D9488");  // teal
    case JunctionCategory::InOrfOutOfFrame:  return QColor("#D97706");  // amber
    case JunctionCategory::DownstreamOrBack: return QColor("#6B7280");  // grey
    }
    return QColor("#6B7280");
}

static QString categoryLabel(JunctionCategory cat)
{
    switch (cat) {
    case JunctionCategory::InOrfInFrame:     return "In ORF + In-frame";
    case JunctionCategory::UpstreamInFrame:  return "Upstream + In-frame";
    case JunctionCategory::InOrfOutOfFrame:  return "In ORF + Out-of-frame";
    case JunctionCategory::DownstreamOrBack: return "Downstream / Backwards";
    }
    return "Unknown";
}

// ────────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────────

JunctionPlotWidget::JunctionPlotWidget(QWidget* parent)
    : QChartView(new QChart, parent)
{
    setRenderHint(QPainter::Antialiasing);
    setRubberBand(QChartView::NoRubberBand);
    setDragMode(QGraphicsView::NoDrag);

    chart()->setTheme(QChart::ChartThemeLight);
    chart()->setAnimationOptions(QChart::NoAnimation);
    chart()->legend()->setVisible(true);
    chart()->legend()->setAlignment(Qt::AlignBottom);
    chart()->setMargins(QMargins(8, 8, 8, 8));

    // Create axes
    m_xAxis = new QValueAxis(this);
    m_xAxis->setTitleText("Position (nt)");
    m_xAxis->setLabelFormat("%d");
    chart()->addAxis(m_xAxis, Qt::AlignBottom);

    m_yAxis = new QValueAxis(this);
    m_yAxis->setTitleText("PPM");
    m_yAxis->setLabelFormat("%.1f");
    chart()->addAxis(m_yAxis, Qt::AlignLeft);
}

// ────────────────────────────────────────────────────────────────────
// Public API
// ────────────────────────────────────────────────────────────────────

void JunctionPlotWidget::setProfile(const GeneJunctionProfile& profile)
{
    m_profile = profile;
    m_annotation = profile.annotation;   // keep for codon lines across collapsed switches
    m_useCollapsed = false;
    m_collapsed.clear();
    m_highlightedPos = -1;

    // Build bar entries from individual sites
    m_bars.clear();
    m_bars.reserve(profile.sites.size());
    for (const auto& site : profile.sites) {
        if (m_inFrameOnly && !site.isInFrame()) continue;

        BarEntry entry;
        entry.position = site.position;
        entry.ppm = site.ppm;
        entry.category = classifyJunction(site);
        entry.rawCount = site.rawCount;
        entry.queryStart = site.queryStart;
        entry.frame = site.frame;
        entry.cdsClass = site.cdsClass;
        m_bars.append(entry);
    }

    buildChart();
}

void JunctionPlotWidget::setCollapsedView(const QVector<CollapsedJunction>& collapsed)
{
    m_collapsed = collapsed;
    m_useCollapsed = true;
    m_highlightedPos = -1;

    m_bars.clear();
    m_bars.reserve(collapsed.size());
    for (const auto& cj : collapsed) {
        if (m_inFrameOnly) {
            // For collapsed view, check dominant frame
            bool inFrame = (cj.dominantFrame == "+0_frame");
            if (!inFrame) continue;
        }

        BarEntry entry;
        entry.position = cj.position;
        entry.ppm = cj.totalPpm;
        entry.cdsClass = cj.cdsClass;
        entry.frame = cj.dominantFrame;
        entry.rawCount = cj.variantCount;
        entry.queryStart = 0;

        // Classify using dominant frame and CDS class
        JunctionSite dummy;
        dummy.cdsClass = cj.cdsClass;
        dummy.frame = cj.dominantFrame;
        entry.category = classifyJunction(dummy);

        m_bars.append(entry);
    }

    buildChart();
}

void JunctionPlotWidget::setFilter(bool inFrameOnly)
{
    if (m_inFrameOnly == inFrameOnly) return;
    m_inFrameOnly = inFrameOnly;

    // Re-apply to existing data
    if (m_useCollapsed) {
        setCollapsedView(m_collapsed);
    } else {
        setProfile(m_profile);
    }
}

void JunctionPlotWidget::highlightPosition(int position)
{
    m_highlightedPos = position;
    updateHighlight();
}

// ────────────────────────────────────────────────────────────────────
// Chart Building
// ────────────────────────────────────────────────────────────────────

void JunctionPlotWidget::buildChart()
{
    // Remove all existing series (but keep axes)
    chart()->removeAllSeries();
    m_highlightLine = nullptr;

    if (m_bars.isEmpty()) {
        m_xAxis->setRange(0, 100);
        m_yAxis->setRange(0, 1);
        m_fullXMin = 0;
        m_fullXMax = 100;
        m_fullYMax = 1;
        return;
    }

    // Determine axis ranges
    qreal xMin = std::numeric_limits<double>::max();
    qreal xMax = std::numeric_limits<double>::lowest();
    qreal yMax = 0.0;

    for (const auto& bar : m_bars) {
        if (bar.position < xMin) xMin = bar.position;
        if (bar.position > xMax) xMax = bar.position;
        if (bar.ppm > yMax) yMax = bar.ppm;
    }

    // Use gene annotation for full range if available
    if (m_profile.annotation.isValid()) {
        xMin = 0;
        xMax = m_profile.annotation.mRNALength;
    } else {
        // Add 5% padding
        qreal pad = (xMax - xMin) * 0.05;
        if (pad < 10) pad = 10;
        xMin = qMax(0.0, xMin - pad);
        xMax += pad;
    }

    // Add 10% y padding
    yMax *= 1.1;
    if (yMax < 1.0) yMax = 1.0;

    m_fullXMin = xMin;
    m_fullXMax = xMax;
    m_fullYMax = yMax;

    buildBars();

    // Configure axes
    m_xAxis->setRange(m_fullXMin, m_fullXMax);
    m_yAxis->setRange(0, m_fullYMax);

    // Create highlight line (initially invisible -- off-screen)
    m_highlightLine = new QLineSeries(this);
    m_highlightLine->setName("");
    QPen hlPen(QColor(220, 38, 38, 160));
    hlPen.setWidth(2);
    hlPen.setStyle(Qt::DashLine);
    m_highlightLine->setPen(hlPen);
    m_highlightLine->append(-1000, 0);
    m_highlightLine->append(-1000, m_fullYMax);
    chart()->addSeries(m_highlightLine);
    m_highlightLine->attachAxis(m_xAxis);
    m_highlightLine->attachAxis(m_yAxis);
    chart()->legend()->markers(m_highlightLine).first()->setVisible(false);

    // ORF boundary markers: green vertical at ATG (orfStart), red at STOP (orfEnd)
    // Matches the original BlastQuery "red bars (start/stop codons)" in the PPM chart.
    m_atgLine  = nullptr;
    m_stopLine = nullptr;
    if (m_annotation.isValid()) {
        auto makeCodonLine = [&](int pos, QColor color, const QString& label) -> QLineSeries* {
            auto* line = new QLineSeries(this);
            line->setName(label);
            QPen pen(color);
            pen.setWidthF(1.5);
            line->setPen(pen);
            line->append(pos, 0);
            line->append(pos, m_fullYMax);
            chart()->addSeries(line);
            line->attachAxis(m_xAxis);
            line->attachAxis(m_yAxis);
            chart()->legend()->markers(line).first()->setVisible(false);
            return line;
        };
        m_atgLine  = makeCodonLine(m_annotation.orfStart, QColor(34, 197, 94, 200),  "ATG");
        m_stopLine = makeCodonLine(m_annotation.orfEnd,   QColor(239, 68,  68, 200),  "STOP");
    }

    updateHighlight();
}

void JunctionPlotWidget::buildBars()
{
    // Calculate a reasonable bar width: roughly 0.3% of x-range, minimum 1 nt
    qreal xRange = m_fullXMax - m_fullXMin;
    qreal barWidth = qMax(1.0, xRange * 0.003);

    // Track which categories are present for legend
    bool hasCategory[4] = {false, false, false, false};

    // Sort bars by position for consistent rendering
    std::sort(m_bars.begin(), m_bars.end(), [](const BarEntry& a, const BarEntry& b) {
        return a.position < b.position;
    });

    // Create one area series per bar (each bar = pair of line series forming a filled rectangle)
    // Group by category for legend entries
    // We create one "dummy" line series per category for the legend, then individual bars
    for (int cat = 0; cat < 4; ++cat) {
        JunctionCategory jc = static_cast<JunctionCategory>(cat);
        QColor color = colorForCategory(jc);

        // Collect bars in this category
        QVector<const BarEntry*> catBars;
        for (const auto& bar : m_bars) {
            if (bar.category == jc) {
                catBars.append(&bar);
            }
        }

        if (catBars.isEmpty()) continue;
        hasCategory[cat] = true;

        // Create each bar as an area series
        bool firstInCategory = true;
        for (const auto* bar : catBars) {
            qreal x1 = bar->position - barWidth / 2.0;
            qreal x2 = bar->position + barWidth / 2.0;
            qreal y = bar->ppm;

            auto* upperLine = new QLineSeries(this);
            upperLine->append(x1, y);
            upperLine->append(x2, y);

            auto* lowerLine = new QLineSeries(this);
            lowerLine->append(x1, 0);
            lowerLine->append(x2, 0);

            auto* area = new QAreaSeries(upperLine, lowerLine);
            area->setColor(color);
            QPen borderPen(color.darker(120));
            borderPen.setWidthF(0.5);
            area->setPen(borderPen);

            if (firstInCategory) {
                area->setName(categoryLabel(jc));
                firstInCategory = false;
            } else {
                area->setName("");
            }

            chart()->addSeries(area);
            area->attachAxis(m_xAxis);
            area->attachAxis(m_yAxis);

            // Hide legend entries for non-first bars in each category
            if (area->name().isEmpty()) {
                auto markers = chart()->legend()->markers(area);
                if (!markers.isEmpty()) {
                    markers.first()->setVisible(false);
                }
            }
        }
    }
}

void JunctionPlotWidget::addBarToChart(const BarEntry& entry, double barWidth)
{
    // Unused -- bars are now built in buildBars()
    Q_UNUSED(entry);
    Q_UNUSED(barWidth);
}

void JunctionPlotWidget::updateHighlight()
{
    if (!m_highlightLine) return;

    if (m_highlightedPos < 0) {
        // Move off-screen
        m_highlightLine->replace(0, QPointF(-1000, 0));
        m_highlightLine->replace(1, QPointF(-1000, m_fullYMax));
    } else {
        m_highlightLine->replace(0, QPointF(m_highlightedPos, 0));
        m_highlightLine->replace(1, QPointF(m_highlightedPos, m_fullYMax));
    }
}

void JunctionPlotWidget::updateCodonLines()
{
    if (m_atgLine) {
        int pos = m_annotation.orfStart;
        m_atgLine->replace(0, QPointF(pos, 0));
        m_atgLine->replace(1, QPointF(pos, m_fullYMax));
    }
    if (m_stopLine) {
        int pos = m_annotation.orfEnd;
        m_stopLine->replace(0, QPointF(pos, 0));
        m_stopLine->replace(1, QPointF(pos, m_fullYMax));
    }
}

void JunctionPlotWidget::setupAxes()
{
    // Axes are already set up in the constructor; this reconfigures after data load
    m_xAxis->setRange(m_fullXMin, m_fullXMax);
    m_yAxis->setRange(0, m_fullYMax);
}

int JunctionPlotWidget::findNearestBar(const QPointF& chartPos) const
{
    if (m_bars.isEmpty()) return -1;

    // Find the bar nearest to the x-coordinate within a tolerance
    qreal xRange = m_fullXMax - m_fullXMin;
    qreal tolerance = xRange * 0.01;  // 1% of visible range
    if (tolerance < 2.0) tolerance = 2.0;

    int bestIdx = -1;
    qreal bestDist = std::numeric_limits<double>::max();

    for (int i = 0; i < m_bars.size(); ++i) {
        qreal dist = std::abs(m_bars[i].position - chartPos.x());
        if (dist < tolerance && dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
        }
    }

    return bestIdx;
}

QString JunctionPlotWidget::tooltipForBar(const BarEntry& entry) const
{
    QString tip;
    tip += QString("Position: %1\n").arg(entry.position);
    tip += QString("PPM: %1\n").arg(entry.ppm, 0, 'f', 2);
    tip += QString("Frame: %1\n").arg(entry.frame);
    tip += QString("CDS: %1\n").arg(entry.cdsClass);

    if (m_useCollapsed) {
        tip += QString("Variants: %1").arg(entry.rawCount);
    } else {
        tip += QString("QueryStart: %1\n").arg(entry.queryStart);
        tip += QString("Raw Count: %1").arg(entry.rawCount);
    }

    return tip;
}

// ────────────────────────────────────────────────────────────────────
// Mouse / Wheel Interaction
// ────────────────────────────────────────────────────────────────────

void JunctionPlotWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        QPointF scenePos = mapToScene(event->pos());
        QPointF chartPos = chart()->mapToValue(scenePos);

        // Check if we clicked on a bar
        int idx = findNearestBar(chartPos);
        if (idx >= 0) {
            emit junctionClicked(m_bars[idx].position);
            event->accept();
            return;
        }

        // Start panning
        m_isPanning = true;
        m_panStartScene = scenePos;
        m_panStartXMin = m_xAxis->min();
        m_panStartXMax = m_xAxis->max();
        setCursor(Qt::ClosedHandCursor);
        event->accept();
        return;
    }

    QChartView::mousePressEvent(event);
}

void JunctionPlotWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton && m_isPanning) {
        m_isPanning = false;
        setCursor(Qt::ArrowCursor);
        event->accept();
        return;
    }

    QChartView::mouseReleaseEvent(event);
}

void JunctionPlotWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (m_isPanning) {
        QPointF scenePos = mapToScene(event->pos());
        QPointF startVal = chart()->mapToValue(m_panStartScene);
        QPointF currentVal = chart()->mapToValue(scenePos);
        qreal dx = startVal.x() - currentVal.x();

        qreal newMin = m_panStartXMin + dx;
        qreal newMax = m_panStartXMax + dx;

        // Clamp to full range
        if (newMin < m_fullXMin) {
            qreal shift = m_fullXMin - newMin;
            newMin += shift;
            newMax += shift;
        }
        if (newMax > m_fullXMax) {
            qreal shift = newMax - m_fullXMax;
            newMin -= shift;
            newMax -= shift;
        }

        m_xAxis->setRange(newMin, newMax);
        emit rangeChanged(newMin, newMax);
        event->accept();
        return;
    }

    // Show tooltip on hover
    QPointF scenePos = mapToScene(event->pos());
    QPointF chartPos = chart()->mapToValue(scenePos);
    int idx = findNearestBar(chartPos);
    if (idx >= 0) {
        QToolTip::showText(event->globalPosition().toPoint(), tooltipForBar(m_bars[idx]), this);
    } else {
        QToolTip::hideText();
    }

    QChartView::mouseMoveEvent(event);
}

void JunctionPlotWidget::wheelEvent(QWheelEvent* event)
{
    // Zoom on x-axis around the mouse position
    QPointF scenePos = mapToScene(event->position().toPoint());
    QPointF chartPos = chart()->mapToValue(scenePos);
    qreal centerX = chartPos.x();

    qreal currentMin = m_xAxis->min();
    qreal currentMax = m_xAxis->max();
    qreal currentRange = currentMax - currentMin;

    // Zoom factor: positive delta = zoom in, negative = zoom out
    qreal factor = (event->angleDelta().y() > 0) ? 0.8 : 1.25;

    qreal newRange = currentRange * factor;

    // Don't zoom out beyond full range
    qreal fullRange = m_fullXMax - m_fullXMin;
    if (newRange > fullRange) newRange = fullRange;

    // Don't zoom in below 20 nucleotides
    if (newRange < 20.0) newRange = 20.0;

    // Keep the mouse position anchored
    qreal relPos = (centerX - currentMin) / currentRange;
    qreal newMin = centerX - relPos * newRange;
    qreal newMax = newMin + newRange;

    // Clamp
    if (newMin < m_fullXMin) {
        newMin = m_fullXMin;
        newMax = newMin + newRange;
    }
    if (newMax > m_fullXMax) {
        newMax = m_fullXMax;
        newMin = newMax - newRange;
    }

    m_xAxis->setRange(newMin, newMax);
    emit rangeChanged(newMin, newMax);

    event->accept();
}

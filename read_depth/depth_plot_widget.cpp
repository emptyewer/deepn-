#include "depth_plot_widget.h"

#include <QAreaSeries>
#include <QGraphicsTextItem>
#include <QLineSeries>
#include <QMouseEvent>
#include <QPen>
#include <QToolTip>
#include <QWheelEvent>

#include <cmath>

DepthPlotWidget::DepthPlotWidget(QWidget* parent)
    : QChartView(new QChart(), parent)
{
    setRenderHint(QPainter::Antialiasing, true);
    setMouseTracking(true);
    setRubberBand(QChartView::NoRubberBand);

    auto* c = chart();
    c->setAnimationOptions(QChart::NoAnimation);
    c->legend()->setVisible(false);
    c->setMargins(QMargins(4, 4, 4, 4));
    c->setBackgroundRoundness(0);

    // Create axes
    m_xAxis = new QValueAxis(this);
    m_xAxis->setTitleText("Position (nt)");
    m_xAxis->setLabelFormat("%d");
    m_xAxis->setGridLineVisible(true);
    m_xAxis->setMinorGridLineVisible(false);
    c->addAxis(m_xAxis, Qt::AlignBottom);

    m_yAxis = new QValueAxis(this);
    m_yAxis->setTitleText("Read Count");
    m_yAxis->setLabelFormat("%d");
    m_yAxis->setGridLineVisible(true);
    m_yAxis->setMinorGridLineVisible(false);
    c->addAxis(m_yAxis, Qt::AlignLeft);

    // Create crosshair lines in the scene
    QPen crossPen(QColor("#9CA3AF"), 1, Qt::DashLine);
    m_crosshairV = new QGraphicsLineItem();
    m_crosshairV->setPen(crossPen);
    m_crosshairV->setVisible(false);
    c->scene()->addItem(m_crosshairV);

    m_crosshairH = new QGraphicsLineItem();
    m_crosshairH->setPen(crossPen);
    m_crosshairH->setVisible(false);
    c->scene()->addItem(m_crosshairH);
}

void DepthPlotWidget::setProfile(const deepn::DepthProfile& profile)
{
    m_profile = profile;
    buildChart();
}

void DepthPlotWidget::setOverlayProfile(const deepn::DepthProfile& secondary)
{
    m_overlayProfile = secondary;
    m_hasOverlay = true;
    buildChart();
}

void DepthPlotWidget::clearOverlay()
{
    m_overlayProfile = deepn::DepthProfile();
    m_hasOverlay = false;
    buildChart();
}

void DepthPlotWidget::setBoundary(const deepn::BoundaryResult& boundary)
{
    m_boundary = boundary;
    m_hasBoundary = true;
    buildChart();
}

void DepthPlotWidget::clearBoundary()
{
    m_boundary = deepn::BoundaryResult();
    m_hasBoundary = false;
    buildChart();
}

void DepthPlotWidget::setJunctionOverlay(const QVector<deepn::JunctionSite>& junctions)
{
    m_junctions = junctions;
    m_hasJunctions = true;
    buildChart();
}

void DepthPlotWidget::clearJunctionOverlay()
{
    m_junctions.clear();
    m_hasJunctions = false;
    buildChart();
}

void DepthPlotWidget::resetZoom()
{
    if (m_xAxis && m_yAxis) {
        m_xAxis->setRange(m_fullXMin, m_fullXMax);
        m_yAxis->setRange(0, m_fullYMax * 1.05);
        m_zoomLevel = 1.0;
        emit rangeChanged(m_fullXMin, m_fullXMax);
    }
}

QChart* DepthPlotWidget::depthChart() const
{
    return chart();
}

void DepthPlotWidget::buildChart()
{
    auto* c = chart();

    // Remove all existing series
    c->removeAllSeries();

    // Remove annotations from the scene (boundary markers, junction markers)
    // Keep crosshair lines, remove everything else we added
    QList<QGraphicsItem*> sceneItems = c->scene()->items();
    for (auto* item : sceneItems) {
        if (item == m_crosshairV || item == m_crosshairH) continue;
        auto* textItem = dynamic_cast<QGraphicsTextItem*>(item);
        if (textItem && textItem->data(0).toString() == "boundary_label") {
            c->scene()->removeItem(textItem);
            delete textItem;
        }
    }

    if (m_profile.points.isEmpty()) {
        m_xAxis->setRange(0, 100);
        m_yAxis->setRange(0, 100);
        return;
    }

    // Determine full data range
    m_fullXMin = m_profile.points.first().position;
    m_fullXMax = m_profile.points.last().position + m_profile.intervalWidth;

    int maxCount = 0;
    for (const auto& pt : m_profile.points) {
        if (pt.count > maxCount) maxCount = pt.count;
    }
    if (m_hasOverlay) {
        for (const auto& pt : m_overlayProfile.points) {
            if (pt.count > maxCount) maxCount = pt.count;
        }
    }
    m_fullYMax = qMax(maxCount, 1);

    // Add primary depth series
    addDepthSeries(m_profile, QColor(PRIMARY_COLOR), 0.7);

    // Add overlay if present
    if (m_hasOverlay && !m_overlayProfile.points.isEmpty()) {
        addDepthSeries(m_overlayProfile, QColor(SECONDARY_COLOR), 0.4);
    }

    // Set axis ranges
    m_xAxis->setRange(m_fullXMin, m_fullXMax);
    m_yAxis->setRange(0, m_fullYMax * 1.05);
    m_zoomLevel = 1.0;

    // Add boundary markers
    if (m_hasBoundary && m_boundary.position > 0) {
        addBoundaryMarker(m_boundary.position,
                          QStringLiteral("3' Boundary: %1 (%2)")
                              .arg(m_boundary.position)
                              .arg(m_boundary.confidenceLabel),
                          QColor(BOUNDARY_COLOR));

        // Add secondary boundaries
        for (const auto& sec : m_boundary.secondary) {
            if (sec.position > 0) {
                addBoundaryMarker(sec.position,
                                  QStringLiteral("2nd: %1 (%2)")
                                      .arg(sec.position)
                                      .arg(sec.confidenceLabel),
                                  QColor("#F97316"));  // orange for secondary
            }
        }
    }

    // Add junction markers
    if (m_hasJunctions && !m_junctions.isEmpty()) {
        addJunctionMarkers(m_junctions);
    }

    emit rangeChanged(m_fullXMin, m_fullXMax);
}

void DepthPlotWidget::addDepthSeries(const deepn::DepthProfile& profile,
                                      const QColor& color, qreal opacity)
{
    auto* c = chart();

    // Upper boundary line (the depth values)
    auto* upperSeries = new QLineSeries(c);
    // Lower boundary (baseline at y=0)
    auto* lowerSeries = new QLineSeries(c);

    for (const auto& pt : profile.points) {
        upperSeries->append(pt.position, pt.count);
        // Also add point at end of interval for step-like appearance
        upperSeries->append(pt.position + profile.intervalWidth, pt.count);
        lowerSeries->append(pt.position, 0);
        lowerSeries->append(pt.position + profile.intervalWidth, 0);
    }

    auto* areaSeries = new QAreaSeries(upperSeries, lowerSeries);
    areaSeries->setName(profile.datasetLabel);

    QPen pen(color, 1.5);
    areaSeries->setPen(pen);

    QColor fillColor = color;
    fillColor.setAlphaF(opacity);
    areaSeries->setBrush(QBrush(fillColor));

    // Border color for the area
    areaSeries->setBorderColor(color);

    c->addSeries(areaSeries);
    areaSeries->attachAxis(m_xAxis);
    areaSeries->attachAxis(m_yAxis);
}

void DepthPlotWidget::addBoundaryMarker(int position, const QString& label,
                                          const QColor& color)
{
    auto* c = chart();

    // Add a vertical line series for the boundary
    auto* lineSeries = new QLineSeries(c);
    lineSeries->append(position, 0);
    lineSeries->append(position, m_fullYMax * 1.05);
    lineSeries->setName(label);

    QPen pen(color, 2, Qt::DashLine);
    lineSeries->setPen(pen);

    c->addSeries(lineSeries);
    lineSeries->attachAxis(m_xAxis);
    lineSeries->attachAxis(m_yAxis);

    // Add text label at the top of the boundary line
    auto* textItem = new QGraphicsTextItem(c);
    textItem->setPlainText(label);
    textItem->setDefaultTextColor(color);
    QFont font;
    font.setPointSize(8);
    font.setBold(true);
    textItem->setFont(font);
    textItem->setData(0, "boundary_label");

    // Position the label near the boundary line
    QPointF chartPos = c->mapToPosition(QPointF(position, m_fullYMax * 0.95), lineSeries);
    textItem->setPos(chartPos.x() + 4, chartPos.y());
}

void DepthPlotWidget::addJunctionMarkers(const QVector<deepn::JunctionSite>& junctions)
{
    auto* c = chart();

    for (const auto& jnc : junctions) {
        auto* lineSeries = new QLineSeries(c);
        lineSeries->append(jnc.position, 0);
        lineSeries->append(jnc.position, m_fullYMax * 0.85);

        QColor color = junctionColor(jnc);
        QPen pen(color, 1.0, Qt::DotLine);
        lineSeries->setPen(pen);
        lineSeries->setName(QStringLiteral("Jnc %1 (%2)")
                                .arg(jnc.position)
                                .arg(jnc.frameLabel()));

        c->addSeries(lineSeries);
        lineSeries->attachAxis(m_xAxis);
        lineSeries->attachAxis(m_yAxis);
    }
}

QColor DepthPlotWidget::junctionColor(const deepn::JunctionSite& site) const
{
    auto category = deepn::classifyJunction(site);
    switch (category) {
    case deepn::JunctionCategory::InOrfInFrame:
        return QColor("#1E40AF");   // deep blue
    case deepn::JunctionCategory::UpstreamInFrame:
        return QColor("#0D9488");   // teal
    case deepn::JunctionCategory::InOrfOutOfFrame:
        return QColor("#D97706");   // amber
    case deepn::JunctionCategory::DownstreamOrBack:
        return QColor("#6B7280");   // grey
    }
    return QColor("#6B7280");
}

int DepthPlotWidget::chartXToPosition(qreal chartX) const
{
    return static_cast<int>(std::round(chartX));
}

int DepthPlotWidget::findDepthAtPosition(int pos) const
{
    // Find the depth point nearest to the given position
    int bestDepth = 0;
    int bestDist = std::numeric_limits<int>::max();

    for (const auto& pt : m_profile.points) {
        // Check if position falls within this interval
        if (pos >= pt.position && pos < pt.position + m_profile.intervalWidth) {
            return pt.count;
        }
        int dist = std::abs(pt.position - pos);
        if (dist < bestDist) {
            bestDist = dist;
            bestDepth = pt.count;
        }
    }
    return bestDepth;
}

void DepthPlotWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton && m_zoomLevel > 1.01) {
        // Start panning when zoomed in
        m_isPanning = true;
        m_lastMousePos = event->position();
        setCursor(Qt::ClosedHandCursor);
        event->accept();
        return;
    }

    if (event->button() == Qt::LeftButton) {
        // Emit position click
        auto* c = chart();
        QPointF chartPos = c->mapToValue(mapToScene(event->position().toPoint()));
        int pos = chartXToPosition(chartPos.x());
        if (pos >= 0) {
            emit positionClicked(pos);
        }
    }

    QChartView::mousePressEvent(event);
}

void DepthPlotWidget::mouseMoveEvent(QMouseEvent* event)
{
    auto* c = chart();

    if (m_isPanning && (event->buttons() & Qt::LeftButton)) {
        // Pan the chart
        QPointF delta = event->position() - m_lastMousePos;
        m_lastMousePos = event->position();

        // Convert pixel delta to axis units
        QRectF plotArea = c->plotArea();
        if (plotArea.width() > 0) {
            qreal xRange = m_xAxis->max() - m_xAxis->min();
            qreal dx = -delta.x() * xRange / plotArea.width();

            qreal newMin = m_xAxis->min() + dx;
            qreal newMax = m_xAxis->max() + dx;

            // Clamp to data range
            if (newMin < m_fullXMin) {
                newMax += (m_fullXMin - newMin);
                newMin = m_fullXMin;
            }
            if (newMax > m_fullXMax) {
                newMin -= (newMax - m_fullXMax);
                newMax = m_fullXMax;
            }
            newMin = qMax(newMin, m_fullXMin);
            newMax = qMin(newMax, m_fullXMax);

            m_xAxis->setRange(newMin, newMax);
            emit rangeChanged(newMin, newMax);
        }
        event->accept();
        return;
    }

    // Update crosshair and emit cursor position
    QPointF scenePos = mapToScene(event->position().toPoint());
    QPointF chartVal = c->mapToValue(scenePos);
    QRectF plotArea = c->plotArea();

    // Check if cursor is within the plot area
    QPointF widgetPlotPos = event->position();
    QRectF mappedPlotArea = QRectF(
        mapFromScene(c->mapToScene(plotArea.topLeft())),
        mapFromScene(c->mapToScene(plotArea.bottomRight()))
    );

    if (mappedPlotArea.contains(widgetPlotPos)) {
        // Update crosshair
        m_crosshairV->setLine(scenePos.x(), c->mapToScene(plotArea.topLeft()).y(),
                               scenePos.x(), c->mapToScene(plotArea.bottomLeft()).y());
        m_crosshairV->setVisible(true);

        m_crosshairH->setLine(c->mapToScene(plotArea.topLeft()).x(), scenePos.y(),
                               c->mapToScene(plotArea.topRight()).x(), scenePos.y());
        m_crosshairH->setVisible(true);

        // Emit cursor position
        int pos = chartXToPosition(chartVal.x());
        int depth = findDepthAtPosition(pos);
        emit cursorMoved(pos, depth);

        // Show tooltip
        QString tip = QStringLiteral("Position: %1 nt\nDepth: %2 reads")
                          .arg(pos).arg(depth);
        QToolTip::showText(event->globalPosition().toPoint(), tip, this);
    } else {
        m_crosshairV->setVisible(false);
        m_crosshairH->setVisible(false);
        QToolTip::hideText();
    }

    QChartView::mouseMoveEvent(event);
}

void DepthPlotWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (m_isPanning) {
        m_isPanning = false;
        if (m_zoomLevel > 1.01) {
            setCursor(Qt::OpenHandCursor);
        } else {
            setCursor(Qt::CrossCursor);
        }
        event->accept();
        return;
    }
    QChartView::mouseReleaseEvent(event);
}

void DepthPlotWidget::wheelEvent(QWheelEvent* event)
{
    // Zoom centered on cursor position
    if (!m_xAxis) return;

    qreal factor = 1.0;
    if (event->angleDelta().y() > 0) {
        factor = 0.8;  // zoom in
    } else {
        factor = 1.25; // zoom out
    }

    qreal currentMin = m_xAxis->min();
    qreal currentMax = m_xAxis->max();
    qreal currentRange = currentMax - currentMin;
    qreal newRange = currentRange * factor;

    // Minimum zoom range: ~200 bp
    if (newRange < 200.0) {
        newRange = 200.0;
    }

    // Maximum zoom range: full data range
    qreal fullRange = m_fullXMax - m_fullXMin;
    if (newRange > fullRange) {
        newRange = fullRange;
    }

    // Center zoom on cursor position
    QPointF chartVal = chart()->mapToValue(mapToScene(event->position().toPoint()));
    qreal cursorX = chartVal.x();

    // Calculate the fraction of the cursor within the current range
    qreal cursorFrac = 0.5;
    if (currentRange > 0) {
        cursorFrac = (cursorX - currentMin) / currentRange;
        cursorFrac = qBound(0.0, cursorFrac, 1.0);
    }

    qreal newMin = cursorX - newRange * cursorFrac;
    qreal newMax = cursorX + newRange * (1.0 - cursorFrac);

    // Clamp to data bounds
    if (newMin < m_fullXMin) {
        newMax += (m_fullXMin - newMin);
        newMin = m_fullXMin;
    }
    if (newMax > m_fullXMax) {
        newMin -= (newMax - m_fullXMax);
        newMax = m_fullXMax;
    }
    newMin = qMax(newMin, m_fullXMin);
    newMax = qMin(newMax, m_fullXMax);

    m_xAxis->setRange(newMin, newMax);
    m_zoomLevel = fullRange / (newMax - newMin);

    // Update cursor style
    if (m_zoomLevel > 1.01) {
        setCursor(Qt::OpenHandCursor);
    } else {
        setCursor(Qt::CrossCursor);
    }

    emit rangeChanged(newMin, newMax);
    event->accept();
}

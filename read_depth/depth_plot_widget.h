#ifndef DEPTH_PLOT_WIDGET_H
#define DEPTH_PLOT_WIDGET_H

#include "data_structures.h"

#include <QChartView>
#include <QGraphicsLineItem>
#include <QValueAxis>

class DepthPlotWidget : public QChartView {
    Q_OBJECT
public:
    explicit DepthPlotWidget(QWidget* parent = nullptr);

    // Set depth data
    void setProfile(const deepn::DepthProfile& profile);
    void setOverlayProfile(const deepn::DepthProfile& secondary);
    void clearOverlay();

    // Set boundary marker
    void setBoundary(const deepn::BoundaryResult& boundary);
    void clearBoundary();

    // Set junction overlay (from MultiQuery++ data)
    void setJunctionOverlay(const QVector<deepn::JunctionSite>& junctions);
    void clearJunctionOverlay();

    // Reset zoom to show full range
    void resetZoom();

    QChart* depthChart() const;

signals:
    void cursorMoved(int position, int depth);
    void rangeChanged(qreal xMin, qreal xMax);
    void positionClicked(int position);

protected:
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private:
    void buildChart();
    void addDepthSeries(const deepn::DepthProfile& profile,
                        const QColor& color, qreal opacity = 1.0);
    void addBoundaryMarker(int position, const QString& label,
                           const QColor& color);
    void addJunctionMarkers(const QVector<deepn::JunctionSite>& junctions);
    QColor junctionColor(const deepn::JunctionSite& site) const;

    int chartXToPosition(qreal chartX) const;
    int findDepthAtPosition(int pos) const;

    deepn::DepthProfile m_profile;
    deepn::DepthProfile m_overlayProfile;
    deepn::BoundaryResult m_boundary;
    QVector<deepn::JunctionSite> m_junctions;
    bool m_hasOverlay = false;
    bool m_hasBoundary = false;
    bool m_hasJunctions = false;

    QValueAxis* m_xAxis = nullptr;
    QValueAxis* m_yAxis = nullptr;

    // Full data range for reset
    qreal m_fullXMin = 0.0;
    qreal m_fullXMax = 1.0;
    qreal m_fullYMax = 1.0;

    // Zoom/pan state
    bool m_isPanning = false;
    QPointF m_lastMousePos;
    qreal m_zoomLevel = 1.0;

    // Crosshair lines
    QGraphicsLineItem* m_crosshairV = nullptr;
    QGraphicsLineItem* m_crosshairH = nullptr;

    // Colors
    static constexpr const char* PRIMARY_COLOR = "#2563EB";
    static constexpr const char* PRIMARY_FILL_COLOR = "#93BBFD";
    static constexpr const char* SECONDARY_COLOR = "#F59E0B";
    static constexpr const char* SECONDARY_FILL_COLOR = "#FDE68A";
    static constexpr const char* BOUNDARY_COLOR = "#DC2626";
};

#endif  // DEPTH_PLOT_WIDGET_H

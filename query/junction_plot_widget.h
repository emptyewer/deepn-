#ifndef JUNCTION_PLOT_WIDGET_H
#define JUNCTION_PLOT_WIDGET_H

#include <data_structures.h>

#include <QChartView>
#include <QValueAxis>

QT_BEGIN_NAMESPACE
class QAreaSeries;
class QChart;
class QLineSeries;
class QScatterSeries;
QT_END_NAMESPACE

struct BarEntry {
    int position = 0;
    double ppm = 0.0;
    deepn::JunctionCategory category;
    int rawCount = 0;
    int queryStart = 0;
    QString frame;
    QString cdsClass;
};

class JunctionPlotWidget : public QChartView
{
    Q_OBJECT

public:
    explicit JunctionPlotWidget(QWidget* parent = nullptr);

    void setProfile(const deepn::GeneJunctionProfile& profile);
    void setCollapsedView(const QVector<deepn::CollapsedJunction>& collapsed);
    void setFilter(bool inFrameOnly);
    void highlightPosition(int position);

    // Color constants for junction categories
    static QColor colorForCategory(deepn::JunctionCategory cat);

signals:
    void junctionClicked(int position);
    void rangeChanged(qreal xMin, qreal xMax);

protected:
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private:
    void buildChart();
    void buildBars();
    void addBarToChart(const BarEntry& entry, double barWidth);
    void updateHighlight();
    void setupAxes();
    int findNearestBar(const QPointF& chartPos) const;
    QString tooltipForBar(const BarEntry& entry) const;

    deepn::GeneJunctionProfile m_profile;
    QVector<deepn::CollapsedJunction> m_collapsed;
    QVector<BarEntry> m_bars;
    bool m_useCollapsed = false;
    bool m_inFrameOnly = false;
    int m_highlightedPos = -1;

    QValueAxis* m_xAxis = nullptr;
    QValueAxis* m_yAxis = nullptr;

    // Full range for reset
    qreal m_fullXMin = 0.0;
    qreal m_fullXMax = 1.0;
    qreal m_fullYMax = 1.0;

    // Interaction state
    bool m_isPanning = false;
    QPointF m_panStartScene;
    qreal m_panStartXMin = 0.0;
    qreal m_panStartXMax = 0.0;

    // Highlight overlay
    QLineSeries* m_highlightLine = nullptr;
};

#endif // JUNCTION_PLOT_WIDGET_H

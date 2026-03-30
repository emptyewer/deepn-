#ifndef DEEPN_SYNC_CONTROLLER_H
#define DEEPN_SYNC_CONTROLLER_H

#include <QObject>
#include <QVector>

QT_BEGIN_NAMESPACE
class QChart;
QT_END_NAMESPACE

namespace deepn {

class SyncController : public QObject {
    Q_OBJECT
public:
    explicit SyncController(QObject* parent = nullptr);

    // Register chart views to synchronize (x-axis only, y-axes independent)
    void addChart(QChart* chart);
    void removeChart(QChart* chart);
    void clear();

    // Set the visible x-range for all charts
    void setXRange(qreal min, qreal max);

    // Enable/disable sync
    void setEnabled(bool enabled);
    bool isEnabled() const;

signals:
    void xRangeChanged(qreal min, qreal max);

private:
    QVector<QChart*> m_charts;
    bool m_enabled = true;
    bool m_syncing = false;  // prevent recursive updates

    void onAxisRangeChanged(qreal min, qreal max);
};

}  // namespace deepn

#endif  // DEEPN_SYNC_CONTROLLER_H

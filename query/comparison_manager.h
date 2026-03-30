#ifndef COMPARISON_MANAGER_H
#define COMPARISON_MANAGER_H

#include <sync_controller.h>

#include <QObject>

class JunctionPlotWidget;

namespace deepn {
class MRNATrackWidget;
}

class ComparisonManager : public QObject
{
    Q_OBJECT

public:
    explicit ComparisonManager(QObject* parent = nullptr);

    void setEnabled(bool enabled);
    bool isEnabled() const;

    void setPrimaryPlot(JunctionPlotWidget* plot);
    void setSecondaryPlot(JunctionPlotWidget* plot);
    void setPrimaryTrack(deepn::MRNATrackWidget* track);
    void setSecondaryTrack(deepn::MRNATrackWidget* track);

    // Synchronize zoom/scroll between primary and secondary charts
    void syncRanges();

private:
    deepn::SyncController m_syncController;
    JunctionPlotWidget* m_primary = nullptr;
    JunctionPlotWidget* m_secondary = nullptr;
    deepn::MRNATrackWidget* m_primaryTrack = nullptr;
    deepn::MRNATrackWidget* m_secondaryTrack = nullptr;
    bool m_enabled = false;

    void onPrimaryRangeChanged(qreal xMin, qreal xMax);
    void onSecondaryRangeChanged(qreal xMin, qreal xMax);
};

#endif // COMPARISON_MANAGER_H

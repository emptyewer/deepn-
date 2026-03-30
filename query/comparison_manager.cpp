#include "comparison_manager.h"
#include "junction_plot_widget.h"

#include <mrna_track_widget.h>

#include <QChart>
#include <QValueAxis>

using namespace deepn;

// ────────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────────

ComparisonManager::ComparisonManager(QObject* parent)
    : QObject(parent)
{
}

// ────────────────────────────────────────────────────────────────────
// Public API
// ────────────────────────────────────────────────────────────────────

void ComparisonManager::setEnabled(bool enabled)
{
    if (m_enabled == enabled) return;
    m_enabled = enabled;

    m_syncController.setEnabled(enabled);

    if (enabled) {
        // Register both charts with the sync controller
        if (m_primary) {
            m_syncController.addChart(m_primary->chart());
        }
        if (m_secondary) {
            m_syncController.addChart(m_secondary->chart());
        }
    } else {
        m_syncController.clear();
    }
}

bool ComparisonManager::isEnabled() const
{
    return m_enabled;
}

void ComparisonManager::setPrimaryPlot(JunctionPlotWidget* plot)
{
    // Disconnect old
    if (m_primary) {
        disconnect(m_primary, &JunctionPlotWidget::rangeChanged,
                   this, &ComparisonManager::onPrimaryRangeChanged);
    }

    m_primary = plot;

    if (m_primary) {
        connect(m_primary, &JunctionPlotWidget::rangeChanged,
                this, &ComparisonManager::onPrimaryRangeChanged);
    }
}

void ComparisonManager::setSecondaryPlot(JunctionPlotWidget* plot)
{
    // Disconnect old
    if (m_secondary) {
        disconnect(m_secondary, &JunctionPlotWidget::rangeChanged,
                   this, &ComparisonManager::onSecondaryRangeChanged);
    }

    m_secondary = plot;

    if (m_secondary) {
        connect(m_secondary, &JunctionPlotWidget::rangeChanged,
                this, &ComparisonManager::onSecondaryRangeChanged);
    }
}

void ComparisonManager::setPrimaryTrack(MRNATrackWidget* track)
{
    m_primaryTrack = track;
}

void ComparisonManager::setSecondaryTrack(MRNATrackWidget* track)
{
    m_secondaryTrack = track;
}

void ComparisonManager::syncRanges()
{
    if (!m_enabled || !m_primary || !m_secondary) return;

    // Take the primary chart's current x-range and apply to the secondary
    QChart* primaryChart = m_primary->chart();
    if (!primaryChart) return;

    auto axes = primaryChart->axes(Qt::Horizontal);
    if (axes.isEmpty()) return;

    auto* xAxis = qobject_cast<QValueAxis*>(axes.first());
    if (!xAxis) return;

    m_syncController.setXRange(xAxis->min(), xAxis->max());
}

// ────────────────────────────────────────────────────────────────────
// Range change propagation
// ────────────────────────────────────────────────────────────────────

void ComparisonManager::onPrimaryRangeChanged(qreal xMin, qreal xMax)
{
    if (!m_enabled) return;

    // Propagate to sync controller, which updates both charts
    m_syncController.setXRange(xMin, xMax);

    // Also update both mRNA tracks
    int iMin = static_cast<int>(xMin);
    int iMax = static_cast<int>(xMax);
    if (m_primaryTrack) m_primaryTrack->setVisibleRange(iMin, iMax);
    if (m_secondaryTrack) m_secondaryTrack->setVisibleRange(iMin, iMax);
}

void ComparisonManager::onSecondaryRangeChanged(qreal xMin, qreal xMax)
{
    if (!m_enabled) return;

    m_syncController.setXRange(xMin, xMax);

    int iMin = static_cast<int>(xMin);
    int iMax = static_cast<int>(xMax);
    if (m_primaryTrack) m_primaryTrack->setVisibleRange(iMin, iMax);
    if (m_secondaryTrack) m_secondaryTrack->setVisibleRange(iMin, iMax);
}

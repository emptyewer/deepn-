#include "sync_controller.h"

#include <QChart>
#include <QDebug>
#include <QValueAxis>

#include <cmath>

namespace deepn {

SyncController::SyncController(QObject* parent)
    : QObject(parent)
{
}

void SyncController::addChart(QChart* chart)
{
    if (!chart || m_charts.contains(chart))
        return;

    m_charts.append(chart);

    // Connect to x-axis range changes on all horizontal value axes
    const auto axes = chart->axes(Qt::Horizontal);
    for (auto* axis : axes) {
        auto* valueAxis = qobject_cast<QValueAxis*>(axis);
        if (valueAxis) {
            connect(valueAxis, &QValueAxis::rangeChanged,
                    this, &SyncController::onAxisRangeChanged);
        }
    }
}

void SyncController::removeChart(QChart* chart)
{
    if (!chart)
        return;

    int idx = m_charts.indexOf(chart);
    if (idx < 0)
        return;

    // Disconnect axis signals
    const auto axes = chart->axes(Qt::Horizontal);
    for (auto* axis : axes) {
        auto* valueAxis = qobject_cast<QValueAxis*>(axis);
        if (valueAxis) {
            disconnect(valueAxis, &QValueAxis::rangeChanged,
                       this, &SyncController::onAxisRangeChanged);
        }
    }

    m_charts.removeAt(idx);
}

void SyncController::clear()
{
    for (auto* chart : m_charts) {
        const auto axes = chart->axes(Qt::Horizontal);
        for (auto* axis : axes) {
            auto* valueAxis = qobject_cast<QValueAxis*>(axis);
            if (valueAxis) {
                disconnect(valueAxis, &QValueAxis::rangeChanged,
                           this, &SyncController::onAxisRangeChanged);
            }
        }
    }
    m_charts.clear();
}

void SyncController::setXRange(qreal min, qreal max)
{
    if (!m_enabled || m_syncing)
        return;

    m_syncing = true;

    for (auto* chart : m_charts) {
        const auto axes = chart->axes(Qt::Horizontal);
        for (auto* axis : axes) {
            auto* valueAxis = qobject_cast<QValueAxis*>(axis);
            if (valueAxis) {
                valueAxis->setRange(min, max);
            }
        }
    }

    m_syncing = false;
    emit xRangeChanged(min, max);
}

void SyncController::setEnabled(bool enabled)
{
    m_enabled = enabled;
}

bool SyncController::isEnabled() const
{
    return m_enabled;
}

void SyncController::onAxisRangeChanged(qreal min, qreal max)
{
    if (!m_enabled || m_syncing)
        return;

    m_syncing = true;

    // Identify the source axis to avoid setting it again
    QValueAxis* sourceAxis = qobject_cast<QValueAxis*>(sender());

    // Propagate the range change to all other charts' x-axes
    for (auto* chart : m_charts) {
        const auto axes = chart->axes(Qt::Horizontal);
        for (auto* axis : axes) {
            auto* valueAxis = qobject_cast<QValueAxis*>(axis);
            if (valueAxis && valueAxis != sourceAxis) {
                // Only update if the range actually differs (avoid spurious signals)
                if (std::abs(valueAxis->min() - min) > 0.01 ||
                    std::abs(valueAxis->max() - max) > 0.01) {
                    valueAxis->setRange(min, max);
                }
            }
        }
    }

    m_syncing = false;
    emit xRangeChanged(min, max);
}

}  // namespace deepn

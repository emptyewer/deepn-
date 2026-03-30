#ifndef DEEPN_EXPORT_ENGINE_H
#define DEEPN_EXPORT_ENGINE_H

#include "data_structures.h"

#include <QSize>
#include <QString>
#include <QVector>

class QWidget;

namespace deepn {

class ExportEngine {
public:
    // CSV export
    static bool exportJunctionCSV(const QString& filePath, const GeneJunctionProfile& profile);
    static bool exportCollapsedCSV(const QString& filePath, const QVector<CollapsedJunction>& collapsed);
    static bool exportDepthCSV(const QString& filePath, const DepthProfile& profile);
    static bool exportBatchSummaryCSV(const QString& filePath,
                                       const QVector<GeneJunctionProfile>& profiles);

    // Figure export (renders a QWidget/QChartView to file)
    static bool exportFigureSVG(const QString& filePath, QWidget* widget, QSize size = {800, 400});
    static bool exportFigurePDF(const QString& filePath, QWidget* widget, QSize size = {800, 400});
    static bool exportFigurePNG(const QString& filePath, QWidget* widget, QSize size = {800, 400}, int dpi = 300);

    static QString lastError();
};

}  // namespace deepn

#endif  // DEEPN_EXPORT_ENGINE_H

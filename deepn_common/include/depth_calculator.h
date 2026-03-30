#ifndef DEEPN_DEPTH_CALCULATOR_H
#define DEEPN_DEPTH_CALCULATOR_H

#include "data_structures.h"

#include <QString>
#include <QVector>

namespace deepn {

class DepthCalculator {
public:
    // Compute depth profile for a gene
    DepthProfile calculate(const QString& dbPath,
                           const QString& gene,
                           const GeneAnnotation& annotation,
                           int intervalWidth = 25,
                           int intervalSpacing = 50) const;

    // Detect 3' boundary from a depth profile
    static BoundaryResult detect3PrimeBoundary(const DepthProfile& profile);

    // Smooth a depth profile (moving average)
    static QVector<DepthPoint> smooth(const QVector<DepthPoint>& points, int windowSize = 3);

    // Downsample for visualization (Largest-Triangle-Three-Buckets)
    static QVector<DepthPoint> downsample(const QVector<DepthPoint>& points, int targetCount);

    QString lastError() const;

private:
    mutable QString m_error;
};

}  // namespace deepn

#endif  // DEEPN_DEPTH_CALCULATOR_H

#ifndef BOUNDARY_DETECTOR_H
#define BOUNDARY_DETECTOR_H

#include "data_structures.h"
#include "depth_calculator.h"

#include <QString>

class BoundaryDetector {
public:
    struct DetailedResult {
        deepn::BoundaryResult primary;
        deepn::InsertExtent insertExtent;
        bool hasInsertExtent = false;
        QString summary;
    };

    // Detect boundary and format results
    static DetailedResult analyze(const deepn::DepthProfile& profile,
                                   int fivePrimeJunction = 0);

    // Format boundary info as rich text for display
    static QString formatBoundaryInfo(const deepn::BoundaryResult& result);

    // Format insert extent as rich text for display
    static QString formatInsertExtent(const deepn::InsertExtent& extent,
                                       const deepn::GeneAnnotation& annotation);
};

#endif  // BOUNDARY_DETECTOR_H

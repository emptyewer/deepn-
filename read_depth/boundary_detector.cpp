#include "boundary_detector.h"

#include <QLocale>

#include <algorithm>
#include <cmath>

BoundaryDetector::DetailedResult BoundaryDetector::analyze(
    const deepn::DepthProfile& profile, int fivePrimeJunction)
{
    DetailedResult result;

    // Run boundary detection from the shared library
    result.primary = deepn::DepthCalculator::detect3PrimeBoundary(profile);

    // Compute insert extent if we have a 5' junction
    if (fivePrimeJunction > 0 && result.primary.position > 0) {
        deepn::InsertExtent extent;
        extent.fivePrimeJunction = fivePrimeJunction;
        extent.threePrimeBoundary = result.primary.position;
        extent.insertLength = extent.threePrimeBoundary - extent.fivePrimeJunction;

        if (extent.insertLength < 0) {
            extent.insertLength = 0;
        }

        // Calculate CDS overlap
        const auto& annot = profile.annotation;
        if (annot.orfStart > 0 && annot.orfEnd > annot.orfStart) {
            int overlapStart = std::max(fivePrimeJunction, annot.orfStart);
            int overlapEnd = std::min(result.primary.position, annot.orfEnd);
            int overlap = std::max(0, overlapEnd - overlapStart);
            int cdsLength = annot.orfEnd - annot.orfStart;

            if (cdsLength > 0) {
                extent.cdsOverlapPercent = (static_cast<double>(overlap) / cdsLength) * 100.0;
            }
        }

        // In-frame is determined by whether the 5' junction position
        // relative to orfStart is divisible by 3
        if (fivePrimeJunction >= profile.annotation.orfStart) {
            int offset = fivePrimeJunction - profile.annotation.orfStart;
            extent.inFrame = (offset % 3 == 0);
        }

        result.insertExtent = extent;
        result.hasInsertExtent = true;
    }

    // Build summary text
    result.summary = formatBoundaryInfo(result.primary);
    if (result.hasInsertExtent) {
        result.summary += "\n" + formatInsertExtent(result.insertExtent, profile.annotation);
    }

    return result;
}

QString BoundaryDetector::formatBoundaryInfo(const deepn::BoundaryResult& result)
{
    if (result.position <= 0) {
        return QStringLiteral("No 3' boundary detected (insufficient coverage data)");
    }

    QLocale locale;
    QString posStr = locale.toString(result.position);

    // Confidence color
    QString confColor;
    if (result.confidenceLabel == "HIGH") {
        confColor = "#059669";  // green
    } else if (result.confidenceLabel == "MEDIUM") {
        confColor = "#D97706";  // amber
    } else {
        confColor = "#DC2626";  // red
    }

    QString html = QStringLiteral(
        "<b>3' Boundary:</b> position <b>%1</b> nt "
        "&nbsp;&nbsp;<span style='color:%2; font-weight:bold;'>(%3 confidence)</span>"
        "<br/>"
        "Depth before: %4 reads &rarr; Depth after: %5 reads"
    ).arg(posStr, confColor, result.confidenceLabel,
          locale.toString(result.depthBefore),
          locale.toString(result.depthAfter));

    // Secondary boundaries
    if (!result.secondary.isEmpty()) {
        html += QStringLiteral("<br/><br/><b>Additional boundaries:</b>");
        for (int i = 0; i < result.secondary.size(); ++i) {
            const auto& sec = result.secondary[i];
            html += QStringLiteral("<br/>&nbsp;&nbsp;%1. position <b>%2</b> nt (%3)")
                        .arg(i + 2)
                        .arg(locale.toString(sec.position))
                        .arg(sec.confidenceLabel);
        }
    }

    return html;
}

QString BoundaryDetector::formatInsertExtent(const deepn::InsertExtent& extent,
                                              const deepn::GeneAnnotation& annotation)
{
    QLocale locale;

    QString frameStr = extent.inFrame
        ? QStringLiteral("<span style='color:#059669; font-weight:bold;'>In-frame</span>")
        : QStringLiteral("<span style='color:#DC2626; font-weight:bold;'>Out-of-frame</span>");

    QString html = QStringLiteral(
        "<b>Insert Reconstruction:</b><br/>"
        "&nbsp;&nbsp;5' junction: position <b>%1</b> nt<br/>"
        "&nbsp;&nbsp;3' boundary: position <b>%2</b> nt<br/>"
        "&nbsp;&nbsp;Insert length: <b>%3</b> nt<br/>"
        "&nbsp;&nbsp;CDS overlap: <b>%4%</b> of ORF (%5-%6)<br/>"
        "&nbsp;&nbsp;Frame: %7"
    ).arg(locale.toString(extent.fivePrimeJunction),
          locale.toString(extent.threePrimeBoundary),
          locale.toString(extent.insertLength),
          QString::number(extent.cdsOverlapPercent, 'f', 1),
          locale.toString(annotation.orfStart),
          locale.toString(annotation.orfEnd),
          frameStr);

    return html;
}

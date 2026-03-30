#ifndef DEEPN_DATA_STRUCTURES_H
#define DEEPN_DATA_STRUCTURES_H

#include <QString>
#include <QVector>

namespace deepn {

// Gene annotation loaded from FASTA reference database
struct GeneAnnotation {
    QString refseq;      // NM_* accession
    QString geneName;    // gene symbol
    QString chromosome;
    int orfStart = 0;    // CDS start position on mRNA
    int orfEnd = 0;      // CDS end position on mRNA
    int mRNALength = 0;  // total mRNA length (= sequence length)
    QString sequence;    // full mRNA nucleotide sequence

    bool isValid() const { return !refseq.isEmpty() && mRNALength > 0; }
};

// A single junction site on a gene (one row in the maps table)
struct JunctionSite {
    int position = 0;       // rstart: nucleotide position on gene reference
    int positionEnd = 0;    // rend: alignment end on gene reference
    double ppm = 0.0;       // abundance (Parts Per Million)
    int queryStart = 0;     // qstart: distance from junction to cDNA match start
    int queryEnd = 0;       // qend: alignment end on query
    QString cdsClass;       // "upstream", "in_orf", "downstream"
    QString frame;          // "+0_frame", "+1_frame", "+2_frame", "backward"
    QString refseq;         // NM_* accession
    QString geneName;       // gene symbol
    int rawCount = 0;       // raw junction count before PPM normalization

    // Convenience: is this an in-frame junction?
    bool isInFrame() const { return frame == "+0_frame"; }

    // Convenience: human-readable frame label
    QString frameLabel() const {
        if (frame == "+0_frame") return "In-frame";
        if (frame == "+1_frame" || frame == "+2_frame") return "Out-of-frame";
        if (frame == "backward") return "Backward";
        return frame;
    }

    // Convenience: human-readable CDS label
    QString cdsLabel() const {
        if (cdsClass == "in_orf") return "In ORF";
        if (cdsClass == "upstream") return "Upstream";
        if (cdsClass == "downstream") return "Downstream";
        return cdsClass;
    }
};

// All junction data for one gene from one dataset
struct GeneJunctionProfile {
    GeneAnnotation annotation;
    QVector<JunctionSite> sites;
    int totalReads = 0;      // total distinct reads in this dataset
    QString datasetLabel;    // e.g., "Selected", "Unselected"
    QString sourceFile;      // SQLite file path
};

// Junctions collapsed by position (sum of PPM across QueryStart variants)
struct CollapsedJunction {
    int position = 0;
    double totalPpm = 0.0;
    int variantCount = 0;        // distinct QueryStart values
    QString dominantFrame;       // most common frame at this position
    QString cdsClass;
    QVector<JunctionSite> variants;
};

// Single depth measurement at an interval position
struct DepthPoint {
    int position = 0;      // start position of the interval on the gene
    int count = 0;         // number of reads covering this interval
    double normalized = 0; // normalized count (e.g., per million)
};

// Complete depth profile for one gene from one dataset
struct DepthProfile {
    GeneAnnotation annotation;
    int intervalWidth = 25;      // W (bp)
    int intervalSpacing = 50;    // D (nt)
    int totalReads = 0;
    QVector<DepthPoint> points;
    QString datasetLabel;
    QString sourceFile;
};

// Detected 3' boundary
struct BoundaryResult {
    int position = 0;
    double confidence = 0.0;        // 0.0 to 1.0
    QString confidenceLabel;        // "HIGH", "MEDIUM", "LOW"
    int depthBefore = 0;            // depth just before boundary
    int depthAfter = 0;             // depth just after boundary
    QVector<BoundaryResult> secondary;  // additional boundaries
};

// Complete insert reconstruction (MultiQuery++ + ReadDepth++ data)
struct InsertExtent {
    int fivePrimeJunction = 0;      // from MultiQuery++
    int threePrimeBoundary = 0;     // from ReadDepth++
    int insertLength = 0;
    double cdsOverlapPercent = 0.0;
    bool inFrame = false;
};

// DESeq2/StatMaker++ result row (for gene navigation)
struct DESeq2Result {
    QString gene;
    double baseMean = 0.0;
    double log2FoldChange = 0.0;
    double lfcSE = 0.0;
    double stat = 0.0;
    double pvalue = 1.0;
    double padj = 1.0;
    QString enrichment;
};

// Junction color classification for visualization
enum class JunctionCategory {
    InOrfInFrame,       // Deep blue - most biologically relevant
    UpstreamInFrame,    // Teal - correct frame but upstream
    InOrfOutOfFrame,    // Amber - within CDS but wrong frame
    DownstreamOrBack    // Grey - past stop codon or backward
};

inline JunctionCategory classifyJunction(const JunctionSite& site) {
    bool inFrame = site.isInFrame();
    if (site.cdsClass == "in_orf" && inFrame)
        return JunctionCategory::InOrfInFrame;
    if (site.cdsClass == "upstream" && inFrame)
        return JunctionCategory::UpstreamInFrame;
    if (site.cdsClass == "in_orf" && !inFrame)
        return JunctionCategory::InOrfOutOfFrame;
    return JunctionCategory::DownstreamOrBack;
}

// Database schema version for migration support
constexpr int MAPS_SCHEMA_VERSION = 2;  // v2 adds rstart/rend columns

}  // namespace deepn

#endif  // DEEPN_DATA_STRUCTURES_H

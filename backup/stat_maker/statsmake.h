#ifndef ANALYZEDEEPN_H
#define ANALYZEDEEPN_H

#include <QString>
#include <QStringList>
#include <QVector>
#include <QMap>
#include <QPair>
#include <vector>
#include <string>
#include "AdvancedNBDispersion.h"

// -------------------------------------------------------------------------
// A structure to store count data for Genes x Conditions x Replicates.
//
// Indices: data[i][b][r] => internally flattened to data[i * (nConds*nReps) + b*nReps + r]
//
// nGenes     = # of genes
// nConds     = # of conditions (typically 2: Baseline and Selected)
// nReplicates= # of replicates for each condition
// data       = flattened 3D array
// -------------------------------------------------------------------------
struct CountArray {
    int nGenes;
    int nConds;
    int nReplicates;
    std::vector<int> data; // Flattened 3D array

    // Indexing operators:
    int &operator()(int gene, int cond, int rep) {
        return data[gene * (nConds * nReplicates) + cond * nReplicates + rep];
    }

    const int &operator()(int gene, int cond, int rep) const {
        return data[gene * (nConds * nReplicates) + cond * nReplicates + rep];
    }
};

// -------------------------------------------------------------------------
// Main data structure similar to "data.deepn" in R.
// It holds:
//   - Vector: CountArray for vector replicates
//   - Bait:   CountArray for bait replicates
//   - multiBait: flag if Bait is multi-dimensional (multiple replicates)
//   - vtr, btr: total read sums for each replicate (for RPM normalization)
//   - omega: overdispersion estimates (Baseline, Selected, baitEffect)
//   - geneNames: list of genes
// -------------------------------------------------------------------------
class DataDeepn {
public:
    CountArray Vector;
    CountArray Bait;
    bool multiBait;

    std::vector<double> vtr;  // Flattened totals for Vector (size = nConds*nReps)
    std::vector<double> btr;  // Flattened totals for Bait   (size = nConds*nReps)

    // Overdispersion [ Baseline, Selected, baitEffect ]
    std::vector<double> omega;

    // Gene names
    std::vector<std::string> geneNames;

    DataDeepn() : multiBait(false) {}
};

// -------------------------------------------------------------------------
// Function declarations
// -------------------------------------------------------------------------

// Reads a set of CSV files (non-selected vs. selected) and assembles them into a CountArray.
// This mirrors readDeepn() in R.
CountArray readDeepn(const QStringList &nFiles,
                     const QStringList &sFiles,
                     std::vector<std::string> &geneNamesOut,
                     bool &genesInitialized);

// Imports vector/bait data from lists of file paths, computing total reads.
// Mirrors import() in R.
DataDeepn importData(const QStringList &vn,
                     const QStringList &vs,
                     const QStringList &bn = QStringList(),
                     const QStringList &bs = QStringList());

// Convert raw counts to RPM. Mirrors rpm(Data) in R.
void rpm(const DataDeepn &data,
         std::vector<double> &rpmVector,
         std::vector<double> &rpmBait);

// Filter out low-abundance genes based on a threshold in RPM, analogous to applyFilter(Data, thresh).
void applyFilter(DataDeepn &data, double threshold, bool base = true);

// Estimate overdispersion using a method-of-moments approach (approx. edgeR).
// Returns a 3-element vector: [Baseline, Selected, baitEffect].
std::vector<double> overdisp(const DataDeepn &data);

// Choose filter thresholds, analogous to chooseFilter(Data, minRPM, maxRPM).
// Returns a vector of (threshold, (#genes passing, overdispSelected)) pairs.
std::vector<QPair<double, QPair<int, double>>>
chooseFilter(const DataDeepn &data, double minRPM, double maxRPM);

// Summarize data in a CSV, partial analog to summary.data.deepn() in R.
void summarizeData(const DataDeepn &data,
                   bool sort,
                   const QString &outfile = QString());

// Main high-level analysis step, analogous to analyzeDeepn(infile, outfile, ...).
void analyzeDeepn(const QString &infile,
                  const QString &outfile = "stat.csv",
                  const QString &msgfile = "messages.txt",
                  bool debug = false,
                  int sortOption = 1);

#endif // ANALYZEDEEPN_H
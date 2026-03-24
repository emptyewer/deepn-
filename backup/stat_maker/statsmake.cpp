#include "statsmake.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>


// -------------------------------------------------------------------------
// parseDeepnCsv: a small utility to parse the CSV files, skipping first 4 lines.
// Returns a vector of (geneName, count).
// R code uses "skip=4, colClasses=c(V1='character', V2='character', V3='numeric')"
//
// This is naive code expecting 3 columns separated by commas.
// -------------------------------------------------------------------------
static std::vector<QPair<QString, double> > parseDeepnCsv(const QString &filename) {
    std::vector<QPair<QString, double> > data;
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        throw std::runtime_error(("Could not open file: " + filename).toStdString());
    }
    QTextStream in(&file);

    // Skip first 4 lines
    for (int i = 0; i < 4; ++i) {
        if (in.atEnd()) break;
        in.readLine();
    }

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        // Split by comma
        QStringList parts = line.split(",");
        if (parts.size() < 3) continue;

        // parts[1] is the gene, parts[2] is numeric count
        QString gene = parts[1].replace(" ", "");  // remove spaces
        double val = parts[2].toDouble();
        data.push_back(qMakePair(gene, val));
    }
    return data;
}

// -------------------------------------------------------------------------
// readDeepn: merges the baseline and selected files for each replicate
// into a single CountArray: [nGenes x 2 conditions x K replicates].
//
// We also do the "divide by minNonZero" logic from R's readDeepn
// to ensure integer counts after scaling.
// -------------------------------------------------------------------------
CountArray readDeepn(const QStringList &nFiles,
                     const QStringList &sFiles,
                     std::vector<std::string> &geneNamesOut,
                     bool &genesInitialized) {
    if (nFiles.size() != sFiles.size()) {
        throw std::runtime_error("Number of 'Non' files must match 'Sel' files!");
    }
    int K = nFiles.size();

    // We'll accumulate gene -> [2 x K] counts in a QMap.
    // That is, for each gene, we have a 2D array [cond=0 or 1][rep=K].
    QMap<QString, std::vector<std::vector<double> > > geneMap;

    for (int k = 0; k < K; ++k) {
        // Parse baseline file
        auto baseData = parseDeepnCsv(nFiles[k]);
        // Parse selected file
        auto selData = parseDeepnCsv(sFiles[k]);

        // Store
        for (auto &p: baseData) {
            QString g = p.first;
            double c = p.second;
            if (!geneMap.contains(g)) {
                geneMap[g] = std::vector<std::vector<double> >(2, std::vector<double>(K, 0.0));
            }
            geneMap[g][0][k] += c; // baseline
        }
        for (auto &p: selData) {
            QString g = p.first;
            double c = p.second;
            if (!geneMap.contains(g)) {
                geneMap[g] = std::vector<std::vector<double> >(2, std::vector<double>(K, 0.0));
            }
            geneMap[g][1][k] += c; // selected
        }
    }

    // Find minNonZero across all counts > 0
    double minNonZero = 1e12;
    for (auto it = geneMap.begin(); it != geneMap.end(); ++it) {
        for (int cond = 0; cond < 2; ++cond) {
            for (int r = 0; r < K; ++r) {
                double val = it.value()[cond][r];
                if (val > 0.0 && val < minNonZero) {
                    minNonZero = val;
                }
            }
        }
    }
    if (minNonZero < 1e-12) {
        minNonZero = 1.0;  // fallback if everything was zero
    }

    // We'll flatten into CountArray
    QStringList geneList = geneMap.keys();
    // If you want them sorted: geneList.sort();
    int nG = geneList.size();

    CountArray out;
    out.nGenes = nG;
    out.nConds = 2; // baseline, selected
    out.nReplicates = K;
    out.data.resize(nG * 2 * K, 0);

    for (int i = 0; i < nG; ++i) {
        const QString &g = geneList[i];
        for (int cond = 0; cond < 2; ++cond) {
            for (int r = 0; r < K; ++r) {
                double val = geneMap[g][cond][r] / minNonZero;
                int ival = (int) std::round(val);
                out(i, cond, r) = ival;
            }
        }
        // If this is the first time we populate gene names globally, push them
        if (!genesInitialized) {
            geneNamesOut.push_back(g.toStdString());
        }
    }

    return out;
}

// -------------------------------------------------------------------------
// importData: mirrors R's import().
// Reads the vector and optional bait files, sets multiBait, etc.
// -------------------------------------------------------------------------
DataDeepn importData(const QStringList &vn,
                     const QStringList &vs,
                     const QStringList &bn,
                     const QStringList &bs) {
    DataDeepn data;
    bool genesInitialized = false;

    // Vector
    data.Vector = readDeepn(vn, vs, data.geneNames, genesInitialized);

    // Bait (optional)
    if (!bn.isEmpty()) {
        data.Bait = readDeepn(bn, bs, data.geneNames, genesInitialized);
        // multiBait if Bait has replicates > 1
        data.multiBait = (data.Bait.nReplicates > 1);
    } else {
        data.Bait.nGenes = 0;
        data.Bait.nConds = 2;
        data.Bait.nReplicates = 0;
        data.Bait.data.clear();
        data.multiBait = false;
    }

    // Compute total read sums for Vector
    data.vtr.resize(data.Vector.nConds * data.Vector.nReplicates, 0.0);
    for (int c = 0; c < data.Vector.nConds; ++c) {
        for (int r = 0; r < data.Vector.nReplicates; ++r) {
            long long sumV = 0;
            for (int i = 0; i < data.Vector.nGenes; ++i) {
                sumV += data.Vector(i, c, r);
            }
            data.vtr[c * data.Vector.nReplicates + r] = (double) sumV;
        }
    }

    // Compute total read sums for Bait
    if (data.Bait.nGenes > 0) {
        data.btr.resize(data.Bait.nConds * data.Bait.nReplicates, 0.0);
        for (int c = 0; c < data.Bait.nConds; ++c) {
            for (int r = 0; r < data.Bait.nReplicates; ++r) {
                long long sumB = 0;
                for (int i = 0; i < data.Bait.nGenes; ++i) {
                    sumB += data.Bait(i, c, r);
                }
                data.btr[c * data.Bait.nReplicates + r] = (double) sumB;
            }
        }
    }

    return data;
}

// -------------------------------------------------------------------------
// rpm: convert raw integer counts to Reads Per Million (RPM).
//  - data: original DataDeepn
//  - rpmVector, rpmBait: outputs of size nGenes*nConds*nReps
// -------------------------------------------------------------------------
void rpm(const DataDeepn &data,
         std::vector<double> &rpmVector,
         std::vector<double> &rpmBait) {
    // Vector
    int vSize = data.Vector.nGenes * data.Vector.nConds * data.Vector.nReplicates;
    rpmVector.resize(vSize, 0.0);

    for (int i = 0; i < data.Vector.nGenes; ++i) {
        for (int c = 0; c < data.Vector.nConds; ++c) {
            for (int r = 0; r < data.Vector.nReplicates; ++r) {
                double denom = data.vtr[c * data.Vector.nReplicates + r];
                auto val = static_cast<double>(data.Vector(i, c, r));
                if (denom > 0) {
                    rpmVector[i * (data.Vector.nConds * data.Vector.nReplicates)
                              + c * data.Vector.nReplicates + r]
                            = (val / denom) * 1.0e6;
                }
            }
        }
    }

    // Bait
    if (data.Bait.nGenes > 0 && data.Bait.nReplicates > 0) {
        int bSize = data.Bait.nGenes * data.Bait.nConds * data.Bait.nReplicates;
        rpmBait.resize(bSize, 0.0);
        for (int i = 0; i < data.Bait.nGenes; ++i) {
            for (int c = 0; c < data.Bait.nConds; ++c) {
                for (int r = 0; r < data.Bait.nReplicates; ++r) {
                    double denom = data.btr[c * data.Bait.nReplicates + r];
                    auto val = static_cast<double>(data.Bait(i, c, r));
                    if (denom > 0) {
                        rpmBait[i * (data.Bait.nConds * data.Bait.nReplicates)
                                + c * data.Bait.nReplicates + r]
                                = (val / denom) * 1.0e6;
                    }
                }
            }
        }
    } else {
        rpmBait.clear();
    }
}

// -------------------------------------------------------------------------
// applyFilter: filter out genes whose RPM is below threshold.
//   - bPass: must exceed threshold in baseline for ALL (if base=TRUE)
//   - sPass: must exceed threshold in selected for ANY
// Then we keep only genes passing both.
// -------------------------------------------------------------------------
void applyFilter(DataDeepn &data, double threshold, bool base) {
    // Compute current RPM
    std::vector<double> rv, rb;
    rpm(data, rv, rb);

    auto getRpmVector = [&](int i, int cond, int rep) {
        return rv[i * (data.Vector.nConds * data.Vector.nReplicates)
                  + cond * data.Vector.nReplicates + rep];
    };
    auto getRpmBait = [&](int i, int cond, int rep) {
        return rb[i * (data.Bait.nConds * data.Bait.nReplicates)
                  + cond * data.Bait.nReplicates + rep];
    };

    std::vector<bool> bPass(data.Vector.nGenes, true);
    std::vector<bool> sPass(data.Vector.nGenes, false);

    for (int i = 0; i < data.Vector.nGenes; ++i) {
        bool passAll = true;
        bool passAny = false;

        if (base) {
            // Baseline => c=0
            // R code: if multiBait => combine mean of vector + all bait baseline
            // We'll do a direct approach:
            double meanV = 0.0;
            if (data.Vector.nReplicates > 0) {
                for (int r = 0; r < data.Vector.nReplicates; ++r) {
                    meanV += getRpmVector(i, 0, r);
                }
                meanV /= data.Vector.nReplicates;
            }
            bool passVector = (meanV > threshold);

            bool passBait = true;
            if (data.Bait.nGenes > 0 && i < data.Bait.nGenes) {
                if (data.multiBait) {
                    for (int r = 0; r < data.Bait.nReplicates; ++r) {
                        double val = getRpmBait(i, 0, r);
                        if (val <= threshold) {
                            passBait = false;
                            break;
                        }
                    }
                } else {
                    // single-bait => 1 replicate
                    double val = getRpmBait(i, 0, 0);
                    passBait = (val > threshold);
                }
            }
            passAll = (passVector && passBait);
        }

        // Selection => must exceed threshold in ANY replicate
        // S = vector selected + bait selected
        bool selectionPass = false;

        // Vector selected
        for (int r = 0; r < data.Vector.nReplicates; ++r) {
            if (getRpmVector(i, 1, r) > threshold) {
                selectionPass = true;
                break;
            }
        }
        // Bait selected
        if (!selectionPass && data.Bait.nGenes > 0 && i < data.Bait.nGenes) {
            if (data.multiBait) {
                for (int r = 0; r < data.Bait.nReplicates; ++r) {
                    if (getRpmBait(i, 1, r) > threshold) {
                        selectionPass = true;
                        break;
                    }
                }
            } else if (data.Bait.nReplicates > 0) {
                if (getRpmBait(i, 1, 0) > threshold) {
                    selectionPass = true;
                }
            }
        }

        passAny = selectionPass;

        bPass[i] = passAll;
        sPass[i] = passAny;
    }

    // Genes that pass both
    std::vector<int> passIndices;
    for (int i = 0; i < data.Vector.nGenes; ++i) {
        if (bPass[i] && sPass[i]) {
            passIndices.push_back(i);
        }
    }

    if ((int) passIndices.size() <= 1) {
        qWarning() << "Fewer than 2 genes pass this filter";
    }

    // Subset the data
    CountArray newVec;
    newVec.nGenes = passIndices.size();
    newVec.nConds = data.Vector.nConds;
    newVec.nReplicates = data.Vector.nReplicates;
    newVec.data.resize(newVec.nGenes * newVec.nConds * newVec.nReplicates, 0);

    CountArray newBait;
    if (data.Bait.nGenes > 0) {
        newBait.nGenes = passIndices.size();
        newBait.nConds = data.Bait.nConds;
        newBait.nReplicates = data.Bait.nReplicates;
        newBait.data.resize(newBait.nGenes * newBait.nConds * newBait.nReplicates, 0);
    }

    std::vector<std::string> newNames(passIndices.size());

    for (int idx = 0; idx < (int) passIndices.size(); ++idx) {
        int oldI = passIndices[idx];
        newNames[idx] = data.geneNames[oldI];

        for (int c = 0; c < newVec.nConds; ++c) {
            for (int r = 0; r < newVec.nReplicates; ++r) {
                newVec(idx, c, r) = data.Vector(oldI, c, r);
            }
        }
        if (data.Bait.nGenes > 0 && oldI < data.Bait.nGenes) {
            for (int c = 0; c < newBait.nConds; ++c) {
                for (int r = 0; r < newBait.nReplicates; ++r) {
                    newBait(idx, c, r) = data.Bait(oldI, c, r);
                }
            }
        }
    }

    data.Vector = newVec;
    data.Bait = newBait;
    data.geneNames = newNames;
}

/**
 * @brief buildSingleConditionData
 *
 * Extracts counts from a given CountArray 'arr' at condition 'cond',
 * building a CountData with a single group = 0,
 * i.e., no difference across replicates.
 *
 * offset[r] can be set to log(librarySize) or 0.0 if you prefer no offset.
 */
CountData buildSingleConditionData(const CountArray &arr, int cond) {
    if (cond < 0 || cond >= arr.nConds) {
        throw std::runtime_error("buildSingleConditionData: invalid cond index");
    }
    // We'll have 'arr.nGenes' genes and 'arr.nReplicates' replicates:
    int G = arr.nGenes;
    int R = arr.nReplicates;

    CountData data(G, R);
    // All replicates in group=0 => single group approach
    for (int r = 0; r < R; ++r) {
        data.group[r] = 0;
        // Optionally compute librarySize as sum of all genes in this condition & replicate
        // or you might have done that earlier. For a simple approach:
        double libSize = 0.0;
        for (int i = 0; i < G; ++i) {
            libSize += arr(i, cond, r);
        }
        if (libSize <= 0.0) {
            libSize = 1.0;
        }
        data.offset[r] = std::log(libSize);
    }

    // Fill in the counts
    for (int i = 0; i < G; ++i) {
        for (int r = 0; r < R; ++r) {
            data.y[i][r] = static_cast<double>(arr(i, cond, r));
        }
    }

    // For logAbundance[i], we can estimate a naive baseline:
    // e.g., the average count across replicates
    for (int i = 0; i < G; ++i) {
        double sumCounts = 0.0;
        for (int r = 0; r < R; ++r) {
            sumCounts += data.y[i][r];
        }
        double meanCount = sumCounts / (double) R;
        if (meanCount < 1.0) meanCount = 1.0;
        data.logAbundance[i] = std::log(meanCount);
        // logFC can remain 0.0 since we have only one group
    }

    return data;
}


/**
 * @brief estimateCommonDispSingle
 *
 * For a single condition in 'arr' (like baseline or selected),
 * build the single-group CountData and run the advanced Cox-Reid approach
 * to get a "common dispersion" for that condition.
 */
double estimateCommonDispSingle(const CountArray &arr, int cond,
                                double alphaStart = 0.1,
                                int maxIter = 100,
                                double tol = 1e-7) {
    // Build a single-condition CountData
    CountData cd = buildSingleConditionData(arr, cond);

    // Fit the dispersion
    double alpha = fitCommonDispersion(cd, alphaStart, maxIter, tol);
    return alpha;
}

// -------------------------------------------------------------------------
// combineBaseline: merges Vector & Bait baseline replicates into one
// condition to estimate "baitEffect" dispersion
// result => CountArray [ nGenes x 1 cond x (nVecReps + nBaitReps) ]
// -------------------------------------------------------------------------
static CountArray combineBaseline(const CountArray &vec, const CountArray &bait) {
    if (vec.nGenes != bait.nGenes) {
        qWarning() << "combineBaseline: mismatch in gene counts. Will proceed with min(#genes).";
    }
    int nG = std::min(vec.nGenes, bait.nGenes);
    int rV = vec.nReplicates;
    int rB = bait.nReplicates;

    CountArray out;
    out.nGenes = nG;
    out.nConds = 1; // single condition
    out.nReplicates = rV + rB;
    out.data.resize(nG * out.nConds * out.nReplicates, 0);

    // copy vector baseline
    for (int i = 0; i < nG; ++i) {
        for (int r = 0; r < rV; ++r) {
            // baseline cond=0 in vec
            out(i, 0, r) = vec(i, 0, r);
        }
    }
    // copy bait baseline
    for (int i = 0; i < nG; ++i) {
        for (int r = 0; r < rB; ++r) {
            out(i, 0, rV + r) = bait(i, 0, r);
        }
    }
    return out;
}

// -------------------------------------------------------------------------
// overdisp: returns [ baseline, selected, baitEffect ] using
// method-of-moments negative binomial approach
// -------------------------------------------------------------------------
std::vector<double> overdisp(const DataDeepn &data) {
    std::vector<double> w(3, 0.0);

    // 1) Baseline for Vector => cond=0
    double dispBase = 0.0;
    if (data.Vector.nGenes > 0 && data.Vector.nReplicates > 0) {
        dispBase = estimateCommonDispSingle(data.Vector, /*cond=*/0, 0.1, 100, 1e-7);
    }

    // 2) Selected for Vector => cond=1
    double dispSel = 0.0;
    if (data.Vector.nGenes > 0 && data.Vector.nReplicates > 0 && data.Vector.nConds > 1) {
        dispSel = estimateCommonDispSingle(data.Vector, /*cond=*/1, 0.1, 100, 1e-7);
    }

    // 3) BaitEffect => combine baseline replicates from Vector + Bait => single condition
    double dispBait = 0.0;
    if (data.Bait.nGenes > 0 && data.Bait.nReplicates > 0) {
        CountArray combined = combineBaseline(data.Vector, data.Bait);
        // Then we fit alpha for that single "condition=0"
        dispBait = estimateCommonDispSingle(combined, 0, 0.1, 100, 1e-7);
    }

    w[0] = dispBase;
    w[1] = dispSel;
    w[2] = dispBait;

    return w;
}

// -------------------------------------------------------------------------
// chooseFilter: tries thresholds from minRPM to maxRPM (log-spaced).
// For each threshold, calls applyFilter on a copy, records (#genes passing, overdispSelected).
// -------------------------------------------------------------------------
std::vector<QPair<double, QPair<int, double>>>
chooseFilter(const DataDeepn &data, double minRPM, double maxRPM) {
    // We'll log-space from minRPM to maxRPM in 19 steps
    int N = 19;
    double logMin = std::log(minRPM);
    double logMax = std::log(maxRPM);
    double step = (logMax - logMin) / (N - 1);

    std::vector<QPair<double, QPair<int, double>>> output;
    output.reserve(N);

    // We'll keep an original copy so we can revert data each iteration
    const DataDeepn &original = data;

    for (int j = 0; j < N; ++j) {
        double t = std::exp(logMin + step * j);

        // copy, filter
        DataDeepn copy = original;
        applyFilter(copy, t, true);
        int gCount = copy.Vector.nGenes;

        // overdisp => selected is index=1
        std::vector<double> w = overdisp(copy);
        double odSel = (w.size() >= 2 ? w[1] : 0.0);

        // store
        output.push_back(qMakePair(t, qMakePair(gCount, odSel)));
    }
    return output;
}

// -------------------------------------------------------------------------
// summarizeData: partial analog to summary.data.deepn() from R.
// - computes average baseline (Vector baseline + optional Bait baseline)
// - computes average Vector selected
// - if single-bait, we do Bait sel + log2 enrichment
// - if multi-bait, we do a naive approach to Bait replicates
// Writes CSV if outfile is not empty, otherwise prints to qDebug.
// -------------------------------------------------------------------------
void summarizeData(const DataDeepn &data, bool sort, const QString &outfile) {
    // Compute RPM
    std::vector<double> rv, rb;
    rpm(data, rv, rb);

    auto getRpmVector = [&](int i, int cond, int rep) {
        return rv[i * (data.Vector.nConds * data.Vector.nReplicates) + cond * data.Vector.nReplicates + rep];
    };
    auto getRpmBait = [&](int i, int cond, int rep) {
        return rb[i * (data.Bait.nConds * data.Bait.nReplicates) + cond * data.Bait.nReplicates + rep];
    };

    // We'll store rows in a small struct for sorting if needed
    struct RowData {
        std::string gene;
        double base;   // average vector baseline
        double vec;    // average vector selected
        double bait1;
        double bait2;
        double enr1;
        double enr2;

        QString toStringMulti() const {
            // Gene,Base,Vec,Bait1,Bait2,Enr1,Enr2
            return QString("%1,%2,%3,%4,%5,%6,%7")
                    .arg(QString::fromStdString(gene))
                    .arg(base)
                    .arg(vec)
                    .arg(bait1)
                    .arg(bait2)
                    .arg(enr1)
                    .arg(enr2);
        }

        QString toStringSingle() const {
            // Gene,Base,Vec,Bait,Enr
            return QString("%1,%2,%3,%4,%5")
                    .arg(QString::fromStdString(gene))
                    .arg(base)
                    .arg(vec)
                    .arg(bait1)  // use bait1 for single-bait
                    .arg(enr1);
        }
    };

    std::vector<RowData> table;
    table.reserve(data.Vector.nGenes);

    for (int i = 0; i < data.Vector.nGenes; ++i) {
        // average vector baseline (cond=0)
        double sumBase = 0.0;
        for (int r = 0; r < data.Vector.nReplicates; ++r) {
            sumBase += getRpmVector(i, 0, r);
        }
        double meanBase = (data.Vector.nReplicates > 0 ? sumBase / data.Vector.nReplicates : 0.0);

        // average vector selected (cond=1)
        double sumSel = 0.0;
        for (int r = 0; r < data.Vector.nReplicates; ++r) {
            sumSel += getRpmVector(i, 1, r);
        }
        double meanSel = (data.Vector.nReplicates > 0 ? sumSel / data.Vector.nReplicates : 0.0);

        RowData row;
        row.gene = data.geneNames[i];
        row.base = meanBase;
        row.vec = meanSel;
        row.bait1 = 0.0;
        row.bait2 = 0.0;
        row.enr1 = 0.0;
        row.enr2 = 0.0;

        // If we have Bait data and i < data.Bait.nGenes
        if (data.Bait.nGenes > 0 && i < data.Bait.nGenes) {
            if (!data.multiBait) {
                // single-bait => 1 replicate
                // average Bait selected
                if (data.Bait.nReplicates > 0) {
                    double sumSelB = 0.0;
                    for (int r = 0; r < data.Bait.nReplicates; ++r) {
                        sumSelB += getRpmBait(i, 1, r);
                    }
                    double meanBaitSel = sumSelB / data.Bait.nReplicates;
                    row.bait1 = meanBaitSel;
                    // log2((bait sel + 0.05)/(vec sel + 0.05))
                    row.enr1 = std::log2((meanBaitSel + 0.05) / (meanSel + 0.05));
                }
            } else {
                // multi-bait => we do a naive approach for replicate 0 vs. replicate 1
                if (data.Bait.nReplicates >= 2) {
                    double bsel0 = getRpmBait(i, 1, 0);
                    double bsel1 = getRpmBait(i, 1, 1);
                    row.bait1 = bsel0;
                    row.bait2 = bsel1;
                    row.enr1 = std::log2((bsel0 + 0.05) / (meanSel + 0.05));
                    row.enr2 = std::log2((bsel1 + 0.05) / (meanSel + 0.05));
                }
                // If we had more than 2 replicates for multiBait,
                // you'd generalize this logic or average them, etc.
            }
        }

        table.push_back(row);
    }

    // Sorting if requested
    if (sort) {
        // If multiBait => sort by enr1 or enr2 (like in the R script).
        // We'll say sortOption=1 => sort by enr1 descending, 2 => enr2 descending
        // but we only have a boolean. We’ll just pick enr1 for demonstration.
        if (data.multiBait) {
            std::sort(table.begin(), table.end(),
                      [](const RowData &a, const RowData &b) {
                          return a.enr1 > b.enr1; // descending
                      });
        } else {
            // single-bait => sort by enr1 descending
            std::sort(table.begin(), table.end(),
                      [](const RowData &a, const RowData &b) {
                          return a.enr1 > b.enr1;
                      });
        }
    }

    // Generate output lines
    QStringList lines;
    if (data.multiBait) {
        lines << "Gene,Base,Vec,Bait1,Bait2,Enr1,Enr2";
        for (auto &rd: table) {
            lines << rd.toStringMulti();
        }
    } else {
        lines << "Gene,Base,Vec,Bait,Enr";
        for (auto &rd: table) {
            lines << rd.toStringSingle();
        }
    }

    if (!outfile.isEmpty()) {
        QFile f(outfile);
        if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream ts(&f);
            for (auto &ln: lines) {
                ts << ln << "\n";
            }
            f.close();
        } else {
            qWarning() << "Cannot open" << outfile << "for writing.";
        }
    } else {
        // Print to debug
        for (auto &ln: lines) {
            qDebug().noquote() << ln;
        }
    }
}

// -------------------------------------------------------------------------
// analyzeDeepn: top-level function that mimics the R version.
//
// 1) Parse a config file with lines like:
//    Vector_Non=someFile.csv
//    Vector_Sel=someFile2.csv
//    Bait1_Non=... Bait1_Sel=...
//    Threshold=5
// 2) importData
// 3) applyFilter
// 4) overdisp
// 5) write messages
// 6) optionally debug => applyFilter(500)
// 7) summarizeData
// -------------------------------------------------------------------------
void analyzeDeepn(const QString &infile,
                  const QString &outfile,
                  const QString &msgfile,
                  bool debug,
                  int sortOption) {
    // Parse config
    QMap<QString, QString> config;
    {
        QFile f(infile);
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qWarning() << "Cannot open config file:" << infile;
            return;
        }
        QTextStream in(&f);
        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            int idx = line.indexOf('=');
            if (idx < 0) continue;
            QString key = line.left(idx).trimmed();
            QString val = line.mid(idx + 1).trimmed();
            config[key] = val;
        }
        f.close();
    }

    // Gather relevant keys
    QStringList vn, vs, bn, bs;
    double thresh = 5.0; // default
    for (auto it = config.begin(); it != config.end(); ++it) {
        const QString &key = it.key();
        const QString &val = it.value();
        if (key.startsWith("Vector_Non")) {
            vn << val;
        } else if (key.startsWith("Vector_Sel")) {
            vs << val;
        } else if (key.contains("Bait") && key.contains("_Non")) {
            bn << val;
        } else if (key.contains("Bait") && key.contains("_Sel")) {
            bs << val;
        } else if (key == "Threshold") {
            thresh = val.toDouble();
        }
    }

    // import
    DataDeepn data = importData(vn, vs, bn, bs);

    // apply filter with threshold
    applyFilter(data, thresh, true);

    // overdisp
    data.omega = overdisp(data);

    // Write a message file
    {
        QFile mf(msgfile);
        if (mf.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&mf);
            out << "Overdispersion estimates\n";
            if (data.omega.size() >= 3) {
                out << "Baseline (vector only):   " << data.omega[0] << "\n";
                out << "Baseline (vector + bait): " << data.omega[2] << "\n";
                out << "Selection:                " << data.omega[1] << "\n";
            }
            out << "-----\n";
            out << "Baseline overdispersions should be around 0 (no overdispersion)\n";
            out << "Overdispersion under selection conditions ideally under 2.\n";
            if (data.omega.size() >= 3) {
                if (data.omega[2] > 2.0 * data.omega[0] && data.omega[2] > 0.15) {
                    out << "\nWARNING: Evidence of secondary bait effects.\n"
                        << "Abundances in presence of bait differ from vector alone\n"
                        << "even in non-selective conditions.\n";
                }
            }
            mf.close();
        }
    }

    // If debug => applyFilter(data, 500)
    if (debug) {
        applyFilter(data, 500.0, true);
    }

    // Summarize
    summarizeData(data, (sortOption != 0), outfile);
}
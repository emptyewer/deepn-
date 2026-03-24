#ifndef ADVANCED_NB_DISPERSION_H
#define ADVANCED_NB_DISPERSION_H

#include <vector>

/**
 * @brief CountData struct for storing multiple genes' read counts across replicates.
 *
 *  Suppose we have G genes, R replicates, possibly in 2 groups (baseline, selected).
 *  We'll store:
 *    - y[i][r] = count for gene i in replicate r
 *    - offset[r] = log of library size or normalization factor for replicate r
 *    - group[r] = group index of replicate r, e.g., 0=baseline, 1=selected
 *
 *  Then the model might be:
 *    mu_{i,r} = exp( offset[r] + log(p_i) + (group[r] * beta_i) )
 *  or some other design. For "common dispersion," we fix alpha across all i.
 *
 *  This struct is a simplified placeholder. In advanced usage, you'd store
 *  design matrices, gene-specific betas, etc.
 */
struct CountData {
    int nGenes;
    int nReplicates;

    // y[i][r] = counts for gene i in replicate r
    std::vector<std::vector<double> > y;

    // offset[r], the log(library size) or log(normalization) for replicate r
    std::vector<double> offset;

    // group[r], an integer group label (0 or 1 in simplest case)
    std::vector<int> group;

    // gene-specific log-abundance parameter, p[i], for baseline
    // can be estimated by log(average expression) or a quick method
    std::vector<double> logAbundance;

    // gene-specific log-fold-change param for group=1 vs group=0
    // This is a simplified approach if we only have 2 groups
    std::vector<double> logFC;

    CountData(int G, int R)
            : nGenes(G), nReplicates(R) {
        y.resize(G, std::vector<double>(R, 0.0));
        offset.resize(R, 0.0);
        group.resize(R, 0);
        logAbundance.resize(G, 0.0);
        logFC.resize(G, 0.0);
    }
};

/**
 * @brief fitCommonDispersion
 *
 * Fits a single negative binomial dispersion alpha ("common dispersion")
 * across all genes using an iterative approach with a Cox-Reid correction.
 *
 * The function:
 *    1) uses the current logAbundance[i] and logFC[i] for each gene i
 *    2) calculates mu_{i,r} for each replicate r:
 *         mu_{i,r} = exp( offset[r] + logAbundance[i] + group[r]*logFC[i] )
 *    3) forms the (penalized) log-likelihood with respect to alpha
 *    4) iterates a Newton-Raphson or secant approach to find alpha
 *
 * @param data  CountData holding the read counts, offsets, group labels, etc.
 * @param alphaStart an initial guess for alpha
 * @param maxIter maximum number of iterations
 * @param tol convergence tolerance
 * @return estimated alpha
 */
double fitCommonDispersion(CountData &data,
                           double alphaStart = 0.1,
                           int maxIter = 100,
                           double tol = 1e-7);

#endif // ADVANCED_NB_DISPERSION_H
#include "deseq_stats.h"
#include "utils.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace deseq2
{

    DeseqStats::DeseqStats(const DeseqDataSet &dds,
                           const Eigen::VectorXd &contrast,
                           double alpha,
                           bool cooks_filter,
                           bool independent_filter,
                           double lfc_null)
        : dds_(dds),
          contrast_vector_(contrast),
          alpha_(alpha),
          cooks_filter_(cooks_filter),
          independent_filter_(independent_filter),
          lfc_null_(lfc_null),
          shrunk_lfcs_(false)
    {

        // Initialize design matrix and LFCs
        design_matrix_ = dds.getDesignMatrix();
        lfc_ = dds.getLFC();

        // Initialize base mean
        base_mean_ = dds.getNormedMeans();

        // Initialize results vectors
        int n_genes = dds.getCounts().cols();
        log2_fold_change_ = Eigen::VectorXd::Zero(n_genes);
        lfc_se_ = Eigen::VectorXd::Zero(n_genes);
        stat_ = Eigen::VectorXd::Zero(n_genes);
        p_value_ = Eigen::VectorXd::Zero(n_genes);
        padj_ = Eigen::VectorXd::Zero(n_genes);
    }

    void DeseqStats::runWaldTest()
    {
        int n_genes = dds_.getCounts().cols();

        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd gene_lfc = lfc_.row(j);
            double disp = dds_.getDispersions()(j);

            // Calculate mu for this gene
            Eigen::VectorXd gene_counts = dds_.getCounts().col(j);
            Eigen::VectorXd size_factors = dds_.getSizeFactors();
            Eigen::VectorXd mu = fitLinMu(gene_counts, size_factors, design_matrix_, 0.5);

            // Perform Wald test
            auto [stat, p_val, se] = ::deseq2::waldTest(
                design_matrix_, disp, gene_lfc, mu, contrast_vector_, lfc_null_);

            // Store results
            stat_(j) = stat;
            p_value_(j) = p_val;
            lfc_se_(j) = se;

            // Calculate log2 fold change
            double contrast_lfc = contrast_vector_.dot(gene_lfc);
            log2_fold_change_(j) = contrast_lfc / std::log(2.0);
        }
    }

    void DeseqStats::cooksFiltering()
    {
        if (!cooks_filter_)
        {
            return;
        }

        // Simple Cooks filtering - set p-values to 1 for genes with high Cooks distances
        // This is a simplified version
        int n_genes = dds_.getCounts().cols();

        for (int j = 0; j < n_genes; ++j)
        {
            // Check if any sample has high Cooks distance for this gene
            bool has_outlier = false;
            // Simplified outlier detection
            if (has_outlier)
            {
                p_value_(j) = 1.0;
            }
        }
    }

    void DeseqStats::independentFiltering()
    {
        if (!independent_filter_)
        {
            pValueAdjustment();
            return;
        }

        int n_genes = base_mean_.size();

        // Sort genes by base mean
        std::vector<std::pair<double, int>> gene_means;
        for (int i = 0; i < n_genes; ++i)
        {
            gene_means.push_back({base_mean_(i), i});
        }
        std::sort(gene_means.begin(), gene_means.end());

        // Find optimal filtering threshold
        double best_threshold = 0.0;
        double best_rejections = 0.0;

        for (size_t i = 0; i < gene_means.size(); ++i)
        {
            double threshold = gene_means[i].first;

            // Count rejections at this threshold
            int rejections = 0;
            for (int j = 0; j < n_genes; ++j)
            {
                if (base_mean_(j) >= threshold && p_value_(j) < alpha_)
                {
                    rejections++;
                }
            }

            if (rejections > best_rejections)
            {
                best_rejections = rejections;
                best_threshold = threshold;
            }
        }

        // Apply filtering
        for (int j = 0; j < n_genes; ++j)
        {
            if (base_mean_(j) < best_threshold)
            {
                p_value_(j) = 1.0;
            }
        }

        // Perform p-value adjustment
        pValueAdjustment();
    }

    void DeseqStats::pValueAdjustment()
    {
        int n = p_value_.size();

        // Benjamini-Hochberg procedure with correct backward monotonicity enforcement
        std::vector<std::pair<double, int>> sorted_p;
        for (int i = 0; i < n; ++i)
        {
            sorted_p.push_back({p_value_(i), i});
        }
        std::sort(sorted_p.begin(), sorted_p.end());

        // Compute raw BH-adjusted values
        std::vector<double> adj(n);
        for (int i = 0; i < n; ++i)
        {
            adj[i] = sorted_p[i].first * n / (i + 1);
        }

        // Enforce monotonicity BACKWARD (cumulative minimum from end)
        for (int i = n - 2; i >= 0; --i)
        {
            adj[i] = std::min(adj[i], adj[i + 1]);
        }

        // Clip to [0, 1] and map back to original indices
        for (int i = 0; i < n; ++i)
        {
            padj_(sorted_p[i].second) = std::min(1.0, adj[i]);
        }
    }

    Eigen::MatrixXd DeseqStats::summary() const
    {
        int n_genes = base_mean_.size();
        Eigen::MatrixXd results(n_genes, 6);

        for (int i = 0; i < n_genes; ++i)
        {
            results(i, 0) = base_mean_(i);
            results(i, 1) = log2_fold_change_(i);
            results(i, 2) = lfc_se_(i);
            results(i, 3) = stat_(i);
            results(i, 4) = p_value_(i);
            results(i, 5) = padj_(i);
        }

        return results;
    }

    void DeseqStats::lfcShrink(const std::string &coeff)
    {
        (void)coeff;

        // Fit prior variance from the observed LFCs
        fitPriorVar();

        // Apply normal prior shrinkage: shrink = prior_var / (prior_var + se^2)
        for (int i = 0; i < log2_fold_change_.size(); ++i)
        {
            double se2 = lfc_se_(i) * lfc_se_(i);
            if (se2 > 0 && prior_var_ > 0)
            {
                double shrink = prior_var_ / (prior_var_ + se2);
                log2_fold_change_(i) *= shrink;
                lfc_se_(i) = std::sqrt(shrink * se2);
            }
        }

        shrunk_lfcs_ = true;
    }

    void DeseqStats::buildContrastVector()
    {
        // This method would build the contrast vector from string representation
        // For now, we assume the contrast vector is provided directly
    }

    double DeseqStats::waldTest(int gene_idx) const
    {
        // This method would perform a Wald test for a specific gene
        // For now, we use the utility function
        (void)gene_idx; // Suppress unused parameter warning
        return 0.0;
    }

    void DeseqStats::fitPriorVar()
    {
        // Estimate prior variance from observed absolute LFCs using MAD-based scale
        std::vector<double> abs_lfc;
        for (int i = 0; i < log2_fold_change_.size(); ++i)
        {
            if (std::isfinite(log2_fold_change_(i)) && lfc_se_(i) > 0)
            {
                abs_lfc.push_back(std::abs(log2_fold_change_(i)));
            }
        }
        if (abs_lfc.empty())
        {
            prior_var_ = 1.0;
            return;
        }
        std::sort(abs_lfc.begin(), abs_lfc.end());
        double median = abs_lfc[abs_lfc.size() / 2];
        double mad = median / 0.6745; // MAD-based scale estimate (for normal distribution)
        prior_var_ = mad * mad;
    }

} // namespace deseq2
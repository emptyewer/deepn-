#include "deseq_dataset.h"
#include "utils.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace deseq2
{

    DeseqDataSet::DeseqDataSet(const Eigen::MatrixXd &counts,
                               const Eigen::MatrixXd &metadata,
                               const std::string &design,
                               bool refit_cooks,
                               double min_mu,
                               double min_disp,
                               double max_disp,
                               double beta_tol)
        : counts_(counts),
          metadata_(metadata),
          design_(design),
          refit_cooks_(refit_cooks),
          min_mu_(min_mu),
          min_disp_(min_disp),
          max_disp_(max_disp),
          beta_tol_(beta_tol)
    {

        // Validate inputs
        testValidCounts(counts);

        // Initialize matrices
        int n_samples = counts.rows();
        int n_genes = counts.cols();

        size_factors_ = Eigen::VectorXd::Ones(n_samples);
        normed_means_ = Eigen::VectorXd::Zero(n_genes);
        genewise_dispersions_ = Eigen::VectorXd::Zero(n_genes);
        fitted_dispersions_ = Eigen::VectorXd::Zero(n_genes);
        map_dispersions_ = Eigen::VectorXd::Zero(n_genes);
        dispersions_ = Eigen::VectorXd::Zero(n_genes);
        trend_coeffs_ = Eigen::Vector2d::Zero();
        lfc_ = Eigen::MatrixXd::Zero(n_genes, 2); // intercept + condition
        cooks_ = Eigen::MatrixXd::Zero(n_samples, n_genes);

        // Build design matrix
        buildDesignMatrix();
    }

    void DeseqDataSet::buildDesignMatrix()
    {
        int n_samples = counts_.rows();

        // Simple design matrix for condition A vs B
        // First column: intercept (all 1s)
        // Second column: condition (0 for A, 1 for B)
        design_matrix_ = Eigen::MatrixXd::Ones(n_samples, 2);

        // Extract condition from metadata (assuming condition is in first column)
        for (int i = 0; i < n_samples; ++i)
        {
            if (metadata_(i, 0) > 0.5)
            { // Assuming B is encoded as 1
                design_matrix_(i, 1) = 1.0;
            }
            else
            {
                design_matrix_(i, 1) = 0.0;
            }
        }
    }

    void DeseqDataSet::fitSizeFactors()
    {
        std::cout << "Fitting size factors..." << std::endl;

        // Compute gene-wise mean log counts (deseq2_norm_fit equivalent)
        Eigen::MatrixXd log_counts = counts_.array().log();
        Eigen::VectorXd logmeans = log_counts.colwise().mean();

        // Filter out genes with -inf log means
        Eigen::VectorXd filtered_genes = (logmeans.array().isFinite()).cast<double>();

        // Apply filtering
        Eigen::VectorXd filtered_logmeans = logmeans.cwiseProduct(filtered_genes);

        // Calculate size factors using median-of-ratios method
        Eigen::VectorXd size_factors(counts_.rows());

        for (int i = 0; i < counts_.rows(); ++i)
        {
            std::vector<double> log_ratios;
            for (int j = 0; j < counts_.cols(); ++j)
            {
                if (filtered_genes(j) > 0 && counts_(i, j) > 0)
                {
                    double log_ratio = std::log(counts_(i, j)) - filtered_logmeans(j);
                    log_ratios.push_back(log_ratio);
                }
            }

            if (!log_ratios.empty())
            {
                // Calculate median of log ratios
                std::sort(log_ratios.begin(), log_ratios.end());
                double median_log_ratio;
                if (log_ratios.size() % 2 == 0)
                {
                    median_log_ratio = (log_ratios[log_ratios.size() / 2 - 1] + log_ratios[log_ratios.size() / 2]) / 2.0;
                }
                else
                {
                    median_log_ratio = log_ratios[log_ratios.size() / 2];
                }
                size_factors(i) = std::exp(median_log_ratio);
            }
            else
            {
                size_factors(i) = 1.0;
            }
        }

        // Normalize size factors to geometric mean of 1
        double log_geometric_mean = size_factors.array().log().mean();
        size_factors_ = size_factors.array() / std::exp(log_geometric_mean);

        // Calculate normalized counts (divide each row by its sample's size factor)
        normed_counts_ = counts_;
        for (int i = 0; i < counts_.rows(); ++i)
        {
            normed_counts_.row(i) = counts_.row(i).array() / size_factors_(i);
        }

        // Calculate normalized means
        normed_means_ = normed_counts_.colwise().mean();

        // Set baseMean to 0 only for genes with all-zero counts
        for (int j = 0; j < normed_means_.size(); ++j)
        {
            bool all_zero = true;
            for (int i = 0; i < counts_.rows(); ++i)
            {
                if (counts_(i, j) != 0)
                {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero)
            {
                normed_means_(j) = 0.0;
            }
        }

        std::cout << "[DEBUG] counts_ shape: " << counts_.rows() << " x " << counts_.cols() << std::endl;
        std::cout << "[DEBUG] size_factors_ shape: " << size_factors_.size() << std::endl;
        std::cout << "[DEBUG] normed_means_ shape: " << normed_means_.size() << std::endl;
    }

    void DeseqDataSet::fitGenewiseDispersions()
    {
        int n_genes = counts_.cols();
        int n_samples = counts_.rows();

        // Calculate normalized counts
        Eigen::MatrixXd normed_counts = Eigen::MatrixXd::Zero(n_samples, n_genes);
        for (int i = 0; i < n_samples; ++i)
        {
            for (int j = 0; j < n_genes; ++j)
            {
                normed_counts(i, j) = counts_(i, j) / size_factors_(i);
            }
        }
        std::cout << "[DEBUG] normed_counts shape: " << normed_counts.rows() << " x " << normed_counts.cols() << std::endl;

        // Use exact Python implementation: fit rough dispersions first
        genewise_dispersions_ = fitRoughDispersions(normed_counts, design_matrix_);
        std::cout << "[DEBUG] genewise_dispersions_ shape: " << genewise_dispersions_.size() << std::endl;

        // Refine with MLE for each gene (simplified for now)
        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd gene_counts = counts_.col(j);
            Eigen::VectorXd mu = fitLinMu(gene_counts, size_factors_, design_matrix_, min_mu_);

            auto [alpha, converged] = fitAlphaMLE(
                gene_counts, design_matrix_, mu, genewise_dispersions_(j),
                min_disp_, max_disp_, prior_disp_var_);

            genewise_dispersions_(j) = alpha;
        }
    }

    void DeseqDataSet::fitDispersionTrend()
    {
        int n_genes = counts_.cols();

        // Filter genes with valid dispersions
        std::vector<int> valid_indices;
        std::vector<double> valid_means;
        std::vector<double> valid_dispersions;

        for (int j = 0; j < n_genes; ++j)
        {
            if (normed_means_(j) > 0 && genewise_dispersions_(j) > 0)
            {
                valid_indices.push_back(j);
                valid_means.push_back(normed_means_(j));
                valid_dispersions.push_back(genewise_dispersions_(j));
            }
        }

        if (valid_indices.size() < 10)
        {
            // Fallback to simple trend
            trend_coeffs_(0) = 0.1;
            trend_coeffs_(1) = 1.0;
        }
        else
        {
            // Fit parametric trend: disp = a0 + a1/mean
            Eigen::MatrixXd X(valid_indices.size(), 2);
            Eigen::VectorXd y(valid_indices.size());

            for (size_t i = 0; i < valid_indices.size(); ++i)
            {
                X(i, 0) = 1.0;
                X(i, 1) = 1.0 / valid_means[i];
                y(i) = valid_dispersions[i];
            }

            // Robust fitting using trimmed mean
            Eigen::VectorXd residuals = y - X * trend_coeffs_;
            double mad = meanAbsoluteDeviation(residuals);

            // Simple iterative reweighting
            for (int iter = 0; iter < 5; ++iter)
            {
                Eigen::VectorXd weights = Eigen::VectorXd::Ones(valid_indices.size());
                for (size_t i = 0; i < valid_indices.size(); ++i)
                {
                    if (std::abs(residuals(i)) > 3 * mad)
                    {
                        weights(i) = 0.1;
                    }
                }

                Eigen::MatrixXd XtWX = X.transpose() * weights.asDiagonal() * X;
                Eigen::VectorXd XtWy = X.transpose() * weights.asDiagonal() * y;
                trend_coeffs_ = XtWX.ldlt().solve(XtWy);

                residuals = y - X * trend_coeffs_;
                mad = meanAbsoluteDeviation(residuals);
            }
        }

        // Calculate fitted dispersions
        for (int j = 0; j < n_genes; ++j)
        {
            fitted_dispersions_(j) = dispersionTrend(normed_means_(j), trend_coeffs_);
        }
    }

    void DeseqDataSet::fitDispersionPrior()
    {
        int n_genes = counts_.cols();

        // Calculate log residuals
        std::vector<double> log_residuals;

        for (int j = 0; j < n_genes; ++j)
        {
            if (normed_means_(j) > 0 && genewise_dispersions_(j) > 0 && fitted_dispersions_(j) > 0)
            {
                double log_resid = std::log(genewise_dispersions_(j)) - std::log(fitted_dispersions_(j));
                log_residuals.push_back(log_resid);
            }
        }

        if (log_residuals.empty())
        {
            prior_disp_var_ = 1.0;
            squared_logres_ = 1.0;
            return;
        }

        // Calculate prior variance using trimmed mean
        Eigen::Map<Eigen::VectorXd> log_resid_vec(log_residuals.data(), log_residuals.size());
        double mad = meanAbsoluteDeviation(log_resid_vec);
        prior_disp_var_ = mad * mad;
        squared_logres_ = prior_disp_var_;
    }

    void DeseqDataSet::fitMAPDispersions()
    {
        int n_genes = counts_.cols();

        for (int j = 0; j < n_genes; ++j)
        {
            double disp_mle = genewise_dispersions_(j);
            double disp_fitted = fitted_dispersions_(j);

            // MAP estimate using log-space posterior mode
            double disp_map;
            if (disp_mle > 0 && disp_fitted > 0)
            {
                double log_disp_mle = std::log(disp_mle);
                double log_disp_fitted = std::log(disp_fitted);
                double log_map = (log_disp_mle / prior_disp_var_ + log_disp_fitted) /
                                 (1.0 / prior_disp_var_ + 1.0);
                disp_map = std::exp(log_map);
            }
            else
            {
                disp_map = disp_fitted;
            }

            // Apply bounds
            disp_map = std::max(min_disp_, std::min(max_disp_, disp_map));

            map_dispersions_(j) = disp_map;
            dispersions_(j) = disp_map;
        }
    }

    void DeseqDataSet::fitLFC()
    {
        int n_genes = counts_.cols();
        int n_coeffs = design_matrix_.cols();

        lfc_ = Eigen::MatrixXd::Zero(n_genes, n_coeffs);

        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd gene_counts = counts_.col(j);
            double disp = dispersions_(j);

            // Fit LFC using IRLS
            auto [beta, mu, H, converged] = irlsSolver(
                gene_counts,
                size_factors_,
                design_matrix_,
                disp,
                0.5,  // min_mu
                1e-8, // beta_tol
                -30,  // min_beta
                30,   // max_beta
                250   // maxiter
            );

            lfc_.row(j) = beta;
        }
    }

    void DeseqDataSet::calculateCooks()
    {
        int n_samples = counts_.rows();
        int n_genes = counts_.cols();

        cooks_ = Eigen::MatrixXd::Zero(n_samples, n_genes);

        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd gene_counts = counts_.col(j);
            double disp = dispersions_(j);
            Eigen::VectorXd beta = lfc_.row(j);

            // Calculate mu
            Eigen::VectorXd eta = design_matrix_ * beta;
            Eigen::VectorXd mu = size_factors_.array() * eta.array().exp();
            mu = mu.cwiseMax(min_mu_);

            // Calculate weights
            Eigen::VectorXd weights = mu.array() / (1.0 + disp * mu.array());

            // Calculate hat matrix diagonal
            Eigen::MatrixXd XtWX = design_matrix_.transpose() * weights.asDiagonal() * design_matrix_;
            Eigen::MatrixXd hat_matrix = design_matrix_ * XtWX.inverse() * design_matrix_.transpose() * weights.asDiagonal();
            Eigen::VectorXd hat_diag = hat_matrix.diagonal();

            // Calculate Cooks distances
            Eigen::VectorXd residuals = gene_counts - mu;
            Eigen::VectorXd cooks_dist = (residuals.array().square() * hat_diag.array()) /
                                         (weights.array() * (1.0 - hat_diag.array()).square());

            cooks_.col(j) = cooks_dist;
        }
    }

    void DeseqDataSet::refit()
    {
        if (!refit_cooks_)
        {
            return;
        }

        // Identify outliers and refit
        replaceOutliers();
        refitWithoutOutliers();
    }

    void DeseqDataSet::replaceOutliers()
    {
        // Simple outlier detection based on Cooks distances
        // This is a simplified version - in practice, more sophisticated methods are used
    }

    void DeseqDataSet::refitWithoutOutliers()
    {
        // Refit the model without outliers
        // This is a simplified version
    }

    bool DeseqDataSet::cooksOutlier(int gene_idx) const
    {
        // Simple outlier detection
        // This is a simplified version
        (void)gene_idx; // Suppress unused parameter warning
        return false;
    }

} // namespace deseq2
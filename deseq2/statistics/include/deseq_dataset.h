#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

namespace deseq2
{

    /**
     * @brief A class to implement dispersion and log fold-change (LFC) estimation.
     *
     * This class implements the DESeq2 pipeline for differential expression analysis.
     * It handles size factor normalization, dispersion estimation, and log fold-change
     * calculation following the DESeq2 methodology.
     */
    class DeseqDataSet
    {
    public:
        /**
         * @brief Constructor for DeseqDataSet
         * @param counts Raw count matrix (samples x genes)
         * @param metadata Sample metadata matrix
         * @param design Design formula string (e.g., "~condition")
         * @param refit_cooks Whether to refit cooks outliers
         * @param min_mu Minimum mean threshold
         * @param min_disp Minimum dispersion threshold
         * @param max_disp Maximum dispersion threshold
         * @param beta_tol Stopping criterion for IRWLS
         */
        DeseqDataSet(const Eigen::MatrixXd &counts,
                     const Eigen::MatrixXd &metadata,
                     const std::string &design = "~condition",
                     bool refit_cooks = true,
                     double min_mu = 0.5,
                     double min_disp = 1e-8,
                     double max_disp = 10.0,
                     double beta_tol = 1e-8);

        /**
         * @brief Fit size factors using median-of-ratios method
         */
        void fitSizeFactors();

        /**
         * @brief Fit genewise dispersions
         */
        void fitGenewiseDispersions();

        /**
         * @brief Fit dispersion trend coefficients
         */
        void fitDispersionTrend();

        /**
         * @brief Fit dispersion prior
         */
        void fitDispersionPrior();

        /**
         * @brief Fit MAP dispersions
         */
        void fitMAPDispersions();

        /**
         * @brief Fit log fold changes
         */
        void fitLFC();

        /**
         * @brief Calculate Cooks distances
         */
        void calculateCooks();

        /**
         * @brief Refit model after outlier removal
         */
        void refit();

        // Getters
        const Eigen::VectorXd &getSizeFactors() const { return size_factors_; }
        const Eigen::VectorXd &getGenewiseDispersions() const { return genewise_dispersions_; }
        const Eigen::VectorXd &getFittedDispersions() const { return fitted_dispersions_; }
        const Eigen::VectorXd &getMAPDispersions() const { return map_dispersions_; }
        const Eigen::VectorXd &getDispersions() const { return dispersions_; }
        const Eigen::MatrixXd &getLFC() const { return lfc_; }
        const Eigen::MatrixXd &getDesignMatrix() const { return design_matrix_; }
        const Eigen::VectorXd &getNormedMeans() const { return normed_means_; }
        const Eigen::MatrixXd &getCounts() const { return counts_; }
        const Eigen::MatrixXd &getNormedCounts() const { return normed_counts_; }

    private:
        // Data matrices
        Eigen::MatrixXd counts_;        // Raw counts (samples x genes)
        Eigen::MatrixXd metadata_;      // Sample metadata
        Eigen::MatrixXd design_matrix_; // Design matrix
        Eigen::VectorXd size_factors_;  // Size factors
        Eigen::VectorXd normed_means_;  // Normalized means
        Eigen::MatrixXd normed_counts_; // Normalized counts

        // Dispersion parameters
        Eigen::VectorXd genewise_dispersions_; // Gene-wise dispersion estimates
        Eigen::VectorXd fitted_dispersions_;   // Fitted dispersions from trend
        Eigen::VectorXd map_dispersions_;      // Maximum a posteriori dispersions
        Eigen::VectorXd dispersions_;          // Final dispersion values
        Eigen::Vector2d trend_coeffs_;         // Dispersion trend coefficients

        // Log fold changes
        Eigen::MatrixXd lfc_; // Log fold changes (genes x coefficients)

        // Cooks distances
        Eigen::MatrixXd cooks_; // Cooks distances

        // Parameters
        std::string design_; // Design formula
        bool refit_cooks_;   // Whether to refit cooks outliers
        double min_mu_;      // Minimum mean threshold
        double min_disp_;    // Minimum dispersion threshold
        double max_disp_;    // Maximum dispersion threshold
        double beta_tol_;    // Beta tolerance for IRWLS

        // Prior parameters
        double prior_disp_var_; // Prior dispersion variance
        double squared_logres_; // Squared log residuals

        // Helper methods
        void buildDesignMatrix();
        void fitMoMDispersions();
        void fitParametricDispersionTrend();
        void replaceOutliers();
        void refitWithoutOutliers();
        bool cooksOutlier(int gene_idx) const;
    };

} // namespace deseq2
#pragma once

#include "deseq_dataset.h"
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace deseq2
{

    /**
     * @brief Statistical tests for differential expression analysis.
     *
     * This class implements p-value estimation for differential gene expression
     * according to the DESeq2 pipeline, including Wald tests and p-value adjustment.
     */
    class DeseqStats
    {
    public:
        /**
         * @brief Constructor for DeseqStats
         * @param dds Fitted DeseqDataSet object
         * @param contrast Contrast vector for testing
         * @param alpha Significance threshold for p-values
         * @param cooks_filter Whether to filter p-values based on cooks outliers
         * @param independent_filter Whether to perform independent filtering
         * @param lfc_null Log fold change under null hypothesis
         */
        DeseqStats(const DeseqDataSet &dds,
                   const Eigen::VectorXd &contrast,
                   double alpha = 0.05,
                   bool cooks_filter = true,
                   bool independent_filter = true,
                   double lfc_null = 0.0);

        /**
         * @brief Run Wald test for differential expression
         */
        void runWaldTest();

        /**
         * @brief Perform Cooks filtering
         */
        void cooksFiltering();

        /**
         * @brief Perform independent filtering
         */
        void independentFiltering();

        /**
         * @brief Perform p-value adjustment
         */
        void pValueAdjustment();

        /**
         * @brief Generate summary results
         * @return Results matrix with statistics
         */
        Eigen::MatrixXd summary() const;

        /**
         * @brief Perform LFC shrinkage
         * @param coeff Coefficient to shrink
         */
        void lfcShrink(const std::string &coeff);

        // Getters
        const Eigen::VectorXd &getBaseMean() const { return base_mean_; }
        const Eigen::VectorXd &getLog2FoldChange() const { return log2_fold_change_; }
        const Eigen::VectorXd &getLfcSE() const { return lfc_se_; }
        const Eigen::VectorXd &getStat() const { return stat_; }
        const Eigen::VectorXd &getPValue() const { return p_value_; }
        const Eigen::VectorXd &getPAdj() const { return padj_; }
        const Eigen::VectorXd &getContrastVector() const { return contrast_vector_; }

    private:
        const DeseqDataSet &dds_;         // Reference to DeseqDataSet
        Eigen::VectorXd contrast_vector_; // Contrast vector
        Eigen::MatrixXd design_matrix_;   // Design matrix
        Eigen::MatrixXd lfc_;             // Log fold changes

        // Results
        Eigen::VectorXd base_mean_;        // Base mean
        Eigen::VectorXd log2_fold_change_; // Log2 fold changes
        Eigen::VectorXd lfc_se_;           // LFC standard errors
        Eigen::VectorXd stat_;             // Wald statistics
        Eigen::VectorXd p_value_;          // P-values
        Eigen::VectorXd padj_;             // Adjusted p-values

        // Parameters
        double alpha_;            // Significance threshold
        bool cooks_filter_;       // Whether to use cooks filtering
        bool independent_filter_; // Whether to use independent filtering
        double lfc_null_;         // Null hypothesis LFC
        bool shrunk_lfcs_;        // Whether LFCs are shrunk
        double prior_var_ = 1.0;  // Prior variance for LFC shrinkage

        // Helper methods
        void buildContrastVector();
        double waldTest(int gene_idx) const;
        void fitPriorVar();
    };

} // namespace deseq2
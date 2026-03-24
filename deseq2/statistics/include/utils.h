#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>

namespace deseq2
{

    /**
     * @brief Load synthetic example data
     * @param modality Data modality ("raw_counts" or "metadata")
     * @param dataset Dataset name ("synthetic")
     * @return Data matrix
     */
    Eigen::MatrixXd loadExampleData(const std::string &modality = "raw_counts",
                                    const std::string &dataset = "synthetic");

    /**
     * @brief Test that count matrix contains valid inputs
     * @param counts Count matrix
     */
    void testValidCounts(const Eigen::MatrixXd &counts);

    /**
     * @brief Calculate dispersion trend
     * @param normed_mean Normalized mean
     * @param coeffs Trend coefficients
     * @return Dispersion trend value
     */
    double dispersionTrend(double normed_mean, const Eigen::Vector2d &coeffs);

    /**
     * @brief Negative binomial negative log-likelihood
     * @param counts Observed counts
     * @param mu Mean parameter
     * @param alpha Dispersion parameter
     * @return Negative log-likelihood
     */
    double nbNLL(const Eigen::VectorXd &counts, const Eigen::VectorXd &mu, double alpha);

    /**
     * @brief IRLS solver for negative binomial GLM
     * @param counts Observed counts
     * @param size_factors Size factors
     * @param design_matrix Design matrix
     * @param disp Dispersion parameter
     * @param min_mu Minimum mean threshold
     * @param beta_tol Beta tolerance
     * @param min_beta Minimum beta threshold
     * @param max_beta Maximum beta threshold
     * @param maxiter Maximum number of iterations
     * @return Tuple of (beta, mu, converged)
     */
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, bool> irlsSolver(
        const Eigen::VectorXd &counts,
        const Eigen::VectorXd &size_factors,
        const Eigen::MatrixXd &design_matrix,
        double disp,
        double min_mu = 0.5,
        double beta_tol = 1e-8,
        double min_beta = -30,
        double max_beta = 30,
        int maxiter = 250);

    /**
     * @brief Fit alpha MLE
     * @param counts Observed counts
     * @param design_matrix Design matrix
     * @param mu Mean estimates
     * @param alpha_hat Initial alpha estimate
     * @param min_disp Minimum dispersion
     * @param max_disp Maximum dispersion
     * @param prior_disp_var Prior dispersion variance
     * @return Tuple of (alpha, converged)
     */
    std::tuple<double, bool> fitAlphaMLE(
        const Eigen::VectorXd &counts,
        const Eigen::MatrixXd &design_matrix,
        const Eigen::VectorXd &mu,
        double alpha_hat,
        double min_disp,
        double max_disp,
        double prior_disp_var = 0.0);

    /**
     * @brief Trimmed mean
     * @param x Input vector
     * @param trim Trim proportion
     * @return Trimmed mean
     */
    double trimmedMean(const Eigen::VectorXd &x, double trim = 0.1);

    /**
     * @brief Fit linear model for mu
     * @param counts Observed counts
     * @param size_factors Size factors
     * @param design_matrix Design matrix
     * @param min_mu Minimum mean threshold
     * @return Fitted mu values
     */
    Eigen::VectorXd fitLinMu(
        const Eigen::VectorXd &counts,
        const Eigen::VectorXd &size_factors,
        const Eigen::MatrixXd &design_matrix,
        double min_mu = 0.5);

    /**
     * @brief Wald test
     * @param design_matrix Design matrix
     * @param disp Dispersion parameter
     * @param lfc Log fold change
     * @param mu Mean estimates
     * @param contrast Contrast vector
     * @param lfc_null Null hypothesis LFC
     * @return Tuple of (statistic, p_value, lfc_se)
     */
    std::tuple<double, double, double> waldTest(
        const Eigen::MatrixXd &design_matrix,
        double disp,
        const Eigen::VectorXd &lfc,
        const Eigen::VectorXd &mu,
        const Eigen::VectorXd &contrast,
        double lfc_null = 0.0);

    /**
     * @brief Fit rough dispersions
     * @param normed_counts Normalized counts
     * @param design_matrix Design matrix
     * @return Rough dispersion estimates
     */
    Eigen::VectorXd fitRoughDispersions(
        const Eigen::MatrixXd &normed_counts,
        const Eigen::MatrixXd &design_matrix);

    /**
     * @brief Robust method of moments dispersion
     * @param normed_counts Normalized counts
     * @param design_matrix Design matrix
     * @return Robust dispersion estimates
     */
    Eigen::VectorXd robustMethodOfMomentsDisp(
        const Eigen::MatrixXd &normed_counts,
        const Eigen::MatrixXd &design_matrix);

    /**
     * @brief Save results to CSV file
     * @param filename Output filename
     * @param results Results matrix
     * @param gene_names Gene names
     */
    void saveResultsToCSV(const std::string &filename,
                          const Eigen::MatrixXd &results,
                          const std::vector<std::string> &gene_names);

    /**
     * @brief Load CSV file
     * @param filename Input filename
     * @return Data matrix
     */
    Eigen::MatrixXd loadCSV(const std::string &filename);

    /**
     * @brief Convert natural log to log2
     * @param x Natural log values
     * @return Log2 values
     */
    Eigen::VectorXd log2(const Eigen::VectorXd &x);

    /**
     * @brief Calculate mean absolute deviation
     * @param x Input vector
     * @return Mean absolute deviation
     */
    double meanAbsoluteDeviation(const Eigen::VectorXd &x);

    /**
     * @brief Load metadata and extract 'condition' column as numeric vector (0 for A, 1 for B)
     * @param filename Input metadata CSV file
     * @return Column vector of encoded condition values
     */
    Eigen::MatrixXd loadMetadataCondition(const std::string &filename);

    /**
     * @brief DESeq2 normalization fit - compute logmeans and filtered_genes
     * @param counts Raw count matrix
     * @return Tuple of (logmeans, filtered_genes)
     */
    std::tuple<Eigen::VectorXd, Eigen::VectorXd> deseq2NormFit(const Eigen::MatrixXd &counts);

    /**
     * @brief DESeq2 normalization transform - compute normalized counts and size factors
     * @param counts Raw count matrix
     * @param logmeans Gene-wise mean log counts
     * @param filtered_genes Genes whose log means are different from -∞
     * @return Tuple of (normalized_counts, size_factors)
     */
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> deseq2NormTransform(
        const Eigen::MatrixXd &counts,
        const Eigen::VectorXd &logmeans,
        const Eigen::VectorXd &filtered_genes);

    /**
     * @brief Convert PPM gene count files to DESeq2-compatible format
     * @param file_paths Vector of file paths to PPM gene count files
     * @param sample_names Vector of sample names (optional, will use filenames if empty)
     * @param total_reads Vector of total reads per sample (optional, will be extracted from files if empty)
     * @return Tuple of (count_matrix, gene_names, sample_names)
     *
     * The input files should have the format:
     * File:,filename.sam
     * , TotalReads , 4599960
     * , TotalHits (count), 4387489
     * Chromosome , GeneName , PPM , NCBI_Acc1, NCBI_Acc2, ...
     * 19_KI270882v1_alt , KIR3DL1 , 0.0, NM_013289, NM_013289, ...
     *
     * The output will be:
     * - count_matrix: Matrix with genes as rows and samples as columns
     * - gene_names: Vector of gene names
     * - sample_names: Vector of sample names
     */
    std::tuple<Eigen::MatrixXd, std::vector<std::string>, std::vector<std::string>>
    convertPpmToDeseq2Format(
        const std::vector<std::string> &file_paths,
        const std::vector<std::string> &sample_names = {},
        const std::vector<int> &total_reads = {});

    /**
     * @brief Extract sample information from PPM file header
     * @param file_path Path to PPM gene count file
     * @return Tuple of (sample_name, total_reads, total_hits)
     */
    std::tuple<std::string, int, int> extractSampleInfo(const std::string &file_path);

    /**
     * @brief Parse PPM gene count data from a single file
     * @param file_path Path to PPM gene count file
     * @param total_reads Total reads for this sample (optional, will be extracted if 0)
     * @return Tuple of (gene_names, ppm_values, raw_counts)
     */
    std::tuple<std::vector<std::string>, std::vector<double>, std::vector<int>>
    parsePpmFile(const std::string &file_path, int total_reads = 0);

} // namespace deseq2
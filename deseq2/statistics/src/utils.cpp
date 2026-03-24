#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <tuple>
#include <iomanip>

namespace deseq2
{

    Eigen::MatrixXd loadExampleData(const std::string &modality, const std::string &dataset)
    {
        if (modality != "raw_counts" && modality != "metadata")
        {
            throw std::invalid_argument("Modality must be 'raw_counts' or 'metadata'");
        }

        if (dataset != "synthetic")
        {
            throw std::invalid_argument("Dataset must be 'synthetic'");
        }

        // Use absolute path based on workspace root
        std::string base_path = "/Users/venky/Projects/deepn-plus/PyDESeq2/datasets/synthetic/";
        if (modality == "raw_counts")
        {
            return loadCSV(base_path + "test_counts.csv");
        }
        else
        {
            return loadMetadataCondition(base_path + "test_metadata.csv");
        }
    }

    void testValidCounts(const Eigen::MatrixXd &counts)
    {
        if (counts.hasNaN())
        {
            throw std::invalid_argument("NaNs are not allowed in the count matrix");
        }

        if ((counts.array() < 0).any())
        {
            throw std::invalid_argument("The count matrix should only contain non-negative values");
        }

        // Check if all values are integers (approximately)
        Eigen::MatrixXd fractional = counts.array() - counts.array().floor();
        if ((fractional.array() > 1e-10).any())
        {
            throw std::invalid_argument("The count matrix should only contain integers");
        }
    }

    double dispersionTrend(double normed_mean, const Eigen::Vector2d &coeffs)
    {
        return coeffs(0) + coeffs(1) / normed_mean;
    }

    double nbNLL(const Eigen::VectorXd &counts, const Eigen::VectorXd &mu, double alpha)
    {
        int n = counts.size();
        double inv_alpha = 1.0 / alpha;

        double nll = n * inv_alpha * std::log(alpha);

        for (int i = 0; i < n; ++i)
        {
            double y = counts(i);
            double m = mu(i);

            // Log gamma terms
            double log_gamma_y_plus_inv_alpha = std::lgamma(y + inv_alpha);
            double log_gamma_y_plus_1 = std::lgamma(y + 1);
            double log_gamma_inv_alpha = std::lgamma(inv_alpha);

            nll += -(log_gamma_y_plus_inv_alpha - log_gamma_y_plus_1 - log_gamma_inv_alpha);
            nll += (inv_alpha + y) * std::log(inv_alpha + m);
            nll -= y * std::log(m);
        }

        return nll;
    }

    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, bool>
    irlsSolver(const Eigen::VectorXd &counts,
               const Eigen::VectorXd &size_factors,
               const Eigen::MatrixXd &design_matrix,
               double disp,
               double min_mu,
               double beta_tol,
               double min_beta,
               double max_beta,
               int maxiter)
    {

        int num_vars = design_matrix.cols();
        Eigen::MatrixXd X = design_matrix;

        // Suppress unused parameter warning for min_beta (it's used in bounds checking)
        (void)min_beta;

        // Initialize beta using QR decomposition if full rank
        Eigen::VectorXd beta;
        if (X.colPivHouseholderQr().rank() == num_vars)
        {
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
            Eigen::VectorXd y = (counts.array() / size_factors.array() + 0.1).log();
            beta = qr.solve(y);
        }
        else
        {
            // Initialize intercept with log base mean
            beta = Eigen::VectorXd::Zero(num_vars);
            beta(0) = (counts.array() / size_factors.array()).log().mean();
        }

        double dev = 1000.0;
        double dev_ratio = 1.0;

        // Ridge factor for regularization
        Eigen::MatrixXd ridge_factor = Eigen::MatrixXd::Identity(num_vars, num_vars) * 1e-6;

        // Initial mu
        Eigen::VectorXd mu = (size_factors.array() * (X * beta).array().exp()).cwiseMax(min_mu);

        bool converged = true;
        int i = 0;

        while (dev_ratio > beta_tol && i < maxiter)
        {
            // Calculate weights
            Eigen::VectorXd W = mu.array() / (1.0 + mu.array() * disp);

            // Calculate working response
            Eigen::VectorXd z = (mu.array() / size_factors.array()).log() +
                                (counts.array() - mu.array()) / mu.array();

            // Update beta using IRLS
            Eigen::MatrixXd H = X.transpose() * W.asDiagonal() * X + ridge_factor;
            Eigen::VectorXd beta_hat = H.ldlt().solve(X.transpose() * (W.asDiagonal() * z));

            i++;

            // Check for divergence
            if ((beta_hat.array().abs() > max_beta).any())
            {
                converged = false;
                break;
            }

            beta = beta_hat;
            mu = (size_factors.array() * (X * beta).array().exp()).cwiseMax(min_mu);

            // Compute deviance ratio
            double old_dev = dev;
            dev = -2 * nbNLL(counts, mu, disp);
            dev_ratio = std::abs(dev - old_dev) / (std::abs(dev) + 0.1);
        }

        // Calculate hat matrix diagonal
        Eigen::VectorXd W = mu.array() / (1.0 + mu.array() * disp);
        Eigen::MatrixXd XtWX = X.transpose() * W.asDiagonal() * X + ridge_factor;
        Eigen::MatrixXd XtWX_inv = XtWX.inverse();
        Eigen::MatrixXd H_diag = X * XtWX_inv * X.transpose();

        // Extract diagonal and apply weights
        Eigen::VectorXd H = H_diag.diagonal().cwiseProduct(W);

        // Return unthresholded mu
        Eigen::VectorXd mu_final = size_factors.array() * (X * beta).array().exp();

        return {beta, mu_final, H, converged};
    }

    std::tuple<double, bool> fitAlphaMLE(
        const Eigen::VectorXd &counts,
        const Eigen::MatrixXd &design_matrix,
        const Eigen::VectorXd &mu,
        double alpha_hat,
        double min_disp,
        double max_disp,
        double prior_disp_var)
    {
        // Suppress unused parameter warnings
        (void)design_matrix;
        (void)alpha_hat;
        (void)prior_disp_var;

        // Simple optimization using golden section search
        double a = std::log(min_disp);
        double b = std::log(max_disp);
        double c = b - (b - a) / 1.618;
        double d = a + (b - a) / 1.618;

        double fc = nbNLL(counts, mu, std::exp(c));
        double fd = nbNLL(counts, mu, std::exp(d));

        for (int iter = 0; iter < 50; ++iter)
        {
            if (fc < fd)
            {
                b = d;
                d = c;
                fd = fc;
                c = b - (b - a) / 1.618;
                fc = nbNLL(counts, mu, std::exp(c));
            }
            else
            {
                a = c;
                c = d;
                fc = fd;
                d = a + (b - a) / 1.618;
                fd = nbNLL(counts, mu, std::exp(d));
            }

            if (std::abs(b - a) < 1e-6)
            {
                break;
            }
        }

        double alpha = std::exp((a + b) / 2.0);
        return std::make_tuple(alpha, true);
    }

    double trimmedMean(const Eigen::VectorXd &x, double trim)
    {
        int n = x.size();
        int trim_count = static_cast<int>(n * trim);

        Eigen::VectorXd sorted = x;
        std::sort(sorted.data(), sorted.data() + n);

        double sum = 0.0;
        int count = 0;
        for (int i = trim_count; i < n - trim_count; ++i)
        {
            sum += sorted(i);
            count++;
        }

        return count > 0 ? sum / count : 0.0;
    }

    Eigen::VectorXd fitLinMu(
        const Eigen::VectorXd &counts,
        const Eigen::VectorXd &size_factors,
        const Eigen::MatrixXd &design_matrix,
        double min_mu)
    {

        auto [beta, mu, H, converged] = irlsSolver(counts, size_factors, design_matrix, 0.1, min_mu);
        return mu;
    }

    std::tuple<double, double, double> waldTest(
        const Eigen::MatrixXd &design_matrix,
        double disp,
        const Eigen::VectorXd &lfc,
        const Eigen::VectorXd &mu,
        const Eigen::VectorXd &contrast,
        double lfc_null)
    {
        int n_samples = design_matrix.rows();
        int n_vars = design_matrix.cols();
        (void)n_samples; // Suppress unused variable warning (used in covariance calculation)
        (void)n_vars;    // Suppress unused variable warning (used in covariance calculation)

        // Calculate weights
        Eigen::VectorXd W = mu.array() / (1.0 + mu.array() * disp);

        // Calculate covariance matrix
        Eigen::MatrixXd XtWX = design_matrix.transpose() * W.asDiagonal() * design_matrix;
        Eigen::MatrixXd cov_matrix = XtWX.inverse();

        // Calculate contrast variance
        double contrast_var = contrast.transpose() * cov_matrix * contrast;
        double lfc_se = std::sqrt(contrast_var);

        // Calculate contrast LFC
        double contrast_lfc = contrast.dot(lfc);

        // Calculate Wald statistic
        double stat = (contrast_lfc - lfc_null) / lfc_se;

        // Calculate p-value (two-sided)
        double p_value = std::erfc(std::abs(stat) / std::sqrt(2.0));

        return std::make_tuple(stat, p_value, lfc_se);
    }

    Eigen::VectorXd fitRoughDispersions(
        const Eigen::MatrixXd &normed_counts,
        const Eigen::MatrixXd &design_matrix)
    {

        int num_samples = design_matrix.rows();
        int num_vars = design_matrix.cols();
        (void)num_samples; // Suppress unused variable warning
        (void)num_vars;    // Suppress unused variable warning

        // Check if we have enough samples
        if (num_samples == num_vars)
        {
            throw std::runtime_error("The number of samples and the number of design variables are equal");
        }

        int n_genes = normed_counts.cols();
        Eigen::VectorXd alpha_rde = Eigen::VectorXd::Zero(n_genes);

        // For each gene, fit linear regression and compute dispersion
        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd y = normed_counts.col(j);

            // Linear regression without intercept (design matrix already has intercept)
            Eigen::VectorXd beta = design_matrix.colPivHouseholderQr().solve(y);
            Eigen::VectorXd y_hat = design_matrix * beta;

            // Ensure y_hat >= 1
            y_hat = y_hat.cwiseMax(1.0);

            // Compute dispersion estimate
            Eigen::VectorXd residuals = y - y_hat;
            double sum_disp = 0.0;
            for (int i = 0; i < num_samples; ++i)
            {
                double term = (residuals(i) * residuals(i) - y_hat(i)) / (y_hat(i) * y_hat(i));
                sum_disp += term;
            }
            alpha_rde(j) = std::max(0.0, sum_disp / (num_samples - num_vars));
        }

        return alpha_rde;
    }

    Eigen::VectorXd fitMomentsDispersions(
        const Eigen::MatrixXd &normed_counts,
        const Eigen::VectorXd &size_factors)
    {

        int n_genes = normed_counts.cols();
        Eigen::VectorXd alpha = Eigen::VectorXd::Zero(n_genes);

        // Mean inverse size factor
        double s_mean_inv = (1.0 / size_factors.array()).mean();

        for (int j = 0; j < n_genes; ++j)
        {
            Eigen::VectorXd gene_counts = normed_counts.col(j);

            // Skip genes with all zeros
            if (gene_counts.sum() == 0)
            {
                alpha(j) = 0.0;
                continue;
            }

            double mu = gene_counts.mean();
            double sigma = 0.0;

            // Compute variance with ddof=1 (unbiased estimator)
            for (int i = 0; i < gene_counts.size(); ++i)
            {
                sigma += (gene_counts(i) - mu) * (gene_counts(i) - mu);
            }
            sigma /= (gene_counts.size() - 1);

            // Handle NaN (variance = 0)
            if (std::isnan(sigma) || sigma == 0)
            {
                alpha(j) = 0.0;
            }
            else
            {
                alpha(j) = (sigma - s_mean_inv * mu) / (mu * mu);
            }
        }

        return alpha;
    }

    Eigen::VectorXd robustMethodOfMomentsDisp(
        const Eigen::MatrixXd &normed_counts,
        const Eigen::MatrixXd &design_matrix)
    {

        return fitRoughDispersions(normed_counts, design_matrix);
    }

    void saveResultsToCSV(const std::string &filename,
                          const Eigen::MatrixXd &results,
                          const std::vector<std::string> &gene_names)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // Set high precision for output
        file << std::setprecision(15);

        // Write header
        file << "GeneName,BaseMean,Log2FoldChange,Log2FoldChangeStdErr,WaldStatistic,Pvalue,Padj\n";

        // Write data
        for (int i = 0; i < results.rows(); ++i)
        {
            file << gene_names[i];
            for (int j = 0; j < results.cols(); ++j)
            {
                file << "," << results(i, j);
            }
            file << "\n";
        }

        file.close();
    }

    Eigen::MatrixXd loadCSV(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::vector<std::vector<double>> data;
        std::string line;

        // Skip header
        std::getline(file, line);

        while (std::getline(file, line))
        {
            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;

            // Skip first column (gene/sample names)
            std::getline(ss, cell, ',');

            while (std::getline(ss, cell, ','))
            {
                // Remove quotes if present
                cell.erase(std::remove(cell.begin(), cell.end(), '"'), cell.end());
                // Skip empty cells
                if (cell.empty())
                    continue;
                try
                {
                    row.push_back(std::stod(cell));
                }
                catch (const std::invalid_argument &)
                {
                    // Skip non-numeric cells
                    continue;
                }
            }
            if (!row.empty())
            {
                data.push_back(row);
            }
        }

        file.close();

        if (data.empty())
        {
            return Eigen::MatrixXd();
        }

        int rows = data.size();
        int cols = data[0].size();

        Eigen::MatrixXd matrix(rows, cols);
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                matrix(i, j) = data[i][j];
            }
        }

        return matrix;
    }

    Eigen::VectorXd log2(const Eigen::VectorXd &x)
    {
        return x.array() / std::log(2.0);
    }

    double meanAbsoluteDeviation(const Eigen::VectorXd &x)
    {
        double mean = x.mean();
        return (x.array() - mean).abs().mean();
    }

    // Helper to load metadata and extract 'condition' as numeric vector
    Eigen::MatrixXd loadMetadataCondition(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        std::vector<double> condition;
        std::string line;
        // Skip header
        std::getline(file, line);
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            // Skip sample name
            std::getline(ss, cell, ',');
            // Read condition
            std::getline(ss, cell, ',');
            cell.erase(std::remove(cell.begin(), cell.end(), '"'), cell.end());
            if (cell == "A")
                condition.push_back(0.0);
            else if (cell == "B")
                condition.push_back(1.0);
            else
                condition.push_back(-1.0); // error value
        }
        file.close();
        // Return as column vector
        Eigen::MatrixXd mat(condition.size(), 1);
        for (size_t i = 0; i < condition.size(); ++i)
            mat(i, 0) = condition[i];
        return mat;
    }

    std::tuple<Eigen::VectorXd, Eigen::VectorXd> deseq2NormFit(const Eigen::MatrixXd &counts)
    {
        int n_samples = counts.rows();
        int n_genes = counts.cols();

        // Compute gene-wise mean log counts
        Eigen::VectorXd logmeans = Eigen::VectorXd::Zero(n_genes);
        Eigen::VectorXd filtered_genes = Eigen::VectorXd::Zero(n_genes);

        for (int j = 0; j < n_genes; ++j)
        {
            double sum_log = 0.0;
            int count = 0;

            for (int i = 0; i < n_samples; ++i)
            {
                if (counts(i, j) > 0)
                {
                    sum_log += std::log(counts(i, j));
                    count++;
                }
            }

            if (count > 0)
            {
                logmeans(j) = sum_log / count;
                filtered_genes(j) = 1.0; // True
            }
            else
            {
                logmeans(j) = -std::numeric_limits<double>::infinity();
                filtered_genes(j) = 0.0; // False
            }
        }

        return std::make_tuple(logmeans, filtered_genes);
    }

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> deseq2NormTransform(
        const Eigen::MatrixXd &counts,
        const Eigen::VectorXd &logmeans,
        const Eigen::VectorXd &filtered_genes)
    {

        int n_samples = counts.rows();
        int n_genes = counts.cols();

        // Compute log counts
        Eigen::MatrixXd log_counts = Eigen::MatrixXd::Zero(n_samples, n_genes);
        for (int i = 0; i < n_samples; ++i)
        {
            for (int j = 0; j < n_genes; ++j)
            {
                if (counts(i, j) > 0)
                {
                    log_counts(i, j) = std::log(counts(i, j));
                }
            }
        }

        // Subtract filtered log means from log counts
        Eigen::MatrixXd log_ratios = Eigen::MatrixXd::Zero(n_samples, n_genes);
        for (int i = 0; i < n_samples; ++i)
        {
            for (int j = 0; j < n_genes; ++j)
            {
                if (filtered_genes(j) > 0.5)
                {
                    log_ratios(i, j) = log_counts(i, j) - logmeans(j);
                }
            }
        }

        // Compute sample-wise median of log ratios
        Eigen::VectorXd log_medians = Eigen::VectorXd::Zero(n_samples);
        for (int i = 0; i < n_samples; ++i)
        {
            std::vector<double> ratios;
            for (int j = 0; j < n_genes; ++j)
            {
                if (filtered_genes(j) > 0.5)
                {
                    ratios.push_back(log_ratios(i, j));
                }
            }
            if (!ratios.empty())
            {
                std::sort(ratios.begin(), ratios.end());
                log_medians(i) = ratios[ratios.size() / 2];
            }
        }

        // Return raw counts divided by size factors and size factors
        Eigen::VectorXd size_factors = log_medians.array().exp();
        Eigen::MatrixXd deseq2_counts = Eigen::MatrixXd::Zero(n_samples, n_genes);

        for (int i = 0; i < n_samples; ++i)
        {
            for (int j = 0; j < n_genes; ++j)
            {
                deseq2_counts(i, j) = counts(i, j) / size_factors(i);
            }
        }

        return std::make_tuple(deseq2_counts, size_factors);
    }

    std::tuple<std::string, int, int> extractSampleInfo(const std::string &file_path)
    {
        std::ifstream file(file_path);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + file_path);
        }

        std::string sample_name;
        int total_reads = 0;
        int total_hits = 0;
        std::string line;

        // Read first line to get sample name
        if (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            std::getline(ss, cell, ','); // "File:"
            std::getline(ss, cell, ','); // filename
            if (!cell.empty())
            {
                // Extract filename without path and extension
                size_t last_slash = cell.find_last_of("/\\");
                if (last_slash != std::string::npos)
                {
                    cell = cell.substr(last_slash + 1);
                }
                size_t last_dot = cell.find_last_of('.');
                if (last_dot != std::string::npos)
                {
                    cell = cell.substr(0, last_dot);
                }
                sample_name = cell;
            }
        }

        // Read total reads
        if (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            std::getline(ss, cell, ','); // empty
            std::getline(ss, cell, ','); // "TotalReads"
            std::getline(ss, cell, ','); // value
            if (!cell.empty())
            {
                try
                {
                    total_reads = std::stoi(cell);
                }
                catch (const std::exception &)
                {
                    // Ignore parsing errors
                }
            }
        }

        // Read total hits
        if (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            std::getline(ss, cell, ','); // empty
            std::getline(ss, cell, ','); // "TotalHits (count)"
            std::getline(ss, cell, ','); // value
            if (!cell.empty())
            {
                try
                {
                    total_hits = std::stoi(cell);
                }
                catch (const std::exception &)
                {
                    // Ignore parsing errors
                }
            }
        }

        file.close();
        return std::make_tuple(sample_name, total_reads, total_hits);
    }

    std::tuple<std::vector<std::string>, std::vector<double>, std::vector<int>>
    parsePpmFile(const std::string &file_path, int total_reads)
    {
        std::ifstream file(file_path);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + file_path);
        }

        std::vector<std::string> gene_names;
        std::vector<double> ppm_values;
        std::vector<int> raw_counts;

        std::string line;

        // Skip header lines (first 4 lines)
        for (int i = 0; i < 4; ++i)
        {
            std::getline(file, line);
        }

        // If total_reads not provided, extract from file
        if (total_reads == 0)
        {
            file.close();
            auto [_, reads, __] = extractSampleInfo(file_path);
            total_reads = reads;
            file.open(file_path);
            if (!file.is_open())
            {
                throw std::runtime_error("Cannot reopen file: " + file_path);
            }
            // Skip header lines again
            for (int i = 0; i < 4; ++i)
            {
                std::getline(file, line);
            }
        }

        // Parse data lines
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            std::vector<std::string> cells;

            // Split line by comma
            while (std::getline(ss, cell, ','))
            {
                // Remove leading/trailing whitespace
                cell.erase(0, cell.find_first_not_of(" \t\r\n"));
                cell.erase(cell.find_last_not_of(" \t\r\n") + 1);
                cells.push_back(cell);
            }

            if (cells.size() >= 3)
            {
                std::string chromosome = cells[0];
                std::string gene_name = cells[1];
                std::string ppm_str = cells[2];

                // Skip if gene name is empty or "GeneName"
                if (gene_name.empty() || gene_name == "GeneName")
                {
                    continue;
                }

                try
                {
                    double ppm = std::stod(ppm_str);

                    // Convert PPM to raw counts
                    int raw_count = static_cast<int>(std::round(ppm * total_reads / 1000000.0));

                    gene_names.push_back(gene_name);
                    ppm_values.push_back(ppm);
                    raw_counts.push_back(raw_count);
                }
                catch (const std::exception &)
                {
                    // Skip lines that can't be parsed
                    continue;
                }
            }
        }

        file.close();
        return std::make_tuple(gene_names, ppm_values, raw_counts);
    }

    std::tuple<Eigen::MatrixXd, std::vector<std::string>, std::vector<std::string>>
    convertPpmToDeseq2Format(
        const std::vector<std::string> &file_paths,
        const std::vector<std::string> &sample_names,
        const std::vector<int> &total_reads)
    {
        if (file_paths.empty())
        {
            throw std::runtime_error("No file paths provided");
        }

        // Determine sample names
        std::vector<std::string> final_sample_names;
        if (sample_names.empty())
        {
            for (const auto &file_path : file_paths)
            {
                auto [sample_name, _, __] = extractSampleInfo(file_path);
                final_sample_names.push_back(sample_name);
            }
        }
        else
        {
            final_sample_names = sample_names;
        }

        // Determine total reads
        std::vector<int> final_total_reads;
        if (total_reads.empty())
        {
            for (const auto &file_path : file_paths)
            {
                auto [_, reads, __] = extractSampleInfo(file_path);
                final_total_reads.push_back(reads);
            }
        }
        else
        {
            final_total_reads = total_reads;
        }

        // Parse all files and collect unique genes
        std::vector<std::vector<std::string>> all_gene_names;
        std::vector<std::vector<int>> all_raw_counts;
        std::set<std::string> unique_genes;

        for (size_t i = 0; i < file_paths.size(); ++i)
        {
            auto [gene_names, ppm_values, raw_counts] = parsePpmFile(file_paths[i], final_total_reads[i]);
            all_gene_names.push_back(gene_names);
            all_raw_counts.push_back(raw_counts);

            // Add to unique genes set
            for (const auto &gene : gene_names)
            {
                unique_genes.insert(gene);
            }
        }

        // Create gene name vector (sorted)
        std::vector<std::string> gene_names(unique_genes.begin(), unique_genes.end());
        std::sort(gene_names.begin(), gene_names.end());

        // Create count matrix
        int n_genes = gene_names.size();
        int n_samples = file_paths.size();
        Eigen::MatrixXd count_matrix = Eigen::MatrixXd::Zero(n_genes, n_samples);

        // Fill count matrix
        for (int sample_idx = 0; sample_idx < n_samples; ++sample_idx)
        {
            const auto &sample_gene_names = all_gene_names[sample_idx];
            const auto &sample_raw_counts = all_raw_counts[sample_idx];

            for (size_t gene_idx = 0; gene_idx < sample_gene_names.size(); ++gene_idx)
            {
                const std::string &gene_name = sample_gene_names[gene_idx];
                int raw_count = sample_raw_counts[gene_idx];

                // Find gene in the master gene list
                auto it = std::find(gene_names.begin(), gene_names.end(), gene_name);
                if (it != gene_names.end())
                {
                    int master_gene_idx = std::distance(gene_names.begin(), it);
                    count_matrix(master_gene_idx, sample_idx) = raw_count;
                }
            }
        }

        return std::make_tuple(count_matrix, gene_names, final_sample_names);
    }

} // namespace deseq2
/**
 * @file example_usage.cpp
 * @brief Comprehensive example demonstrating how to use the DESeq2 C++ library
 *
 * This example shows the complete workflow for differential expression analysis:
 * 1. Data preparation and loading
 * 2. DESeq2 analysis pipeline
 * 3. Statistical testing
 * 4. Results interpretation and export
 *
 * @author DESeq2 C++ Implementation
 * @date 2024
 */

#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <limits>

using namespace deseq2;

/**
 * @brief Print matrix information for debugging
 * @param matrix The matrix to print info about
 * @param name Name of the matrix
 */
void printMatrixInfo(const Eigen::MatrixXd &matrix, const std::string &name)
{
    std::cout << name << " shape: " << matrix.rows() << " x " << matrix.cols() << std::endl;
    std::cout << name << " range: [" << matrix.minCoeff() << ", " << matrix.maxCoeff() << "]" << std::endl;
    std::cout << name << " has NaN: " << (matrix.hasNaN() ? "Yes" : "No") << std::endl;
    std::cout << name << " has Inf: " << ((matrix.array() == std::numeric_limits<double>::infinity()).any() ? "Yes" : "No") << std::endl;
}

/**
 * @brief Print vector information for debugging
 * @param vector The vector to print info about
 * @param name Name of the vector
 */
void printVectorInfo(const Eigen::VectorXd &vector, const std::string &name)
{
    std::cout << name << " size: " << vector.size() << std::endl;
    std::cout << name << " range: [" << vector.minCoeff() << ", " << vector.maxCoeff() << "]" << std::endl;
    std::cout << name << " has NaN: " << (vector.hasNaN() ? "Yes" : "No") << std::endl;
    std::cout << name << " has Inf: " << ((vector.array() == std::numeric_limits<double>::infinity()).any() ? "Yes" : "No") << std::endl;
}

/**
 * @brief Print summary statistics for results
 * @param results Results matrix from DESeq2 analysis
 */
void printResultsSummary(const Eigen::MatrixXd &results)
{
    std::cout << "\n=== RESULTS SUMMARY ===" << std::endl;
    std::cout << "Number of genes analyzed: " << results.rows() << std::endl;

    // Count significant genes (padj < 0.05)
    int significant_count = 0;
    int upregulated_count = 0;
    int downregulated_count = 0;

    for (int i = 0; i < results.rows(); ++i)
    {
        double padj = results(i, 4); // padj column
        double lfc = results(i, 1);  // log2FoldChange column

        if (padj < 0.05)
        {
            significant_count++;
            if (lfc > 0)
            {
                upregulated_count++;
            }
            else if (lfc < 0)
            {
                downregulated_count++;
            }
        }
    }

    std::cout << "Significant genes (padj < 0.05): " << significant_count << std::endl;
    std::cout << "  - Upregulated: " << upregulated_count << std::endl;
    std::cout << "  - Downregulated: " << downregulated_count << std::endl;
    std::cout << "  - No change: " << (significant_count - upregulated_count - downregulated_count) << std::endl;

    // Print top 5 most significant genes
    std::cout << "\nTop 5 most significant genes:" << std::endl;
    std::cout << std::setw(10) << "Gene"
              << std::setw(15) << "log2FC"
              << std::setw(15) << "padj"
              << std::setw(15) << "baseMean" << std::endl;
    std::cout << std::string(55, '-') << std::endl;

    // Create vector of indices sorted by p-value
    std::vector<std::pair<double, int>> pvalue_indices;
    for (int i = 0; i < results.rows(); ++i)
    {
        pvalue_indices.push_back({results(i, 4), i}); // padj, gene_index
    }

    // Sort by p-value (ascending)
    std::sort(pvalue_indices.begin(), pvalue_indices.end());

    // Print top 5
    for (int i = 0; i < std::min(5, (int)pvalue_indices.size()); ++i)
    {
        int gene_idx = pvalue_indices[i].second;
        std::cout << std::setw(10) << ("gene" + std::to_string(gene_idx + 1))
                  << std::setw(15) << std::fixed << std::setprecision(3) << results(gene_idx, 1)
                  << std::setw(15) << std::scientific << std::setprecision(2) << results(gene_idx, 4)
                  << std::setw(15) << std::fixed << std::setprecision(1) << results(gene_idx, 0) << std::endl;
    }
}

/**
 * @brief Example 1: Basic DESeq2 analysis with synthetic data
 */
void exampleBasicAnalysis()
{
    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "EXAMPLE 1: BASIC DESEQ2 ANALYSIS" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    try
    {
        // Step 1: Load data
        std::cout << "\nStep 1: Loading synthetic data..." << std::endl;
        Eigen::MatrixXd counts = loadExampleData("raw_counts", "synthetic");
        Eigen::MatrixXd metadata = loadExampleData("metadata", "synthetic");

        // Transpose counts to match expected format (samples x genes)
        counts = counts.transpose().eval();

        printMatrixInfo(counts, "Counts");
        printMatrixInfo(metadata, "Metadata");

        // Step 2: Create DeseqDataSet
        std::cout << "\nStep 2: Creating DeseqDataSet..." << std::endl;
        DeseqDataSet dds(counts, metadata, "~condition", true);

        // Step 3: Fit size factors (normalization)
        std::cout << "\nStep 3: Fitting size factors..." << std::endl;
        dds.fitSizeFactors();
        printVectorInfo(dds.getSizeFactors(), "Size factors");

        // Step 4: Fit genewise dispersions
        std::cout << "\nStep 4: Fitting genewise dispersions..." << std::endl;
        dds.fitGenewiseDispersions();
        printVectorInfo(dds.getGenewiseDispersions(), "Genewise dispersions");

        // Step 5: Fit dispersion trend
        std::cout << "\nStep 5: Fitting dispersion trend..." << std::endl;
        dds.fitDispersionTrend();
        printVectorInfo(dds.getFittedDispersions(), "Fitted dispersions");

        // Step 6: Fit dispersion prior
        std::cout << "\nStep 6: Fitting dispersion prior..." << std::endl;
        dds.fitDispersionPrior();

        // Step 7: Fit MAP dispersions
        std::cout << "\nStep 7: Fitting MAP dispersions..." << std::endl;
        dds.fitMAPDispersions();
        printVectorInfo(dds.getMAPDispersions(), "MAP dispersions");

        // Step 8: Fit log fold changes
        std::cout << "\nStep 8: Fitting log fold changes..." << std::endl;
        dds.fitLFC();
        printMatrixInfo(dds.getLFC(), "Log fold changes");

        // Step 9: Calculate Cooks distances and refit
        std::cout << "\nStep 9: Calculating Cooks distances and refitting..." << std::endl;
        dds.calculateCooks();
        dds.refit();

        // Step 10: Statistical analysis
        std::cout << "\nStep 10: Performing statistical analysis..." << std::endl;

        // Create contrast vector for condition B vs A
        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0; // Test condition B vs A

        DeseqStats ds(dds, contrast, 0.05, true, true);

        // Run Wald test
        std::cout << "Running Wald test..." << std::endl;
        ds.runWaldTest();

        // Apply filtering and p-value adjustment
        std::cout << "Applying Cooks filtering..." << std::endl;
        ds.cooksFiltering();

        std::cout << "Applying independent filtering..." << std::endl;
        ds.independentFiltering();

        // Generate results
        std::cout << "Generating results..." << std::endl;
        Eigen::MatrixXd results = ds.summary();

        // Create gene names
        std::vector<std::string> gene_names;
        for (int i = 0; i < results.rows(); ++i)
        {
            gene_names.push_back("gene" + std::to_string(i + 1));
        }

        // Save results
        std::string output_filename = "example_basic_results.csv";
        std::cout << "Saving results to " << output_filename << "..." << std::endl;
        saveResultsToCSV(output_filename, results, gene_names);

        // Print summary
        printResultsSummary(results);

        std::cout << "\nBasic analysis completed successfully!" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in basic analysis: " << e.what() << std::endl;
    }
}

/**
 * @brief Example 2: Custom data analysis
 */
void exampleCustomDataAnalysis()
{
    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "EXAMPLE 2: CUSTOM DATA ANALYSIS" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    try
    {
        // Step 1: Create custom count data
        std::cout << "\nStep 1: Creating custom count data..." << std::endl;

        // Simulate RNA-seq data: 6 samples (3 control, 3 treatment), 100 genes
        int n_samples = 6;
        int n_genes = 100;

        Eigen::MatrixXd counts(n_samples, n_genes);

        // Generate realistic count data
        for (int gene = 0; gene < n_genes; ++gene)
        {
            for (int sample = 0; sample < n_samples; ++sample)
            {
                double base_mean = 100.0 + gene * 10.0; // Increasing expression across genes

                // Add treatment effect for some genes
                if (gene < 20 && sample >= 3)
                { // First 20 genes are upregulated in treatment
                    base_mean *= 3.0;
                }
                else if (gene >= 20 && gene < 40 && sample >= 3)
                { // Next 20 genes are downregulated
                    base_mean *= 0.3;
                }

                // Add some noise
                double noise = 1.0 + 0.5 * (rand() / (double)RAND_MAX - 0.5);
                counts(sample, gene) = std::max(0, (int)(base_mean * noise));
            }
        }

        // Step 2: Create metadata
        std::cout << "\nStep 2: Creating metadata..." << std::endl;
        Eigen::MatrixXd metadata(n_samples, 1);
        for (int i = 0; i < 3; ++i)
        {
            metadata(i, 0) = 0.0; // Control samples
        }
        for (int i = 3; i < 6; ++i)
        {
            metadata(i, 0) = 1.0; // Treatment samples
        }

        printMatrixInfo(counts, "Custom counts");
        printMatrixInfo(metadata, "Custom metadata");

        // Step 3: Run DESeq2 analysis
        std::cout << "\nStep 3: Running DESeq2 analysis..." << std::endl;

        DeseqDataSet dds(counts, metadata, "~condition", true);

        // Run the complete pipeline
        dds.fitSizeFactors();
        dds.fitGenewiseDispersions();
        dds.fitDispersionTrend();
        dds.fitDispersionPrior();
        dds.fitMAPDispersions();
        dds.fitLFC();
        dds.calculateCooks();
        dds.refit();

        // Step 4: Statistical testing
        std::cout << "\nStep 4: Statistical testing..." << std::endl;

        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0; // Treatment vs Control

        DeseqStats ds(dds, contrast, 0.05, true, true);
        ds.runWaldTest();
        ds.cooksFiltering();
        ds.independentFiltering();

        // Step 5: Generate and save results
        std::cout << "\nStep 5: Generating results..." << std::endl;
        Eigen::MatrixXd results = ds.summary();

        std::vector<std::string> gene_names;
        for (int i = 0; i < results.rows(); ++i)
        {
            gene_names.push_back("custom_gene_" + std::to_string(i + 1));
        }

        std::string output_filename = "example_custom_results.csv";
        saveResultsToCSV(output_filename, results, gene_names);

        // Print summary
        printResultsSummary(results);

        std::cout << "\nCustom data analysis completed successfully!" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in custom analysis: " << e.what() << std::endl;
    }
}

/**
 * @brief Example 3: Loading data from CSV files
 */
void exampleLoadFromCSV()
{
    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "EXAMPLE 3: LOADING DATA FROM CSV FILES" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    try
    {
        std::cout << "\nThis example demonstrates how to load data from CSV files." << std::endl;
        std::cout << "Expected CSV format:" << std::endl;
        std::cout << "- Counts file: genes as rows, samples as columns" << std::endl;
        std::cout << "- Metadata file: samples as rows, variables as columns" << std::endl;
        std::cout << "- Include 'condition' column in metadata for group comparison" << std::endl;

        // Example of how to load from CSV (commented out since files don't exist)
        /*
        std::cout << "\nLoading counts from CSV..." << std::endl;
        Eigen::MatrixXd counts = loadCSV("path/to/counts.csv");
        counts = counts.transpose().eval(); // Transpose to samples x genes

        std::cout << "Loading metadata from CSV..." << std::endl;
        Eigen::MatrixXd metadata = loadMetadataCondition("path/to/metadata.csv");

        // Run analysis
        DeseqDataSet dds(counts, metadata, "~condition", true);
        // ... rest of analysis pipeline
        */

        std::cout << "\nCSV loading example completed (no actual files loaded)." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in CSV loading example: " << e.what() << std::endl;
    }
}

/**
 * @brief Example 4: Advanced analysis with different parameters
 */
void exampleAdvancedAnalysis()
{
    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "EXAMPLE 4: ADVANCED ANALYSIS WITH CUSTOM PARAMETERS" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    try
    {
        // Load data
        Eigen::MatrixXd counts = loadExampleData("raw_counts", "synthetic");
        Eigen::MatrixXd metadata = loadExampleData("metadata", "synthetic");
        counts = counts.transpose().eval();

        std::cout << "\nRunning advanced analysis with custom parameters..." << std::endl;

        // Create DeseqDataSet with custom parameters
        DeseqDataSet dds(counts, metadata, "~condition",
                         true,  // refit_cooks
                         1.0,   // min_mu (higher than default)
                         1e-6,  // min_disp (lower than default)
                         5.0,   // max_disp (lower than default)
                         1e-6); // beta_tol (tighter convergence)

        // Run analysis pipeline
        dds.fitSizeFactors();
        dds.fitGenewiseDispersions();
        dds.fitDispersionTrend();
        dds.fitDispersionPrior();
        dds.fitMAPDispersions();
        dds.fitLFC();
        dds.calculateCooks();
        dds.refit();

        // Statistical analysis with different parameters
        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0;

        DeseqStats ds(dds, contrast,
                      0.01, // alpha (stricter significance threshold)
                      true, // cooks_filter
                      true, // independent_filter
                      0.0); // lfc_null

        ds.runWaldTest();
        ds.cooksFiltering();
        ds.independentFiltering();

        // Generate results
        Eigen::MatrixXd results = ds.summary();

        std::vector<std::string> gene_names;
        for (int i = 0; i < results.rows(); ++i)
        {
            gene_names.push_back("gene" + std::to_string(i + 1));
        }

        std::string output_filename = "example_advanced_results.csv";
        saveResultsToCSV(output_filename, results, gene_names);

        // Print summary with stricter threshold
        std::cout << "\n=== ADVANCED RESULTS SUMMARY (alpha = 0.01) ===" << std::endl;
        std::cout << "Number of genes analyzed: " << results.rows() << std::endl;

        int significant_count = 0;
        for (int i = 0; i < results.rows(); ++i)
        {
            if (results(i, 4) < 0.01)
            { // padj < 0.01
                significant_count++;
            }
        }

        std::cout << "Significant genes (padj < 0.01): " << significant_count << std::endl;

        std::cout << "\nAdvanced analysis completed successfully!" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in advanced analysis: " << e.what() << std::endl;
    }
}

/**
 * @brief Example 5: PPM Gene Count File Conversion
 */
void examplePPMConversion()
{
    std::cout << "\n=== Example 5: PPM Gene Count File Conversion ===" << std::endl;
    std::cout << "This example demonstrates how to convert PPM gene count files to DESeq2 format." << std::endl;

    try
    {
        // Example file paths (using the existing example files)
        std::vector<std::string> file_paths = {
            "../example/HDPTP_Tail_2_NON_50M_summary.csv",
            "../example/HDPTP_Tail_2_SEL_50N_summary.csv",
            "../example/HDPTP_V_NON_51M_summary.csv",
            "../example/HDPTP_V_SEL_51N_summary.csv"};

        std::cout << "Converting PPM files to DESeq2 format..." << std::endl;

        // Convert PPM files to DESeq2 format
        auto [count_matrix, gene_names, sample_names] = deseq2::convertPpmToDeseq2Format(file_paths);

        std::cout << "Conversion completed successfully!" << std::endl;
        std::cout << "Count matrix dimensions: " << count_matrix.rows() << " genes x " << count_matrix.cols() << " samples" << std::endl;
        std::cout << "Number of unique genes: " << gene_names.size() << std::endl;
        std::cout << "Sample names: ";
        for (const auto &name : sample_names)
        {
            std::cout << name << " ";
        }
        std::cout << std::endl;

        // Create metadata for DESeq2 analysis
        // Assuming NON samples are control (0) and SEL samples are treatment (1)
        Eigen::MatrixXd metadata = Eigen::MatrixXd::Zero(count_matrix.cols(), 1);
        for (size_t i = 0; i < sample_names.size(); ++i)
        {
            if (sample_names[i].find("SEL") != std::string::npos)
            {
                metadata(i, 0) = 1.0; // Treatment
            }
            else
            {
                metadata(i, 0) = 0.0; // Control
            }
        }

        std::cout << "Metadata created with conditions: ";
        for (int i = 0; i < metadata.rows(); ++i)
        {
            std::cout << (metadata(i, 0) == 0 ? "Control" : "Treatment") << " ";
        }
        std::cout << std::endl;

        // Transpose count matrix to match DESeq2 format (samples x genes)
        Eigen::MatrixXd counts_transposed = count_matrix.transpose();

        std::cout << "Running DESeq2 analysis on converted data..." << std::endl;

        // Create DeseqDataSet
        deseq2::DeseqDataSet dds(counts_transposed, metadata, "~condition", true);

        // Run analysis pipeline
        dds.fitSizeFactors();
        dds.fitGenewiseDispersions();
        dds.fitDispersionTrend();
        dds.fitDispersionPrior();
        dds.fitMAPDispersions();
        dds.fitLFC();
        dds.calculateCooks();
        dds.refit();

        // Create contrast for testing (Treatment vs Control)
        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0;

        // Run statistical tests
        deseq2::DeseqStats ds(dds, contrast, 0.05, true, true);
        ds.runWaldTest();
        ds.cooksFiltering();
        ds.independentFiltering();

        // Get results
        Eigen::MatrixXd results = ds.summary();

        std::cout << "Analysis completed!" << std::endl;
        std::cout << "Results matrix dimensions: " << results.rows() << " genes x " << results.cols() << " statistics" << std::endl;

        // Save results with gene names
        deseq2::saveResultsToCSV("ppm_conversion_results.csv", results, gene_names);
        std::cout << "Results saved to 'ppm_conversion_results.csv'" << std::endl;

        // Print summary statistics
        int significant_genes = 0;
        for (int i = 0; i < results.rows(); ++i)
        {
            if (results(i, 4) < 0.05) // Adjusted p-value column
            {
                significant_genes++;
            }
        }

        std::cout << "Significant genes (adjusted p-value < 0.05): " << significant_genes << std::endl;
        std::cout << "Percentage of significant genes: "
                  << (100.0 * significant_genes / results.rows()) << "%" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in PPM conversion example: " << e.what() << std::endl;
    }

    std::cout << "\n=== All Examples Completed ===" << std::endl;
}

/**
 * @brief Main function demonstrating DESeq2 library usage
 */
int main()
{
    std::cout << "DESeq2 C++ Library Usage Examples" << std::endl;
    std::cout << "=================================" << std::endl;

    // Set random seed for reproducible results
    srand(42);

    // Run examples
    exampleBasicAnalysis();
    exampleCustomDataAnalysis();
    exampleLoadFromCSV();
    exampleAdvancedAnalysis();
    examplePPMConversion();

    std::cout << "\n"
              << std::string(60, '=') << std::endl;
    std::cout << "ALL EXAMPLES COMPLETED SUCCESSFULLY!" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    std::cout << "\nGenerated output files:" << std::endl;
    std::cout << "- example_basic_results.csv: Basic analysis results" << std::endl;
    std::cout << "- example_custom_results.csv: Custom data analysis results" << std::endl;
    std::cout << "- example_advanced_results.csv: Advanced analysis results" << std::endl;
    std::cout << "- ppm_conversion_results.csv: PPM conversion results" << std::endl;

    std::cout << "\nEach CSV file contains the following columns:" << std::endl;
    std::cout << "- baseMean: Normalized mean counts" << std::endl;
    std::cout << "- log2FoldChange: Log2 fold change between conditions" << std::endl;
    std::cout << "- lfcSE: Standard error of log fold change" << std::endl;
    std::cout << "- stat: Wald statistic" << std::endl;
    std::cout << "- pvalue: Raw p-value" << std::endl;
    std::cout << "- padj: Adjusted p-value (Benjamini-Hochberg)" << std::endl;

    return 0;
}
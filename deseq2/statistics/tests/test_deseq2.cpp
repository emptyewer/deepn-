#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <string>

using namespace deseq2;

int main()
{
    try
    {
        std::cout << "Loading synthetic data..." << std::endl;

        // Load data
        Eigen::MatrixXd counts = loadExampleData("raw_counts", "synthetic");
        Eigen::MatrixXd metadata = loadExampleData("metadata", "synthetic");

        // Transpose counts to match expected format (samples x genes)
        counts = counts.transpose().eval();

        std::cout << "Counts matrix shape: " << counts.rows() << " x " << counts.cols() << std::endl;
        std::cout << "Metadata matrix shape: " << metadata.rows() << " x " << metadata.cols() << std::endl;

        // Create DeseqDataSet
        std::cout << "Creating DeseqDataSet..." << std::endl;
        DeseqDataSet dds(counts, metadata, "~condition", true);

        // Step 1: Fit size factors
        std::cout << "Fitting size factors..." << std::endl;
        dds.fitSizeFactors();

        // Step 2: Fit genewise dispersions
        std::cout << "Fitting genewise dispersions..." << std::endl;
        dds.fitGenewiseDispersions();

        // Step 3: Fit dispersion trend
        std::cout << "Fitting dispersion trend..." << std::endl;
        dds.fitDispersionTrend();

        // Step 4: Fit dispersion prior
        std::cout << "Fitting dispersion prior..." << std::endl;
        dds.fitDispersionPrior();

        // Step 5: Fit MAP dispersions
        std::cout << "Fitting MAP dispersions..." << std::endl;
        dds.fitMAPDispersions();

        // Step 6: Fit log fold changes
        std::cout << "Fitting log fold changes..." << std::endl;
        dds.fitLFC();

        // Step 7: Calculate Cooks distances and refit
        std::cout << "Calculating Cooks distances..." << std::endl;
        dds.calculateCooks();
        dds.refit();

        // Step 8: Statistical analysis
        std::cout << "Performing statistical analysis..." << std::endl;

        // Create contrast vector for condition B vs A
        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0; // Test condition B vs A

        DeseqStats ds(dds, contrast, 0.05, true, true);

        // Run Wald test
        std::cout << "Running Wald test..." << std::endl;
        ds.runWaldTest();

        // Cooks filtering
        std::cout << "Applying Cooks filtering..." << std::endl;
        ds.cooksFiltering();

        // Independent filtering and p-value adjustment
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
        std::string output_filename = "cpp_results.csv";
        std::cout << "Saving results to " << output_filename << "..." << std::endl;
        saveResultsToCSV(output_filename, results, gene_names);

        // Print summary statistics
        std::cout << "\nResults summary:" << std::endl;
        std::cout << "Number of genes: " << results.rows() << std::endl;
        std::cout << "Number of significant genes (padj < 0.05): "
                  << (ds.getPAdj().array() < 0.05).count() << std::endl;

        // Compare with Python results
        std::cout << "\nComparing with Python results..." << std::endl;
        Eigen::MatrixXd python_results = loadCSV("/Users/venky/Projects/deepn-plus/PyDESeq2/output_files/synthetic_example/results.csv");

        if (python_results.rows() == results.rows() && python_results.cols() == results.cols())
        {
            double max_diff = 0.0;
            for (int i = 0; i < results.rows(); ++i)
            {
                for (int j = 0; j < results.cols(); ++j)
                {
                    double diff = std::abs(results(i, j) - python_results(i, j));
                    max_diff = std::max(max_diff, diff);
                }
            }

            std::cout << "Maximum difference between C++ and Python results: " << max_diff << std::endl;

            if (max_diff < 1e-6)
            {
                std::cout << "SUCCESS: C++ and Python results match exactly!" << std::endl;
            }
            else
            {
                std::cout << "WARNING: Results differ by more than 1e-6" << std::endl;
            }
        }
        else
        {
            std::cout << "ERROR: Result matrices have different dimensions" << std::endl;
            std::cout << "C++ results: " << results.rows() << " x " << results.cols() << std::endl;
            std::cout << "Python results: " << python_results.rows() << " x " << python_results.cols() << std::endl;
        }

        std::cout << "\nTest completed successfully!" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
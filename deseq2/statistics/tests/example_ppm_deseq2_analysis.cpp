#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

/**
 * @brief Complete DESeq2 analysis example using PPM gene count files
 *
 * This example demonstrates the full workflow:
 * 1. Convert PPM files to DESeq2 format
 * 2. Create metadata for experimental design
 * 3. Run complete DESeq2 analysis pipeline
 * 4. Perform statistical testing
 * 5. Generate and save results
 */
int main()
{
    std::cout << "============================================================" << std::endl;
    std::cout << "COMPLETE DESEQ2 ANALYSIS WITH PPM CONVERSION" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "This example demonstrates a complete DESeq2 analysis workflow" << std::endl;
    std::cout << "using PPM gene count files from the examples folder." << std::endl;
    std::cout << std::endl;

    try
    {
        // Step 1: Define input files
        std::cout << "Step 1: Loading PPM gene count files..." << std::endl;
        std::vector<std::string> file_paths = {
            "../example/HDPTP_Tail_2_NON_50M_summary.csv",
            "../example/HDPTP_Tail_2_SEL_50N_summary.csv",
            "../example/HDPTP_V_NON_51M_summary.csv",
            "../example/HDPTP_V_SEL_51N_summary.csv"};

        std::cout << "Input files:" << std::endl;
        for (size_t i = 0; i < file_paths.size(); ++i)
        {
            std::cout << "  " << (i + 1) << ". " << file_paths[i] << std::endl;
        }
        std::cout << std::endl;

        // Step 2: Convert PPM files to DESeq2 format
        std::cout << "Step 2: Converting PPM files to DESeq2 format..." << std::endl;
        auto [count_matrix, gene_names, sample_names] = deseq2::convertPpmToDeseq2Format(file_paths);

        std::cout << "Conversion completed successfully!" << std::endl;
        std::cout << "  - Count matrix: " << count_matrix.rows() << " genes × " << count_matrix.cols() << " samples" << std::endl;
        std::cout << "  - Unique genes: " << gene_names.size() << std::endl;
        std::cout << "  - Samples: ";
        for (const auto &name : sample_names)
        {
            std::cout << name << " ";
        }
        std::cout << std::endl
                  << std::endl;

        // Step 3: Create experimental metadata
        std::cout << "Step 3: Creating experimental metadata..." << std::endl;

        // Define experimental design
        // Based on sample names, we can identify different conditions:
        // - HDPTP_Tail_2_NON_50M: Control (NON = non-selected)
        // - HDPTP_Tail_2_SEL_50N: Treatment (SEL = selected)
        // - HDPTP_V_NON_51M: Control (NON = non-selected)
        // - HDPTP_V_SEL_51N: Treatment (SEL = selected)

        Eigen::MatrixXd metadata = Eigen::MatrixXd::Zero(count_matrix.cols(), 2);

        // Column 1: Condition (0 = Control, 1 = Treatment)
        // Column 2: Replicate group (0 = HDPTP_Tail, 1 = HDPTP_V)
        for (size_t i = 0; i < sample_names.size(); ++i)
        {
            const std::string &name = sample_names[i];

            // Determine condition (NON = Control, SEL = Treatment)
            if (name.find("NON") != std::string::npos)
            {
                metadata(i, 0) = 0.0; // Control
            }
            else if (name.find("SEL") != std::string::npos)
            {
                metadata(i, 0) = 1.0; // Treatment
            }

            // Determine replicate group
            if (name.find("HDPTP_Tail") != std::string::npos)
            {
                metadata(i, 1) = 0.0; // HDPTP_Tail group
            }
            else if (name.find("HDPTP_V") != std::string::npos)
            {
                metadata(i, 1) = 1.0; // HDPTP_V group
            }
        }

        std::cout << "Experimental design:" << std::endl;
        for (size_t i = 0; i < sample_names.size(); ++i)
        {
            std::string condition = (metadata(i, 0) == 0) ? "Control" : "Treatment";
            std::string group = (metadata(i, 1) == 0) ? "HDPTP_Tail" : "HDPTP_V";
            std::cout << "  " << sample_names[i] << ": " << condition << " (" << group << ")" << std::endl;
        }
        std::cout << std::endl;

        // Step 4: Prepare data for DESeq2 analysis
        std::cout << "Step 4: Preparing data for DESeq2 analysis..." << std::endl;

        // Transpose count matrix to match DESeq2 format (samples × genes)
        Eigen::MatrixXd counts_transposed = count_matrix.transpose();

        std::cout << "  - Transposed count matrix: " << counts_transposed.rows() << " samples × " << counts_transposed.cols() << " genes" << std::endl;
        std::cout << "  - Metadata matrix: " << metadata.rows() << " samples × " << metadata.cols() << " variables" << std::endl;

        // Validate data
        deseq2::testValidCounts(counts_transposed);
        std::cout << "  - Data validation passed" << std::endl;
        std::cout << std::endl;

        // Step 5: Run DESeq2 analysis pipeline
        std::cout << "Step 5: Running DESeq2 analysis pipeline..." << std::endl;

        // Create DeseqDataSet with condition as the main factor
        std::cout << "  - Creating DeseqDataSet..." << std::endl;
        deseq2::DeseqDataSet dds(counts_transposed, metadata, "~condition", true);

        // Fit size factors
        std::cout << "  - Fitting size factors..." << std::endl;
        dds.fitSizeFactors();

        // Fit dispersions
        std::cout << "  - Fitting genewise dispersions..." << std::endl;
        dds.fitGenewiseDispersions();

        std::cout << "  - Fitting dispersion trend..." << std::endl;
        dds.fitDispersionTrend();

        std::cout << "  - Fitting dispersion prior..." << std::endl;
        dds.fitDispersionPrior();

        std::cout << "  - Fitting MAP dispersions..." << std::endl;
        dds.fitMAPDispersions();

        // Fit log fold changes
        std::cout << "  - Fitting log fold changes..." << std::endl;
        dds.fitLFC();

        // Handle outliers
        std::cout << "  - Calculating Cooks distances and refitting..." << std::endl;
        dds.calculateCooks();
        dds.refit();

        std::cout << "DESeq2 analysis pipeline completed!" << std::endl;
        std::cout << std::endl;

        // Step 6: Statistical testing
        std::cout << "Step 6: Performing statistical testing..." << std::endl;

        // Create contrast for Treatment vs Control
        Eigen::VectorXd contrast(2);
        contrast << 0.0, 1.0; // Test Treatment vs Control

        std::cout << "  - Testing Treatment vs Control..." << std::endl;

        // Create DeseqStats object
        deseq2::DeseqStats ds(dds, contrast, 0.05, true, true);

        // Run statistical tests
        std::cout << "  - Running Wald test..." << std::endl;
        ds.runWaldTest();

        std::cout << "  - Applying Cooks filtering..." << std::endl;
        ds.cooksFiltering();

        std::cout << "  - Applying independent filtering..." << std::endl;
        ds.independentFiltering();

        // Get results
        Eigen::MatrixXd results = ds.summary();

        std::cout << "Statistical testing completed!" << std::endl;
        std::cout << std::endl;

        // Step 7: Results analysis and summary
        std::cout << "Step 7: Analyzing results..." << std::endl;

        std::cout << "Results matrix dimensions: " << results.rows() << " genes × " << results.cols() << " statistics" << std::endl;

        // Count significant genes
        int significant_genes = 0;
        int upregulated = 0;
        int downregulated = 0;

        for (int i = 0; i < results.rows(); ++i)
        {
            if (results(i, 4) < 0.05) // Adjusted p-value column
            {
                significant_genes++;
                if (results(i, 1) > 0) // Log2 fold change column
                {
                    upregulated++;
                }
                else
                {
                    downregulated++;
                }
            }
        }

        std::cout << "Differential expression summary:" << std::endl;
        std::cout << "  - Total genes analyzed: " << results.rows() << std::endl;
        std::cout << "  - Significant genes (padj < 0.05): " << significant_genes << std::endl;
        std::cout << "  - Upregulated: " << upregulated << std::endl;
        std::cout << "  - Downregulated: " << downregulated << std::endl;
        std::cout << "  - Percentage significant: " << std::fixed << std::setprecision(2)
                  << (100.0 * significant_genes / results.rows()) << "%" << std::endl;
        std::cout << std::endl;

        // Step 8: Find top differentially expressed genes
        std::cout << "Step 8: Top differentially expressed genes..." << std::endl;

        // Create vector of gene indices sorted by adjusted p-value
        std::vector<std::pair<double, int>> gene_pvalues;
        for (int i = 0; i < results.rows(); ++i)
        {
            gene_pvalues.push_back({results(i, 4), i}); // Adjusted p-value
        }

        // Sort by p-value (ascending)
        std::sort(gene_pvalues.begin(), gene_pvalues.end());

        std::cout << "Top 10 most significant genes:" << std::endl;
        std::cout << std::setw(20) << "Gene" << std::setw(15) << "log2FC"
                  << std::setw(15) << "padj" << std::setw(15) << "baseMean" << std::endl;
        std::cout << std::string(65, '-') << std::endl;

        for (int i = 0; i < std::min(10, static_cast<int>(gene_pvalues.size())); ++i)
        {
            int gene_idx = gene_pvalues[i].second;
            std::string gene_name = gene_names[gene_idx];
            double log2fc = results(gene_idx, 1);
            double padj = results(gene_idx, 4);
            double base_mean = results(gene_idx, 0);

            std::cout << std::setw(20) << gene_name
                      << std::setw(15) << std::fixed << std::setprecision(3) << log2fc
                      << std::setw(15) << std::scientific << std::setprecision(2) << padj
                      << std::setw(15) << std::fixed << std::setprecision(1) << base_mean << std::endl;
        }
        std::cout << std::endl;

        // Step 9: Save results
        std::cout << "Step 9: Saving results..." << std::endl;

        // Save main results
        deseq2::saveResultsToCSV("ppm_deseq2_results.csv", results, gene_names);
        std::cout << "  - Main results saved to 'ppm_deseq2_results.csv'" << std::endl;

        // Save significant genes only
        std::vector<std::string> significant_gene_names;
        Eigen::MatrixXd significant_results(significant_genes, results.cols());
        int sig_idx = 0;

        for (int i = 0; i < results.rows(); ++i)
        {
            if (results(i, 4) < 0.05)
            {
                significant_gene_names.push_back(gene_names[i]);
                significant_results.row(sig_idx) = results.row(i);
                sig_idx++;
            }
        }

        deseq2::saveResultsToCSV("ppm_deseq2_significant.csv", significant_results, significant_gene_names);
        std::cout << "  - Significant genes saved to 'ppm_deseq2_significant.csv'" << std::endl;

        // Save metadata
        std::ofstream metadata_file("ppm_deseq2_metadata.csv");
        if (metadata_file.is_open())
        {
            metadata_file << "sample,condition,group" << std::endl;
            for (size_t i = 0; i < sample_names.size(); ++i)
            {
                std::string condition = (metadata(i, 0) == 0) ? "Control" : "Treatment";
                std::string group = (metadata(i, 1) == 0) ? "HDPTP_Tail" : "HDPTP_V";
                metadata_file << sample_names[i] << "," << condition << "," << group << std::endl;
            }
            metadata_file.close();
            std::cout << "  - Metadata saved to 'ppm_deseq2_metadata.csv'" << std::endl;
        }

        // Save converted count matrix
        std::ofstream count_file("ppm_deseq2_counts.csv");
        if (count_file.is_open())
        {
            count_file << "Gene";
            for (const auto &sample : sample_names)
            {
                count_file << "," << sample;
            }
            count_file << std::endl;

            for (int i = 0; i < count_matrix.rows(); ++i)
            {
                count_file << gene_names[i];
                for (int j = 0; j < count_matrix.cols(); ++j)
                {
                    count_file << "," << count_matrix(i, j);
                }
                count_file << std::endl;
            }
            count_file.close();
            std::cout << "  - Count matrix saved to 'ppm_deseq2_counts.csv'" << std::endl;
        }

        std::cout << std::endl;

        // Step 10: Summary statistics
        std::cout << "Step 10: Summary statistics..." << std::endl;

        // Calculate summary statistics
        double total_counts = count_matrix.sum();
        double mean_counts_per_gene = total_counts / count_matrix.rows();
        double mean_counts_per_sample = total_counts / count_matrix.cols();

        // Count genes with zero expression
        int zero_expression_genes = 0;
        for (int i = 0; i < count_matrix.rows(); ++i)
        {
            if (count_matrix.row(i).sum() == 0)
            {
                zero_expression_genes++;
            }
        }

        std::cout << "Data summary:" << std::endl;
        std::cout << "  - Total counts: " << std::fixed << std::setprecision(0) << total_counts << std::endl;
        std::cout << "  - Mean counts per gene: " << std::fixed << std::setprecision(1) << mean_counts_per_gene << std::endl;
        std::cout << "  - Mean counts per sample: " << std::fixed << std::setprecision(0) << mean_counts_per_sample << std::endl;
        std::cout << "  - Genes with zero expression: " << zero_expression_genes
                  << " (" << std::fixed << std::setprecision(1) << (100.0 * zero_expression_genes / count_matrix.rows()) << "%)" << std::endl;

        // Size factors summary
        const Eigen::VectorXd &size_factors = dds.getSizeFactors();
        std::cout << "  - Size factors range: [" << size_factors.minCoeff() << ", " << size_factors.maxCoeff() << "]" << std::endl;

        // Dispersion summary
        const Eigen::VectorXd &dispersions = dds.getDispersions();
        std::cout << "  - Dispersion range: [" << dispersions.minCoeff() << ", " << dispersions.maxCoeff() << "]" << std::endl;

        std::cout << std::endl;

        std::cout << "============================================================" << std::endl;
        std::cout << "ANALYSIS COMPLETED SUCCESSFULLY!" << std::endl;
        std::cout << "============================================================" << std::endl;
        std::cout << std::endl;
        std::cout << "Generated output files:" << std::endl;
        std::cout << "  - ppm_deseq2_results.csv: Complete DESeq2 results" << std::endl;
        std::cout << "  - ppm_deseq2_significant.csv: Significant genes only" << std::endl;
        std::cout << "  - ppm_deseq2_metadata.csv: Sample metadata" << std::endl;
        std::cout << "  - ppm_deseq2_counts.csv: Converted count matrix" << std::endl;
        std::cout << std::endl;
        std::cout << "Results columns:" << std::endl;
        std::cout << "  - baseMean: Normalized mean counts" << std::endl;
        std::cout << "  - log2FoldChange: Log2 fold change (Treatment vs Control)" << std::endl;
        std::cout << "  - lfcSE: Standard error of log fold change" << std::endl;
        std::cout << "  - stat: Wald statistic" << std::endl;
        std::cout << "  - pvalue: Raw p-value" << std::endl;
        std::cout << "  - padj: Adjusted p-value (Benjamini-Hochberg)" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
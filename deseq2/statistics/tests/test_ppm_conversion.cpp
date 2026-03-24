#include "utils.h"
#include <iostream>
#include <vector>
#include <string>

/**
 * @brief Test program for PPM gene count file conversion
 *
 * This program demonstrates how to convert PPM gene count files to DESeq2 format.
 * It uses the existing example files to show the conversion process.
 */
int main()
{
    std::cout << "=== PPM Gene Count File Conversion Test ===" << std::endl;
    std::cout << "This program converts PPM gene count files to DESeq2-compatible format." << std::endl;

    try
    {
        // Example file paths (using the existing example files)
        std::vector<std::string> file_paths = {
            "../example/HDPTP_Tail_2_NON_50M_summary.csv",
            "../example/HDPTP_Tail_2_SEL_50N_summary.csv",
            "../example/HDPTP_V_NON_51M_summary.csv",
            "../example/HDPTP_V_SEL_51N_summary.csv"};

        std::cout << "\nInput files:" << std::endl;
        for (const auto &file_path : file_paths)
        {
            std::cout << "  - " << file_path << std::endl;
        }

        std::cout << "\nStep 1: Extracting sample information..." << std::endl;
        for (const auto &file_path : file_paths)
        {
            auto [sample_name, total_reads, total_hits] = deseq2::extractSampleInfo(file_path);
            std::cout << "  Sample: " << sample_name
                      << ", Total Reads: " << total_reads
                      << ", Total Hits: " << total_hits << std::endl;
        }

        std::cout << "\nStep 2: Converting PPM files to DESeq2 format..." << std::endl;

        // Convert PPM files to DESeq2 format
        auto [count_matrix, gene_names, sample_names] = deseq2::convertPpmToDeseq2Format(file_paths);

        std::cout << "Conversion completed successfully!" << std::endl;
        std::cout << "Count matrix dimensions: " << count_matrix.rows() << " genes x " << count_matrix.cols() << " samples" << std::endl;
        std::cout << "Number of unique genes: " << gene_names.size() << std::endl;

        std::cout << "\nSample names:" << std::endl;
        for (size_t i = 0; i < sample_names.size(); ++i)
        {
            std::cout << "  Sample " << (i + 1) << ": " << sample_names[i] << std::endl;
        }

        std::cout << "\nFirst 10 genes:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(10), gene_names.size()); ++i)
        {
            std::cout << "  Gene " << (i + 1) << ": " << gene_names[i] << std::endl;
        }

        std::cout << "\nCount matrix preview (first 5 genes x 4 samples):" << std::endl;
        for (int i = 0; i < std::min(5, static_cast<int>(count_matrix.rows())); ++i)
        {
            std::cout << "  " << gene_names[i] << ": ";
            for (int j = 0; j < count_matrix.cols(); ++j)
            {
                std::cout << count_matrix(i, j) << " ";
            }
            std::cout << std::endl;
        }

        // Calculate some statistics
        double total_counts = count_matrix.sum();
        double mean_counts_per_gene = total_counts / count_matrix.rows();
        double mean_counts_per_sample = total_counts / count_matrix.cols();

        std::cout << "\nStatistics:" << std::endl;
        std::cout << "  Total counts across all genes and samples: " << total_counts << std::endl;
        std::cout << "  Mean counts per gene: " << mean_counts_per_gene << std::endl;
        std::cout << "  Mean counts per sample: " << mean_counts_per_sample << std::endl;

        // Count genes with zero expression
        int zero_expression_genes = 0;
        for (int i = 0; i < count_matrix.rows(); ++i)
        {
            if (count_matrix.row(i).sum() == 0)
            {
                zero_expression_genes++;
            }
        }

        std::cout << "  Genes with zero expression: " << zero_expression_genes
                  << " (" << (100.0 * zero_expression_genes / count_matrix.rows()) << "%)" << std::endl;

        // Save the converted data
        std::cout << "\nStep 3: Saving converted data..." << std::endl;

        // Save count matrix
        std::ofstream count_file("converted_counts.csv");
        if (count_file.is_open())
        {
            // Write header
            count_file << "Gene";
            for (const auto &sample : sample_names)
            {
                count_file << "," << sample;
            }
            count_file << std::endl;

            // Write data
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
            std::cout << "  Count matrix saved to 'converted_counts.csv'" << std::endl;
        }

        // Save metadata
        std::ofstream metadata_file("converted_metadata.csv");
        if (metadata_file.is_open())
        {
            metadata_file << "sample,condition" << std::endl;
            for (size_t i = 0; i < sample_names.size(); ++i)
            {
                std::string condition = (sample_names[i].find("SEL") != std::string::npos) ? "Treatment" : "Control";
                metadata_file << sample_names[i] << "," << condition << std::endl;
            }
            metadata_file.close();
            std::cout << "  Metadata saved to 'converted_metadata.csv'" << std::endl;
        }

        std::cout << "\n=== Conversion Test Completed Successfully ===" << std::endl;
        std::cout << "\nThe converted data is now ready for DESeq2 analysis!" << std::endl;
        std::cout << "You can use the count matrix and metadata with the DeseqDataSet class." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
# PPM Gene Count File Conversion Guide

This guide explains how to convert PPM (Parts Per Million) gene count files to DESeq2-compatible format using the C++ DESeq2 library.

## Overview

The PPM conversion functionality allows you to convert gene count files in PPM format to the raw count format required by DESeq2 for differential expression analysis. This is particularly useful when working with data from gene counting pipelines that output results in PPM format.

## Input File Format

Your PPM files should have the following format:

```csv
File:,HDPTP_Tail_2_NON_50M.sam
 , TotalReads , 4599960
 , TotalHits (count), 4387489
Chromosome , GeneName , PPM , NCBI_Acc
19_KI270882v1_alt , KIR3DL1 , 0.0,NM_013289,NM_013289,NM_013289,...
19_KI270882v1_alt , KIR2DL1 , 0.0,NM_014218,NM_014218,NM_014218,...
19_KI270882v1_alt , KIR2DL5B , 33.5043575038,NM_001018081,NM_001018081,...
```

### File Structure

1. **Header Information** (first 4 lines):
   - Line 1: File name information
   - Line 2: Total reads count
   - Line 3: Total hits count
   - Line 4: Column headers

2. **Data Lines**:
   - Column 1: Chromosome
   - Column 2: Gene name
   - Column 3: PPM value
   - Column 4+: NCBI accession numbers (optional)

## Conversion Process

The conversion process involves:

1. **PPM to Raw Counts**: Converting PPM values to raw read counts using the formula:
   ```
   Raw Count = PPM × Total Reads / 1,000,000
   ```

2. **Data Aggregation**: Combining multiple files into a single count matrix

3. **Format Standardization**: Creating DESeq2-compatible output format

## API Functions

### Main Conversion Function

```cpp
std::tuple<Eigen::MatrixXd, std::vector<std::string>, std::vector<std::string>> 
convertPpmToDeseq2Format(
    const std::vector<std::string> &file_paths,
    const std::vector<std::string> &sample_names = {},
    const std::vector<int> &total_reads = {}
);
```

**Parameters:**
- `file_paths`: Vector of paths to PPM files
- `sample_names`: Optional vector of sample names (uses filenames if empty)
- `total_reads`: Optional vector of total reads per sample (extracted from files if empty)

**Returns:**
- `count_matrix`: Matrix with genes as rows and samples as columns
- `gene_names`: Vector of gene names
- `sample_names`: Vector of sample names

### Helper Functions

#### Extract Sample Information

```cpp
std::tuple<std::string, int, int> extractSampleInfo(const std::string &file_path);
```

Extracts sample name, total reads, and total hits from a PPM file header.

#### Parse Single PPM File

```cpp
std::tuple<std::vector<std::string>, std::vector<double>, std::vector<int>> 
parsePpmFile(const std::string &file_path, int total_reads = 0);
```

Parses a single PPM file and returns gene names, PPM values, and raw counts.

## Usage Examples

### Basic Usage

```cpp
#include "utils.h"
#include <vector>
#include <string>

int main()
{
    // Define file paths
    std::vector<std::string> file_paths = {
        "sample1_ppm.csv",
        "sample2_ppm.csv",
        "sample3_ppm.csv"
    };

    // Convert PPM files to DESeq2 format
    auto [count_matrix, gene_names, sample_names] = 
        deseq2::convertPpmToDeseq2Format(file_paths);

    // Use the converted data for DESeq2 analysis
    // count_matrix is ready for DeseqDataSet constructor
}
```

### Advanced Usage with Custom Names

```cpp
#include "utils.h"
#include <vector>
#include <string>

int main()
{
    std::vector<std::string> file_paths = {
        "control_rep1.csv",
        "control_rep2.csv",
        "treatment_rep1.csv",
        "treatment_rep2.csv"
    };

    // Custom sample names
    std::vector<std::string> sample_names = {
        "Control_1", "Control_2", "Treatment_1", "Treatment_2"
    };

    // Custom total reads (optional)
    std::vector<int> total_reads = {5000000, 4800000, 5200000, 5100000};

    // Convert with custom parameters
    auto [count_matrix, gene_names, final_sample_names] = 
        deseq2::convertPpmToDeseq2Format(file_paths, sample_names, total_reads);
}
```

### Integration with DESeq2 Analysis

```cpp
#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"

int main()
{
    // Convert PPM files
    std::vector<std::string> file_paths = {
        "control1.csv", "control2.csv",
        "treatment1.csv", "treatment2.csv"
    };
    
    auto [count_matrix, gene_names, sample_names] = 
        deseq2::convertPpmToDeseq2Format(file_paths);

    // Create metadata
    Eigen::MatrixXd metadata(4, 1);
    metadata << 0, 0, 1, 1; // Control, Control, Treatment, Treatment

    // Transpose for DESeq2 format (samples x genes)
    Eigen::MatrixXd counts_transposed = count_matrix.transpose();

    // Run DESeq2 analysis
    deseq2::DeseqDataSet dds(counts_transposed, metadata, "~condition");
    dds.fitSizeFactors();
    dds.fitGenewiseDispersions();
    dds.fitDispersionTrend();
    dds.fitDispersionPrior();
    dds.fitMAPDispersions();
    dds.fitLFC();

    // Statistical testing
    Eigen::VectorXd contrast(2);
    contrast << 0.0, 1.0; // Treatment vs Control
    
    deseq2::DeseqStats ds(dds, contrast, 0.05);
    ds.runWaldTest();
    ds.cooksFiltering();
    ds.independentFiltering();
    
    Eigen::MatrixXd results = ds.summary();
    
    // Save results with gene names
    deseq2::saveResultsToCSV("results.csv", results, gene_names);
}
```

## Test Program

A test program is included to demonstrate the conversion functionality:

```bash
# Build the test program
cd deseq2
mkdir build && cd build
cmake ..
make test_ppm_conversion

# Run the test
./test_ppm_conversion
```

The test program will:
1. Convert example PPM files
2. Display conversion statistics
3. Save converted data to CSV files
4. Show how to integrate with DESeq2 analysis

## Output Files

The conversion process generates:

1. **converted_counts.csv**: Count matrix with genes as rows and samples as columns
2. **converted_metadata.csv**: Sample metadata with conditions

## Error Handling

The conversion functions include comprehensive error handling:

- **File not found**: Throws `std::runtime_error` with file path
- **Invalid PPM values**: Skips lines that can't be parsed
- **Missing total reads**: Extracts from file header automatically
- **Empty files**: Returns empty matrices with appropriate warnings

## Performance Considerations

- **Memory usage**: Large files are processed line by line to minimize memory usage
- **Gene deduplication**: Uses efficient set-based deduplication
- **Sorting**: Gene names are sorted for consistent output

## Troubleshooting

### Common Issues

1. **File format errors**: Ensure your PPM files follow the expected format
2. **Missing total reads**: The converter will extract from file headers automatically
3. **Gene name conflicts**: Duplicate gene names are handled automatically
4. **Memory issues**: For very large files, consider processing in batches

### Debug Information

The test program provides detailed output including:
- Sample information extraction
- Conversion statistics
- Data preview
- Error messages

## Integration with Existing Workflows

The converted data can be easily integrated with:

- **R DESeq2**: Save as CSV and import into R
- **Python**: Use pandas to read the CSV files
- **Other C++ analysis**: Direct use with the DESeq2 C++ library

## Example Data

The library includes example PPM files in the `example/` directory:
- `HDPTP_Tail_2_NON_50M_summary.csv`
- `HDPTP_Tail_2_SEL_50N_summary.csv`
- `HDPTP_V_NON_51M_summary.csv`
- `HDPTP_V_SEL_51N_summary.csv`

These files demonstrate the expected format and can be used for testing the conversion functionality. 
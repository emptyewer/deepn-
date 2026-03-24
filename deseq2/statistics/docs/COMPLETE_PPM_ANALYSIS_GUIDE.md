# Complete PPM DESeq2 Analysis Guide

This guide demonstrates how to perform a complete DESeq2 differential expression analysis using PPM (Parts Per Million) gene count files from the examples folder.

## Overview

The `example_ppm_deseq2_analysis.cpp` program provides a complete workflow that:

1. **Converts PPM files** to DESeq2-compatible format
2. **Creates experimental metadata** based on sample names
3. **Runs the full DESeq2 analysis pipeline**
4. **Performs statistical testing** for differential expression
5. **Generates comprehensive results** and summary statistics
6. **Saves multiple output files** for further analysis

## Quick Start

### Building the Example

```bash
cd deseq2
mkdir build && cd build
cmake ..
make example_ppm_deseq2_analysis
```

### Running the Analysis

```bash
./example_ppm_deseq2_analysis
```

## Input Data

The example uses the following PPM files from the `example/` folder:

- `HDPTP_Tail_2_NON_50M_summary.csv` - Control sample (HDPTP_Tail, non-selected)
- `HDPTP_Tail_2_SEL_50N_summary.csv` - Treatment sample (HDPTP_Tail, selected)
- `HDPTP_V_NON_51M_summary.csv` - Control sample (HDPTP_V, non-selected)
- `HDPTP_V_SEL_51N_summary.csv` - Treatment sample (HDPTP_V, selected)

### Experimental Design

The analysis automatically interprets the experimental design based on sample names:

- **Condition**: NON (Control) vs SEL (Treatment)
- **Group**: HDPTP_Tail vs HDPTP_V (replicate groups)

## Analysis Workflow

### Step 1: PPM File Conversion

The program converts PPM values to raw counts using the formula:

```
Raw Count = PPM × Total Reads / 1,000,000
```

**Output**: Count matrix with 19,117 genes × 4 samples

### Step 2: Experimental Metadata Creation

Automatically creates metadata based on sample names:

- **Column 1**: Condition (0 = Control, 1 = Treatment)
- **Column 2**: Replicate group (0 = HDPTP_Tail, 1 = HDPTP_V)

### Step 3: DESeq2 Analysis Pipeline

Runs the complete DESeq2 workflow:

1. **Size Factor Normalization** - Median-of-ratios method
2. **Dispersion Estimation** - Gene-wise, trend, and MAP dispersions
3. **Log Fold Change Estimation** - Using negative binomial GLM
4. **Outlier Detection** - Cook's distance filtering and refitting

### Step 4: Statistical Testing

Performs differential expression analysis:

- **Wald Test** - Treatment vs Control comparison
- **Cooks Filtering** - Removes outlier samples
- **Independent Filtering** - Removes low-count genes

### Step 5: Results Analysis

Provides comprehensive results summary:

- **19,117 genes analyzed**
- **1,763 significant genes** (9.22%)
- **1,621 upregulated** in treatment
- **142 downregulated** in treatment

## Output Files

The analysis generates four output files:

### 1. `ppm_deseq2_results.csv` - Complete Results

Contains all genes with the following columns:

- **GeneName**: Gene identifier
- **BaseMean**: Normalized mean counts
- **Log2FoldChange**: Log2 fold change (Treatment vs Control)
- **Log2FoldChangeStdErr**: Standard error of log fold change
- **WaldStatistic**: Wald test statistic
- **Pvalue**: Raw p-value
- **Padj**: Adjusted p-value (Benjamini-Hochberg)

### 2. `ppm_deseq2_significant.csv` - Significant Genes Only

Contains only genes with adjusted p-value < 0.05 (1,763 genes)

### 3. `ppm_deseq2_metadata.csv` - Sample Information

Sample metadata with conditions and groups:

```csv
sample,condition,group
HDPTP_Tail_2_NON_50M,Control,HDPTP_Tail
HDPTP_Tail_2_SEL_50N,Treatment,HDPTP_Tail
HDPTP_V_NON_51M,Control,HDPTP_V
HDPTP_V_SEL_51N,Treatment,HDPTP_V
```

### 4. `ppm_deseq2_counts.csv` - Converted Count Matrix

Raw count matrix with genes as rows and samples as columns.

## Key Results

### Top Differentially Expressed Genes

The analysis identified several highly significant genes:

1. **RWDD1** - log2FC: 14.02, padj: 4.17e-185
2. **DLG3** - log2FC: 10.53, padj: 6.95e-117
3. **HGF** - log2FC: 11.68, padj: 1.06e-112
4. **RLIM** - log2FC: 9.64, padj: 4.00e-99
5. **CYB5R4** - log2FC: 9.49, padj: 1.60e-93

### Data Quality Metrics

- **Total counts**: 12,556,453
- **Mean counts per gene**: 656.8
- **Mean counts per sample**: 3,139,113
- **Genes with zero expression**: 7,386 (38.6%)
- **Size factors range**: [0.2, 4.5]
- **Dispersion range**: [0.0, 0.3]

## Customization Options

### Modifying Input Files

To use your own PPM files, modify the `file_paths` vector:

```cpp
std::vector<std::string> file_paths = {
    "path/to/your/control1.csv",
    "path/to/your/treatment1.csv",
    "path/to/your/control2.csv",
    "path/to/your/treatment2.csv"
};
```

### Adjusting Statistical Parameters

Modify the statistical testing parameters:

```cpp
// Change significance threshold
deseq2::DeseqStats ds(dds, contrast, 0.01, true, true); // alpha = 0.01

// Modify contrast vector for different comparisons
Eigen::VectorXd contrast(2);
contrast << 0.0, 1.0; // Treatment vs Control
```

### Custom Experimental Design

For different experimental designs, modify the metadata creation:

```cpp
// Add more experimental factors
Eigen::MatrixXd metadata = Eigen::MatrixXd::Zero(count_matrix.cols(), 3);
// Column 1: Condition
// Column 2: Replicate group
// Column 3: Additional factor
```

## Integration with Other Tools

### R Integration

Load results into R for further analysis:

```r
# Load results
results <- read.csv("ppm_deseq2_results.csv")

# Create volcano plot
library(ggplot2)
ggplot(results, aes(x = Log2FoldChange, y = -log10(Padj))) +
  geom_point(aes(color = Padj < 0.05)) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
```

### Python Integration

Use pandas for data analysis:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
results = pd.read_csv("ppm_deseq2_results.csv")

# Filter significant genes
significant = results[results['Padj'] < 0.05]

# Create volcano plot
plt.scatter(results['Log2FoldChange'], -np.log10(results['Padj']))
plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(Adjusted P-value)')
plt.title('Volcano Plot')
plt.show()
```

## Troubleshooting

### Common Issues

1. **File not found errors**: Ensure PPM files are in the correct location
2. **Memory issues**: For very large datasets, consider processing in batches
3. **Convergence warnings**: May occur with low-count genes, usually handled automatically

### Performance Optimization

- **Large datasets**: The program processes files efficiently but may take time for very large gene sets
- **Memory usage**: Approximately 1GB for 20,000 genes × 4 samples
- **Processing time**: ~30 seconds for the example dataset

## Advanced Usage

### Batch Processing

For multiple experimental conditions:

```cpp
// Process multiple file sets
for (const auto &file_set : file_sets) {
    auto [count_matrix, gene_names, sample_names] =
        deseq2::convertPpmToDeseq2Format(file_set);
    // Run analysis for each set
}
```

### Quality Control

Add quality control steps:

```cpp
// Check for low-quality samples
for (int i = 0; i < counts_transposed.rows(); ++i) {
    double total_counts = counts_transposed.row(i).sum();
    if (total_counts < 1000000) {
        std::cout << "Warning: Low counts in sample " << i << std::endl;
    }
}
```

## Example Output Summary

The complete analysis demonstrates:

- **Robust PPM conversion** with proper count estimation
- **Comprehensive DESeq2 pipeline** execution
- **Statistical rigor** with multiple testing corrections
- **Rich results** with detailed gene-level statistics
- **Multiple output formats** for downstream analysis

This example serves as a template for analyzing PPM gene count data with DESeq2, providing a complete workflow from raw data to interpretable results.

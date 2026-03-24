# C++ DESeq2 Implementation

This is a C++ implementation of the DESeq2 algorithm for differential expression analysis, based on the Python PyDESeq2 library.

## Overview

The implementation includes:

- **DeseqDataSet**: Core class for dispersion and log fold-change estimation
- **DeseqStats**: Statistical tests for differential expression analysis
- **Utils**: Utility functions for data loading, mathematical operations, and file I/O

## Dependencies

- **Eigen3**: Linear algebra library
- **Boost**: Filesystem and system libraries
- **C++17**: Modern C++ standard

## Building

### Prerequisites

Install the required dependencies:

```bash
# On macOS with Homebrew
brew install eigen boost

# On Ubuntu/Debian
sudo apt-get install libeigen3-dev libboost-all-dev

# On CentOS/RHEL
sudo yum install eigen3-devel boost-devel
```

### Build Instructions

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake ..

# Build
make -j$(nproc)

# Run the test
./test_deseq2
```

## Usage

The main test program (`test_deseq2`) replicates the Python step-by-step workflow:

1. **Data Loading**: Loads synthetic count and metadata from PyDESeq2
2. **Size Factor Fitting**: Normalizes data using median-of-ratios method
3. **Dispersion Estimation**: Fits gene-wise dispersions and trend
4. **Log Fold Change**: Estimates log fold changes using IRLS
5. **Statistical Testing**: Performs Wald tests and p-value adjustment
6. **Results Comparison**: Compares output with Python results

## File Structure

```
deseq2/
├── include/
│   ├── deseq_dataset.h    # DeseqDataSet class header
│   ├── deseq_stats.h      # DeseqStats class header
│   └── utils.h           # Utility functions header
├── src/
│   ├── deseq_dataset.cpp  # DeseqDataSet implementation
│   ├── deseq_stats.cpp    # DeseqStats implementation
│   ├── utils.cpp         # Utility functions implementation
│   └── test_deseq2.cpp   # Main test program
├── CMakeLists.txt        # Build configuration
└── README.md            # This file
```

## Algorithm Details

### Size Factor Normalization

Uses the median-of-ratios method to normalize for sequencing depth differences between samples.

### Dispersion Estimation

1. **Gene-wise dispersions**: Maximum likelihood estimation for each gene
2. **Dispersion trend**: Parametric fit of dispersion vs. mean relationship
3. **MAP dispersions**: Maximum a posteriori estimates using empirical Bayes

### Statistical Testing

1. **Wald test**: Tests for differential expression using negative binomial GLM
2. **P-value adjustment**: Benjamini-Hochberg procedure for multiple testing
3. **Independent filtering**: Optimizes power by filtering low-count genes

## Output

The program generates:

- `cpp_results.csv`: Results matrix with columns:
  - `baseMean`: Normalized mean counts
  - `log2FoldChange`: Log2 fold change between conditions
  - `lfcSE`: Standard error of log fold change
  - `stat`: Wald statistic
  - `pvalue`: Raw p-value
  - `padj`: Adjusted p-value

## Comparison with Python

The C++ implementation is designed to produce results that match the Python PyDESeq2 library exactly. The test program automatically compares outputs and reports any differences.

## Notes

- This is a simplified implementation focused on the core DESeq2 algorithm
- Some advanced features (e.g., complex design matrices, LFC shrinkage) are simplified
- The implementation follows the same mathematical framework as the original DESeq2 R package

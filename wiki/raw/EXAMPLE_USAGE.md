# DESeq2 C++ Library - Usage Examples

This document provides comprehensive examples and step-by-step instructions for using the DESeq2 C++ library in your projects.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Building the Examples](#building-the-examples)
3. [Example Overview](#example-overview)
4. [Step-by-Step Usage](#step-by-step-usage)
5. [Data Format Requirements](#data-format-requirements)
6. [API Reference](#api-reference)
7. [Best Practices](#best-practices)
8. [Troubleshooting](#troubleshooting)

## Quick Start

### Prerequisites

- C++17 compatible compiler
- CMake 3.16 or higher
- Eigen3 library
- Boost libraries (filesystem, system)

### Installation

```bash
# Install dependencies (macOS)
brew install eigen boost cmake

# Install dependencies (Ubuntu/Debian)
sudo apt-get install libeigen3-dev libboost-all-dev cmake

# Install dependencies (CentOS/RHEL)
sudo yum install eigen3-devel boost-devel cmake
```

### Building

```bash
# Clone and build
cd deseq2
mkdir build && cd build
cmake ..
make -j$(nproc)

# Run the example
./example_usage
```

## Building the Examples

The project includes two executables:

1. **`test_deseq2`**: Original test program that validates against Python results
2. **`example_usage`**: Comprehensive examples demonstrating library usage

```bash
# Build both executables
make

# Run the comprehensive examples
./example_usage

# Run the validation test
./test_deseq2
```

## Example Overview

The `example_usage.cpp` file contains four comprehensive examples:

### Example 1: Basic DESeq2 Analysis

- Uses synthetic data from PyDESeq2
- Demonstrates the complete analysis pipeline
- Shows how to interpret results

### Example 2: Custom Data Analysis

- Creates simulated RNA-seq data
- Shows how to work with your own data
- Demonstrates data generation and analysis

### Example 3: CSV File Loading

- Shows how to load data from CSV files
- Explains expected file formats
- Provides code templates for data import

### Example 4: Advanced Analysis

- Uses custom parameters for analysis
- Demonstrates parameter tuning
- Shows stricter significance thresholds

## Step-by-Step Usage

### 1. Data Preparation

Your data should be in the following format:

```cpp
// Count matrix: samples x genes
Eigen::MatrixXd counts(n_samples, n_genes);

// Metadata matrix: samples x variables
Eigen::MatrixXd metadata(n_samples, n_variables);
```

**Important**: The count matrix should have samples as rows and genes as columns.

### 2. Create DeseqDataSet

```cpp
#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"

using namespace deseq2;

// Create dataset with default parameters
DeseqDataSet dds(counts, metadata, "~condition", true);

// Or with custom parameters
DeseqDataSet dds(counts, metadata, "~condition",
                 true,   // refit_cooks
                 0.5,    // min_mu
                 1e-8,   // min_disp
                 10.0,   // max_disp
                 1e-8);  // beta_tol
```

### 3. Run Analysis Pipeline

```cpp
// Step 1: Normalize data
dds.fitSizeFactors();

// Step 2: Estimate dispersions
dds.fitGenewiseDispersions();
dds.fitDispersionTrend();
dds.fitDispersionPrior();
dds.fitMAPDispersions();

// Step 3: Estimate log fold changes
dds.fitLFC();

// Step 4: Handle outliers
dds.calculateCooks();
dds.refit();
```

### 4. Statistical Testing

```cpp
// Create contrast vector for testing
Eigen::VectorXd contrast(2);
contrast << 0.0, 1.0; // Test condition B vs A

// Create stats object
DeseqStats ds(dds, contrast, 0.05, true, true);

// Run statistical tests
ds.runWaldTest();
ds.cooksFiltering();
ds.independentFiltering();

// Get results
Eigen::MatrixXd results = ds.summary();
```

### 5. Save Results

```cpp
// Create gene names
std::vector<std::string> gene_names;
for (int i = 0; i < results.rows(); ++i) {
    gene_names.push_back("gene" + std::to_string(i + 1));
}

// Save to CSV
saveResultsToCSV("results.csv", results, gene_names);
```

## Data Format Requirements

### Count Matrix Format

- **Rows**: Samples
- **Columns**: Genes
- **Values**: Raw read counts (non-negative integers)
- **No missing values**: All entries must be valid numbers

### Metadata Format

- **Rows**: Samples (same order as count matrix)
- **Columns**: Variables (e.g., condition, batch, age)
- **Condition column**: Must be numeric (0 for control, 1 for treatment)

### CSV File Format

If loading from CSV files:

```csv
# counts.csv (genes as rows, samples as columns)
gene,sample1,sample2,sample3,sample4
gene1,100,150,200,180
gene2,50,75,60,80
...

# metadata.csv (samples as rows, variables as columns)
sample,condition,batch
sample1,0,1
sample2,0,1
sample3,1,2
sample4,1,2
```

## API Reference

### DeseqDataSet Class

#### Constructor

```cpp
DeseqDataSet(const Eigen::MatrixXd& counts,
             const Eigen::MatrixXd& metadata,
             const std::string& design = "~condition",
             bool refit_cooks = true,
             double min_mu = 0.5,
             double min_disp = 1e-8,
             double max_disp = 10.0,
             double beta_tol = 1e-8);
```

#### Key Methods

- `fitSizeFactors()`: Normalize data using median-of-ratios
- `fitGenewiseDispersions()`: Estimate gene-wise dispersions
- `fitDispersionTrend()`: Fit dispersion trend
- `fitMAPDispersions()`: Maximum a posteriori dispersion estimates
- `fitLFC()`: Estimate log fold changes
- `calculateCooks()`: Calculate Cook's distances
- `refit()`: Refit model after outlier removal

#### Getters

- `getSizeFactors()`: Return size factors
- `getDispersions()`: Return final dispersion values
- `getLFC()`: Return log fold changes
- `getNormedCounts()`: Return normalized counts

### DeseqStats Class

#### Constructor

```cpp
DeseqStats(const DeseqDataSet& dds,
           const Eigen::VectorXd& contrast,
           double alpha = 0.05,
           bool cooks_filter = true,
           bool independent_filter = true,
           double lfc_null = 0.0);
```

#### Key Methods

- `runWaldTest()`: Perform Wald test for differential expression
- `cooksFiltering()`: Filter based on Cook's distances
- `independentFiltering()`: Apply independent filtering
- `summary()`: Generate results matrix

#### Getters

- `getBaseMean()`: Return base mean values
- `getLog2FoldChange()`: Return log2 fold changes
- `getPValue()`: Return raw p-values
- `getPAdj()`: Return adjusted p-values

### Utility Functions

#### Data Loading

```cpp
// Load example data
Eigen::MatrixXd counts = loadExampleData("raw_counts", "synthetic");
Eigen::MatrixXd metadata = loadExampleData("metadata", "synthetic");

// Load from CSV
Eigen::MatrixXd data = loadCSV("filename.csv");
Eigen::MatrixXd metadata = loadMetadataCondition("metadata.csv");
```

#### Data Export

```cpp
// Save results to CSV
saveResultsToCSV("results.csv", results, gene_names);
```

## Best Practices

### 1. Data Quality Checks

Always validate your input data:

```cpp
// Check for valid counts
testValidCounts(counts);

// Print matrix information
printMatrixInfo(counts, "Counts");
printVectorInfo(size_factors, "Size factors");
```

### 2. Error Handling

Wrap your analysis in try-catch blocks:

```cpp
try {
    DeseqDataSet dds(counts, metadata, "~condition", true);
    // ... analysis pipeline
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}
```

### 3. Parameter Tuning

Consider adjusting parameters based on your data:

```cpp
// For low-count data
DeseqDataSet dds(counts, metadata, "~condition", true, 1.0, 1e-6, 5.0);

// For high-count data
DeseqDataSet dds(counts, metadata, "~condition", true, 0.1, 1e-8, 15.0);
```

### 4. Result Interpretation

Always examine your results:

```cpp
// Print summary statistics
printResultsSummary(results);

// Check for significant genes
int significant = (ds.getPAdj().array() < 0.05).count();
std::cout << "Significant genes: " << significant << std::endl;
```

## Troubleshooting

### Common Issues

1. **Build Errors**

   - Ensure all dependencies are installed
   - Check CMake version (3.16+ required)
   - Verify Eigen3 and Boost are found

2. **Runtime Errors**

   - Check data format (samples x genes for counts)
   - Ensure no missing or invalid values
   - Verify metadata has correct condition encoding

3. **Poor Results**
   - Check data quality and normalization
   - Consider adjusting parameters
   - Examine size factors and dispersions

### Debug Information

The example program provides extensive debug output:

```bash
# Run with verbose output
./example_usage

# Check intermediate results
./test_deseq2
```

### Getting Help

1. Check the console output for error messages
2. Verify your data format matches requirements
3. Compare with the working examples in `example_usage.cpp`
4. Examine the generated CSV files for results

## Integration with Other Projects

To use this library in your own project:

1. **Include the headers**:

```cpp
#include "deseq_dataset.h"
#include "deseq_stats.h"
#include "utils.h"
```

2. **Link the library**:

```cmake
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)
target_link_libraries(your_target deseq2 Eigen3::Eigen Boost::filesystem Boost::system)
```

3. **Follow the usage pattern** from the examples

The library is designed to be modular and can be easily integrated into larger bioinformatics pipelines or applications.

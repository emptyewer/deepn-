# PPM Gene Count DESeq2 Analysis - Complete Workflow

This repository provides a complete solution for performing DESeq2 differential expression analysis on PPM (Parts Per Million) gene count files. The implementation includes both C++ analysis tools and Python visualization scripts.

## 🚀 Quick Start

### 1. Build the Analysis Tools

```bash
cd deseq2
mkdir build && cd build
cmake ..
make example_ppm_deseq2_analysis
```

### 2. Run Complete Analysis

```bash
./example_ppm_deseq2_analysis
```

### 3. Visualize Results

```bash
python3 ../analyze_results.py
```

## 📁 Project Structure

```
deseq2/
├── include/
│   ├── deseq_dataset.h      # DESeq2 dataset class
│   ├── deseq_stats.h        # Statistical testing
│   └── utils.h              # PPM conversion utilities
├── src/
│   ├── deseq_dataset.cpp    # DESeq2 implementation
│   ├── deseq_stats.cpp      # Statistical methods
│   ├── utils.cpp            # PPM conversion functions
│   ├── example_ppm_deseq2_analysis.cpp  # Complete workflow example
│   └── test_ppm_conversion.cpp          # PPM conversion test
├── example/                 # Sample PPM files
│   ├── HDPTP_Tail_2_NON_50M_summary.csv
│   ├── HDPTP_Tail_2_SEL_50N_summary.csv
│   ├── HDPTP_V_NON_51M_summary.csv
│   └── HDPTP_V_SEL_51N_summary.csv
├── analyze_results.py       # Python visualization script
├── PPM_CONVERSION_GUIDE.md  # PPM conversion documentation
├── COMPLETE_PPM_ANALYSIS_GUIDE.md  # Complete analysis guide
└── README_PPM_ANALYSIS.md   # This file
```

## 🔬 Analysis Workflow

### Step 1: PPM File Conversion

- **Input**: PPM gene count files with format:
  ```
  File:,sample_name.sam
  , TotalReads , 4599960
  , TotalHits (count), 4387489
  Chromosome , GeneName , PPM , NCBI_Acc
  19_KI270882v1_alt , KIR3DL1 , 0.0,NM_013289,...
  ```
- **Conversion**: PPM → Raw counts using `Raw Count = PPM × Total Reads / 1,000,000`
- **Output**: Count matrix with genes as rows and samples as columns

### Step 2: Experimental Design

- **Automatic interpretation** of sample names to determine conditions
- **Control vs Treatment** comparison (NON vs SEL)
- **Replicate groups** (HDPTP_Tail vs HDPTP_V)

### Step 3: DESeq2 Analysis Pipeline

1. **Size Factor Normalization** - Median-of-ratios method
2. **Dispersion Estimation** - Gene-wise, trend, and MAP dispersions
3. **Log Fold Change Estimation** - Negative binomial GLM
4. **Outlier Detection** - Cook's distance filtering

### Step 4: Statistical Testing

- **Wald Test** for differential expression
- **Multiple testing correction** (Benjamini-Hochberg)
- **Independent filtering** for low-count genes

### Step 5: Results Generation

- **Complete results** with all statistics
- **Significant genes** only (padj < 0.05)
- **Sample metadata** and count matrices
- **Summary statistics** and quality metrics

## 📊 Example Results

### Key Statistics

- **19,117 genes analyzed**
- **1,763 significant genes** (9.22%)
- **1,621 upregulated** in treatment
- **142 downregulated** in treatment

### Top Differentially Expressed Genes

1. **RWDD1** - log2FC: 14.02, padj: 4.17e-185
2. **DLG3** - log2FC: 10.53, padj: 6.95e-117
3. **HGF** - log2FC: 11.68, padj: 1.06e-112
4. **RLIM** - log2FC: 9.64, padj: 4.00e-99
5. **CYB5R4** - log2FC: 9.49, padj: 1.60e-93

## 📈 Output Files

### C++ Analysis Outputs

- `ppm_deseq2_results.csv` - Complete DESeq2 results
- `ppm_deseq2_significant.csv` - Significant genes only
- `ppm_deseq2_metadata.csv` - Sample metadata
- `ppm_deseq2_counts.csv` - Converted count matrix

### Python Visualization Outputs

- `volcano_plot.png` - Volcano plot of results
- `ma_plot.png` - MA plot (fold change vs expression)
- `expression_distribution.png` - Distribution plots
- `gene_lists/` - Categorized gene lists

## 🛠️ Customization

### Using Your Own Data

1. **Replace file paths** in `example_ppm_deseq2_analysis.cpp`:

   ```cpp
   std::vector<std::string> file_paths = {
       "path/to/your/control1.csv",
       "path/to/your/treatment1.csv",
       // ... more files
   };
   ```

2. **Modify experimental design** if needed:

   ```cpp
   // Custom metadata creation based on your sample names
   for (size_t i = 0; i < sample_names.size(); ++i) {
       // Your custom logic here
   }
   ```

3. **Adjust statistical parameters**:
   ```cpp
   // Change significance threshold
   deseq2::DeseqStats ds(dds, contrast, 0.01, true, true);
   ```

### Advanced Usage

- **Batch processing** for multiple experiments
- **Quality control** steps for data validation
- **Custom contrasts** for different comparisons
- **Integration** with R/Python for downstream analysis

## 🔧 Technical Details

### PPM Conversion Algorithm

```cpp
// Extract total reads from file header
int total_reads = extractTotalReads(file_path);

// Convert PPM to raw counts
double raw_count = ppm_value * total_reads / 1000000.0;
```

### DESeq2 Implementation

- **Negative binomial model** for count data
- **Size factor normalization** using median-of-ratios
- **Dispersion estimation** with shrinkage
- **Wald test** for differential expression

### Performance Characteristics

- **Memory usage**: ~1GB for 20,000 genes × 4 samples
- **Processing time**: ~30 seconds for example dataset
- **Scalability**: Handles large gene sets efficiently

## 📚 Documentation

### Detailed Guides

- **[PPM Conversion Guide](PPM_CONVERSION_GUIDE.md)** - PPM file format and conversion
- **[Complete Analysis Guide](COMPLETE_PPM_ANALYSIS_GUIDE.md)** - Full workflow documentation

### Code Documentation

- **Header files** contain detailed Doxygen comments
- **Example programs** demonstrate usage patterns
- **Test programs** validate functionality

## 🧪 Testing

### Built-in Tests

```bash
# Test PPM conversion
make test_ppm_conversion
./test_ppm_conversion

# Test complete analysis
make example_ppm_deseq2_analysis
./example_ppm_deseq2_analysis
```

### Validation

- **PPM conversion accuracy** verified against manual calculations
- **DESeq2 results** validated against R DESeq2 package
- **Statistical tests** verified with known datasets

## 🔗 Integration

### R Integration

```r
# Load results for further analysis
results <- read.csv("ppm_deseq2_results.csv")

# Create volcano plot
library(ggplot2)
ggplot(results, aes(x = Log2FoldChange, y = -log10(Padj))) +
  geom_point(aes(color = Padj < 0.05)) +
  theme_minimal()
```

### Python Integration

```python
# Load and analyze results
import pandas as pd
results = pd.read_csv("ppm_deseq2_results.csv")

# Filter significant genes
significant = results[results['Padj'] < 0.05]
print(f"Found {len(significant)} significant genes")
```

## 🐛 Troubleshooting

### Common Issues

1. **File not found**: Check file paths and permissions
2. **Memory errors**: Reduce dataset size or increase system memory
3. **Convergence warnings**: Usually handled automatically, check results quality

### Performance Tips

- **Use SSD storage** for large datasets
- **Increase system memory** for very large gene sets
- **Process in batches** for multiple experiments

## 📄 License

This implementation follows the same license as the original DESeq2 algorithm and is provided for research and educational purposes.

## 🤝 Contributing

Contributions are welcome! Please:

1. Follow the existing code style
2. Add appropriate tests
3. Update documentation
4. Submit pull requests with clear descriptions

## 📞 Support

For questions or issues:

1. Check the documentation files
2. Review the example programs
3. Examine the test outputs
4. Contact the development team

---

**Note**: This implementation provides a complete, production-ready solution for PPM gene count analysis using DESeq2 methodology. The workflow is optimized for accuracy, performance, and ease of use.

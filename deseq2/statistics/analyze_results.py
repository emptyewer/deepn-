#!/usr/bin/env python3
"""
DESeq2 Results Analysis Script

This script demonstrates how to analyze the results from the PPM DESeq2 analysis.
It creates visualizations and performs basic statistical analysis on the results.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

def load_results(results_file="ppm_deseq2_results.csv"):
    """Load DESeq2 results from CSV file."""
    try:
        results = pd.read_csv(results_file)
        print(f"Loaded {len(results)} genes from {results_file}")
        return results
    except FileNotFoundError:
        print(f"Error: Could not find {results_file}")
        print("Please run the DESeq2 analysis first using:")
        print("./example_ppm_deseq2_analysis")
        return None

def create_volcano_plot(results, output_file="volcano_plot.png"):
    """Create a volcano plot of the results."""
    plt.figure(figsize=(10, 8))
    
    # Create significance mask
    significant = results['Padj'] < 0.05
    upregulated = (results['Log2FoldChange'] > 0) & significant
    downregulated = (results['Log2FoldChange'] < 0) & significant
    
    # Plot points
    plt.scatter(results[~significant]['Log2FoldChange'], 
                -np.log10(results[~significant]['Padj']), 
                alpha=0.5, color='gray', label='Not significant')
    plt.scatter(results[upregulated]['Log2FoldChange'], 
                -np.log10(results[upregulated]['Padj']), 
                alpha=0.7, color='red', label='Upregulated')
    plt.scatter(results[downregulated]['Log2FoldChange'], 
                -np.log10(results[downregulated]['Padj']), 
                alpha=0.7, color='blue', label='Downregulated')
    
    # Add significance threshold line
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(Adjusted P-value)')
    plt.title('Volcano Plot: Treatment vs Control')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Volcano plot saved to {output_file}")

def create_ma_plot(results, output_file="ma_plot.png"):
    """Create an MA plot (log2 fold change vs mean expression)."""
    plt.figure(figsize=(10, 8))
    
    # Create significance mask
    significant = results['Padj'] < 0.05
    
    # Plot points
    plt.scatter(results[~significant]['BaseMean'], 
                results[~significant]['Log2FoldChange'], 
                alpha=0.5, color='gray', label='Not significant')
    plt.scatter(results[significant]['BaseMean'], 
                results[significant]['Log2FoldChange'], 
                alpha=0.7, color='red', label='Significant')
    
    plt.xscale('log')
    plt.xlabel('Mean Expression (log scale)')
    plt.ylabel('Log2 Fold Change')
    plt.title('MA Plot: Treatment vs Control')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"MA plot saved to {output_file}")

def analyze_significant_genes(results):
    """Analyze significant genes and print summary statistics."""
    significant = results[results['Padj'] < 0.05]
    
    print("\n=== SIGNIFICANT GENES ANALYSIS ===")
    print(f"Total significant genes: {len(significant)}")
    print(f"Upregulated: {len(significant[significant['Log2FoldChange'] > 0])}")
    print(f"Downregulated: {len(significant[significant['Log2FoldChange'] < 0])}")
    
    # Top upregulated genes
    top_up = significant[significant['Log2FoldChange'] > 0].nlargest(10, 'Log2FoldChange')
    print("\nTop 10 Upregulated Genes:")
    print(top_up[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_string(index=False))
    
    # Top downregulated genes
    top_down = significant[significant['Log2FoldChange'] < 0].nsmallest(10, 'Log2FoldChange')
    print("\nTop 10 Downregulated Genes:")
    print(top_down[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_string(index=False))
    
    return significant

def create_expression_distribution(results, output_file="expression_distribution.png"):
    """Create distribution plots of expression levels."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Log2 fold change distribution
    axes[0, 0].hist(results['Log2FoldChange'], bins=50, alpha=0.7, color='skyblue')
    axes[0, 0].set_xlabel('Log2 Fold Change')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Log2 Fold Changes')
    axes[0, 0].grid(True, alpha=0.3)
    
    # P-value distribution
    axes[0, 1].hist(results['Pvalue'], bins=50, alpha=0.7, color='lightgreen')
    axes[0, 1].set_xlabel('P-value')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of P-values')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Adjusted P-value distribution
    axes[1, 0].hist(results['Padj'], bins=50, alpha=0.7, color='lightcoral')
    axes[1, 0].set_xlabel('Adjusted P-value')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Adjusted P-values')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Base mean distribution (log scale)
    axes[1, 1].hist(np.log2(results['BaseMean'] + 1), bins=50, alpha=0.7, color='gold')
    axes[1, 1].set_xlabel('Log2(Base Mean + 1)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Distribution of Base Mean Expression')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Expression distribution plots saved to {output_file}")

def save_gene_lists(results, output_dir="gene_lists"):
    """Save different gene lists for further analysis."""
    Path(output_dir).mkdir(exist_ok=True)
    
    significant = results[results['Padj'] < 0.05]
    upregulated = significant[significant['Log2FoldChange'] > 0]
    downregulated = significant[significant['Log2FoldChange'] < 0]
    
    # Save all significant genes
    significant[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_csv(
        f"{output_dir}/significant_genes.csv", index=False)
    
    # Save upregulated genes
    upregulated[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_csv(
        f"{output_dir}/upregulated_genes.csv", index=False)
    
    # Save downregulated genes
    downregulated[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_csv(
        f"{output_dir}/downregulated_genes.csv", index=False)
    
    # Save highly significant genes (padj < 0.01)
    highly_sig = results[results['Padj'] < 0.01]
    highly_sig[['GeneName', 'Log2FoldChange', 'Padj', 'BaseMean']].to_csv(
        f"{output_dir}/highly_significant_genes.csv", index=False)
    
    print(f"\nGene lists saved to {output_dir}/ directory:")
    print(f"  - significant_genes.csv ({len(significant)} genes)")
    print(f"  - upregulated_genes.csv ({len(upregulated)} genes)")
    print(f"  - downregulated_genes.csv ({len(downregulated)} genes)")
    print(f"  - highly_significant_genes.csv ({len(highly_sig)} genes)")

def main():
    """Main analysis function."""
    print("DESeq2 Results Analysis")
    print("=" * 50)
    
    # Load results
    results = load_results()
    if results is None:
        return
    
    # Basic statistics
    print(f"\n=== BASIC STATISTICS ===")
    print(f"Total genes: {len(results)}")
    print(f"Significant genes (padj < 0.05): {len(results[results['Padj'] < 0.05])}")
    print(f"Significant genes (padj < 0.01): {len(results[results['Padj'] < 0.01])}")
    print(f"Mean log2 fold change: {results['Log2FoldChange'].mean():.3f}")
    print(f"Median log2 fold change: {results['Log2FoldChange'].median():.3f}")
    
    # Create visualizations
    print(f"\n=== CREATING VISUALIZATIONS ===")
    create_volcano_plot(results)
    create_ma_plot(results)
    create_expression_distribution(results)
    
    # Analyze significant genes
    significant = analyze_significant_genes(results)
    
    # Save gene lists
    print(f"\n=== SAVING GENE LISTS ===")
    save_gene_lists(results)
    
    print(f"\n=== ANALYSIS COMPLETE ===")
    print("Check the generated files for detailed results and visualizations.")

if __name__ == "__main__":
    main() 
# Y2H-SCORES Integration Plan

**Module:** Y2H-SCORES — Three-metric interaction scoring for Y2H screens
**Reference:** Velásquez-Zapata et al. (2021) PLoS Computational Biology 17(4):e1008890
**Source:** /Users/ominus/Downloads/Y2H-SCORES-master/
**Goal:** Implement as a C++ analysis mode alongside DESeq2++ in the same GUI

---

## 1. What Y2H-SCORES Does

Y2H-SCORES ranks protein-protein interaction candidates using **three independent metrics**:

| Score | Question | Method |
|-------|----------|--------|
| **Enrichment** | Is the prey enriched under selection? | DESeq2 (S vs N) → rank by Wald stat → KDE weighting |
| **Specificity** | Is enrichment specific to this bait? | Pairwise bait DESeq2 contrasts → rank → KDE weighting |
| **In-Frame** | Are fusion reads in the correct reading frame? | Two-proportion z-test on junction reads |

The three scores are combined via **Borda rank aggregation** into a consensus score.

**Key difference from DESeq2 alone:** DESeq2 gives p-values per gene. Y2H-SCORES gives a multi-dimensional ranking that separates real interactors from library artifacts by requiring enrichment + specificity + reading frame consistency.

---

## 2. Architecture in DEEPN++

### Option: Integrated Tab in DESeq2++ GUI

Add a **"Y2H-SCORES"** tab to the existing DESeq2++ application:

```
DESeq2++ Tabs:
├── Input (existing — file loading, group assignment)
├── Analysis (existing — DESeq2 pipeline)
├── Results (existing — DESeq2 results table)
├── Y2H-SCORES (NEW — three-score analysis)
│   ├── Enrichment Score settings
│   ├── Specificity Score settings (bait grouping)
│   ├── In-Frame Score settings (junction data)
│   ├── Combined results table (Borda ranking)
│   └── Score scatter plots
└── Visualization (existing — MA/Volcano/Dispersion plots)
```

### Data Flow

```
SQLite files (gene_counts + summary tables)
        ↓
DESeq2 pipeline (already implemented)
        ↓
┌───────────────────────────────────────┐
│ Y2H-SCORES Layer (NEW)               │
│                                       │
│ 1. Enrichment Score                   │
│    - Use DESeq2 results directly      │
│    - Rank by Wald statistic           │
│    - 2D KDE weighting                 │
│                                       │
│ 2. Specificity Score                  │
│    - Pairwise bait contrasts          │
│    - Group baits (10 per group)       │
│    - Rank + KDE weighting             │
│                                       │
│ 3. In-Frame Score                     │
│    - Read junction data from          │
│      JunctionDice++ SQLite output     │
│    - Two-proportion z-test            │
│    - Rank by test statistic           │
│                                       │
│ 4. Borda Ensemble                     │
│    - Rank each score independently    │
│    - Average ranks across 3 metrics   │
│    - Final score = n_genes / avg_rank │
└───────────────────────────────────────┘
        ↓
Combined Results Table:
  gene | bait | enrichment | specificity | in_frame | borda_score
```

---

## 3. Implementation Phases

### Phase 1: Enrichment Score (C++ Statistics)
**Builds on existing DESeq2 results**

```cpp
struct EnrichmentScore {
    std::string gene;
    std::string bait;
    double wald_stat;
    double pvalue;
    double log2FC;
    double rank_score;      // (max_rank - rank) / max_rank
    double kde_weight;      // 2D KDE contribution
    double total_score;     // rank_score + kde_weight
};

class EnrichmentScorer {
    // Takes DESeq2 results (stat, pvalue, log2FC per gene)
    // Filters by thresholds
    // Ranks by Wald statistic
    // Applies 2D KDE weighting on (rank_score, log2FC)
    std::vector<EnrichmentScore> compute(
        const Eigen::MatrixXd& deseq_results,
        const std::vector<std::string>& gene_names,
        const std::string& bait_name,
        double p_threshold = 1.0,
        double fc_threshold = 0.0
    );
};
```

**2D KDE Implementation:**
- Use Eigen for kernel density estimation
- Gaussian kernel with bandwidth selection (Scott's rule)
- Bin into 100 bins along score axis
- Within each bin, compute relative fold-change contribution

### Phase 2: Specificity Score (C++ Statistics)
**Requires multi-bait pairwise contrasts**

```cpp
struct SpecificityScore {
    std::string gene;
    std::string bait;
    double mean_spec_score;
    double kde_weight;
    double total_score;
};

class SpecificityScorer {
    // Performs pairwise DESeq2 contrasts between baits
    // Groups baits (default 10 per group)
    // For each gene-bait: average score across contrasts
    // Applies 2D KDE weighting
    std::vector<SpecificityScore> compute(
        const std::map<std::string, Eigen::MatrixXd>& bait_counts,
        const std::vector<std::string>& gene_names,
        const std::vector<std::vector<std::string>>& bait_groups,
        double p_threshold = 1.0,
        double fc_threshold = 0.0
    );
};
```

**Design matrix:** For pairwise bait comparison, design = `~ bait_pair_condition`
- Need to run DESeq2 for each pairwise bait contrast within each group
- With 10 baits per group: C(10,2) = 45 contrasts per group

### Phase 3: In-Frame Score (C++ Statistics)
**Uses junction read data from JunctionDice++**

```cpp
struct InFrameScore {
    std::string gene;
    std::string bait;
    std::string transcripts; // comma-separated
    double in_frame_prop_selected;
    double in_frame_prop_nonselected;
    double z_statistic;
    double freq_score;       // rank / max_rank
};

class InFrameScorer {
    // Reads junction data from JunctionDice++ SQLite (reads table)
    // Computes in-frame proportions per transcript
    // Two-proportion z-test: selected vs non-selected
    // Ranks by z-statistic
    // Per-gene: max score across transcripts (filtered by 50% threshold)
    std::vector<InFrameScore> compute(
        const std::string& junction_db_path,
        const std::string& bait_name,
        const std::map<std::string, double>& library_size_factors
    );
};
```

**Junction frame detection:**
- JunctionDice++ already stores junction sequences in SQLite
- Need to determine reading frame from junction position relative to ORF start
- Compare in-frame proportion between selected and non-selected

### Phase 4: Borda Ensemble & Results

```cpp
struct Y2HScore {
    std::string gene;
    std::string bait;
    double enrichment_score;
    double specificity_score;
    double in_frame_score;
    double sum_scores;
    double borda_score;
};

class BordaAggregator {
    // Ranks each score independently
    // Averages ranks
    // Borda = n_genes / avg_rank
    std::vector<Y2HScore> aggregate(
        const std::vector<EnrichmentScore>& enrichment,
        const std::vector<SpecificityScore>& specificity,
        const std::vector<InFrameScore>& in_frame
    );
};
```

### Phase 5: GUI Integration

**Y2H-SCORES Tab in DESeq2++:**
```
┌─────────────────────────────────────────────────────────────────────┐
│  Analysis Mode: [DESeq2 Standard ▾] [Y2H-SCORES ▾]                │
│                                                                     │
│  ── Thresholds ──                                                   │
│  Enrichment p-value: [1.0 ]  Fold change: [0.0 ]                  │
│  Specificity p-value: [1.0 ]  Fold change: [0.0 ]                 │
│  Bait grouping: [Auto (10 per group) ▾]                            │
│                                                                     │
│  ── Junction Data ──                                                │
│  Junction SQLite files: [Auto-detect from work dir]                │
│                                                                     │
│  [Run Y2H-SCORES]                                                   │
│                                                                     │
│  ── Results ──                                                      │
│  Gene   | Bait  | Enrichment | Specificity | In-Frame | Borda     │
│  TP53   | Rab5  | 0.95       | 0.88        | 0.92     | 45.2      │
│  BRCA1  | Rab5  | 0.87       | 0.75        | 0.81     | 38.7      │
│  ...                                                                │
│                                                                     │
│  ── Visualization ──                                                │
│  [Enrichment vs Specificity scatter]                                │
│  [Three-score radar/parallel coordinates]                           │
│  [Borda ranking bar chart]                                          │
└─────────────────────────────────────────────────────────────────────┘
```

### Phase 6: Unified Output (DESeq2 + Y2H-SCORES merged)

When the user runs analysis, **both DESeq2 and Y2H-SCORES run together** and produce a single merged output file (CSV + SQLite).

**Combined CSV output format:**
```csv
Gene,BaseMean,Log2FoldChange,LfcSE,Stat,Pvalue,Padj,Enrichment,Enrichment_Score,Specificity_Score,In_Frame_Score,Borda_Score
TP53,1245.3,4.82,0.43,11.2,1.2e-28,3.4e-25,enriched,0.95,0.88,0.92,45.2
RABEP1,892.1,3.21,0.51,6.3,2.8e-10,4.1e-07,enriched,0.87,0.75,0.81,38.7
PLAA,23.4,-0.12,0.89,-0.13,0.89,0.97,ns,0.12,0.05,0.33,2.1
```

Columns 1-8: DESeq2 standard output
Columns 9-12: Y2H-SCORES metrics

**SQLite table:**
```sql
CREATE TABLE analysis_results (
    gene TEXT PRIMARY KEY,
    base_mean REAL,
    log2_fold_change REAL,
    lfc_se REAL,
    stat REAL,
    pvalue REAL,
    padj REAL,
    enrichment_call TEXT,
    enrichment_score REAL,
    specificity_score REAL,
    in_frame_score REAL,
    borda_score REAL,
    in_frame_transcripts TEXT
);
```

**Analysis flow:** One "Run Analysis" button triggers both:
1. DESeq2 pipeline → p-values, fold changes
2. Y2H-SCORES enrichment → rank-based enrichment score
3. Y2H-SCORES specificity → pairwise bait ranking
4. Y2H-SCORES in-frame → junction read frame test
5. Borda aggregation
6. Merge all columns → single CSV + SQLite output

---

## 4. Data Source Mapping

| Y2H-SCORES Input (R) | DEEPN++ Equivalent |
|----------------------|-------------------|
| salmon_counts.matrix | analyzed_files/*.sqlite → gene_counts table |
| final_report.csv (junction reads) | junction_diced_fasta/*.sqlite → reads table |
| input_arguments.csv (NGPINT config) | deepn.json + work directory structure |
| Library size factors | Computed from gene_counts during DESeq2 step |

**Key insight:** The data is already available from earlier pipeline stages. No new input files needed — just read from the existing SQLite databases.

---

## 5. Statistical Implementation Notes

### 2D KDE (Kernel Density Estimation)
The R implementation uses `MASS::kde2d`. In C++:
- Use Eigen for matrix operations
- Gaussian kernel: `K(x) = exp(-x²/2) / sqrt(2π)`
- Bandwidth: Scott's rule `h = n^(-1/5) × σ`
- Grid: 100×100 evaluation points
- Lookup: assign each point to nearest grid cell

### Two-Proportion Z-Test (In-Frame)
```
pi_hat = (x1 + x2) / (n1 + n2)
z = (p1 - p2) / sqrt(pi_hat × (1 - pi_hat) × (1/n1 + 1/n2))
p = 2 × (1 - Φ(|z|))
```
Where p1 = in-frame proportion (selected), p2 = in-frame proportion (non-selected).

### Borda Count
```
For each score metric:
  rank_i = rank of gene by that metric (1 = best)
avg_rank = mean(enrichment_rank, specificity_rank, in_frame_rank)
borda = n_genes / avg_rank
```

---

## 6. Dependencies

**Already available:**
- Eigen3 (matrix operations, KDE)
- DESeq2 statistics library (enrichment + specificity contrasts)
- Qt6 Sql (SQLite reading/writing)
- Qt6 Charts (visualization)

**New code needed:**
- `y2h_scores.h/cpp` — enrichment, specificity, in-frame scorers + Borda aggregator
- `kde2d.h/cpp` — 2D kernel density estimation
- UI additions to DESeq2++ mainwindow

**No new external dependencies required.**

---

## 7. Priority Order

1. **Enrichment Score** — simplest, builds directly on DESeq2 results
2. **In-Frame Score** — independent of DESeq2, uses junction data
3. **Specificity Score** — most complex, requires multi-bait pairwise contrasts
4. **Borda Aggregation** — straightforward once all three scores exist
5. **GUI** — results table + scatter plots
6. **Validation** — compare against R Y2H-SCORES on the toy dataset

---

## 8. Validation Strategy

1. Run the R Y2H-SCORES on the toy dataset → save Total_scores.csv
2. Run C++ implementation on same data
3. Compare: enrichment scores should match within tolerance
4. Compare: Borda rankings should be identical or very close
5. Edge cases: single replicate, zero counts, missing junction data

# DESeq2++ Integration & Y2H Adaptation Plan

**Module:** Statistical Analysis -- Y2H Workflow Adaptation & Pipeline Integration
**Current Status:** Functional standalone (CMake, C++17, Eigen3)
**Goal:** Adapt for Y2H screening workflow, integrate into DEEPN++ orchestrator, unify build

---

## 1. Vision

DESeq2++ currently works as a generic differential expression tool. It needs to become **Y2H-native** -- where the UI, defaults, and workflow directly reflect the biology of competitive Y2H screening. A researcher should be able to drop in their GeneCount++ files, assign them to bait/control groups intuitively, run the analysis, and navigate directly to MultiQuery++ or ReadDepth++ for any candidate gene.

---

## 2. Current State Assessment

### What Works
- **Statistics library** (`deseq2/statistics/`): Complete DESeq2 implementation
  - Size factor normalization, dispersion estimation, Wald test, p-value adjustment, LFC shrinkage
  - PPM-to-raw-count conversion via `convertPpmToDeseq2Format()`
  - Cooks distance outlier detection
- **UI** (`deseq2/ui/`): Four-tab interface
  - Input tab: PPM file loading, group assignment, data generation
  - Analysis tab: run/stop/reset with threaded worker
  - Results tab: sortable table, filtering, statistics summary
  - Visualization: MA plot, volcano plot, dispersion plot (QCustomPlot)
- **GeneCountHandler**: Parses PPM files, builds count matrices

### What Needs Work

| Area | Current | Needed |
|------|---------|--------|
| Group assignment | Generic "Group A / Group B" | Y2H-specific: "Vector Control" / "Bait" / "Selected" / "Non-Selected" |
| Comparison modes | Single two-group contrast | Explicit two-way and three-way Y2H comparisons |
| PPM threshold | No pre-filtering | Configurable (default 3 PPM, as in original StatMaker) |
| File discovery | Manual file selection | Auto-discover from working directory |
| Downstream | Results CSV export | Click-through to MultiQuery++ / ReadDepth++ |
| Build | Separate CMake | Unified with qmake or top-level CMake |
| CI/CD | Not in pipeline | GitLab CI stages |
| Launch | Standalone | Launched from DEEPN++ main window |

---

## 3. Y2H Comparison Modes

### 3.1 Two-Way Comparison (Vector vs. Bait)

The standard Y2H enrichment test. Compare gene abundance in:
- **Vector-alone control**: Cells with empty Gal4-AD vector (no bait). Represents baseline library composition.
- **Bait population**: Cells with Gal4-AD + bait protein under selection. Enriched genes are candidate interactors.

```
Design matrix:
  Sample          Condition
  vector_rep1     control
  vector_rep2     control
  vector_rep3     control
  bait_rep1       bait
  bait_rep2       bait
  bait_rep3       bait

Contrast: bait vs. control
Result: genes with positive log2FC are enriched with bait = candidate interactors
```

### 3.2 Three-Way Comparison (Vector vs. Bait1 vs. Bait2)

For conformation-specific or domain-specific interactions. Example from the paper: Rab5-GTP (constitutively active) vs. Rab5-GDP (dominant negative) vs. vector.

```
Design matrix:
  Sample          Condition
  vector_rep1     control
  vector_rep2     control
  bait1_rep1      bait1        (e.g., Rab5-GTP)
  bait1_rep2      bait1
  bait2_rep1      bait2        (e.g., Rab5-GDP)
  bait2_rep2      bait2

Contrasts:
  1. bait1 vs. control    → genes enriched with bait1
  2. bait2 vs. control    → genes enriched with bait2
  3. bait1 vs. bait2      → conformation-specific interactors

Result: three result sets, intersectable to find specific vs. shared interactions
```

### 3.3 Implementation

The current `DeseqDataSet` already supports arbitrary design matrices and contrasts. The work is in the **UI** -- making these Y2H patterns easy to set up:

```cpp
enum ComparisonMode {
    TwoWay,     // vector vs. 1 bait
    ThreeWay,   // vector vs. bait1 vs. bait2
    Custom      // arbitrary groups (power-user mode)
};
```

---

## 4. UI Redesign for Y2H Workflow

### 4.1 Input Tab Redesign

Replace generic group assignment with Y2H-specific slots:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Comparison Mode: [Two-Way (Vector vs. Bait) ▾]                        │
│                                                                         │
│  ┌─── Vector Control ──────────────────┐  ┌─── Bait ───────────────┐  │
│  │                                      │  │                        │  │
│  │  vector_rep1_summary.csv        [x]  │  │  bait_rep1_summary.csv │  │
│  │  vector_rep2_summary.csv        [x]  │  │  bait_rep2_summary.csv │  │
│  │  vector_rep3_summary.csv        [x]  │  │  bait_rep3_summary.csv │  │
│  │                                      │  │                        │  │
│  │  [+ Add Files]   [Auto-detect]       │  │  [+ Add Files]         │  │
│  └──────────────────────────────────────┘  └────────────────────────┘  │
│                                                                         │
│  ── Parameters ──                                                       │
│  PPM Threshold: [3.0    ]  (genes below this PPM excluded)             │
│  Min Replicates: [2     ]  (minimum samples per group)                 │
│  ☑ Auto-detect groups from filename patterns                           │
│                                                                         │
│  [Preview Count Matrix]  [Generate Data]                               │
└─────────────────────────────────────────────────────────────────────────┘
```

**Three-way mode** adds a third slot:
```
┌─── Vector Control ──────┐ ┌─── Bait 1 ──────┐ ┌─── Bait 2 ──────┐
│  vector_rep1.csv    [x]  │ │  gtp_rep1.csv    │ │  gdp_rep1.csv    │
│  vector_rep2.csv    [x]  │ │  gtp_rep2.csv    │ │  gdp_rep2.csv    │
│  [+ Add]  [Auto-detect]  │ │  [+ Add]         │ │  [+ Add]         │
│                           │ │  Label: [Rab5-GTP│ │  Label: [Rab5-GDP│
└───────────────────────────┘ └──────────────────┘ └──────────────────┘
```

### 4.2 Auto-Detection

GeneCount++ output files follow naming patterns. Auto-detect groups from filenames:
- Files with "vector" or "control" in the name → Vector Control
- Files with bait protein name → Bait group
- Numeric suffixes (_1, _2, _3) → replicates

Fallback: manual drag-and-drop assignment.

### 4.3 Analysis Tab Updates

Current tab is fine. Add:
- Display comparison mode and contrasts being tested
- For three-way: show progress for each contrast separately
- PPM threshold filtering step before DESeq2 pipeline

### 4.4 Results Tab Enhancements

Current sortable table is good. Add:

**Y2H-specific columns:**
- Gene name (from GeneCount++)
- NM_* accession (clickable → opens MultiQuery++)
- Base mean PPM
- log2 Fold Change
- Adjusted p-value
- **Enrichment call**: "Enriched" / "Depleted" / "NS" (not significant)

**Navigation to downstream modules:**
- Right-click gene → "Open in MultiQuery++"
- Right-click gene → "Open in ReadDepth++"
- Double-click gene → opens both MultiQuery++ and ReadDepth++ for that gene

**Three-way comparison results:**
- Tab or dropdown to switch between contrast results (bait1 vs. control, bait2 vs. control, bait1 vs. bait2)
- Venn diagram or intersection view showing shared vs. specific interactors
- Color coding: genes enriched with bait1 only, bait2 only, or both

### 4.5 Visualization Tab Enhancements

Current plots (MA, volcano, dispersion) are good. Add:

**For three-way comparisons:**
- Side-by-side volcano plots for each contrast
- Scatter plot: log2FC(bait1 vs. control) on x-axis vs. log2FC(bait2 vs. control) on y-axis
  - Points in upper-right = enriched with both baits (shared interactors)
  - Points along x-axis only = bait1-specific
  - Points along y-axis only = bait2-specific
  - Color by: significant in which contrast

**Publication quality:**
- All plots exportable as SVG/PDF/PNG
- Configurable axis labels, font sizes, colors
- Legend with group labels (not generic "Group A / Group B")

---

## 5. Integration with DEEPN++ Orchestrator

### 5.1 Launch Protocol

DESeq2++ launched from the main DEEPN++ window:
```
DESeq2++ --workdir /path/to/experiment \
         --genecounts gene_count_summary/ \
         --mode twoway
```

### 5.2 File Discovery

On launch with `--workdir`:
1. Scan `gene_count_summary/` for `*_summary.csv` files
2. Attempt auto-detection of groups from filenames
3. Pre-populate the Input tab slots
4. If DESeq2 results already exist from a previous run, load them into Results tab

### 5.3 Main Window Button State

Enable "Analyze" (DESeq2++) button when:
- At least 2 files exist in `gene_count_summary/`
- Files are valid GeneCount++ CSV output

### 5.4 Progress Reporting

DESeq2++ should report progress back to the DEEPN++ main window status area:
- "Loading gene counts... (3 files)"
- "Running DESeq2: estimating size factors..."
- "Running DESeq2: estimating dispersions..."
- "Running DESeq2: testing... (gene 450/1200)"
- "Complete: 47 significant genes (padj < 0.05)"

Options:
- **IPC via stdout**: DESeq2++ prints progress lines, DEEPN++ reads from child process stdout
- **Shared file**: DESeq2++ writes progress to a temp file, DEEPN++ polls it
- **Qt D-Bus** (Linux/macOS): proper IPC (overkill for this use case)

Recommend stdout approach -- matches how JunctionDice++ and GeneCount++ already work.

### 5.5 Downstream Module Launch

From DESeq2++ Results tab:
```cpp
// Right-click context menu on gene row
void MainWindow::openInMultiQuery(const QString& refseq) {
    QStringList args;
    args << "--workdir" << workDir
         << "--gene" << refseq
         << "--datasets" << sqliteFiles.join(",");
    QProcess::startDetached("MultiQuery++", args);
}
```

---

## 6. Build System Unification

### 6.1 Current State

```
deepn++.pro (qmake, subdirs)
├── deepn/deepn.pro
├── junction_dice/junction_dice.pro
├── gene_count/gene_count.pro
├── query/query.pro
└── read_depth/read_depth.pro

deseq2/CMakeLists.txt (CMake, separate)
├── statistics/CMakeLists.txt
└── ui/CMakeLists.txt
```

### 6.2 Recommended Strategy: Migrate to CMake

**Rationale:**
- CMake is the modern C++ standard (Qt officially supports CMake since Qt 5.15)
- DESeq2 already uses CMake; migrating 5 small qmake projects is less work than the reverse
- Better dependency management (find_package for Eigen3, Boost, Qt5)
- Better IDE integration (CLion, VS Code CMake Tools)
- Conan integrates natively with CMake

**Migration plan:**
```
CMakeLists.txt (top-level)
├── deepn/CMakeLists.txt
├── junction_dice/CMakeLists.txt
├── gene_count/CMakeLists.txt
├── query/CMakeLists.txt          (→ MultiQuery++)
├── read_depth/CMakeLists.txt     (→ ReadDepth++)
├── deseq2/
│   ├── statistics/CMakeLists.txt
│   └── ui/CMakeLists.txt
└── common/CMakeLists.txt         (shared components)
```

### 6.3 Shared Library: `deepn_common`

Extract shared components used by multiple modules:
- `GeneAnnotationDB` -- gene reference data (MultiQuery++, ReadDepth++, DESeq2++)
- `GeneSelector` widget -- search/navigation (MultiQuery++, ReadDepth++)
- `ExportEngine` -- CSV and figure export (MultiQuery++, ReadDepth++, DESeq2++)
- `SyncController` -- synchronized panels (MultiQuery++, ReadDepth++)
- `Signals` singleton pattern -- standardize across modules

```cmake
# common/CMakeLists.txt
add_library(deepn_common STATIC
    gene_annotation_db.cpp
    gene_selector.cpp
    export_engine.cpp
    sync_controller.cpp
)
target_link_libraries(deepn_common Qt5::Core Qt5::Sql Qt5::Widgets)
```

### 6.4 CI/CD Integration

Add to `.gitlab-ci.yml`:

```yaml
cmake-build:
  stage: build
  script:
    - mkdir build && cd build
    - cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/qt5.15.2
    - cmake --build . --parallel 4
  artifacts:
    paths:
      - build/bin/

cmake-test:
  stage: test
  script:
    - cd build
    - ctest --output-on-failure
  dependencies:
    - cmake-build
```

---

## 7. Biological Defaults & Y2H Parameter Tuning

### 7.1 PPM Threshold Pre-Filter

Original StatMaker used a PPM threshold of 3 (genes below 3 PPM in all samples excluded before analysis). Implement as a pre-processing step in the UI:

```cpp
// Before building count matrix, filter genes below PPM threshold
void GeneCountHandler::applyPpmThreshold(double threshold) {
    // Remove genes where max PPM across ALL samples < threshold
    // This reduces noise from extremely rare plasmids
}
```

Default: 3 PPM. Configurable in the Input tab.

### 7.2 Overdispersion Awareness

The paper documented:
- Non-selective conditions: phi = 0.06 (near-Poisson)
- Selective conditions: phi = 3.75 (heavy overdispersion)

The DESeq2 negative binomial model handles this naturally, but the UI should surface it:
- After dispersion estimation, report the estimated overdispersion
- Flag if dispersion is unusually high or low relative to Y2H expectations
- In the dispersion plot, optionally overlay reference lines for expected Y2H ranges

### 7.3 Replicate Recommendations

The paper found that 750mL culture volumes with 3+ replicates gave reproducible results. The UI should:
- Warn if < 3 replicates per group
- Warn if groups have unbalanced replicate counts
- Not block analysis (some experiments have constraints) but flag the limitation

### 7.4 Result Interpretation Guidance

For Y2H-naive users, provide context in the Results tab:
- "47 genes significantly enriched (padj < 0.05, log2FC > 1)"
- "These are candidate interactors for your bait protein"
- "Use MultiQuery++ to inspect individual candidates"
- Tooltip on enrichment column: "Positive log2FC = enriched with bait = candidate interactor"

---

## 8. Output Format for Downstream Modules

### 8.1 Results CSV Format

```csv
gene,refseq,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,enrichment
TBC1D15,NM_146001,1245.3,4.82,0.43,11.2,1.2e-28,3.4e-25,enriched
RABEP1,NM_007278,892.1,3.21,0.51,6.3,2.8e-10,4.1e-07,enriched
PLAA,NM_001031689,23.4,-0.12,0.89,-0.13,0.89,0.97,ns
```

### 8.2 Navigation Integration

Results CSV is read by MultiQuery++ and ReadDepth++ for:
- Ranked gene list (prev/next navigation)
- Filtering (padj < threshold, |log2FC| > threshold)
- Sorting (by padj, log2FC, baseMean)
- The `refseq` column is the key for querying depth databases

---

## 9. Implementation Phases

### Phase 1: Y2H Input UI
**Goal:** Replace generic group assignment with Y2H-specific slots

- [ ] Add `ComparisonMode` enum and selector (Two-Way / Three-Way / Custom)
- [ ] Redesign Input tab: Vector Control / Bait slots with drag-and-drop
- [ ] Three-Way mode: add third slot with custom labels
- [ ] PPM threshold parameter with default 3.0
- [ ] Auto-detect groups from filename patterns
- [ ] Replicate count warnings
- [ ] Preserve existing GeneCountHandler logic (it works)

**Dependencies:** None (UI-only changes to existing module)

### Phase 2: Three-Way Comparison Support
**Goal:** Run multiple contrasts from a single analysis

- [ ] Build design matrix for three-group comparison
- [ ] Run Wald test for each pairwise contrast
- [ ] Store results per contrast in the data model
- [ ] Results tab: contrast selector dropdown
- [ ] Three-way scatter plot (log2FC bait1 vs. log2FC bait2)
- [ ] Venn/intersection view for shared vs. specific interactors

**Dependencies:** Phase 1, existing DeseqDataSet supports arbitrary contrasts

### Phase 3: DEEPN++ Orchestrator Integration
**Goal:** Launch from main window, auto-discover files, report progress

- [ ] CLI argument parsing (`--workdir`, `--genecounts`, `--mode`)
- [ ] Auto-discover GeneCount++ output files in working directory
- [ ] Pre-populate Input tab from discovered files
- [ ] Stdout progress reporting for DEEPN++ main window
- [ ] Button state in DEEPN++ (enable "Analyze" when gene counts exist)
- [ ] Load previous results if DESeq2 output CSV already exists

**Dependencies:** Phase 1

### Phase 4: Downstream Navigation
**Goal:** Click-through from DESeq2++ results to MultiQuery++ / ReadDepth++

- [ ] Right-click context menu on Results table rows
- [ ] "Open in MultiQuery++" action (launches with --gene argument)
- [ ] "Open in ReadDepth++" action
- [ ] Double-click opens both
- [ ] Results CSV format standardized with `refseq` column for downstream lookup
- [ ] Gene name column linked to NM_* accession

**Dependencies:** Phase 1, MultiQuery++ Phase 2 (launch protocol), ReadDepth++ Phase 2

### Phase 5: Build System Migration
**Goal:** Unified CMake build for all modules

- [ ] Create top-level CMakeLists.txt
- [ ] Migrate `deepn/` from qmake to CMake
- [ ] Migrate `junction_dice/` from qmake to CMake
- [ ] Migrate `gene_count/` from qmake to CMake
- [ ] Migrate `query/` (MultiQuery++) from qmake to CMake
- [ ] Migrate `read_depth/` (ReadDepth++) from qmake to CMake
- [ ] Create `common/` shared library (GeneAnnotationDB, GeneSelector, ExportEngine)
- [ ] Integrate DESeq2 CMake into top-level build
- [ ] Update CI/CD pipeline for CMake-based builds
- [ ] Add test stage to CI (DESeq2 tests + any new tests)
- [ ] Remove old .pro files after migration verified

**Dependencies:** All modules functional, CI/CD access

### Phase 6: Visualization Enhancements
**Goal:** Publication-quality Y2H-specific plots

- [ ] Three-way scatter plot (bait1 log2FC vs. bait2 log2FC)
- [ ] Color coding by significance in each contrast
- [ ] Intersection/Venn display for shared vs. specific interactors
- [ ] Y2H-specific axis labels and legends (not generic "Group A/B")
- [ ] Overdispersion reference lines on dispersion plot
- [ ] SVG/PDF/PNG export with publication formatting
- [ ] Configurable font sizes, colors, dimensions

**Dependencies:** Phase 2

---

## 10. Testing Strategy

### 10.1 Unit Tests (Statistics Library)

Existing test files in `deseq2/statistics/tests/`:
- `test_example.cpp` -- basic pipeline usage
- `test_ppm.cpp` -- PPM analysis
- `test_deseq2.cpp` -- unit tests
- `test_convert.cpp` -- format conversion

Add:
- Three-way comparison test (3 groups, verify all 3 contrasts)
- PPM threshold filtering test
- Edge cases: single replicate, zero-count genes, all genes significant

### 10.2 Integration Tests

- End-to-end: load GeneCount++ CSV files → run DESeq2 → verify results CSV format
- Verify results CSV is readable by MultiQuery++ data loader
- Verify gene navigation from DESeq2++ launches MultiQuery++ with correct arguments

### 10.3 Validation Against Original StatMaker

If test data from the original DEEPN pipeline is available:
- Run same dataset through both original StatMaker (R/JAGS) and DESeq2++ (C++)
- Compare: do they identify the same top candidates?
- Note: different statistical methods (Bayesian vs. frequentist) will give different p-values, but the ranked order should be similar

---

## 11. Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| Eigen3 version incompatibility across platforms | Build failure on some platforms | Pin Eigen3 version in CMake, bundle as submodule if needed |
| Three-way comparison design matrix edge cases | Incorrect contrasts, wrong p-values | Validate against R DESeq2 with same input data |
| CMake migration breaks existing builds | CI/CD regression | Keep qmake files until CMake is verified on all 3 platforms |
| Qt5 static linking with Eigen3/Boost on Windows | Link errors with MXE | Test Windows build early in Phase 5, may need CMake toolchain file adjustments |
| Auto-detection misassigns files to groups | Wrong statistical comparison | Always show preview of group assignment before analysis; require user confirmation |

---

## 12. Open Questions

1. **Should the build migration (Phase 5) happen first?** It would simplify everything else, but it's high-risk. Could do a parallel track: keep qmake working while building CMake alongside.
2. **R validation**: Is original DEEPN test data available for comparison? Would help validate the three-way comparison implementation.
3. **IPC mechanism**: Stdout-based progress (simple, matches existing pattern) vs. shared memory / D-Bus (more robust)? Recommend stdout for consistency.
4. **Custom comparison mode**: Should we support arbitrary multi-group designs (e.g., 4+ conditions), or is Two-Way / Three-Way sufficient for all Y2H use cases?

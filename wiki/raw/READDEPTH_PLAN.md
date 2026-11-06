# ReadDepth++ Implementation Plan

**Module:** Read Depth Visualization & 3' Boundary Detection
**Replaces:** Original ReadDepth module
**Status:** Skeleton only (MainWindow shell)
**Build:** qmake, Qt5 C++20, QCustomPlot / Qt Charts

---

## 1. Vision

ReadDepth++ completes the prey plasmid reconstruction story. MultiQuery++ reveals where the 5' junction is (where Gal4-AD meets cDNA). ReadDepth++ reveals where the 3' end is -- the point where read coverage drops to zero. Together, they define the full extent of the prey insert without physical plasmid isolation.

The module provides an **interactive depth plot** with zoom/pan, cursor tracking, automatic 3' boundary detection, CDS annotation overlay, and publication-quality export. Multiple datasets can be compared to show how depth profiles change under selection.

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  ReadDepth++ Application                                        │
│                                                                 │
│  ┌──────────────┐  ┌──────────────────────┐  ┌──────────────┐  │
│  │ GeneSelector │  │  DepthCalculator     │  │ GeneAnnotation│  │
│  │ (query bar)  │──│  (interval counting) │──│ (ref data)   │  │
│  └──────┬───────┘  └─────────┬────────────┘  └──────┬───────┘  │
│         │                    │                       │          │
│         v                    v                       v          │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  DepthPlotWidget                                        │    │
│  │  ┌───────────────────────────────────────────────────┐  │    │
│  │  │  Interactive Depth Chart (zoom/pan/cursor)        │  │    │
│  │  │  + CDS annotation track                           │  │    │
│  │  │  + 3' boundary marker                             │  │    │
│  │  │  + Optional MultiQuery++ junction overlay         │  │    │
│  │  └───────────────────────────────────────────────────┘  │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────┐  │
│  │ BoundaryDetect│  │ ExportEngine │  │ BatchRunner          │  │
│  │ (auto 3' find)│  │ (CSV/SVG/PDF)│  │ (top N genes)        │  │
│  └──────────────┘  └──────────────┘  └──────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Core Algorithm: Interval-Based Depth Counting

### 3.1 How It Works

For a queried gene, ReadDepth counts how many sequencing reads contain specific subsequences of the gene's cDNA. This measures coverage at each position.

```
Gene cDNA sequence (e.g., 3000 nt):
ATGCCTGAACTTGGC...                                              ...TAAGGCTAA

Divide into intervals of W bp (default W=25), spaced D bp apart (default D=50):
Interval 1: position 1-25    → "ATGCCTGAACTTGGCAAACTGGTTT"
Interval 2: position 51-75   → "GGCATTCCTAAGGTTCCTGAACTGG"
Interval 3: position 101-125 → "TCAGGTTTCCAAGGTCCTTTGAACC"
...

For each interval:
  Count = number of reads in SQLite `reads` table whose `sequence` contains the interval substring

Plot: position (x) vs count (y)
```

### 3.2 Optimized Algorithm

The naive approach (substring search every read for every interval) is O(R * I * L) where R=reads, I=intervals, L=read length. For large datasets this is slow.

**Optimized approach using the maps table:**

Instead of searching raw read sequences, use the mapping data:
1. Query `maps` table for all reads mapping to the target gene
2. For each mapped read, we know `qstart` and `qend` (the portion of the read that mapped)
3. The read covers positions `qstart` to `qend` on the gene
4. For each interval, count how many reads have coverage spanning that interval

```sql
-- For interval at position P with width W:
SELECT COUNT(DISTINCT read) FROM maps
WHERE gene = :gene AND qstart <= :P AND qend >= :P + :W;
```

This is O(I * log(R)) with proper indexes -- much faster.

**Fallback for unmapped regions:** If the maps table doesn't cover all positions (e.g., reads that matched but weren't in the maps table), fall back to substring search on the reads table for specific intervals.

### 3.3 Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| Interval width (W) | 25 bp | 10-100 bp | Size of each counting window |
| Interval spacing (D) | 50 nt | 25-500 nt | Distance between interval starts |
| Minimum depth | 0 | 0-100 | Filter: don't plot intervals below this count |

Original ReadDepth allowed adjustment in 50 nt increments. We'll allow finer control with a slider + direct entry.

---

## 4. Automatic 3' Boundary Detection

### 4.1 Algorithm

The 3' boundary is where read depth drops from "has coverage" to "no coverage." This is a step function detection problem.

```
Depth profile:
  ████████████████████████████████████░░░░░░░░░░░░░░░░░░░░░░░
  ^                                  ^                        ^
  5' junction (from MultiQuery++)    3' boundary (we detect)  gene end

Algorithm:
1. Smooth the depth profile (moving average, window = 3 intervals)
2. Compute the gradient (difference between consecutive intervals)
3. Find the position with the largest negative gradient
4. Refine: walk backward from the drop to find the last interval with depth > threshold
5. Report this as the estimated 3' boundary (± interval spacing resolution)
```

### 4.2 Confidence Scoring

Not all drop-offs are clean. Score the boundary detection:
- **High confidence**: Sharp drop (>80% decrease in 1-2 intervals), consistent across datasets
- **Medium confidence**: Gradual taper (50-80% decrease over 3-5 intervals)
- **Low confidence**: Noisy profile, multiple drop-offs, or very low overall depth

Display the confidence alongside the boundary position.

### 4.3 Multiple Inserts

Some genes may have prey plasmids with different 3' boundaries (different library clones). The algorithm should detect all major drop-off points, not just the first one. Report as: "Primary 3' boundary at position X (confidence: high), secondary at position Y (confidence: medium)."

---

## 5. UI Layout Specification

### 5.1 Main View

```
┌─────────────────────────────────────────────────────────────────────────┐
│  [Gene: NM_146001 ▾]  [TBC1D15]  [Dataset: experiment1.sqlite ▾]      │
│  [◀ Prev]  [Next ▶]  [Batch: Top 50 ▾]  Interval: [25]bp  Gap: [50]nt│
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Read Depth                                                             │
│                                                                         │
│  Count ▲                                                                │
│   800  ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓                  │
│   600  ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓                  │
│   400  ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓                │
│   200  ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░              │
│     0  ┼──────────────────────────────────────────────────────── pos ▶  │
│        0     500    1000   1500   2000   2500   3000                    │
│                                            │                            │
│  ┌─ Gene Annotation ──────────────────── ──┼──────────────────────┐     │
│  │ ░░░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░░░░░░░░░░░░░░░░░░░░ │     │
│  │ 5'UTR         CDS                3'UTR                        │     │
│  │      ▲ATG                    ▲STOP                            │     │
│  └───────────────────────────────────────────────────────────────┘     │
│                                                                         │
│  3' Boundary: position 2,650 (confidence: HIGH)   ▼ detected here      │
│  Insert extent: 678 → 2,650 (1,972 nt)                                 │
│                                                                         │
│  Cursor: position 1,245  |  Depth: 612 reads  |  Zoom: 100%           │
├─────────────────────────────────────────────────────────────────────────┤
│  [Export CSV]  [Export Figure (SVG/PDF/PNG)]  [Batch Report]           │
└─────────────────────────────────────────────────────────────────────────┘
```

### 5.2 Interactive Controls

**Zoom:**
- Mouse wheel zooms in/out on x-axis (centered on cursor position)
- Zoom range: full gene view (100%) to single-interval view
- Zoom indicator in status bar

**Pan:**
- Click-drag to pan along x-axis when zoomed in
- Scroll bar along bottom for quick navigation

**Cursor:**
- Crosshair cursor tracks mouse position
- Status bar shows: position, depth count, interval index
- Click to place a persistent marker (for boundary annotation)

**Interval controls:**
- Spinner for interval width (W): 10-100 bp, step 5
- Spinner for interval spacing (D): 25-500 nt, step 25
- Changes trigger immediate recalculation and replot

### 5.3 Dataset Comparison Mode

```
┌─────────────────────────────────────────────────────────────────────────┐
│  [Gene: NM_146001]  [Compare Mode ✓]  [Overlay ○ / Split ●]           │
├─────────────────────────────────────┬───────────────────────────────────┤
│  Unselected Population              │  Selected Population              │
│                                     │                                   │
│  Count ▲                            │  Count ▲                          │
│   50   ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓   │  800   ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │
│   25   ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓   │  400   ┤  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │
│    0   ┼───────────────────── pos    │    0   ┼───────────────── pos    │
│                                     │                                   │
│  ░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░░░    │  ░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░░  │
│       ← scroll/zoom linked →       │       ← scroll/zoom linked →      │
│                                     │                                   │
│  3' Boundary: 2,650 (HIGH)         │  3' Boundary: 2,650 (HIGH)        │
└─────────────────────────────────────┴───────────────────────────────────┘
```

Two comparison modes:
- **Split** (default): Side-by-side synchronized panels (same as MultiQuery++)
- **Overlay**: Both depth profiles on same axes with transparency. Selected in bold color, unselected in semi-transparent. Useful for seeing enrichment directly.

---

## 6. Integration with MultiQuery++

### 6.1 Linked Context

When both MultiQuery++ and ReadDepth++ are open for the same gene:
- MultiQuery++ shows the 5' junction positions
- ReadDepth++ shows the 3' boundary
- Together they define the complete insert: `[5' junction position] → [3' boundary]`

### 6.2 Junction Overlay on Depth Plot

Optionally overlay MultiQuery++ junction positions on the depth plot:
- Vertical dashed lines at each junction position
- Colored by the same scheme (in-frame = blue, OOF = amber, etc.)
- Helps researchers see the relationship between junction sites and coverage

### 6.3 Insert Extent Summary

When both datasets are available, display:
```
Insert Reconstruction:
  5' junction: position 678 (from MultiQuery++)
  3' boundary: position 2,650 (from ReadDepth++)
  Insert length: 1,972 nt
  Covers: 68% of CDS (542-2,180)
  In-frame: YES
```

---

## 7. Class Design

```
MainWindow
├── GeneSelector              -- search bar, prev/next, sort/filter
├── DepthPlotWidget           -- QCustomPlot with depth area chart
│   ├── DepthChartLayer       -- the actual depth bars/area fill
│   ├── AnnotationTrack       -- CDS/UTR gene structure below plot
│   ├── BoundaryMarker        -- vertical line at detected 3' boundary
│   ├── JunctionOverlay       -- optional vertical lines from MultiQuery++
│   └── CursorTracker         -- crosshair + status bar updates
├── DepthCalculator           -- computes depth profile from SQLite
│   ├── IntervalCounter       -- interval-based counting algorithm
│   └── DepthProfile          -- result: vector of (position, count) pairs
├── BoundaryDetector          -- automatic 3' boundary detection
│   ├── Smoother              -- moving average filter
│   └── GradientAnalyzer      -- finds largest drop-off
├── GeneAnnotationDB          -- shared with MultiQuery++ (gene ref data)
├── ComparisonManager         -- split/overlay mode management
│   └── SyncController        -- links zoom/scroll
├── BatchRunner               -- process top N genes
│   └── BatchWorker           -- QThread background processing
└── ExportEngine              -- CSV + figure export
    ├── CsvExporter           -- depth profile data
    └── FigureExporter        -- SVG/PDF/PNG
```

---

## 8. Data Structures

```cpp
// Single depth measurement at an interval position
struct DepthPoint {
    int position;       // start position of the interval on the gene
    int count;          // number of reads covering this interval
    double normalized;  // optional: normalized count (e.g., per million)
};

// Complete depth profile for one gene from one dataset
struct DepthProfile {
    QString refseq;                 // NM_* accession
    QString geneName;               // gene symbol
    int intervalWidth;              // W (bp)
    int intervalSpacing;            // D (nt)
    int mRNALength;                 // total gene length
    int orfStart;                   // CDS start
    int orfEnd;                     // CDS end
    int totalReads;                 // total reads for this gene
    QVector<DepthPoint> points;     // the depth profile
    QString datasetLabel;
    QString sourceFile;
};

// Detected 3' boundary
struct BoundaryResult {
    int position;                   // estimated 3' boundary position
    double confidence;              // 0.0 to 1.0
    QString confidenceLabel;        // "HIGH", "MEDIUM", "LOW"
    int dropMagnitude;              // depth just before vs after boundary
    QVector<BoundaryResult> secondary;  // additional boundaries if multiple inserts
};

// Complete insert reconstruction (with MultiQuery++ data)
struct InsertExtent {
    int fivePrimeJunction;          // from MultiQuery++ (0 if unknown)
    int threePrimeBoundary;         // from ReadDepth++
    int insertLength;               // 3' - 5'
    double cdsOverlapPercent;       // how much of CDS the insert covers
    bool inFrame;                   // based on 5' junction frame classification
};
```

---

## 9. Export Capabilities

### 9.1 Data Export
- **Depth CSV**: position, count, normalized_count for each interval
- **Boundary CSV**: gene, 3'_position, confidence, insert_length
- **Batch CSV**: summary across all queried genes

### 9.2 Figure Export
- **SVG/PDF/PNG**: Depth plot with CDS annotation, boundary marker, legend
- DPI control for PNG (300 default)
- Configurable dimensions (default: 8" x 4" for publication single-column)
- Optional: include junction overlay from MultiQuery++ in exported figure

---

## 10. Integration with DEEPN++ Orchestrator

### 10.1 Launch Protocol
```
ReadDepth++ --workdir /path/to/experiment \
            --datasets file1.sqlite,file2.sqlite \
            --generef data/gene_dictionary/ \
            --gene NM_146001 \
            --junctions multiquery_results.csv  # optional
```

### 10.2 File Discovery
On launch, discover:
- `analyzed_files/*.sqlite` -- depth databases
- `gene_count_summary/*.csv` -- for gene lookup
- DESeq2++ results (for ranked navigation)
- MultiQuery++ junction data (for overlay, if available)

### 10.3 Button State in Main Window
Enable "Read Depth" button when:
- At least one `.sqlite` file exists in `analyzed_files/`
- The SQLite database has a populated `reads` table

---

## 11. Implementation Phases

### Phase 1: Core Depth Algorithm
**Goal:** Calculate and return depth profiles from SQLite data

- [ ] Implement `IntervalCounter` -- divide gene into intervals, count reads per interval
- [ ] Implement maps-table-based counting (fast path using qstart/qend)
- [ ] Implement sequence-based counting (fallback using reads.sequence substring search)
- [ ] Implement `DepthProfile` data structure
- [ ] Share `GeneAnnotationDB` with MultiQuery++ (same SQLite gene reference)
- [ ] Unit tests: known depth profiles from test data

**Dependencies:** Gene reference data in SQLite format (shared with MultiQuery++)

### Phase 2: Interactive Depth Plot
**Goal:** Display depth profile with zoom, pan, cursor tracking

- [ ] Implement `DepthPlotWidget` with QCustomPlot
- [ ] Area fill chart (filled bars or continuous area) for depth profile
- [ ] X-axis: nucleotide position; Y-axis: read count
- [ ] Mouse wheel zoom (centered on cursor)
- [ ] Click-drag pan when zoomed in
- [ ] Cursor crosshair with status bar (position, depth, zoom level)
- [ ] Gene annotation track below plot (CDS/UTR regions, start/stop markers)
- [ ] Interval parameter controls (width spinner, spacing spinner) with live recalculation

**Dependencies:** Phase 1 complete

### Phase 3: 3' Boundary Detection
**Goal:** Automatically find and mark the 3' boundary

- [ ] Implement `Smoother` -- moving average filter on depth profile
- [ ] Implement `GradientAnalyzer` -- find largest negative gradient
- [ ] Implement boundary refinement (walk back to last depth > threshold)
- [ ] Confidence scoring (sharp vs. gradual vs. noisy)
- [ ] Multiple boundary detection (for genes with multiple prey clones)
- [ ] Visual boundary marker on the plot (vertical line + label)
- [ ] Insert extent display: "3' boundary at position X (confidence: Y)"

**Dependencies:** Phase 2 complete

### Phase 4: Dataset Comparison
**Goal:** Compare depth profiles between selected and unselected populations

- [ ] Split panel mode with synchronized zoom/scroll
- [ ] Overlay mode with transparency
- [ ] Toggle between split and overlay
- [ ] Independent y-axes (depth scales differ between datasets)
- [ ] Side-by-side boundary comparison

**Dependencies:** Phase 2 complete

### Phase 5: MultiQuery++ Integration & Navigation
**Goal:** Link with junction data, navigate through candidates

- [ ] Junction overlay on depth plot (vertical dashed lines from MultiQuery++ data)
- [ ] Insert extent summary (5' junction + 3' boundary = full insert)
- [ ] Gene navigation: prev/next through DESeq2++ ranked list
- [ ] Search by NM_* accession or gene name
- [ ] Filter by p-value, log2FC

**Dependencies:** MultiQuery++ Phase 1 (shared GeneAnnotationDB), DESeq2++ results format

### Phase 6: Export & Batch Mode
**Goal:** Publication-quality output and multi-gene processing

- [ ] CSV export: depth profile, boundary results, insert extents
- [ ] Figure export: SVG/PDF/PNG with annotation, boundary marker, legend
- [ ] Batch runner: process top N genes, generate per-gene figures + summary CSV
- [ ] Batch progress reporting via QThread
- [ ] Export dialog with format/DPI/size controls

**Dependencies:** All previous phases

### Phase 7: Integration & Polish
**Goal:** Launch from DEEPN++, cross-platform testing

- [ ] CLI argument parsing for orchestrator launch
- [ ] Auto-discovery of files in working directory
- [ ] Button state integration in DEEPN++ main window
- [ ] Performance optimization for large datasets (cache depth profiles)
- [ ] Cross-platform testing

**Dependencies:** All previous phases

---

## 12. Shared Components with MultiQuery++

These components should be implemented once and shared:
- `GeneAnnotationDB` -- gene reference data loading
- `GeneSelector` -- search bar, prev/next navigation, sort/filter
- `ExportEngine` -- CSV and figure export
- `BatchRunner` -- background batch processing pattern
- `SyncController` -- synchronized zoom/scroll for split panels

Consider extracting these into a shared library (`deepn_common/`) used by both MultiQuery++ and ReadDepth++.

---

## 13. Open Questions

1. **Maps-based vs. sequence-based counting**: The maps table has qstart/qend for read positions on the gene. Is this sufficient for depth counting, or do we need the actual sequence substring search for accuracy?
2. **Shared library**: Should MultiQuery++ and ReadDepth++ share a common library, or is code duplication across two small modules acceptable?
3. **Combined module**: Should ReadDepth be a tab/panel within MultiQuery++ rather than a separate application? The user selected "Interactive with zoom/pan/cursor" as a standalone module, but integration was also an option.
4. **Performance**: For genes with >100K mapped reads, is the interval counting fast enough for interactive parameter changes? May need to pre-compute at multiple resolutions.

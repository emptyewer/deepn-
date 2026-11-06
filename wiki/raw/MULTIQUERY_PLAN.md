# MultiQuery++ Implementation Plan

**Module:** Junction Analysis & Visualization Dashboard
**Replaces:** Original BlastQuery module
**Status:** Skeleton only (MainWindow + DockWidget scaffold)
**Build:** qmake, Qt5 C++20, Qt Charts / QCustomPlot

---

## 1. Vision

MultiQuery++ is the researcher's primary inspection tool after DESeq2++ identifies candidate interactors. For each gene, it answers: **where do the junction fusions land, how abundant are they, and are they biologically meaningful?**

The UI is a **linked dashboard** (not tabs) where selecting a junction in the plot highlights the table row, and selecting a table row pans the plot to that position. A detail panel shows sequence context and metadata. Multiple datasets are compared via **synchronized split panels** with linked scroll/zoom.

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  MultiQuery++ Application                                       │
│                                                                 │
│  ┌──────────────┐  ┌─────────────────────┐  ┌───────────────┐  │
│  │ GeneSelector │  │ JunctionDataModel   │  │ GeneAnnotation│  │
│  │ (query bar)  │──│ (loads SQLite data)  │──│ (ref sequences│  │
│  └──────┬───────┘  └─────────┬───────────┘  │  ORF coords)  │  │
│         │                    │               └───────┬───────┘  │
│         v                    v                       v          │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  DashboardView                                          │    │
│  │  ┌─────────────────────────┬──────────────────────┐     │    │
│  │  │  JunctionPlotWidget     │  GeneDetailPanel     │     │    │
│  │  │  (bar chart + mRNA      │  (sequence, domain,  │     │    │
│  │  │   track + CDS annot)    │   frame, metadata)   │     │    │
│  │  ├─────────────────────────┤                      │     │    │
│  │  │  JunctionTableWidget    │                      │     │    │
│  │  │  (sortable results)     │                      │     │    │
│  │  └─────────────────────────┴──────────────────────┘     │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────┐  │
│  │ BatchRunner  │  │ ExportEngine │  │ ComparisonManager    │  │
│  │ (top N genes)│  │ (CSV/SVG/PDF)│  │ (sync split panels)  │  │
│  └──────────────┘  └──────────────┘  └──────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Data Model

### 3.1 Input Data Sources

**Primary: SQLite depth databases** (`analyzed_files/*.sqlite`)
Created by JunctionDice++ during dicing. Schema:
```sql
-- Reads table: raw sequencing reads with extracted cDNA sequences
CREATE TABLE reads (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    read TEXT,       -- read name (@SRR...)
    sequence TEXT    -- extracted cDNA sequence
);
CREATE INDEX reads_idx ON reads (read);

-- Maps table: BLAT/BLAST mapping results per read
CREATE TABLE maps (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    read TEXT,       -- read name (FK to reads)
    gene TEXT,       -- gene name
    qstart INTEGER,  -- query start position (= QueryStart in original)
    qend INTEGER,    -- query end position
    refseq TEXT,     -- NM_* accession number
    frame TEXT,      -- reading frame classification
    location TEXT    -- CDS classification (upstream/in_orf/downstream)
);
CREATE INDEX maps_idx ON maps (read, gene, refseq);
```

**Secondary: Gene reference data** (`data/`)
- Gene lists (`.prn` files): NM accession, gene name, chromosome, ORF start/end, full sequence
- Exon dictionaries (`.p` files): Python pickle format -- need C++ reader or convert to SQLite/JSON at build time
- BLAST databases: for sequence lookup

**Tertiary: Gene count summaries** (`gene_count_summary/*.csv`)
- Gene name, PPM values per sample
- Used for quick lookup of which genes to query

### 3.2 Core Data Structures

```cpp
// A single junction site on a gene
struct JunctionSite {
    int position;           // nucleotide position on mRNA
    double ppm;             // abundance (Parts Per Million)
    int queryStart;         // distance from junction tag to cDNA match
    QString cdsClass;       // "upstream", "in_orf", "downstream"
    QString frame;          // "in_frame", "out_of_frame", "backwards", "intron"
    QString refseq;         // NM_* accession
    QString geneName;       // gene symbol
    int rawCount;           // raw junction count (before PPM normalization)
};

// All junction data for one gene from one dataset
struct GeneJunctionProfile {
    QString refseq;                 // NM_* accession
    QString geneName;               // gene symbol
    int orfStart;                   // CDS start position on mRNA
    int orfEnd;                     // CDS end position on mRNA
    int mRNALength;                 // total mRNA length
    QString sequence;               // full mRNA sequence (for display)
    QVector<JunctionSite> sites;    // all junction sites
    int totalReads;                 // total reads in this dataset
    QString datasetLabel;           // e.g., "Selected", "Unselected"
    QString sourceFile;             // SQLite file path
};

// Filtered view: collapsed by position
struct CollapsedJunction {
    int position;
    double totalPpm;                // sum of PPM across QueryStart variants
    int variantCount;               // number of distinct QueryStart values
    QString dominantFrame;          // most common frame at this position
    QString cdsClass;
    QVector<JunctionSite> variants; // individual entries before collapsing
};
```

### 3.3 Gene Annotation Model

The exon dictionaries (`.p` files) are Python pickles. For C++ access:

**Option A (recommended):** Pre-convert `.p` files to SQLite at build/install time using a Python script. Ship the SQLite version.
```sql
CREATE TABLE genes (
    refseq TEXT PRIMARY KEY,
    gene_name TEXT,
    chromosome TEXT,
    orf_start INTEGER,
    orf_end INTEGER,
    mrna_length INTEGER,
    sequence TEXT
);
```

**Option B:** Convert to JSON at build time, load at runtime.

**Option C:** Use pybind11/embedded Python (rejected -- adds runtime dependency).

---

## 4. UI Layout Specification

### 4.1 Main Dashboard (Single Gene View)

```
┌─────────────────────────────────────────────────────────────────────────┐
│  [Gene: NM_146001 ▾]  [TBC1D15]  [Dataset: experiment1.sqlite ▾]  [⚙] │
│  [◀ Prev]  [Next ▶]  [Batch: Top 50 ▾]                                │
├───────────────────────────────────────────────────┬─────────────────────┤
│                                                   │                     │
│  Junction Abundance Plot                          │  Gene Detail        │
│                                                   │                     │
│  PPM ▲                                            │  Gene: TBC1D15      │
│  500 ┤                                            │  RefSeq: NM_146001  │
│  400 ┤  ██                                        │  mRNA: 3,421 nt     │
│  300 ┤  ██                                        │  CDS: 542-2,180     │
│  200 ┤  ██  ██                                    │                     │
│  100 ┤  ██  ██  ██      ██                        │  ── Selected ──     │
│    0 ┼──██──██──██──────██──────── position ▶     │  Position: 678      │
│      500    1000   1500   2000                    │  PPM: 423.5         │
│                                                   │  Frame: In-frame    │
│  ┌─ mRNA Track ─────────────────────────────┐     │  CDS: In ORF        │
│  │ ░░░░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░░░░ │     │  QueryStart: 1      │
│  │ 5'UTR        CDS              3'UTR      │     │                     │
│  │       ▲ATG                 ▲STOP         │     │  Sequence context:  │
│  └──────────────────────────────────────────┘     │  ...ATGCCTGAA...    │
│                                                   │                     │
├───────────────────────────────────────────────────┤                     │
│                                                   │                     │
│  Junction Results Table                           │                     │
│  ┌──────┬───────┬────────┬─────────┬──────────┐   │                     │
│  │ Pos  │  PPM  │ QStart │   CDS   │  Frame   │   │                     │
│  ├──────┼───────┼────────┼─────────┼──────────┤   │                     │
│  │► 678 │ 423.5 │    1   │  In ORF │ In-frame │   │                     │
│  │  912 │ 287.2 │    1   │  In ORF │ In-frame │   │                     │
│  │ 1456 │  89.1 │    1   │  In ORF │ Out-of-f │   │                     │
│  │ 2301 │  12.4 │    3   │  Down   │ Out-of-f │   │                     │
│  └──────┴───────┴────────┴─────────┴──────────┘   │                     │
│  [Collapse by Position ☐]  [Show In-Frame Only ☐] │                     │
│                                                   │                     │
├───────────────────────────────────────────────────┴─────────────────────┤
│  [Export CSV]  [Export Figure (SVG/PDF/PNG)]  [Export Batch Report]     │
└─────────────────────────────────────────────────────────────────────────┘
```

### 4.2 Dataset Comparison (Synchronized Split Panels)

```
┌─────────────────────────────────────────────────────────────────────────┐
│  [Gene: NM_146001 ▾]  [TBC1D15]  [Compare Mode ✓]                     │
├──────────────────────────────────┬──────────────────────────────────────┤
│  Unselected Population           │  Selected Population                 │
│                                  │                                      │
│  PPM ▲                           │  PPM ▲                               │
│   50 ┤  ██  ██  ██  ██  ██  ██  │  500 ┤  ██                           │
│   25 ┤  ██  ██  ██  ██  ██  ██  │  400 ┤  ██                           │
│    0 ┼──────────────────── pos ▶ │  200 ┤  ██  ██                       │
│                                  │    0 ┼──██──██──────── pos ▶         │
│  ░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░   │  ░░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░      │
│       ← scroll/zoom linked →    │       ← scroll/zoom linked →         │
├──────────────────────────────────┼──────────────────────────────────────┤
│  Pos │  PPM │ Frame              │  Pos │  PPM │ Frame                  │
│  678 │  8.2 │ In-fr              │  678 │ 423  │ In-fr   ▲ enriched    │
│  912 │  5.1 │ In-fr              │  912 │ 287  │ In-fr   ▲ enriched    │
│ 1456 │  4.8 │ OOF                │ 1456 │  89  │ OOF     ▲ enriched    │
└──────────────────────────────────┴──────────────────────────────────────┘
```

**Synchronization behavior:**
- Zoom level shared between panels (same nucleotide range visible)
- Scroll position linked (pan one, both pan)
- Y-axes independent (PPM scales will differ dramatically between selected/unselected)
- Table selection linked (click row in left, highlights same position in right)
- Enrichment markers in the comparison table (▲/▼ arrows with fold-change)

### 4.3 Gene Navigation Bar

```
┌─────────────────────────────────────────────────────────────────┐
│  Search: [NM_146001    🔍]  Gene: TBC1D15                      │
│  [◀ Prev] [Next ▶]  Showing 3 of 847 candidates                │
│  Sort by: [Adjusted p-value ▾]  Filter: [padj < 0.05 ▾]       │
│  Batch: [Generate Top 50 Report]                                │
└─────────────────────────────────────────────────────────────────┘
```

- **Search**: NM_* accession or gene name (autocomplete from gene count data)
- **Prev/Next**: Navigate through ranked candidate list (from DESeq2++ results)
- **Sort/Filter**: Reorder candidates by p-value, log2FC, base mean
- **Batch**: Trigger batch export for top N genes

---

## 5. Junction Plot Specification

### 5.1 Bar Chart

- **X-axis**: Nucleotide position along the mRNA (1 to mRNA length)
- **Y-axis**: Junction abundance in PPM
- **Bar width**: Proportional to zoom level; minimum 3px at any zoom
- **Bar colors** (biological significance):

| Category | Color | Meaning |
|----------|-------|---------|
| In ORF + In-frame | Deep blue (#1E40AF) | True in-frame fusion -- most biologically relevant |
| Upstream + In-frame | Teal (#0D9488) | Correct frame but upstream of CDS start -- may still produce functional protein |
| In ORF + Out-of-frame | Amber (#D97706) | Within CDS but wrong frame -- usually non-functional |
| Downstream / Backwards | Grey (#6B7280) | Past stop codon or reverse orientation |

- **Hover tooltip**: Position, PPM, Frame, CDS, QueryStart, raw count
- **Click**: Selects junction, highlights table row, populates detail panel

### 5.2 mRNA Track (Below Plot)

A horizontal annotation track showing gene structure:
- **5' UTR**: Light grey region
- **CDS**: Dark region with start codon (ATG) and stop codon marked with vertical lines
- **3' UTR**: Light grey region
- **Scale bar**: Nucleotide positions at regular intervals
- Synchronized with the junction plot x-axis (zoom/pan linked)

### 5.3 Optional Sequence Display

Below the mRNA track, optionally show the actual nucleotide sequence:
- CDS nucleotides in **black monospace**
- UTR nucleotides in **grey monospace**
- Only visible at high zoom levels (< 200bp visible)
- Codon boundaries marked with thin separators

---

## 6. Batch Analysis Mode

### 6.1 Batch Runner

Input: DESeq2++ results CSV (gene list ranked by adjusted p-value)
Parameters:
- Top N genes (default: 50)
- p-value cutoff (default: 0.05)
- log2FC cutoff (default: 1.0)
- Datasets to include

### 6.2 Batch Output

For each gene, generate:
1. Junction plot (SVG/PDF)
2. Results table (CSV)
3. Summary statistics (total junctions, in-frame percentage, dominant position)

Aggregate report:
- Summary CSV: gene, total_junctions, max_ppm, in_frame_count, dominant_position, dominant_frame
- HTML report with embedded SVG plots (optional stretch goal)

---

## 7. Export Capabilities

### 7.1 Data Export
- **CSV**: Full junction results table with all columns
- **Filtered CSV**: Collapsed-by-position view
- **Batch CSV**: Summary across all queried genes

### 7.2 Figure Export
- **SVG**: Vector format for publication (editable in Illustrator/Inkscape)
- **PDF**: Direct publication-ready format
- **PNG**: Raster at configurable DPI (default 300 for publication, 150 for screen)

Figure export includes:
- Junction plot with color legend
- mRNA track with CDS annotation
- Title bar with gene name, accession, dataset label
- Axis labels and scale

QCustomPlot supports all three formats natively via `savePdf()`, `savePng()`, `saveSvg()` (with QSvgGenerator).

---

## 8. Integration with DEEPN++ Orchestrator

### 8.1 Launch Protocol

MultiQuery++ is launched from the DEEPN++ main window as a child process:
```
MultiQuery++ --workdir /path/to/experiment \
             --datasets file1.sqlite,file2.sqlite \
             --genecounts gene_count_summary/ \
             --generef data/gene_dictionary/ \
             --results deseq2_results.csv
```

### 8.2 File Discovery

On launch, automatically discover:
- `analyzed_files/*.sqlite` -- depth databases from JunctionDice++
- `gene_count_summary/*.csv` -- gene counts from GeneCount++
- DESeq2++ results CSV (if present, enables ranked navigation)

### 8.3 Button State in Main Window

Enable "Multi Query" button when:
- At least one `.sqlite` file exists in `analyzed_files/`
- At least one `.csv` file exists in `gene_count_summary/`

---

## 9. Class Design

```
MainWindow
├── GeneSelector          -- search bar, prev/next, sort/filter
├── DashboardView         -- QSplitter-based layout manager
│   ├── JunctionPlotWidget    -- QCustomPlot with bar chart + mRNA track
│   ├── JunctionTableWidget   -- QTableView with sortable model
│   └── GeneDetailPanel       -- QWidget with labels and sequence view
├── ComparisonManager     -- manages synchronized split panel mode
│   ├── SyncController        -- links zoom/scroll between panels
│   └── EnrichmentCalculator  -- computes fold-change between datasets
├── JunctionDataModel     -- loads SQLite, computes PPM, classifies junctions
│   ├── SqliteLoader          -- reads depth database
│   └── JunctionClassifier    -- assigns CDS/frame based on gene annotation
├── GeneAnnotationDB      -- loads gene reference data (SQLite format)
├── BatchRunner           -- processes top N genes, collects results
│   └── BatchWorker           -- QThread for background batch processing
└── ExportEngine          -- CSV writer + figure exporter
    ├── CsvExporter
    └── FigureExporter        -- SVG/PDF/PNG via QCustomPlot
```

---

## 10. Implementation Phases

### Phase 1: Data Layer & Core Model (Foundation)
**Goal:** Load junction data from SQLite, classify junctions, load gene annotations

- [ ] Convert exon dictionary `.p` files to SQLite (Python script + ship SQLite)
- [ ] Implement `GeneAnnotationDB` -- load gene reference data
- [ ] Implement `SqliteLoader` -- read depth database, join reads + maps tables
- [ ] Implement `JunctionClassifier` -- compute PPM, classify CDS/frame per junction
- [ ] Implement `JunctionDataModel` -- expose `GeneJunctionProfile` for a queried gene
- [ ] Unit tests: load test SQLite, verify junction classification matches expected values

**Dependencies:** Gene reference data in SQLite format, test depth database

### Phase 2: Single-Gene Dashboard (Core UI)
**Goal:** Working dashboard for one gene from one dataset

- [ ] Implement `JunctionPlotWidget` -- bar chart with color-coded bars, x-axis = position, y-axis = PPM
- [ ] Implement mRNA track below the plot (CDS/UTR annotation, start/stop markers)
- [ ] Implement `JunctionTableWidget` -- sortable table (Position, PPM, QStart, CDS, Frame)
- [ ] Implement `GeneDetailPanel` -- metadata, sequence context for selected junction
- [ ] Wire bidirectional linking: click bar → highlight row → populate detail (and reverse)
- [ ] Implement `GeneSelector` -- NM_* search with autocomplete, gene name display
- [ ] Layout with QSplitter (plot+table left, detail right)

**Dependencies:** Phase 1 complete

### Phase 3: Dataset Comparison (Split Panels)
**Goal:** Compare selected vs. unselected (or bait1 vs. bait2) side by side

- [ ] Implement `ComparisonManager` -- create side-by-side DashboardView instances
- [ ] Implement `SyncController` -- link x-axis zoom/scroll between panels
- [ ] Independent y-axes (PPM scales differ between datasets)
- [ ] Linked table selection (click position in one panel, highlights same in other)
- [ ] Enrichment markers in table (fold-change arrows, ▲/▼)
- [ ] Dataset selector dropdown for each panel

**Dependencies:** Phase 2 complete

### Phase 4: Navigation & Filtering
**Goal:** Fast iteration through candidate genes

- [ ] Load DESeq2++ results CSV for ranked gene list
- [ ] Prev/Next navigation through candidates
- [ ] Sort by adjusted p-value, log2FC, base mean
- [ ] Filter by p-value threshold, log2FC cutoff
- [ ] "Collapse by Position" toggle for filtered results view
- [ ] "In-Frame Only" filter toggle
- [ ] Keyboard shortcuts: arrow keys for prev/next, Enter to select from search

**Dependencies:** Phase 2 complete, DESeq2++ results CSV format

### Phase 5: Export & Batch Mode
**Goal:** Publication-quality outputs and automated multi-gene analysis

- [ ] Implement `CsvExporter` -- results table, filtered table, comparison table
- [ ] Implement `FigureExporter` -- SVG, PDF, PNG at configurable DPI
- [ ] Figure includes: plot, mRNA track, legend, title, axes, dataset label
- [ ] Implement `BatchRunner` + `BatchWorker` -- process top N genes on background thread
- [ ] Batch output: per-gene CSV + figure, summary CSV, progress reporting
- [ ] Export dialog with format selection and DPI/size controls

**Dependencies:** Phases 2-4 complete

### Phase 6: Integration & Polish
**Goal:** Launch from DEEPN++ orchestrator, cross-platform testing

- [ ] CLI argument parsing for orchestrator launch
- [ ] Auto-discovery of files in working directory
- [ ] Button state integration in DEEPN++ main window
- [ ] Sequence display at high zoom levels (codon-level view)
- [ ] Hover tooltips on plot bars
- [ ] Cross-platform testing (macOS, Linux, Windows)
- [ ] Performance profiling with large datasets

**Dependencies:** All previous phases

---

## 11. Open Questions

1. **Exon dictionary format**: Should we pre-convert `.p` files to SQLite during build, or ship a conversion tool that runs on first launch?
2. **Plotting library**: QCustomPlot (already used by DESeq2++) vs. Qt Charts (already imported in query.pro). Recommend QCustomPlot for consistency and better export support.
3. **Gene sequence source**: Can we extract mRNA sequences from the BLAST database files, or do we need a separate sequence file?
4. **Memory**: For batch mode with 50+ genes, should we reuse one plot widget and re-render, or create plot images off-screen?

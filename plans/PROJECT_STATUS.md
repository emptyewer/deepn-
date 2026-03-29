# DEEPN++ Project Status

**Date:** 2026-03-26
**Branch:** master (pending: GeneCount++ multi-file crash fix)
**Last Activity:** 2026-03-24 (CMake build system + DESeq2 module merge)

---

## Project Overview

DEEPN++ is a C++ rewrite of the **DEEPN** (Dynamic Enrichment for Evaluation of Protein Networks) pipeline, originally published in *Cell Reports* (Piper, Pashkova, Peterson, Krishnamani, Breheny, Stamnes; 2016; [PMC5594928](https://pmc.ncbi.nlm.nih.gov/articles/PMC5594928/)). The original user guide is at [github.com/emptyewer/DEEPN](https://github.com/emptyewer/DEEPN/releases).

DEEPN is a batch-processing method for **Yeast Two-Hybrid (Y2H)** screens. Instead of picking individual colonies, it uses competitive growth and high-throughput sequencing to identify protein-protein interactions in parallel:

1. **Bait** proteins (Gal4-DBD fusions) are mated with a library of **prey** proteins (Gal4-AD fused to cDNA fragments), generating millions of diploid yeast strains.
2. Under **selective growth** (SD-His +/- 3AT), plasmids encoding proteins that interact with the bait **enrich** in the population, while non-interactors **deplete**.
3. By sequencing the population before and after selection, interactors are identified by their enrichment -- even extremely rare ones (e.g., 0.002% of population).
4. **Differential comparisons** between baits (e.g., GTP-locked vs. GDP-locked Rab5) reveal conformation-specific interactions.

---

## Architectural Shift: Original vs. Rewrite

The rewrite changes the fundamental input stage. The original DEEPN required **pre-mapped .sam files** from Tophat2/Bowtie2, splitting reads into mapped (for GeneCount) and unmapped (for JunctionMake). DEEPN++ takes **raw FASTQ.gz files** and handles mapping internally via BLAT/BLAST, combining the junction finding and mapping into a single JunctionDice++ step.

```
 ORIGINAL DEEPN (2016)                    DEEPN++ (Current Rewrite)
 ──────────────────────                    ────────────────────────

 .fastq files                             .fastq.gz files
      |                                        |
      v                                        v
 Tophat2 / Bowtie2                        JunctionDice++
 (external, genome-wide)                  (built-in, junction-first)
      |                                        |
      ├──> mapped .sam                         ├──> .jdice.fasta
      │         |                              │    (cDNA portions only)
      │         v                              │         |
      │    GeneCount                           │    BLAT / BLASTN / megaBLAST
      │    (counts per gene)                   │    (against gene DB, not genome)
      │                                        │         |
      └──> unmapped .sam                       │         v
                |                              │    mapped_files/
                v                              │         |
           JunctionMake                        │         v
           (find junctions)                    │    GeneCount++
                |                              │    (counts per gene)
                v                              │
           BLASTN                              └──> analyzed_files/
           (map junction reads)                     (SQLite depth DBs)
                |                                        |
                v                                        v
           BlastQuery                             DESeq2++
           (junction plots)                       (enrichment statistics)
                |                                        |
                v                                        v
           ReadDepth                              MultiQuery++ (skeleton)
           (3' border mapping)                    ReadDepth++  (skeleton)
                |
                v
           StatMaker
           (R/JAGS Bayesian)
```

**Key implications of this architectural change:**
- No Tophat2/Bowtie2 dependency (major portability improvement)
- Mapping is against a **gene database** (RefSeq mRNAs) rather than the full genome
- Junction finding and mapping happen in a **single pass** per file
- The depth SQLite database is created during dicing (not as a separate step)
- BLAT/BLASTN/megaBLAST are bundled with the application

---

## Pipeline Architecture

```
                         DEEPN++ Pipeline
                    (Competitive Y2H Screen Analysis)

 Sequencing of Y2H population (before/after selective growth)
                              |
                              v
                    Raw FASTQ.gz files
                              |
                              v
 ┌─────────────────────────────────────────────────────────────┐
 │  JunctionDice++                                             │
 │                                                             │
 │  The junction sequence is the Gal4-AD / cDNA fusion point.  │
 │  Reads spanning this junction contain the prey identity.    │
 │                                                             │
 │  1. Scan FASTQ reads for the junction sequence              │
 │     (forward + reverse complement matching)                 │
 │  2. Extract the cDNA portion of junction-spanning reads     │
 │  3. Write extracted sequences to FASTA                      │
 │  4. Map to reference gene database via BLAT/BLASTN/megaBLAST│
 │  5. Create SQLite depth database for downstream ReadDepth++ │
 └──────────────────────────┬──────────────────────────────────┘
                            |
                            v
                   Mapping output files
            (*.blat.txt / *.blastn.txt / *.megablast.txt)
                            |
                            v
 ┌─────────────────────────────────────────────────────────────┐
 │  GeneCount++                                                │
 │                                                             │
 │  Quantifies plasmid abundance = proxy for interaction       │
 │  strength. Each mapped read is a vote for that gene's       │
 │  prey plasmid being present in the population.              │
 │                                                             │
 │  1. Parse BLAT/BLAST output, filter by quality              │
 │     (identity >97%, bitscore >50)                           │
 │  2. Group hits by read, track best-scoring genes            │
 │  3. Count reads per gene -> PPM normalization               │
 │     PPM = (gene_count / total_reads) x 1,000,000           │
 │  4. Store in SQLite, export CSV/XLSX                        │
 └──────────────────────────┬──────────────────────────────────┘
                            |
                            v
                  PPM gene count files
            (per-sample: gene -> abundance)
                            |
                            v
 ┌─────────────────────────────────────────────────────────────┐
 │  DESeq2++                                                   │
 │                                                             │
 │  Statistical comparison: bait vs. vector-alone controls,    │
 │  or bait1 vs. bait2 (e.g., Rab5-GTP vs. Rab5-GDP).        │
 │  Genes significantly enriched with bait = candidate         │
 │  interactors.                                               │
 │                                                             │
 │  1. Library-size normalization (median-of-ratios)           │
 │  2. Negative binomial GLM, dispersion estimation            │
 │  3. Wald test per gene, p-value adjustment (BH)            │
 │  4. Log fold change shrinkage                               │
 │  5. Visualization: MA plot, volcano plot, dispersion plot   │
 └──────────────────────────┬──────────────────────────────────┘
                            |
                            v
              Ranked candidate interactors
         (gene, log2FC, adjusted p-value, base mean)
                            |
              ┌─────────────┴─────────────┐
              v                           v
 ┌────────────────────────┐  ┌────────────────────────┐
 │  MultiQuery++          │  │  ReadDepth++           │
 │  (NOT IMPLEMENTED)     │  │  (NOT IMPLEMENTED)     │
 │                        │  │                        │
 │  See detailed spec     │  │  See detailed spec     │
 │  below from original   │  │  below from original   │
 │  BlastQuery module.    │  │  ReadDepth module.      │
 └────────────────────────┘  └────────────────────────┘
```

---

## Module-by-Module Status

### 1. DEEPN (Main Orchestrator) -- FUNCTIONAL

| Aspect | Status |
|--------|--------|
| UI | Complete -- database selector, folder picker, junction/exclude inputs, action buttons |
| Database loading | Complete -- reads `deepn.json` for reference genomes |
| Directory monitoring | Complete -- QFileSystemWatcher tracks 6 subdirectories |
| Button state logic | Complete -- enables/disables based on available files |
| Child process launching | Complete -- spawns JunctionDice++, GeneCount++, etc. |
| Progress display | Complete -- DQTextEdit custom widget for real-time status |

**Working directory layout managed:**
- `fastq/` -- raw sequencing reads from Y2H population
- `junction_diced_fasta/` -- extracted cDNA sequences from junction-spanning reads
- `query_files/` -- BLAST queries (*.bq)
- `mapped_files/` -- alignment output (*.megablast.txt, *.blat.txt, *.blastn.txt)
- `gene_count_summary/` -- per-sample gene abundance in PPM (*.csv)
- `analyzed_files/` -- SQLite databases and downstream results (*.rd, *.sqlite)

**Comparison with original:** The original DEEPN expected `mapped_sam_files/` and `unmapped_sam_files/` directories (pre-mapped by Tophat2). DEEPN++ starts from raw FASTQ and creates its own directory structure.

### 2. JunctionDice++ -- FUNCTIONAL (differences from original noted)

Replaces the original JunctionMake module. Scans sequencing reads for the Gal4-AD/cDNA fusion junction, extracts the cDNA portion, and maps it to identify which prey gene is present.

| Aspect | Status |
|--------|--------|
| FASTQ.gz reading | Complete -- zstr/zlib decompression |
| Junction detection | Complete -- single-pass regex with N-tolerant matching + reverse complement |
| Sequence translation | Complete -- DNA to protein (for reading frame validation) |
| Repeat filtering | Complete -- regex for 20+ repeat bases |
| FASTA output | Complete -- writes `.jdice.fasta` |
| Mapping dispatch | Complete -- launches BLAT, BLASTN, or megaBLAST against reference gene DB |
| Depth DB creation | Complete -- SQLite with reads/maps tables (WAL mode) for downstream ReadDepth++ |
| Statistics tracking | Complete -- DiceStat + MapStat structs, JSON summary |
| Threading | Complete -- JDWorker on QThread per file |

**Differences from original JunctionMake:**

| Feature | Original JunctionMake | JunctionDice++ |
|---------|----------------------|----------------|
| Input | Unmapped .sam files | Raw FASTQ.gz files |
| Junction search | **Three-tier offset search**: primary (positions 9-23 of 50bp seq), secondary (offset -4), tertiary (offset -8). Each extracts a 15bp tag. Falls back progressively. | **Single regex pass**: takes rightmost `matchLength` chars of junction, builds N-tolerant pattern with wildcard last 3 chars |
| Non-matching reads | Discarded (only junction reads output) | **Also written to FASTA** (lines 260-264 in jdworker.cpp) -- all reads go to output regardless of junction match |
| Read limit | None | **Hard-coded 1M read cap** (`jdworker.cpp:286`: `if (stat->dstat.totalReads > 1000000) break; // TODO: REmove this line`) |
| Junction sequences | Library-specific presets: mouse cDNA (Clontech pGADT7) and yeast genomic (Phil James pGAD-C1/C2/C3) with specific 50bp sequences | Configurable via UI input field |
| Mapping | Separate BLASTN step after junction extraction | Integrated -- BLAT/BLASTN/megaBLAST dispatched directly after dicing |
| Output format | `.junctions.txt` + `.blast.txt` + `.p` (Python dict for BlastQuery) | `.jdice.fasta` + mapping output + `.sqlite` depth DB |

**Action items for JunctionDice++:**
- **Remove the 1M read limit** (`jdworker.cpp:286`) -- this is a leftover debug cap
- **Consider implementing three-tier search** to match original's robustness against cloning mismatches
- **Filter non-junction reads from output** -- currently all reads are written, which inflates the FASTA and wastes mapping time
- `static QRegularExpression` on line 200-203 means all workers in the process share one compiled pattern; this is fine as long as all workers use the same junction sequence (they do, currently)

### 3. GeneCount++ -- FUNCTIONAL (bug fix pending)

Replaces the original GeneCount module. Counts how many sequencing reads map to each gene, which is a direct proxy for that prey plasmid's abundance in the Y2H population.

| Aspect | Status |
|--------|--------|
| BLAT/BLAST parsing | Complete -- MapHit class parses tab-delimited output |
| Quality filtering | Complete -- identity >97%, bitscore >50 |
| Read grouping | Complete -- ReadHits groups multi-hit reads, tracks best bitscore |
| Reading frame tracking | Complete -- computes frame (0/1/2 fwd/rev) and location (upstream/in_orf/downstream) |
| SQLite storage | Complete -- in-memory + on-disk databases, batch inserts with index |
| PPM normalization | Complete -- (raw_count / total_reads) x 1,000,000 |
| CSV/XLSX output | Complete -- QSimpleXlsxWriter |
| Threading | Complete -- GCWorker per file |
| Multi-file parallel processing | **Bug fix pending** -- hardcoded SQLite connection names caused crashes. Fix: unique per-worker connection names via QUuid + moved DB setup to worker thread. |

**Differences from original GeneCount:**

| Feature | Original GeneCount | GeneCount++ |
|---------|-------------------|-------------|
| Input | Mapped .sam files (from Tophat2) | BLAT/BLAST tab-delimited output (from JunctionDice++) |
| Output | `_summary.csv` (chromosome, gene, PPM, NM accessions) + `_ChrGene.csv` (chromosome mapping) | SQLite database + CSV/XLSX |
| Metrics | TotalReads and TotalReads(PPM) separately tracked | Single read count |
| Chromosome info | Outputs per-chromosome gene assignment | Not tracked (maps to gene DB, not genome) |

**Note:** The distinction between TotalReads (all mapped reads) and TotalReads(PPM) (reads to known exons) from the original is less relevant in DEEPN++ since mapping is against a gene database rather than the full genome.

### 4. DESeq2++ (Statistics + UI) -- FUNCTIONAL

**Added 2026-03-24.** Replaces the original R/JAGS-based StatMaker program.

#### Statistics Library (`deseq2/statistics/`)

| Component | Status |
|-----------|--------|
| DeseqDataSet | Complete -- size factors, genewise/fitted/MAP/final dispersions, LFC, Cooks distance, refit |
| DeseqStats | Complete -- Wald test, Cooks filtering, independent filtering, p-value adjustment, LFC shrinkage |
| Utils | Complete -- fitAlphaMLE, IRLS solver, rough dispersions, PPM-to-DESeq2 conversion, CSV I/O |
| Tests | Present -- 4 test files (example usage, PPM analysis, unit tests, conversion tests) |

**Dependencies:** Eigen3 (linear algebra), Boost (filesystem)

#### UI Application (`deseq2/ui/`)

| Component | Status |
|-----------|--------|
| Input tab | Complete -- PPM file loading, group assignment (bait vs. control), data generation |
| Analysis tab | Complete -- run/stop/reset with progress bar, AnalysisWorker on QThread |
| Results tab | Complete -- sortable table, p-value/log2FC filtering, statistics summary |
| Visualization tab | Complete -- MA plot, volcano plot, dispersion plot via QCustomPlot |
| GeneCountHandler | Complete -- parses PPM files from GeneCount++ output, auto/manual group assignment, count matrix generation |
| Export | Complete -- CSV export of results |

**Differences from original StatMaker:**

| Feature | Original StatMaker | DESeq2++ |
|---------|-------------------|----------|
| Platform | Mac OS X only | Cross-platform (Qt5) |
| Statistical method | Bayesian MCMC via R + JAGS; posterior probability of enrichment | Frequentist: negative binomial GLM, Wald test, BH-adjusted p-values |
| Dependencies | R, JAGS, Bioconductor, deepn R package | Eigen3, Boost (bundled, no external runtime deps) |
| Comparison modes | Explicit **two-way** (vector vs. 1 bait) and **three-way** (vector vs. bait1 vs. bait2) with drag-and-drop file assignment | General design matrix (supports arbitrary contrasts) |
| PPM threshold | Configurable (default 3 PPM); genes below threshold excluded pre-analysis | Filtering via independent filtering step in DeseqStats |
| File input | Drag-and-drop GeneCount CSVs into Vector/Bait/Selected/Non-Selected slots | Tab-based PPM file loading with group assignment |
| Output | `statmaker_output.csv` with p-values 0-1 | CSV with base mean, log2FC, SE, Wald stat, p-value, adjusted p-value |
| Visualization | None (external tools like GraphPad Prism) | Built-in MA plot, volcano plot, dispersion plot (QCustomPlot) |

**Note:** The original StatMaker's drag-and-drop interface for assigning files to Vector/Bait categories was very user-friendly for the Y2H workflow. DESeq2++'s group assignment serves the same purpose but could be enhanced with explicit "Vector Control" / "Bait" / "Selected" / "Non-Selected" labeling to match the Y2H mental model.

### 5. MultiQuery++ -- SKELETON ONLY

**Must implement the original BlastQuery module's functionality.** This is the primary tool researchers use to inspect individual candidate genes after statistical ranking.

#### Original BlastQuery Specification (from User Guide)

**Input:** `.p` files (Python dictionaries) from JunctionMake's `blast_results_query/` folder. In DEEPN++, the equivalent data source is the SQLite depth databases created by JunctionDice++.

**Query:** User enters an NCBI reference number (NM_*** identifier) found in GeneCount summary CSV files. For yeast, hybrid nomenclature: `YDR182W_CDC1`.

**Four-tab interface:**

**Tab 1: Results**
| Column | Description |
|--------|-------------|
| Position | First nucleotide position of the identified insert in the matched gene |
| #Junctions | Abundance count in PPM for that junction site |
| QueryStart | From BLASTN `q.start`; nucleotides between junction tag and matched cDNA position. Usually 1, but cloning artifacts increase this, affecting reading frame calculation |
| CDS | Classification: "Upstream of coding region", "Within coding region (In ORF)", or "Downstream of coding region" |
| Frame | Reading frame: "In-frame" (correct), "Intron" (genomic library, intron-interrupted), or "Backwards" (reverse orientation, QueryStart flagged as 1000) |

Position + QueryStart values jointly determine CDS location and reading frame.

**Tab 2: Filtered Results**
- Collapses junctions with same Position but different QueryStart values
- Shows only Positions present in the leftmost (first) dataset
- Purpose: compare unselected library vs. selected population at matching positions

**Tab 3: Plot**
- **Y-axis:** Junction abundance (PPM)
- **X-axis:** Position along mRNA/gene sequence (nucleotide coordinates)
- **Bar colors:**
  - **Dark blue:** Within CDS and in correct reading frame
  - **Cyan:** Upstream of CDS start but correct translational reading frame
  - **Grey:** Downstream of CDS stop, or out-of-frame
  - **Red bars:** Translation start (ATG) and stop codon positions
- **Below plot:** mRNA/gene sequence text with black text for protein coding sequence and grey text for untranslated regions
- **Multiple datasets:** Compare unselected vs. selected populations side-by-side

**Tab 4: Export**
- "Save CSV" button exports Results and Filtered Results tables
- Compatible with Excel and GraphPad Prism

#### Current MultiQuery++ State

| Aspect | Status |
|--------|--------|
| UI framework | Minimal -- MainWindow + stackable DockWidget scaffold |
| Dataset selection dropdown | Not implemented |
| NM_*** query input | Not implemented |
| Results table (Position, #Junctions, QueryStart, CDS, Frame) | Not implemented |
| Filtered Results table | Not implemented |
| Junction position plot with color-coded bars | Not implemented |
| mRNA sequence display with CDS/UTR coloring | Not implemented |
| CSV export | Not implemented |
| Qt Charts | Imported but unused |

### 6. ReadDepth++ -- SKELETON ONLY

**Must implement the original ReadDepth module's functionality.** This module determines the 3' end of cDNA inserts, allowing reconstruction of the prey plasmid.

#### Original ReadDepth Specification (from User Guide)

**Purpose:** Identify the 3' boundary of cDNA inserts by measuring read depth across the gene's sequence in fixed intervals.

**Method:**
1. For a queried gene, find all mapped reads corresponding to that gene
2. Divide the gene's cDNA sequence into **25bp intervals**
3. Count how many reads contain each 25bp interval
4. Plot read depth (y-axis) vs. position along cDNA (x-axis)

**Parameters:**
- Interval distance adjustable via up/down arrows in **50 nucleotide increments**
- Direct numerical entry allowed (minimum 50 nt intervals)

**Output:** Graph showing read coverage along the full cDNA length. For enriched prey plasmids:
- Read depth is high from the 5' junction site through the insert
- Read depth drops to zero at the **3' end of the cDNA insert**
- This drop-off reveals the 3' boundary, completing the prey plasmid reconstruction (5' junction from BlastQuery + 3' border from ReadDepth)

**Data source:** Mapped .sam files (original) or the SQLite depth databases from JunctionDice++ (DEEPN++)

#### Current ReadDepth++ State

| Aspect | Status |
|--------|--------|
| UI framework | Minimal -- MainWindow shell |
| Gene query input | Not implemented |
| 25bp interval counting | Not implemented |
| Interval distance adjustment (50 nt increments) | Not implemented |
| Read depth vs. position plot | Not implemented |
| 3' boundary detection | Not implemented |
| Qt Charts | Imported but unused |

---

## Biological Context & Key Concepts

### The Junction Sequence

The junction sequence is the 50bp DNA boundary between the **Gal4 Activation Domain (AD)** coding sequence and the **cDNA insert** in the prey plasmid. The last 3 nucleotides of the junction define the reading frame for the downstream cDNA.

**Library-specific junction sequences (from original user guide):**

| Library | Junction Sequence (50bp) |
|---------|-------------------------|
| Mouse cDNA (Clontech pGADT7 Mate/Plate) | `AATTCCACCCAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGG` |
| Yeast Genomic (Phil James pGAD-C1/C2/C3) | `ATGATGAAGATACCCCACCAAACCCAAAAAAAGAGATCGAATTCCCGGGG` |
| hORFeome (from deepn++ config) | `CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC` |

**Requirements for custom junction sequences:**
- Must be exactly **50 nucleotides**
- Must be in **UPPERCASE**
- Positioned immediately upstream of the cDNA/fragment fusion site
- Last 3 nucleotides must define a complete codon for the preceding reading frame

### Three-Tier Junction Search (Original, not yet in DEEPN++)

The original JunctionMake used a three-tier fallback search to compensate for cloning errors:

1. **Primary search:** Extract 15bp tag from positions 9-23 of the 50bp junction. Search reads for this tag. Extract everything downstream as the "Downstream Reading Frame."
2. **Secondary search (offset -4):** If primary fails, search 15bp from 4 nucleotides upstream. Same downstream reading frame.
3. **Tertiary search (offset -8):** If secondary fails, search 15bp from 8 nucleotides upstream. Same downstream reading frame.

All three tiers extract the same reading frame -- the offsets handle base substitutions/insertions near the junction.

JunctionDice++ currently uses a **single regex pass** with N-tolerant matching on the rightmost `matchLength` characters. This is simpler but may miss junctions that the three-tier approach would catch.

### PPM (Parts Per Million) Normalization

PPM = (gene_read_count / total_reads_in_sample) x 1,000,000

This normalizes for sequencing depth. A gene at 100 PPM in the unselected population that rises to 10,000 PPM after selection is a strong interaction candidate.

The original distinguished **TotalReads** (all mapped reads in the .sam file) from **TotalReads(PPM)** (reads mapping to known exons), since genome-wide mapping produces reads in non-exonic regions. DEEPN++ maps against a gene database so this distinction is less relevant.

### Overdispersion in Y2H Count Data

Under non-selective conditions, count variability is near-Poisson (phi=0.06). Under selection, overdispersion jumps to phi=3.75 because interacting plasmids amplify at different rates. This is why the negative binomial model (DESeq2) is appropriate.

### Reading Frame and CDS Classification

Not all junction fusions produce functional proteins. The original BlastQuery classified each junction by:
- **CDS location:** Upstream / In ORF / Downstream of coding sequence
- **Frame:** In-frame / Out-of-frame / Backwards (reverse orientation) / Intron (genomic libraries)
- **QueryStart:** Distance from junction tag to matched cDNA position (usually 1; larger values indicate cloning artifacts that shift the reading frame)

The paper discovered that even out-of-frame fusions can produce signal (e.g., RabEP1/Rabaptin via low-level frameshifting), so frame classification is informative but not strictly filtering.

### Yeast Library Considerations

Yeast genomic libraries (pGAD-C1, C2, C3) use three separate vectors with insertions/deletions at the cloning site to sample all three reading frames. The original JunctionMake used a consensus junction sequence and adjusted tags accordingly. cDNA libraries (mouse/human) have directional cloning so "Backwards" inserts are rare; genomic libraries frequently show reverse-orientation inserts and intron interruptions.

---

## Build System Status

### qmake (Primary -- Original Modules)

The root `deepn++.pro` builds: deepn, junction_dice, gene_count, query, read_depth.

- **macOS:** Targets macOS 14.0+, Qt5 static, C++20
- **Linux:** Docker-based (Ubuntu 20.04 + Qt 5.15.2 static)
- **Windows:** Docker-based (MXE MinGW cross-compilation)

### CMake (DESeq2 Module)

Separate CMake build for `deseq2/statistics` and `deseq2/ui`:
- C++17 standard
- Dependencies: Qt5 (static), Eigen3, Boost, QCustomPlot, QXlsx
- Builds `deseq2_statistics` static library + `DESeq2++` UI executable

### CI/CD (GitLab)

4-stage pipeline: base -> docker -> build -> publish

| Platform | Method | Status |
|----------|--------|--------|
| macOS | Dedicated runner, qmake | Configured |
| Linux | Docker container, qmake | Configured |
| Windows | Docker container, MXE qmake | Configured |
| Release | 7z archives per platform | Configured |

**Note:** CI only covers the qmake modules. The CMake-based DESeq2 module is not yet integrated into the CI pipeline.

### Submodules (6)

| Submodule | Purpose | Status |
|-----------|---------|--------|
| samtools | BAM/SAM file handling (htslib) | Configured |
| statgen (libStatGen) | Statistical genetics library | Configured |
| ncbi/zlib | NCBI's zlib copy | Configured |
| zlib | Compression library | Built (confirmed in zlib-build.log) |
| blat | BLAT sequence alignment tool | Precompiled macOS binary present |
| xlsx (QSimpleXlsxWriter) | Excel output | Configured |

### Conan Dependency

- `seqan/2.4.0` -- sequence analysis library (used by JunctionDice++)

---

## Known Issues

### Critical (Affects Correctness)

1. **JunctionDice++ writes non-junction reads to output** -- `jdworker.cpp:260-264`: reads that do NOT match the junction pattern are still written to the FASTA output. This means non-junction reads get mapped and counted, inflating gene counts with noise. The original JunctionMake only output junction-spanning reads.

2. **JunctionDice++ has a hard-coded 1M read limit** -- `jdworker.cpp:286`: `if (stat->dstat.totalReads > 1000000) break; // TODO: REmove this line`. This debug cap will silently truncate any real dataset.

### Build Issues

3. **NCBI BLAST build incomplete** -- `ncbi-build.log` shows missing `Makefile.mk`. The `makeblastdb` and `makembindex` tools are not being compiled. BLAT works as a fallback.

4. **DESeq2 not in CI** -- The CMake-based DESeq2 module has no corresponding stages in `.gitlab-ci.yml`.

5. **Build system fragmentation** -- Two parallel build systems (qmake for original modules, CMake for DESeq2). No unified top-level build.

6. **Qt5 static path hardcoded** -- CMake references `/Users/venky/Softwares/qt5-static`.

### Runtime Issues

7. **GeneCount++ multi-file crash** -- Hardcoded SQLite connection names caused concurrent workers to overwrite each other's connections. **Fix implemented, pending build verification.**

### Incomplete Modules

8. **MultiQuery++** -- UI skeleton only. See detailed spec above from original BlastQuery module.

9. **ReadDepth++** -- UI skeleton only. See detailed spec above from original ReadDepth module.

### Code Quality

10. **PythonQt integration abandoned** -- Dead code in `backup/pythonqt/`.

11. **No automated tests in CI** -- Test files exist for DESeq2 but no test stage in the pipeline.

### Design Gaps vs. Original

12. **No three-tier junction search** -- JunctionDice++ uses single-pass regex vs. original's three-tier offset strategy. May miss junctions with cloning mismatches.

13. **No chromosome-level output** -- GeneCount++ doesn't output `_ChrGene.csv` (chromosome mapping). Less relevant since mapping is against gene DB, but may be useful for some analyses.

14. **No explicit Y2H comparison UI** -- Original StatMaker had drag-and-drop slots for "Vector Control" / "Bait" / "Selected" / "Non-Selected". DESeq2++ uses generic group assignment which is more flexible but less intuitive for the Y2H workflow.

---

## Evolution from Original DEEPN

| Component | Original (2016) | DEEPN++ (Current) | Parity |
|-----------|-----------------|-------------------|--------|
| Language | Python/R scripts | C++20 / C++17 | Rewritten |
| GUI | Per-module GUIs | Qt5 desktop app suite | Improved |
| Input | Pre-mapped .sam files (Tophat2 required) | Raw FASTQ.gz (self-contained) | **Improved** |
| Junction search | Three-tier offset (15bp tag, offsets -4, -8) | Single regex with N-tolerant matching | **Simplified** (see issue #12) |
| Junction output | Junction reads only | All reads (junction + non-junction) | **Bug** (see issue #1) |
| Read limit | None | 1M debug cap left in | **Bug** (see issue #2) |
| Gene counting | Mapped .sam -> PPM + chromosome info | BLAST output -> PPM via SQLite | Comparable |
| Statistics | StatMaker (R + JAGS, Bayesian MCMC) | DESeq2++ (Eigen3, frequentist NB-GLM) | **Improved** (no R/JAGS deps) |
| StatMaker UI | Drag-and-drop Vector/Bait/Selected slots | Generic group assignment | Different approach |
| Junction viz | BlastQuery (4-tab: Results, Filtered, Plot, Export) | MultiQuery++ (skeleton) | **Not implemented** |
| Read depth | ReadDepth (25bp interval counting, adjustable) | ReadDepth++ (skeleton) | **Not implemented** |
| Supported genomes | mm10, hg38, sacCer3 | mm10, hg38, sacCer3 | Same |
| Library junctions | Mouse cDNA (Clontech), Yeast (pGAD-C1/C2/C3) presets | Configurable via UI | More flexible |
| Platforms | Mac OS X 10.10+, Windows 7+ (32-bit) | macOS 14+, Linux, Windows (64-bit) | **Improved** |
| Distribution | GitHub release | GitLab CI -> 7z per platform | **Improved** |
| External deps | Tophat2, Bowtie2, Samtools, R, JAGS, Bioconductor | BLAT/BLAST (bundled) | **Greatly reduced** |

---

## Development Timeline

| Date | Milestone |
|------|-----------|
| 2016-09 | Original DEEPN paper published (Piper et al., Cell Reports) |
| 2020-09 | DEEPN++ project initialized, Docker CI/CD setup |
| 2020-09 to 2021-03 | Docker pipeline iterations (Linux, Windows builds) |
| 2021-03 | Docker user permissions fix |
| 2022-03 | Python integration exploration (pybind11 -> PythonQt) |
| 2023-02 | samtools/statgen integration, GeneCount++ working |
| 2024-04 | Project rename from deepn-plus to deepn++ |
| 2026-03-24 | **Major update:** CMake build system, DESeq2++ module, code quality rules, project restructuring |
| 2026-03-26 | GeneCount++ multi-file crash diagnosed and fixed; project status assessment |

---

## Reference Data

Pre-bundled reference genomes and gene databases:

| Database | Organism | Use Case | Files |
|----------|----------|----------|-------|
| hg38 | Human | hORFeome Y2H libraries | GeneList.prn, exonDict.p, BLAST DBs |
| mm10 | Mouse | Mouse cDNA Y2H library (as in the paper) | GeneList.prn, exonDict.p, BLAST DBs |
| sacCer3 | Yeast (S. cerevisiae) | Yeast interaction screens | GeneList.prn, exonDict.p, BLAST DBs |

**Database contents (from user guide):**
- Only annotated mRNAs with NM_* nomenclature (no microRNAs, lncRNAs, or theoretical splice variants)
- Yeast genes: protein coding sequence plus 100bp untranslated regions on each end
- Yeast nomenclature: hybrid format combining SGD systematic names with GenBank common names (e.g., `YDR182W_CDC1`)
- Built-in BLASTN databases for all three organisms

Configuration in `deepn.json` and `inputs/genomes.yml`.

---

## Technology Stack

| Layer | Technology |
|-------|-----------|
| Language | C++20 (qmake modules), C++17 (DESeq2) |
| GUI Framework | Qt 5.15+ (static builds) |
| Build Systems | qmake, CMake, Conan |
| Linear Algebra | Eigen3 |
| Plotting | QCustomPlot (DESeq2), Qt Charts (Query/ReadDepth) |
| Spreadsheet | QSimpleXlsxWriter, QXlsx |
| Database | SQLite (in-memory + on-disk, WAL mode) |
| Sequence Analysis | SeqAn 2.4.0, BLAT, NCBI BLAST+ |
| Compression | zlib / zstr |
| BAM Handling | samtools / htslib |
| CI/CD | GitLab CI, Docker |
| Platforms | macOS 14+, Ubuntu 20.04, Windows (MinGW) |

---

## What's Next (Suggested Priorities)

### Critical Fixes
1. **Remove the 1M read debug cap** in `jdworker.cpp:286` -- silently truncates real data
2. **Stop writing non-junction reads** in `jdworker.cpp:260-264` -- inflates gene counts with noise

### Pipeline Completion
3. **Implement MultiQuery++** -- Four-tab interface matching original BlastQuery: Results table (Position, #Junctions/PPM, QueryStart, CDS, Frame), Filtered Results, color-coded junction position plot (dark blue=in CDS+in-frame, cyan=upstream+correct frame, grey=downstream/out-of-frame, red=start/stop codons), mRNA sequence display, CSV export
4. **Implement ReadDepth++** -- 25bp interval read depth counting with adjustable interval distance (50 nt increments), plot read depth vs. position to identify 3' insert boundaries

### Enhancements
5. **Consider three-tier junction search** -- to match original's robustness against cloning mismatches
6. **Add Y2H-specific labels to DESeq2++ UI** -- "Vector Control" / "Bait" / "Selected" / "Non-Selected" to match the mental model from original StatMaker

### Build & Infrastructure
7. **Unify build system** -- single CMake or qmake build for all modules
8. **Fix NCBI BLAST build** -- resolve missing Makefile
9. **Add DESeq2 to CI pipeline**
10. **Add test stage to CI**

### Cleanup
11. **Remove dead PythonQt code** in `backup/pythonqt/`

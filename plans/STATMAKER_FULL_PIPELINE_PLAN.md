# StatMaker Full Pipeline Plan

**Status: COMPLETE as of 2026-04-01**
All 5 phases implemented. 29/29 tests pass. Additional correctness fixes applied per paper review (see plans/Y2H_CORRECTNESS_FIXES.md).

---

**Original Goal:** Run the full DESeq2 + Y2H-SCORES pipeline in one StatMaker analysis pass, write merged results to a dedicated SQLite database under `analyzed_files/`, and fix data-integrity issues such as corrupted gene names in the results table.

---

## 1. Desired Outcome

StatMaker should become the single analysis GUI for:
- DESeq2 differential enrichment
- Y2H-SCORES enrichment
- Y2H-SCORES specificity
- Y2H-SCORES in-frame scoring
- Borda aggregation

One analysis run should produce:
- an in-memory results model for the UI
- a merged CSV export
- a merged SQLite database at:
  - `<workdir>/analyzed_files/statmaker_results.sqlite`

If `<workdir>/analyzed_files/` does not exist, StatMaker should create it automatically.

---

## 2. Current Gaps

### 2.1 Full pipeline is not yet running

Current worker behavior in [`deseq2/ui/src/main_window.cpp`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp):
- runs DESeq2
- runs Y2H enrichment scoring
- runs Borda aggregation with enrichment-only input

Missing from the actual GUI pipeline:
- specificity scoring
- in-frame scoring from junction SQLite inputs
- merged full-pipeline SQLite schema
- auto-targeting output into `analyzed_files/`

### 2.2 Output path is wrong for the final product

Current SQLite output is written near the first input file as `deseq2_results.sqlite`.

That is not the right long-term layout because:
- it mixes analysis output with input files
- it uses a DESeq2-specific name for a StatMaker result
- it is inconsistent with the DEEPN++ working-directory model

### 2.3 Gene-name corruption is an active correctness risk

The results table sometimes shows binary-looking garbage in gene names. This must be treated as a correctness issue before any downstream export is trusted.

Probable root cause:
- `AnalysisWorker::runAnalysis()` emits a pointer to a stack-local `AnalysisResults`
- that pointer crosses a thread boundary
- the receiving slot copies from memory that may already be invalid

This is high priority because it can silently corrupt:
- the UI
- CSV exports
- SQLite exports
- downstream navigation into MultiQuery++ / ReadDepth++

---

## 3. Target Pipeline

### 3.1 Input Discovery

On launch with a working directory:
1. Ensure `<workdir>/analyzed_files/` exists.
2. Discover GeneCount outputs from the configured input directory.
3. Discover relevant junction/depth SQLite files for in-frame scoring.
4. Load previous `statmaker_results.sqlite` if present.

### 3.2 Analysis Run

One StatMaker run should execute:
1. DESeq2 normalization, dispersion fitting, LFC fitting, and statistical tests
2. Enrichment scoring from DESeq2 results
3. Specificity scoring from pairwise bait contrasts
4. In-frame scoring from junction/depth SQLite data
5. Borda aggregation across all available Y2H metrics
6. Merge all columns into one canonical result model
7. Write merged output to CSV and `analyzed_files/statmaker_results.sqlite`

### 3.3 Canonical Output Location

Recommended layout:
```text
<workdir>/
├── gene_count_summary/
├── analyzed_files/
│   ├── sample_a.sqlite
│   ├── sample_b.sqlite
│   └── statmaker_results.sqlite
```

The analysis database should be distinct from per-sample depth databases.

---

## 4. SQLite Output Plan

### 4.1 Database Name

Use:
- `statmaker_results.sqlite`

Do not keep using:
- `deseq2_results.sqlite`

### 4.2 Schema

Recommended primary table:

```sql
CREATE TABLE analysis_results (
    gene TEXT PRIMARY KEY,
    bait TEXT,
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
    in_frame_transcripts TEXT,
    active_contrast_label TEXT,
    created_at TEXT
);
```

Recommended supporting tables:
- `analysis_summary`
- `contrasts`
- `run_parameters`
- `input_files`

### 4.3 Write Policy

Write behavior should be:
- create `analyzed_files/` if missing
- open `statmaker_results.sqlite`
- replace prior result tables for the current run
- wrap inserts in a transaction
- include run metadata so results are reproducible

---

## 5. Gene Name Integrity Plan

### 5.1 Immediate Fix

Replace cross-thread pointer passing with one of:
- queued signal carrying `AnalysisResults` by value
- `std::shared_ptr<AnalysisResults>` with explicit ownership
- heap allocation owned and freed by the receiver

The current stack-pointer approach should be removed first.

### 5.2 Validation Layer

Before results are shown or persisted:
- verify `geneNames.size() == results.rows()`
- reject or log names with embedded NULs or non-printable control bytes
- verify `QString::fromStdString()` does not produce malformed display text

### 5.3 Export Consistency Tests

Add a fixture that validates the same gene names through:
- worker output
- UI table
- CSV export
- SQLite export

---

## 6. Product Rename Plan

All user-visible naming should change from DESeq2++ to StatMaker.

Rename scope:
- window title
- application name in menus and dialogs
- bundle/executable naming
- output file names
- progress text
- documentation and plans

Non-goal for the first pass:
- renaming the source directory `deseq2/`

That repository-level rename can happen later once build and packaging references are stabilized.

---

## 7. Implementation Phases

### Phase 1: Correctness First

- Fix `AnalysisResults` lifetime across threads
- Verify gene-name integrity in UI and exports
- Add regression coverage for corrupted-name scenarios

### Phase 2: Output Layout

- Resolve working directory consistently
- autocreate `<workdir>/analyzed_files/`
- write merged results to `statmaker_results.sqlite`
- keep sample depth databases untouched

### Phase 3: Full Y2H-SCORES Wiring

- wire specificity scoring into the worker pipeline
- wire in-frame scoring into the worker pipeline
- pass actual UI thresholds into the scorers
- merge full Y2H output into `AnalysisResults`

### Phase 4: Unified StatMaker UI

- rename visible DESeq2 labels to StatMaker
- expose merged results columns where useful
- show which Y2H metrics were actually computed in the run

### Phase 5: Persistence and Reload

- load `statmaker_results.sqlite` on startup when present
- retain active contrast metadata
- make downstream tools consume the new result file name

---

## 8. Acceptance Criteria

- One StatMaker run executes DESeq2 and all intended Y2H-SCORES stages
- `<workdir>/analyzed_files/` is created automatically if absent
- merged results are written to `<workdir>/analyzed_files/statmaker_results.sqlite`
- gene names are stable and readable in the UI, CSV, and SQLite
- user-visible naming consistently says StatMaker
- logs clearly state which Y2H metrics ran and which were skipped

---

## 9. Dependencies

- final decision on which SQLite tables downstream tools will read
- junction/depth SQLite discovery rules for in-frame scoring
- a safe ownership model for cross-thread analysis results

---

## 10. Recommended Order

1. Fix gene-name corruption and result ownership.
2. Move output into `analyzed_files/statmaker_results.sqlite`.
3. Wire the missing Y2H specificity and in-frame stages.
4. Rename the product surface to StatMaker.
5. Update downstream readers to prefer the new SQLite file.

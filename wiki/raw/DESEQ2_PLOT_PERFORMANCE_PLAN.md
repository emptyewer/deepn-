# StatMaker Plot Performance Plan

**Area:** `deseq2/ui` visualization path for the unified StatMaker GUI  
**Primary Goal:** Make plot interactions feel responsive on large result sets without regressing export quality or three-way comparison support.  
**Secondary Goal:** Clarify the current `Y2H-SCORES` execution path and keep it from becoming an accidental source of UI latency.

---

## 1. Current State

The active plotting implementation is in [`deseq2/ui/src/main_window.cpp`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp). The current application should be treated as **StatMaker** at the product level, even though the source tree still uses the `deseq2/` path. The UI no longer uses the archived R StatMaker plotting code during normal execution; plot rendering is currently handled with `Qt Charts`.

Current behavior:
- `refreshPlot()` rebuilds a plot from scratch every time the plot type changes or the visualization tab is opened.
- Each plot allocates a new `QChart`, new series objects, new axes, and repopulates all points one by one.
- MA, volcano, dispersion, and three-way scatter all iterate the full result matrix on the UI thread.
- The results table is also fully repopulated and auto-resized after analysis, which can compound the perception that plotting is slow.

Relevant hotspots:
- [`deseq2/ui/src/main_window.cpp#L2027`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L2027)
- [`deseq2/ui/src/main_window.cpp#L2045`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L2045)
- [`deseq2/ui/src/main_window.cpp#L2108`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L2108)
- [`deseq2/ui/src/main_window.cpp#L2216`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L2216)
- [`deseq2/ui/src/main_window.cpp#L2293`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L2293)
- [`deseq2/ui/src/main_window.cpp#L1642`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L1642)

---

## 2. Root Causes

### 2.1 Full chart reconstruction on every refresh
Every plot call creates a new `QChart` and replacement series instead of reusing a persistent chart/series model. This guarantees repeated allocation, layout, legend rebuilds, and axis attachment work.

### 2.2 Point-by-point append on the UI thread
Large result sets are appended with repeated `series->append(x, y)` calls. That is one of the slowest ways to feed `Qt Charts` scatter series at scale.

### 2.3 No downsampling or viewport-aware rendering
The plots attempt to render every gene even though the widget is typically only a few hundred pixels wide. For tens of thousands of rows, much of that work is visually redundant.

### 2.4 Expensive refresh triggers
`refreshPlot()` is called:
- after analysis completes
- when switching to the visualization tab
- when switching plot type

This means the user can repeatedly pay the full rebuild cost even when underlying data has not changed.

### 2.5 Results table work masks plotting cost
The table path currently:
- creates `QTableWidgetItem` objects for every cell
- calls `resizeColumnsToContents()`
- enables sorting and sorts the full table

That makes the application feel slower overall and can be mistaken for plot slowness.

---

## 3. Y2H-SCORES Status

`Y2H-SCORES` is currently being run, but only partially.

Observed behavior:
- The analysis worker computes enrichment scores after DESeq2 results are produced.
- It also runs Borda aggregation using enrichment-only input.
- Those outputs are stored in `AnalysisResults` and written to SQLite / CSV export.

Evidence:
- [`deseq2/ui/src/main_window.cpp#L258`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L258)
- [`deseq2/ui/src/main_window.cpp#L264`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L264)
- [`deseq2/ui/src/main_window.cpp#L1530`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L1530)
- [`deseq2/ui/src/main_window.cpp#L1957`](/Volumes/Projects/deepn++/deseq2/ui/src/main_window.cpp#L1957)

Current limitations:
- `SpecificityScorer` is present in the statistics library but does not appear to be called.
- `InFrameScorer` is present but does not appear to be called.
- The exposed Y2H settings widgets are only partially wired; the enrichment call currently uses hardcoded thresholds (`1.0`, `0.0`) instead of the UI values.
- Y2H scores are persisted, but they are not added to the visible results table columns.

Conclusion:
- `Y2H-SCORES` is active enough to add compute time after DESeq2.
- It is not yet fully integrated as a complete Y2H ranking pipeline.
- Plot performance work should assume the long-term product is **StatMaker**, not a DESeq2-only tool.

---

## 4. Performance Plan

### Phase 1: Measure Before Changing

Add lightweight timing around:
- `updateResultsTable()`
- `refreshPlot()`
- each individual plot builder
- `writeResultsToSqlite()`
- Y2H enrichment + Borda steps

Implementation notes:
- Use `QElapsedTimer` with progress-log output.
- Record row counts and visible point counts in logs.
- Capture timings for datasets around 1k, 10k, and 50k genes.

Success criteria:
- We know whether the dominant cost is table creation, chart creation, point insertion, or post-analysis Y2H work.

### Phase 2: Remove Avoidable Rebuilds

Refactor visualization state so the chart view owns persistent chart objects or persistent series per plot type.

Changes:
- Keep one chart instance alive per plot type, or keep one chart and replace only data series.
- Track a dirty flag keyed by:
  - active contrast
  - plot type
  - p-value threshold
  - fold-change threshold
  - point size
- Do not redraw when the visualization tab is reopened unless data or display settings changed.

Expected impact:
- Eliminates repeated chart/legend/axis setup for common navigation flows.

### Phase 3: Batch Series Population

Replace repeated `append(x, y)` loops with batched point creation.

Changes:
- Build `QVector<QPointF>` buffers first.
- Use bulk replacement APIs where available instead of one-point-at-a-time insertion.
- Precompute transformed coordinates:
  - `log10(baseMean)`
  - `-log10(pvalue)`
  - contrast-pair scatter coordinates

Expected impact:
- Significant reduction in UI-thread overhead for large plots.

### Phase 4: Add Adaptive Downsampling

Do not render every gene when the number of points is much larger than the plot resolution.

Strategy:
- Keep full-resolution data for export.
- Render a downsampled view for interactive display.
- Start with a simple policy:
  - if points <= 5,000: render all
  - if points > 5,000: bin by x-axis pixel bucket and keep representative points
- Preserve all significant hits, even when non-significant points are downsampled.

Expected impact:
- Keeps volcano and MA plots responsive on large analyses while preserving the biologically important outliers.

### Phase 5: Move Plot Preparation Off the UI Thread

Separate data preparation from widget mutation.

Changes:
- Build plot-ready point buffers in a worker task after analysis completes.
- Cache per-plot datasets inside `AnalysisResults` or a dedicated plot-cache structure.
- Limit UI-thread work to attaching prepared buffers to series and updating axes.

Expected impact:
- Better responsiveness during contrast switches and plot-type changes.

### Phase 6: Fix Table-Side Performance Noise

The table path should be optimized in the same pass because it competes for the same user-perceived latency budget.

Changes:
- Disable sorting while populating.
- Avoid `resizeColumnsToContents()` on every refresh for large tables.
- Consider replacing `QTableWidget` with `QTableView + QAbstractTableModel`.
- Consider lazy pagination or a significance-only default view for very large result sets.

Expected impact:
- Reduces the chance that chart work is blamed for table bottlenecks.

### Phase 7: Make Y2H-SCORES Explicit and Cheap

Treat Y2H scoring as a separately measurable post-analysis stage.

Changes:
- Measure enrichment and Borda timing separately.
- Wire enrichment thresholds from the UI instead of hardcoded values.
- Make `SpecificityScorer` and `InFrameScorer` opt-in until the data plumbing is complete.
- Do not recompute Y2H scores on purely visual actions.

Expected impact:
- Prevents confusion about whether plotting is slow versus post-analysis ranking being slow.

---

## 5. Recommended Implementation Order

1. Add timing instrumentation.
2. Stop unnecessary plot refreshes.
3. Batch point population with cached transformed coordinates.
4. Disable heavy table behaviors during population.
5. Add adaptive downsampling for interactive plots.
6. Move plot-data preparation off the UI thread.
7. Clean up Y2H-SCORES wiring and expose its timing in logs.

This order gives the fastest path to visible improvement with minimal architectural risk.

---

## 6. Acceptance Criteria

For a representative large result set:
- Switching to the Visualization tab should feel near-instant when nothing changed.
- Switching plot types should complete in under 300 ms for cached/downsampled interactive plots.
- First render after analysis should complete in under 1 second for typical large datasets.
- Export should still use full-resolution data.
- Timing logs should clearly separate:
  - DESeq2 analysis time
  - Y2H scoring time
  - results table update time
  - plot preparation time
  - plot attach/render time

---

## 7. Risks

- `Qt Charts` may still be a poor fit for very large scatter plots even after batching.
- If performance remains inadequate after Phases 1 to 4, the fallback should be to replace `Qt Charts` scatter rendering with a lighter custom plotting path or a more performant plotting library.
- Downsampling must preserve significant hits and extreme outliers, or users will lose trust in the plots.

---

## 8. Definition Of Done

This work is done when:
- plot latency is measured and reduced with reproducible numbers
- interactive redraws avoid full recomputation
- large datasets remain usable in the visualization tab
- `Y2H-SCORES` execution status is explicit in logs and code
- the code path clearly separates analysis, ranking, table population, and plotting costs

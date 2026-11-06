# DEEPN++ project wiki

DEEPN++ is a native C++/Qt project for yeast two-hybrid screening analysis. Its
documents cover StatMaker statistical analysis, Y2H-SCORES integration,
MultiQuery++ junction inspection, ReadDepth++ visualization, and DESeq2/PPM
analysis workflows. This wiki preserves the existing plans and guides as source
documents; their descriptions of status and relative paths reflect their
original context. The repository [README](../README.md) remains its entry point.

## Search

Collection: `deepn-wiki`. The orca registers the collection against the shared
checkout's `wiki/` after merge. These commands apply after registration:

```bash
qmd search -c deepn-wiki "StatMaker" --json -n 5
qmd get "qmd://deepn-wiki/raw/STATMAKER_FULL_PIPELINE_PLAN.md:1:40"
```

Start with a literal search, then read bounded sections including frontmatter
where present. Check each document's status and date before using it as current
guidance.

## Structure

| Location | Purpose |
| --- | --- |
| `README.md` | Wiki home, search instructions and source index |
| `raw/` | Original documents, preserved without content changes |
| `log/` | Dated activity records |
| Topic pages (`<topic>.md`) | Maintained knowledge and decisions as they are recorded |

## Raw-document index

| File | What it is |
| --- | --- |
| [DEEPN.pdf](raw/DEEPN.pdf) | Existing DEEPN PDF reference |
| [DESEQ2_INTEGRATION_PLAN.md](raw/DESEQ2_INTEGRATION_PLAN.md) | StatMaker integration and Y2H adaptation plan |
| [DESEQ2_PLOT_PERFORMANCE_PLAN.md](raw/DESEQ2_PLOT_PERFORMANCE_PLAN.md) | StatMaker plot performance plan |
| [MULTIQUERY_PLAN.md](raw/MULTIQUERY_PLAN.md) | MultiQuery++ implementation plan |
| [PROJECT_STATUS.md](raw/PROJECT_STATUS.md) | Existing DEEPN++ project status document |
| [READDEPTH_PLAN.md](raw/READDEPTH_PLAN.md) | ReadDepth++ implementation plan |
| [STATMAKER_FULL_PIPELINE_PLAN.md](raw/STATMAKER_FULL_PIPELINE_PLAN.md) | StatMaker full pipeline plan |
| [Y2H_SCORES_INTEGRATION_PLAN.md](raw/Y2H_SCORES_INTEGRATION_PLAN.md) | Y2H-SCORES integration plan |
| [README_MACOS_BUNDLE.md](raw/README_MACOS_BUNDLE.md) | macOS application bundle guide |
| [COMPLETE_PPM_ANALYSIS_GUIDE.md](raw/COMPLETE_PPM_ANALYSIS_GUIDE.md) | Complete PPM/DESeq2 analysis workflow |
| [EXAMPLE_USAGE.md](raw/EXAMPLE_USAGE.md) | DESeq2 C++ library usage examples |
| [PPM_CONVERSION_GUIDE.md](raw/PPM_CONVERSION_GUIDE.md) | PPM gene count conversion guide |
| [README_PPM_ANALYSIS.md](raw/README_PPM_ANALYSIS.md) | PPM analysis overview and links to its companion guides |

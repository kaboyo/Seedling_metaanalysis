# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Seed/seedling 16S microbiome **meta-analysis** across published studies. Two parallel
implementations of the same config-driven pipeline: a dependency-light Python workflow and a
richer R/Quarto notebook. Both read the same `analysis_config.json`.

## Commands

```bash
# Python pipeline (full run → analysis_results/)
python3 microbiome_analysis.py --config analysis_config.json
python3 microbiome_analysis.py --config analysis_config.smoke.json   # 120 samples, 99 perms

# Tests — use system python3, NOT .venv/bin/python (venv has no pytest)
python3 -m pytest tests/ -q                                          # 19 tests, ~30s

# R/Quarto notebook (the main analysis document)
quarto render phyloseq_analysis.qmd -P config:analysis_config.json
quarto render phyloseq_analysis.qmd -P config:analysis_config.smoke.json

# R phylogenetic/DESeq2 add-on, invoked as a subprocess by the Python pipeline
Rscript --vanilla phyloseq_analysis.R --config analysis_config.json \
    [--filtered-features analysis_results/qc/filtered_feature_ids.csv]

# THE MASTER SCRIPT — whole analysis in one document (pipeline + figures + audit figures)
quarto render seedling_analysis_master.qmd -P config:analysis_config.json
quarto render seedling_analysis_master.qmd -P config:analysis_config.smoke.json   # fast pass
quarto render seedling_analysis_master.qmd -P run_validation:true                 # also rerun the audit

# All 42 analysis figures alone (seconds — reads persisted outputs only)
quarto render figures_all.qmd -P config:analysis_config.json

# Independent manuscript audit: recompute headline stats, then regenerate fig1-fig4
python3 validation/verify_manuscript.py --config analysis_config.json
quarto render manuscript_figures.qmd -P config:analysis_config.json
```

Environment: Python 3.14 (`/usr/local/bin/python3`), R 4.3.2, Quarto 1.7.32.

## Architecture

**Config is the single source of truth.** `analysis_config.json` (full) and
`analysis_config.smoke.json` (fast) differ only in `permutations`, `io_threads`, `max_samples`,
`indicator_permutations`, and `outputs.directory`. Never hardcode paths, thresholds, or grouping
variables into scripts — add them to the config and read them via `get_config_value()` (R) or the
config dict (Python). Paths in the config are resolved relative to the config file's directory.

### There is exactly one `.qmd`: `seedling_analysis_master.qmd`

[seedling_analysis_master.qmd](seedling_analysis_master.qmd) is the whole project in one file —
19,381 lines, 325 chunks, assembled from all 26 notebooks that used to sit around the repo.

| Parts | Source | Executes |
| --- | --- | --- |
| 1 | `phyloseq_analysis.qmd` — the pipeline; writes what Part 2 reads | **yes** |
| 2 | `figures_all.qmd` — all 42 analysis figures | **yes** |
| 3 | `manuscript_figures.qmd` — audit Fig 1–4, plus an optional `run_validation` step | **yes** |
| 4 | `seedling_meta_analysis.qmd` — an earlier, simpler version of Part 1 | no |
| 5–8 | the 19 `joined-scripts/` notebooks (project, beta-diversity, diversity-analysis, reports) | no |
| 9 | `rphylopic`, `shinyphaser`, `presentation` | no |

**Parts 4–9 carry `eval = FALSE` on every chunk, and must keep it.** They are archival code from the
strawberry/EXC projects and earlier iterations: they reference `Data/` paths, `alpha_df`,
`all_dunn_results`, `TreatmentCode.x` and `attach()` on objects that do not exist in this repo.
Removing `eval = FALSE` fails the render immediately. To revive any of it, lift the chunk out and
adapt it to this schema (`run_accession`, `Developmental_Stage`, `Datesets/`) first.

The master is **generated**, not hand-written — by `merge_master.py` then `append_archive.py`
(preserved in `tasks/plans/`). Originals are in `tasks/plans/archived-qmd-2026-09-12/`. Prefer
regenerating over hand-editing 19,381 lines.

Three things about the merged file that look like mistakes and are not:

- **Helpers are defined more than once** (`save_nature_figure`, `theme_nature_base`, `get_palette`,
  `resolve_path`, `figure_size`, `safe_file_label`, `theme_pipeline_classic`, `theme_pipeline_faceted`).
  R takes the most recent definition, so each part runs against the helpers it was written for.
  Collapsing them into one definition will silently change figures — don't, without re-checking each.
- **Chunk labels carry prefixes** — `fx-`/`ms-` in Parts 2–3, and `sma-`, `pj1-`…`rp5-`, `aux1-`… in
  Parts 4–9. Quarto refuses to render a document with duplicate labels, and `setup`/`session`/`echo`
  collided repeatedly across 26 sources.
- **Part 3 is gated on `MS_READY`**, set by `ms-setup` from whether `<outputs.directory>/validation/`
  exists. As a standalone notebook a missing directory was a hard `stop()`; inside the master that
  killed Parts 1 and 2 as well and made any smoke run impossible.

Four underlying execution paths, all config-driven:

1. [microbiome_analysis.py](microbiome_analysis.py) — ~1700 lines, flat module of pure functions,
   orchestrated by `run_pipeline()`. Deliberately limited to pandas/numpy/scipy/matplotlib
   (see [requirements-microbiome.txt](requirements-microbiome.txt)): PCoA, NMDS, PERMANOVA,
   PERMDISP, Bray-Curtis, Jaccard, Hedges' g, and random-effects meta-analysis are all
   hand-implemented rather than pulled from scikit-bio. Keep it that way — scikit-bio, biom, and
   statsmodels are not installed.
2. [phyloseq_analysis.qmd](phyloseq_analysis.qmd) — 76-chunk Quarto notebook, the primary
   manuscript-oriented analysis (phyloseq, DESeq2, vegan, picante, decontam, indicspecies,
   GUniFrac). Structured around Simonin et al. 2022 (*New Phytologist*) meta-analysis logic.
3. [phyloseq_analysis.R](phyloseq_analysis.R) — narrow R sidecar for Faith's PD and DESeq2 DA.
   `run_phyloseq_analysis()` in the Python pipeline shells out to it and parses a
   `__R_STATUS__ {json}` line from stdout; `integrate_r_outputs()` then merges
   `diversity/faith_pd.csv` back into the alpha table. If `Rscript` is missing the step is
   **skipped, not failed** — preserve that degradation.
4. **The validation path** — [validation/verify_manuscript.py](validation/verify_manuscript.py) +
   [manuscript_figures.qmd](manuscript_figures.qmd). This one is deliberately **not** wired into the
   shared helpers. Read the warning below before refactoring it.

### Figures have one home: `figures_all.qmd`

[figures_all.qmd](figures_all.qmd) generates **all 42 analysis figures** and is the only place to
add or change one. It reads the pipeline's *persisted* outputs — `setup/ps_clean.rds` (1,886 × 2,000,
429 KB), `cache/*_distance.rds`, `diversity/pcoa_*.csv`, `diversity/umap_data.csv`, `tables/*.csv` —
so it renders in about a minute and never opens the 609 MB count table. It computes no new
statistics: alpha diversity comes from `estimate_richness()` on the filtered object, and everything
else is a redraw of a stored result.

Every figure is wrapped in `fig_section(name, needs = <input paths>)`, which records a **skip** and
continues when an input is absent rather than aborting the render, then writes
`tables/figure_manifest.csv` accounting for all 44 attempts. Check that manifest after any render.

All figures land in `analysis_results/figures/` — one directory, no exceptions. `phyloseq_analysis.qmd`
historically scattered them across `figures/`, `diversity/` and `setup/`, and wrote data files
(`pairwise_*.rds`, `umap_data.csv`) into `figures/`; both are cleaned up. `diversity/umap_data.csv` is
now the canonical UMAP coordinate path.

Two figures skip legitimately: `pcoa_bray_genus_by_timepoint` (the config lists `Timepoint` as a
secondary group but `Datesets/metadata.csv` has no such column) and `pcoa_bray_genus_by_geography`
(the pipeline never exports that ordination table). Fixing the second belongs in the pipeline's
variable loop, not here — a figures document should not invent an ordination.

Gotcha when editing it: `as.data.frame(sample_data(ps))` keeps the phyloseq class and every dplyr
verb rejects it — use `as(sample_data(ps), "data.frame")` then `as_tibble()`. And the by-variable
PCoA exports store the grouping column generically as `group`, not under the variable's own name,
while the overall exports use real column names; `pcoa_plot()` accepts either.

### The validation path is independent on purpose

`validation/verify_manuscript.py` reimplements PERMANOVA, PERMDISP, Bray-Curtis, Shannon and the
filtering cascade from the raw inputs in numpy/scipy. It duplicates logic that already exists in
`microbiome_analysis.py`, and that duplication is the point — its docstring says it plainly:

> *"Deliberately NOT a re-run of phyloseq_analysis.qmd: the statistics are reimplemented from the raw
> inputs so that agreement is evidence rather than tautology."*

It is the audit that caught, among other things, that Table 4 was computed on ASV-level top-100 raw
counts at n=1,553 rather than genus-level relative abundances at n=1,883 as the caption claimed. If
you refactor it to call the pipeline's helpers, it will agree with the pipeline by construction —
including wherever both are wrong — and `analysis_results/validation/discrepancy_audit.csv` stops
meaning anything. **Do not deduplicate these statistics.** Shared config, paths and entry-point
conventions are fine and already done; shared statistical code is not.

Its one optional dependency is `statsmodels`, used at a single `smf.ols` call for the alpha
depth-model table. It is **not installed**; the import is guarded and that one table is skipped,
following the same degrade-don't-fail rule as the missing-`Rscript` path. `pip3 install statsmodels`
if you need `TableS2_alpha_depth_models.csv` regenerated.

`manuscript_figures.qmd` recomputes Fig 1 and Fig 2a from raw data but reads Fig 2b/2c, 3 and 4 from
the audit's tables of record in `analysis_results/validation/supplementary/` — so those panels show
published values rather than a fresh permutation run. Its `validate` chunk reconciles every
recomputed quantity against those tables; all six checks currently pass, with Kruskal-Wallis H at
36.354 recomputed vs 36.357 in TableS1.

### Memory discipline (non-negotiable)

`Datesets/asv-counts.csv` is **609 MB** (2,456 samples × 129,070 ASVs) and `feature-table.tsv` is
1.3 GB. Nothing may densely materialize these:

- Python uses a **two-pass streaming** design: `first_pass_counts()` computes per-sample and
  per-feature statistics row by row, `select_feature_indices()` picks the retained features, and
  `second_pass_collect()` reads only those columns. Distance matrices are blocked
  (`bray_curtis_distance(block_size=32)`, `jaccard_distance(block_size=128)`).
- R passes `nrows = analysis.max_samples` to `data.table::fread()` for smoke runs and selects only
  matching ASV columns when a filtered-feature list is supplied.
- Beta diversity runs on a separate `ps_beta` object restricted to `filters.beta_top_taxa` (100)
  genera; distances are cached to `analysis_results/cache/*.rds` and reused across Bray-Curtis,
  Jaccard, NMDS, PCoA, PERMANOVA, and PERMDISP.
- `memory.use_mori` / `memory.mirai_workers` in the config gate an optional `mori::share()` +
  `mirai` shared-memory backend; it must stay optional.

### Data and analysis conventions

- Sample key is **`run_accession`** throughout. `load_metadata()` tries UTF-8 then falls back to
  `latin-1` (author names contain accented characters) and drops `Unnamed:` columns.
- Grouping: `Developmental_Stage` (Seeds vs Seedlings) is primary; `study_alias`,
  `Scientific_Name`, `Geography`, `Timepoint` are secondary. Because this is a meta-analysis,
  `study_alias` is a confounder, not just another factor — per-study effect sizes
  (`per_study_effect_sizes()`) feed `random_effects_meta()` for forest plots.
- Taxonomy filtering keeps Bacteria only and drops Chloroplast (order), Mitochondria (family), and
  Cyanobacteria/Cyanobacteriota — configured under `filters`, applied by `build_taxonomy_mask()`.
- Rank columns in the Python code are `["Phyla", "Family", "Genus"]` (note: `Phyla`, not `Phylum`).
- Rarefaction is a sensitivity analysis only; primary compositional results are non-rarefied.

## Gotchas

- **Bioconductor masks dplyr.** With phyloseq/DESeq2 loaded, `BiocGenerics` masks `count`,
  `rename`, `filter`, `select`, `intersect`, `setdiff`, `union`, `slice`, `arrange`. Always write
  `dplyr::count()` etc. — especially inside `purrr::map()` lambdas, where the masking surfaces more
  reliably than in `|>` chains.
- **`Datesets/` is spelled that way on purpose.** Don't "fix" it; the configs and both pipelines
  reference it.
- Large inputs and `analysis_results*/` are tracked in git and only partially ignored
  ([.gitignore](.gitignore) covers just `.env`, `node_modules/`, `dist/`, `.DS_Store`). Don't add
  more multi-hundred-MB artifacts, and watch for accidental `* (1).csv` / `* (1).tiff` duplicates —
  several already exist in `Datesets/` and `analysis_results/`.
- R style here: native pipe `|>`, `.by` argument, `join_by()`, `purrr::map()` + `list_rbind()`.
- Figures are TIFF at publication settings via `theme_nature_publication()`; the palette comes from
  [wildagain_joint_pellet.R](wildagain_joint_pellet.R) through `get_pellet()`, with an
  RColorBrewer fallback when that file is absent.

## Layout

| Path | Role |
| --- | --- |
| [Datesets/](Datesets/) | Canonical inputs: `asv-counts.csv`, `taxa.csv`, `metadata.csv`, `tree.nexus` |
| [analysis_results/](analysis_results/) | Full-run outputs: `qc/`, `diversity/`, `taxonomy/`, `differential_abundance/`, `figures/`, `tables/`, `manuscript/`, `cache/` |
| `analysis_results/manuscript/` | **Off-limits.** Manuscript documents of record — do not edit, move, or delete anything here without being asked explicitly. |
| `analysis_results/validation/` | Audit inputs, tables of record (`supplementary/`), figures, `discrepancy_audit.csv`, `revision_memo.md`. Written by the validation path; `README.txt` there explains the four findings. |
| `analysis_results_smoke/` | Smoke-run outputs, same layout |
| [joined-scripts/](joined-scripts/) | Merged notebooks from the earlier EXC / strawberry projects — canonical variants under `project/`, `beta-diversity/`, `diversity-analysis/`; the rest are **archival** and may still use old variable names (see [joined-scripts/README.md](joined-scripts/README.md)) |
| [tests/](tests/) | `unittest`-style tests in `MicrobiomeAnalysisTests`, run under pytest |
| [plan.md](plan.md) | Living ExecPlan: implementation status, figure/table plan, methods prose |
| `*.slurm`, `merge_fastp.sh`, `script_multiqc.sh`, `qiime_processor.py` | Upstream HPC read-processing and QIIME2 stage that produced `Datesets/` |

`.qza`/`.biom`/`.rmarkdown` files and `phyloseq_analysis 2.rmarkdown` / `3.rmarkdown` are
intermediates and snapshots — read them for context, but edit `.qmd` and `.py` sources.

The ~30 top-level dot-directories (`.roo`, `.windsurf`, `.continue`, `.goose`, …) are empty
scaffolding from other agent tools; only `.agents/skills/` has content, and it holds installed
skills, not project rules. There are no Cursor or Copilot instruction files.

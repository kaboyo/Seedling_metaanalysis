# ExecPlan: Generalist / specialist analysis with a dendrogram figure

> This ExecPlan is a living document maintained in accordance with ~/.claude/PLANS.md.
> Last updated: 2026-10-06 — EXECUTED. The user chose method C (both analyses, two rings).
> Where the text below differs from what was built, the Outcomes section is authoritative.

## Purpose

Classify each bacterial genus as a **generalist** (found at both developmental stages in
proportion), a **seed specialist**, a **seedling specialist**, or **too rare to call**, and
show the result on a phylogenetic dendrogram so that it can be seen whether specialisation
follows phylogeny. One new table of record and one new figure.

## Context and Orientation

- `seedling_analysis_master.qmd` is the whole project in one file. **Part 1** (the pipeline)
  computes statistics and writes tables. **Part 2** (the figures) only redraws persisted
  results; it must not compute new inferential statistics. So the classification goes in
  Part 1 and the dendrogram in Part 2.
- The filtered dataset is `analysis_results/setup/ps_clean.rds`: 1,886 samples, 2,823 ASVs,
  283 genera, with the phylogenetic tree attached. Part 1 holds it as `ps`; Part 2 reloads it.
- Part 2 chunk `fx-fig-tree-genus` (master line ~6417) already draws three tree figures and
  defines the helpers this figure reuses: `TREE_FILE`, `TREE_EXCLUDE_GENERA`,
  `TREE_PHYLUM_COLOURS`, `TREE_TAXON_COLOURS`, `tree_stage_rings()`. `make_tree_figures.R`
  runs that chunk alone.
- Sample key `run_accession`; stage column `Developmental_Stage` (config
  `analysis.primary_group`), levels in `differential_abundance.groups` (`Seeds`, `Seedlings`).
- Two facts about the data that constrain the method:
  - **98.8% of ASVs occur in a single study**, so ASV-level specialisation across studies is
    an artefact of study. The analysis is therefore done at **genus** level.
  - Three of the eight studies sampled seedlings only (`PRJNA649549`, `PRJNA810761`,
    `PRJNA811248`). A genus seen only there would be called a seedling specialist merely
    because nobody sequenced the seeds. The analysis is restricted to the **five studies
    that sampled both stages** (1,690 samples).
- R packages needed are all installed: `vegan` (has `clamtest`), `ggtree`, `ggtreeExtra`,
  `ggnewscale`, `ape`, `phangorn`. No new dependency.

## Plan of Work

**Method (recommended; alternatives in the Decision Log).** The multinomial species
classification method, CLAM (Chazdon et al. 2011, *Ecology* 92:1332), as implemented in
`vegan::clamtest()`. It compares a genus's summed counts in two habitats, corrects for
unequal sampling, and assigns one of four classes. Parameters come from the config, never
hardcoded: `specialization` = 2/3 (supermajority rule), `alpha` = 0.005, `coverage_limit` = 10
— vegan's defaults, which are those recommended by Chazdon et al. for single-taxon calls.

1. **Config.** Add `clam_specialization`, `clam_alpha`, `clam_coverage_limit` under
   `taxonomy_analysis` in both `analysis_config.json` and `analysis_config.smoke.json`.
2. **Part 1, new chunk `generalist-specialist`** placed after `indicator-species`. A function
   `perform_generalist_specialist(ps, group_col, study_col, tables_dir, ...)` that:
   - sums ASV counts to genus;
   - keeps the studies in which both stage levels occur; if fewer than one such study or
     only one stage is present (the smoke run), prints a message and returns
     `status = "skipped"` — degrade, do not fail;
   - runs `clamtest()` on the pooled samples, with the stage factor levels fixed so that the
     class names never flip between runs;
   - runs `clamtest()` again within each study, and for each genus counts the studies where
     it was classifiable (`n_studies_tested`) and where the class matched the pooled class
     (`n_studies_agree`). This is the check against the largest study (1,352 of the 1,690
     samples) driving the pooled call;
   - writes `tables/generalist_specialist_genus.csv` (one row per genus: Phylum, Genus,
     reads per stage, pooled class, `n_studies_tested`, `n_studies_agree`, one class column
     per study) and `tables/generalist_specialist_summary.csv` (class counts, pooled and per
     study).
3. **Part 2, new `fig_section("tree_generalist_specialist")`** inside chunk
   `fx-fig-tree-genus`, so it shares the tree helpers and `make_tree_figures.R` draws it
   unchanged. It needs the genus table, `ps_clean.rds` and the tree. Design:
   - the pipeline's own tree pruned to one tip per classified genus (that genus's most
     abundant ASV), midpoint-rooted, fan layout — same construction as `tree_genus_by_stage`;
   - tip points coloured by phylum (existing phylum colours);
   - **ring 1: class** — seed specialist navy, seedling specialist vermillion (the reserved
     stage colours), generalist a colour used by no other legend in the figure;
     too-rare genera are left off the tree;
   - **ring 2: agreement** — in how many studies the genus received the same class;
   - **rings 3–4:** mean relative abundance in seeds and in seedlings (`tree_stage_rings()`);
   - tip labels in italics for the most abundant genera only, as many as stay legible.
   Caption added to the caption registry at the top of the master.
4. **Launcher `make_generalist_specialist.R`** — same pattern as `make_tree_figures.R`: runs
   the Part 2 setup chunks to load `ps`, then the Part 1 `generalist-specialist` chunk, so
   the table can be produced in seconds without the multi-hour full render. The figure is
   then drawn with `Rscript make_tree_figures.R`.
5. **Docs.** `make_tree_figures.R` header and the figure count in `CLAUDE.md`.

Nothing under `analysis_results/manuscript/` or `validation/` is touched. Rollback is
`git checkout` of the edited files plus deleting the two new tables, the figure and the
launcher.

## Milestones

**M1 — classification table.** Config keys, Part 1 chunk, launcher.
Run `Rscript make_generalist_specialist.R analysis_config.json`.
Verify: both CSVs exist; pooled class counts are 58 generalists, 153 seed specialists,
40 seedling specialists, 5 too rare (256 genera) — the values from the feasibility probe;
`Rscript make_generalist_specialist.R analysis_config.smoke.json` skips with a message
rather than erroring.

**M2 — dendrogram.** Figure section and caption.
Run `Rscript make_tree_figures.R analysis_config.json`.
Verify: `analysis_results/figures/tree_generalist_specialist.tiff` is written; inspected by
eye — every ring matches its legend, no colour means two things, labels legible; the three
existing tree figures are byte-for-byte unaffected in content (same tables written).

**M3 — docs and regression.** Update `CLAUDE.md` and launcher headers;
`python3 -m pytest tests/ -q` still reports 19 passed.

## Validation and Acceptance

- `analysis_results/tables/generalist_specialist_genus.csv` has one row per genus and a
  `class` column with exactly the four CLAM classes.
- `analysis_results/figures/tree_generalist_specialist.tiff` opens and shows a fan
  dendrogram with a class ring whose colours match its legend.
- The full master render is **not** run as part of this plan (it takes hours); the new
  chunks are exercised through the two launchers, which evaluate the master's own code.

## Interfaces and Dependencies

- No new packages.
- New config keys: `taxonomy_analysis.clam_specialization`, `.clam_alpha`,
  `.clam_coverage_limit`.
- New function: `perform_generalist_specialist()` in Part 1.
- New outputs: two CSVs in `tables/`, one TIFF in `figures/`.
- New file: `make_generalist_specialist.R`.

## Progress

- [x] 2026-10-06 Method approved: C, both analyses as two rings
- [x] 2026-10-06 M1 classification table
- [x] 2026-10-06 M2 dendrogram
- [x] 2026-10-06 M3 regression (20 Python tests pass); `CLAUDE.md` not updated, see Outcomes

## Surprises and Discoveries

- 2026-10-06 (probe): `clamtest()` lists its two specialist classes in the order the groups
  appear in the data. The class *names* are always right; only their order varies. The code
  reads the totals by name (`Total_<stage>`), so nothing depends on the order.
- 2026-10-06: 158 of the 250 classified genera were testable in one study only, and 48 got
  the pooled class in no single study. Drawn plainly, the tree showed a whole clade of
  "seed specialists" (Acidobacteriota, Planctomycetota, Verrucomicrobiota) that is the 12
  seed samples of `PRJNA735565`. Only 46 genera hold their class in two or more studies.
- 2026-10-06: `ggtreeExtra::geom_fruit()` fails if its data shares a column name with the
  data attached to the tree; and `ggnewscale` renames `fill`, so `override.aes = list(fill =)`
  on another legend is silently ignored.
- 2026-10-06 (probe): pooled over all eight studies the result is 61 / 152 / 70 / 0; dropping
  the three seedling-only studies removes 30 "seedling specialists", which were artefacts of
  those studies having no seed samples.

## Decision Log

- 2026-10-06 — **Genus level, not ASV.** 98.8% of ASVs occur in one study.
- 2026-10-06 — **Both-stage studies only.** See the second discovery above.
- 2026-10-06 — **Method: C chosen by the user** (both analyses, two rings).
- 2026-10-06 — **Agreement shown as opacity, not a third ring.** A stage tile is solid where
  the class holds in two or more studies and pale otherwise; this keeps the figure at the two
  rings asked for while not hiding the single-study calls.
- 2026-10-06 — **Host breadth: Levins' B against a quasiswap null**, each host reduced to one
  average sample of `host_breadth_depth` reads so that sample number does not set breadth.
  `EcolUtils` (the usual implementation) is not installed; `vegan::nullmodel` does the same.
- 2026-10-06 — **Config keys live in a new `generalist_specialist` section**, not under
  `taxonomy_analysis` as first written.
- 2026-10-06 — Options as put to the user:
  - A. *Stage specialisation by CLAM* (recommended, written up above). Matches the
    manuscript's primary contrast, seeds vs seedlings.
  - B. *Host breadth by Levins' niche breadth* across the 16 host plant species, with a
    permutation null. Answers a different question (does a genus occur on many hosts?).
    Weak here: half of the genera (141 of 283) occur on a single host species and host is
    largely confounded with study.
  - C. *Both*, as two rings on the one dendrogram.

## Limitations to state in the methods

CLAM treats counts as individuals. Sequence reads are not independent individuals, so its
significance is overstated; the per-study agreement ring is the guard against over-reading
a pooled call. Study remains a confounder that pooling does not remove.

## Outcomes and Retrospective

Built:

- `analysis_config.json`, `analysis_config.smoke.json`: new `generalist_specialist` section.
- `seedling_analysis_master.qmd`: Part 1 chunk `generalist-specialist`; a fourth
  `fig_section` in Part 2 chunk `fx-fig-tree-genus`; one caption.
- `make_generalist_specialist.R` (new); `make_tree_figures.R` header.
- `analysis_results/tables/generalist_specialist_genus.csv` (282 genera) and
  `generalist_specialist_summary.csv`; `analysis_results/figures/tree_generalist_specialist.tiff`.

Results (full config): stage, CLAM, five both-stage studies pooled — 58 generalists,
152 seed specialists, 40 seedling specialists, 5 too rare. Host breadth, 16 host species,
999 null tables — 6 generalists, 25 specialists, 151 neither, 100 too rare.

Verified: both launchers run clean on the full config; the three existing tree tables are
byte-identical before and after; 20 Python tests pass; figure inspected by eye.

Added after the first delivery, at the user's request (2026-10-06):

- A host-context heatmap ring: one column per host plant genus (12), each cell the mean
  relative abundance of the bacterial genus on that host, sqrt colour scale, absent cells
  in a faint grey so a column can be traced to its name in the opening of the tree. The
  values are written to `tables/tree_generalist_specialist.csv`.
- Larger fonts, and the canvas resized to 7.2 x 8.3 in (double-column width) with the
  legends beneath the tree; the heatmap needs the full page width.

Not done:

- The full master render was not run, so `tables/figure_captions.*` and
  `figure_manifest.csv` do not list the new figure yet.
- The in-function skip path (data with one stage only) was not exercised: there is no smoke
  `ps_clean.rds` to run it on.
- `CLAUDE.md` was left alone. Its figure counts (42 / 44) were already out of date before
  this work and need a wider refresh than this plan covers.

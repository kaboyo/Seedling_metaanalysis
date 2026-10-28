# Plan: Figure 2 (composition with beta diversity) and the Cell Reports colours throughout

## Context

Two requests, to be done together because the first figure should be born in the new colours:

1. **Composition and beta diversity in one figure.** The class-level Sankey and the PCoAs
   belong together. The PCoAs are to be coloured by host species, by cultivar and by plant
   family, and shown separately for seeds and for seedlings. Your answers: keep the
   **class-level Sankey**; a **separate ordination for each stage**; number it **Figure 2**,
   which makes five main figures.
2. **Cell Reports colours everywhere.** Every figure in `analysis_results/by_analysis/` is
   to use the colours of `CELL-REPORTS-D-26-05893.pdf`.

## Part A — Figure 2: composition and beta diversity

Target: one page, 183 mm wide, about 235 mm tall.

| Panel | Content | Source |
| --- | --- | --- |
| A | Class composition in seeds and seedlings of each plant family (Sankey) | existing `fig_alluvial_class_by_plant_family` |
| B | PCoA coloured by host species: seeds above, seedlings below | new |
| C | PCoA coloured by cultivar: seeds above, seedlings below | new |
| D | PCoA coloured by plant family: seeds above, seedlings below | new |

B–D are one new figure, `pcoa_css_by_stage.tiff`: three columns (species, cultivar, family),
two rows (seeds, seedlings), a legend under each column.

**The ordinations are computed in Part 1 and only drawn in Part 2**, as the project requires.

- New Part 1 chunk `pcoa-by-stage`, straight after `beta-diversity-css`. It reads the
  CSS-normalised genus-level Bray–Curtis distances that chunk already caches
  (`cache/bray_css_genus.rds`), takes the seed samples and the seedling samples apart, and
  runs `cmdscale(k = 2, eig = TRUE)` on each, exactly as `css_pcoa_plot()` does for the
  existing PCoAs. It writes `diversity/pcoa_css_genus_by_stage.csv`: sample, stage, two
  coordinates, the variance each axis explains, host species, cultivar, plant family.
- New Part 2 chunk `fx-fig-pcoa-by-stage`, in the beta-diversity section. It draws the six
  panels from that table. Axis labels carry each stage's own variance explained.
- Samples with no cultivar recorded are drawn grey as "not reported", so the three columns
  show the same samples.

**What this does not include:** PERMANOVA within each stage. The panels are ordinations
only; the legend cites the existing overall tests (Table S4a). Say if you want per-stage
tests added, which would be a new analysis and a new methods sentence.

**Numbering** (config `main_figures`, merged by the existing chunk `fx-fig-main-figures`):

| New | Content | Was |
| --- | --- | --- |
| Figure 1 | Study overview | Figure 1 |
| Figure 2 | Sankey (A) + PCoAs by stage (B–D) | `Composition_Sankey` and `Beta_diversity` (both retired) |
| Figure 3 | Host phylogeny and phylosymbiosis | Figure 2 |
| Figure 4 | Stage effects within plant families and cultivars | Figure 3 |
| Figure 5 | Host generalists and specialists | Figure 4 |

## Part B — Cell Reports colours throughout

**The palette.** I sampled the dominant colours of the paper's figure pages. These are
approximations read off the PDF, not the authors' own colour codes; the first step of the
work is to sample them again at full resolution and fix the final values.

| Role | Colours (approximate) |
| --- | --- |
| Categorical, in order of use | plum `#845C64`, teal `#0C9494`, orange `#DC7C2C`, slate blue `#547CAC`, sage `#7C9C7C`, sand `#F4CC8C`, purple `#84648C`, coral `#D46C74`, dark green `#3C5C44`, lavender grey `#B4BCCC`, pale yellow `#F4E49C`, steel `#4C5C74` |
| Seeds / seedlings | orange `#DC7C2C` / teal `#0C9494` |
| Sequential (bars, heat) | pale sage to dark green, as on the paper's phylogeny bars |
| Diverging (Hedges' g, log odds) | orange – white – teal, so the stage heatmaps keep "seed side / seedling side" |
| Host generalist / specialist | orange `#EC9C64` / slate blue `#547CAC`, as on the paper's host-breadth bars |

**Where colours are set today**, and what changes:

- One new `palette` section in `analysis_config.json` and `analysis_config.smoke.json` holds
  the roles above. Nothing is hard-coded in the script.
- Part 1 chunk `color-palettes` and Part 2 chunk `fx-palettes` build their base palettes
  from it: `NATURE_THICK`, `STAGE_PELLET`, `distinct_palette()`, `get_palette()`,
  `PIPELINE_BAR_FILL`.
- The local colour constants are re-pointed to palette roles instead of fixed hex codes:
  `HOST_CLASS_COLOURS`, `OXYGEN_COLOURS`, `BREADTH_CLASS_COLOURS`, `TRAIT_SOURCE_COLOURS`,
  `HOST_FAMILY_COLOURS`, `TREE_TAXON_COLOURS`, `STAGE_HEAT_FAMILY_COLOURS`, the palette inside
  `build_taxa_strips()`, the country colours of the overview, and the `BrBG`, `BuPu`
  and seed-to-seedling gradients of the heatmaps.
- Rules kept from earlier corrections (`tasks/lessons.md`): seeds and seedlings keep their
  two colours for nothing else; categorical colours are solid and clearly distinct; no colour
  means two things in one figure.

**Redrawing.** Colours live in the plotting code, so every figure is redrawn; no statistic
changes.

| Figures | How they are redrawn |
| --- | --- |
| Part 2 figures (overview, alpha, trees, heatmaps, Sankey, volcanoes, network, and so on; 34 figures) | run all Part 2 chunks from the master, as verified today |
| Taxa strips (9) | `make_composition.R` |
| CSS PCoAs (8) | their drawing loop moves from Part 1 to the new Part 2 PCoA chunk, reading the cached distances, so they are redrawn with the rest |
| Geographic maps (2) and DESeq2 volcano (1) | Part 1 chunks; I will try them through a launcher on the saved data. If one cannot be rebuilt without the raw 609 MB table, it keeps its old colours until the next full render, and I will say which |

## Files to change

- `seedling_analysis_master.qmd` — Part 1: `color-palettes`, new chunk `pcoa-by-stage`,
  the palette in `build_taxa_strips()`, the CSS PCoA drawing loop. Part 2: `fx-palettes`,
  new chunk `fx-fig-pcoa-by-stage`, the colour constants listed above, captions.
- `analysis_config.json`, `analysis_config.smoke.json` — `palette`; `main_figures`
  (five figures); the new figures in `analysis_folders`.
- `make_composition.R` or a new `make_beta_figures.R` — to run the two new chunks.
- `analysis_results/by_analysis/00_main_figures/figure_legends.md` — five legends and the
  table of where the manuscript's current figures go.
- `analysis_results/figure_legends_host_analyses.md` — colour words in the legends
  ("magenta", "navy", "vermillion") changed to the new colours.
- `CLAUDE.md` — the palette's home; `tasks/plans/figure2-and-cell-palette.md` — this plan
  with progress, per the project's rule.

Not touched: `analysis_results/manuscript/` (the manuscript and the supplementary workbook
hold no figures), the validation path and its audit figures, all tables.

## Order of work

1. Sample the paper's colours at full resolution; fix the palette; check that the categorical
   colours are distinguishable, including for colour-blind readers.
2. Put the palette in the config and wire the two base-palette chunks to it. Redraw two test
   figures (overview, stage heatmap) and look at them before going further.
3. Re-point the local colour constants; redraw everything; look at every main figure.
4. Add the two PCoA chunks; draw `pcoa_css_by_stage.tiff`; build Figure 2; renumber.
5. Legends, captions, project notes.

## Verification

- All Part 2 chunks run from the master with no error (the runner used today), and the
  figure manifest shows no failure.
- The per-stage PCoA chunk is checked against the existing one: run on all samples it must
  reproduce the axes of the current PCoA (14.8% and 8.6%).
- No table changes: checksums of `tables/`, `diversity/` and `differential_abundance/`
  before and after are identical, apart from the one new coordinates table.
- A search of the master for the old hex codes of `NATURE_THICK`, the magenta and the navy /
  vermillion pair returns nothing outside the archival Parts 4–9.
- Every main figure and the Sankey are opened and checked by eye: legends complete, no two
  categories sharing a colour, text readable against the new fills.
- All chunks parse; `python3 -m pytest tests/ -q` passes (20 tests).
- The full master render (hours) is not run.

## Things to know before approving

- **Five main figures**, not four, as you chose.
- **Separate ordinations cannot be compared between the two rows.** A sample's position in
  the seeds panel and in the seedlings panel are on different axes. The legend will say so.
- **Twenty cultivars and sixteen host species need more colours than the paper's palette
  has** (about twelve distinct ones). For those two legends the palette is extended with
  lighter and darker steps of the same hues; neighbouring shades will be harder to tell apart
  than they are now.
- **Cell Reports' own figures use other colours as well** (bright ones for host orders, for
  example). I am taking the muted palette of its overview, bar and tree figures as "the"
  palette.
- **A caution.** `CELL-REPORTS-D-26-05893.pdf` is a submitted manuscript, and a copy sits
  in your `Documents/Review/cell/` folder. If you received it as a reviewer, it is
  confidential until published. A colour scheme and a panel arrangement are not its
  intellectual content, but your methods text and figure plan now follow it closely; that is
  worth a thought before submission, and citing it once it is published.
- Nothing is committed.

## Progress (2026-10-12)

- [x] Palette sampled from the PDF pages and set in the config (`palette`).
- [x] Both base-palette chunks and every local colour constant read it; ColorBrewer and
      viridis scales replaced. Host generalist / specialist became purple / dark slate, not
      orange / slate blue as first planned: orange is the seed colour.
- [x] Part 1 chunk `pcoa-by-stage`, Part 2 chunk `fx-fig-pcoa-by-stage`, launcher
      `make_beta_figures.R`. On all samples the method reproduces 14.8% and 8.6%.
- [x] Figure 2 built; five main figures in the config.
- [x] All 31 Part 2 chunks run; 35 figures written, 5 switched off. Strips and the eight
      Part 1 PCoAs redrawn (PERMANOVA cache reused, distances unchanged).
- [x] 2026-10-13: Figure 2 panel A is now `composition_dendrogram` (chunk
  `fx-fig-composition-dendrogram`): circular, seeds and seedlings clustered apart and joined at
  the root, species circle, stage and family rings, class and genus bars. The class Sankey
  moved to Figure S2 (15 supplementary figures). Legends shortened everywhere; the network
  (S15) redrawn by module; Figure S11 panel B labelled with the bar's own value.
- [ ] Not redrawn: `volcano_deseq2`, `fig1_geographic_distribution`,
      `geographic_distribution_map`. They need Part 1 objects that are not saved; they take
      the new colours at the next full render.

## Surprises

- The plan said the Part 1 PCoA loop would move to Part 2. It was simpler and safer to run
  the Part 1 chunk as it is: its PERMANOVA cache was reused and nothing was recomputed.
- Running the two map chunks on the saved filtered data gave a different study table from
  the pipeline's (the pipeline draws them at an earlier filtering step). The table was put
  back and the redrawn maps discarded.

# Microbiome Analysis Plan for Seedling Meta-analysis

## Summary
- Implement a reproducible comparative-ecology workflow for the microbiome data under `Datesets/`.
- Use `Developmental_Stage` as the primary grouping variable and treat `study_alias`, `Scientific_Name`, and `Geography` as secondary explanatory variables.
- Prefer `Datesets/tree.nexus` as the canonical tree input. If only `Datesets/tree.qza` exists, derive `tree.nexus` from the embedded Newick tree first.

## Inputs
- `Datesets/asv-counts.csv`
- `Datesets/taxa.csv`
- `Datesets/metadata.csv`
- `Datesets/tree.nexus`

## Current State
- `asv-counts.csv` is large at roughly 609 MB with 2,456 sample rows and 129,070 ASV columns.
- `metadata.csv` loads with `latin-1` fallback and currently has 2,457 metadata rows, one more than the ASV table.
- `taxa.csv` contains many unassigned features, so assigned-only summaries are needed alongside full summaries.
- The local environment has `pandas`, `numpy`, `matplotlib`, and `seaborn`, but does not currently have `qiime`, `scikit-bio`, `biom`, `scipy`, `statsmodels`, or `sklearn`.

## Implementation Status
- Added `microbiome_analysis.py` as the main workflow entrypoint.
- Added `tests/test_microbiome_analysis.py` for metadata loading, Nexus tree parsing, and a tiny end-to-end run.
- Added `analysis_config.json` with `Datesets/tree.nexus` as the canonical tree input.
- Converted `Datesets/tree.qza` to `Datesets/tree.nexus`.
- Added an `outputs.generate_plots` toggle so tests and smoke runs can skip Matplotlib-heavy figure generation.
- Verified the workflow on a 120-sample smoke run written to `analysis_results_smoke/`.
- Confirmed smoke-run outputs for QC, alpha diversity, distance matrices, ordination, PERMANOVA, PERMDISP, taxonomy summaries, and tree-status reporting.
- Patched distance-matrix calculations to avoid divide warnings on zero-denominator comparisons.
- Updated `phyloseq_analysis.qmd` with an optional Mori/mirai memory backend: `mori::share()` is used when installed to keep large R objects in OS-level shared memory for worker reuse, and `mirai::daemons(memory = ...)` can be enabled for byte-capped task queues.
- Added in-render beta-distance caching and configured beta-diversity analyses to use the top `beta_top_taxa` taxa, so Bray-Curtis/Jaccard/ordination steps reuse one distance object instead of recomputing large matrices.
- Patched the Quarto ASV loader so smoke runs with `analysis.max_samples` read only the needed sample rows, and filtered-feature runs read only matching ASV columns when a filtered feature list is provided.
- Applied the local `nature-paper-workflow` to the Quarto analysis by adding a manuscript package layer with a figure plan, table plan, claim-evidence register, FAIR data checklist, and code availability record under `analysis_results/manuscript/` when the notebook is rendered.

## Analysis Workflow

### 1. Environment and Configuration
- Keep the first implementation self-contained in Python with `pandas` and `numpy` fallbacks where higher-level microbiome packages are missing.
- Store workflow settings in `analysis_config.json`, including:
  - input paths
  - filtering thresholds
  - primary and secondary grouping variables
  - permutation count
  - output directory
  - a `generate_plots` switch for smoke tests and non-graphical runs

### 2. Data Loading and Validation
- Load metadata with UTF-8 first and `latin-1` fallback.
- Drop unnamed metadata columns and require `run_accession`.
- Validate:
  - sample IDs in metadata vs count table
  - unique sample IDs
  - tree tip overlap with retained features
  - taxonomy key overlap with ASV features

### 3. QC and Filtering
- Stream the ASV matrix instead of loading the full dense table into memory.
- In the Quarto workflow, pass `nrows = analysis.max_samples` to `data.table::fread()` for smoke runs so the 609 MB count table is not fully materialized before subsetting.
- When a filtered-feature list is provided, select only the first sample-ID column plus matching ASV feature columns during `fread()` rather than reading all features first.
- Compute per-sample:
  - library size
  - observed ASVs
  - Shannon diversity
  - Simpson diversity
  - Pielou evenness
  - unassigned fraction
- Compute per-feature prevalence and total abundance.
- Filter with config-driven thresholds for:
  - minimum library size
  - minimum feature prevalence
  - minimum total count
  - maximum retained feature count for downstream feature-space analyses

### 4. Diversity Analysis
- Write alpha-diversity outputs for all retained samples.
- Build taxon-collapsed abundance tables at:
  - phylum
  - family
  - genus
- Use the genus-level relative abundance table as the default beta-diversity input.
- Limit beta-diversity to the top configured taxa to keep distance computation practical in the current environment.
- In the R/Quarto workflow, build a dedicated `ps_beta` object from the cleaned phyloseq object and restrict it to `filters.beta_top_taxa` before Bray-Curtis, Jaccard, NMDS, PCoA, PERMANOVA, PERMDISP, and UMAP.
- Cache each beta-distance object in memory and as compact RDS under `analysis_results/cache/`; write full CSV distance matrices only when `outputs.write_distance_matrices` is explicitly enabled.
- Compute:
  - Bray-Curtis distance
  - Jaccard distance
  - PCoA from Bray-Curtis
- Run:
  - PERMANOVA
  - PERMDISP
  for the primary and secondary grouping variables.

### 5. Phylogenetic Handling
- Prefer `.nexus` for tree parsing.
- Accept `.qza` only as an upstream source for conversion to `.nexus`.
- Until `scikit-bio` is installed, record Faith's PD and UniFrac as skipped with explicit status outputs rather than silently failing.

### 6. Reporting Outputs
- Write outputs under `analysis_results/` with these subdirectories:
  - `qc/`
  - `diversity/`
  - `taxonomy/`
  - `derived/`
- Required tables:
  - sample QC summary
  - feature filtering summary
  - alpha-diversity table
  - Bray-Curtis distance matrix
  - Jaccard distance matrix
  - PCoA coordinates
  - PERMANOVA results
  - PERMDISP results
  - phylum/family/genus relative abundance tables
  - phylogenetic status table
- Required figures:
  - library-size histogram
  - alpha-diversity by primary group
  - Bray-Curtis PCoA
  - unassigned-fraction by primary group
  - top-genus stacked summary by primary group

### 7. Nature-Style Manuscript Package
- Treat manuscript claims as output-backed statements, not narrative assertions. Each major claim must point to a concrete table, figure, or status output.
- Write manuscript planning outputs under `analysis_results/manuscript/`:
  - `nature_figure_plan.csv`
  - `nature_table_plan.csv`
  - `claim_evidence_register.csv`
  - `fair_data_availability_checklist.csv`
  - `code_availability_record.csv`
- Use the following main-figure spine unless the rendered results argue otherwise:
  - Figure 1: dataset scope and sequencing depth
  - Figure 2: dominant taxonomic composition
  - Figure 3: alpha diversity by developmental stage
  - Figure 4: Bray-Curtis beta-diversity ordination
  - Figure 5: core and flexible microbiome structure
  - Figure 6: differential or indicator taxa
- Keep statistical interpretations conditional on exported evidence: PERMANOVA must be read with PERMDISP, alpha-diversity claims must include exact test outputs, and core-microbiome claims must state prevalence and abundance thresholds.
- Before journal submission, replace local-only data paths with accessions, DOIs, repository links, software versions, and a commit/tag.

### 8. Memory Backend and Async Execution
- Keep `memory.enabled: true` in the JSON configs so the notebook uses memory-aware defaults.
- If the `mori` package is installed, call `mori::share()` on large count, metadata, and taxonomy objects after loading. This preserves a single OS-level shared-memory backing for objects sent to local worker processes.
- If `mirai_workers` is set above 1 and `mirai` is installed, start local daemons with `daemons(memory = memory.mirai_memory_mb)`; leave `mirai_workers: 1` for the current serial notebook path.
- Treat Mori as an optimization layer, not a hard dependency. The notebook must continue without `mori` or `mirai`, using top-taxa reduction and cached distance objects as the baseline memory control.
- Prefer compact `dist` RDS outputs over full dense CSV distance matrices in production runs because converting `dist` objects with `as.matrix()` can temporarily double memory use.

## Validation
- Metadata must load successfully with fallback encoding.
- The metadata/count-table row mismatch must be surfaced in QC outputs.
- Retained features must be recorded after filtering.
- Beta-diversity matrices must be square and sample-aligned.
- PCoA output must contain one row per analyzed sample.
- Taxonomy-collapsed tables must preserve per-sample totals after normalization.
- Tree overlap with retained features must be reported before any phylogenetic metric is attempted.

## Assumptions
- The first implementation targets comparative ecology rather than differential-abundance modeling or formal effect-size meta-analysis.
- Group comparisons are one-factor tests in this first pass.
- Missing phylogenetic packages do not block the rest of the workflow.
- `plan.md` should be kept current as implementation decisions change.

---

## Bioinformatic and Statistical Analysis (R/Phyloseq Pipeline)

### Overview

All downstream bioinformatic and statistical analyses were performed in R (v4.4+) using the Quarto reproducible reporting framework (`phyloseq_analysis.qmd`). The pipeline is config-driven via a JSON configuration file (`analysis_config.json`) that specifies input paths, filtering thresholds, grouping variables, and output directories. The primary grouping variable was `Developmental_Stage` (Seeds vs. Seedlings), with `study_alias`, `Scientific_Name`, and `Geography` treated as secondary explanatory variables.

### Data Import and Phyloseq Object Construction

Amplicon sequence variant (ASV) count tables, sample metadata, taxonomy assignments, and a rooted phylogenetic tree were imported and merged into a single phyloseq object (McMurdie and Holmes, 2013). The count matrix was transposed so that features (ASVs) were represented as rows and samples as columns (`taxa_are_rows = TRUE`). Taxonomy assignments were matched to ASV identifiers via the `OTUID` column, and the phylogenetic tree (Nexus format) was read using the ape package (Paradis and Schliep, 2019).

### Taxonomy Filtering and Quality Control

Prior to any diversity or composition analysis, non-target sequences were removed in the following order: (i) non-bacterial kingdoms were excluded, retaining only features classified as Bacteria; (ii) chloroplast sequences (Order: Chloroplast) were removed; (iii) mitochondrial sequences (Family: Mitochondria) were removed; and (iv) remaining Cyanobacteria/Cyanobacteriota sequences were excluded. Soil samples were removed from the dataset to focus analyses on seed and seedling developmental stages. The list of retained feature identifiers was saved for reproducibility.

### Decontamination

Putative contaminant ASVs were identified using prevalence-based methods implemented in the decontam package (Davis et al., 2018). Where negative control samples were available in the metadata, the `isContaminant()` function was applied with `method = "prevalence"`. Identified contaminants were removed from the phyloseq object prior to downstream analyses.

### Data Cleaning

Following decontamination, negative control and blank samples were removed based on sample identifier patterns. Samples with fewer than 500 total reads (configurable via `min_library_size`) were excluded. Zero-sum taxa were pruned, and taxonomy labels were cleaned by removing uncultured designations (e.g., `d__uncultured`, `s__metagenome`) and taxonomic rank prefixes (e.g., `d__`, `p__`, `c__`). ASVs were renamed with sequential numbering (e.g., `ASV 1 : original_id`) for readability.

### Taxonomic Composition

Phylum-level taxonomic composition was assessed by aggregating ASV counts to the phylum rank. Relative abundances were computed per sample, and the top 12 phyla by mean relative abundance were retained for visualization; remaining phyla were collapsed into an "Other" category. Horizontal stacked bar charts were generated to display mean relative abundance per study, faceted by plant species (grouped by plant family) and developmental stage (Seeds vs. Seedlings). Plant species were ordered by plant family membership to enable visual comparison of related host taxa.

### Heatmap Visualization

The top 25 taxa by mean relative abundance were visualized using a log10-transformed relative abundance heatmap, implemented with both ggplot2 tile-based and pheatmap clustered heatmap representations. Samples were grouped by the primary grouping variable to facilitate visual identification of taxon-group associations.

### Core Microbiome Identification

Core microbiome membership was determined using the abundance-occupancy framework implemented in the microbiome package (Lahti and Shetty, 2017). ASVs were classified as "core" if they exceeded configurable prevalence (default: 50%) and mean relative abundance (default: 0.1%) thresholds. Global core taxa were identified across all samples, and group-specific core analyses were performed for each developmental stage using `microbiome::core_members()`. Abundance-occupancy scatter plots were generated to visualize the relationship between mean relative abundance and prevalence, with core and flexible taxa distinguished by color. A stacked bar chart summarized the proportion of total relative abundance attributed to core versus flexible taxa.

### Indicator Species Analysis

Indicator taxa associated with each developmental stage were identified using the indicator value (IndVal.g) method implemented in the indicspecies package (De Cáceres and Legendre, 2009). The `multipatt()` function was applied to relative abundance data with 999 permutations. Significance was assessed at *P* < 0.05. A complementary analysis using the phi coefficient of association (`func = "r.g"`) was also performed. Significant indicator taxa were annotated with genus-level taxonomy and visualized as horizontal bar charts ordered by indicator value.

### Alpha Diversity

Alpha diversity metrics were calculated on the cleaned phyloseq object using `phyloseq::estimate_richness()`, including Observed ASV richness, Shannon diversity index, Simpson diversity index, and Fisher's alpha. Pielou's evenness was derived as the ratio of the Shannon index to the natural logarithm of observed richness. Faith's phylogenetic diversity (PD) was computed using `picante::pd()` where a phylogenetic tree was available (Faith, 1992).

Data normality was assessed using the Shapiro-Wilk test. Statistical differences in alpha diversity among developmental stages were evaluated using the Kruskal-Wallis rank-sum test, followed by pairwise Wilcoxon rank-sum tests with Bonferroni correction for multiple comparisons. Diversity distributions were visualized as boxplots overlaid with jittered points (boxplot alpha = 0.3, jitter alpha = 0.8). Additionally, a faceted richness plot was generated to display observed ASV richness per plant species, ordered by median richness, with median values indicated by vertical lines and a log10-scaled x-axis.

### Rarefied Alpha Diversity

To account for uneven sequencing depth, alpha diversity metrics (Observed, Shannon, Simpson) were recomputed on data rarefied to the minimum library size (or a configurable depth, default: 500 reads) using `vegan::rrarefy()`. Kruskal-Wallis tests were repeated on rarefied data to confirm robustness of alpha diversity patterns. Rarefaction curves were generated at evenly spaced sequencing depths using `vegan::rarefy()` to assess sampling sufficiency.

### Beta Diversity

Community dissimilarity was assessed using Bray-Curtis and Jaccard distance matrices computed via `phyloseq::distance()`. Where a phylogenetic tree was available, weighted and unweighted UniFrac distances were additionally computed using `phyloseq::UniFrac()` (Lozupone and Knight, 2005).

**Principal Coordinates Analysis (PCoA).** Ordination was performed on each distance matrix using classical multidimensional scaling (`stats::cmdscale()` with `eig = TRUE`). The proportion of variance explained by each axis was reported. PCoA plots were annotated with 95% confidence ellipses (`ggplot2::stat_ellipse()`), with points colored by developmental stage and, in supplementary plots, additionally shaped by developmental stage and colored by plant species.

**Non-metric Multidimensional Scaling (NMDS).** NMDS ordination was performed on the Bray-Curtis distance matrix using `vegan::metaMDS()` with *k* = 2 dimensions, 50 random starts, and a maximum of 100 iterations. The final stress value was reported, and a Shepard diagram was generated to assess the goodness of ordination fit. NMDS plots included 95% confidence ellipses around group centroids.

**UMAP.** As a complementary dimensionality reduction approach, Uniform Manifold Approximation and Projection (UMAP) was applied to centered log-ratio (CLR)-transformed abundance data using the umap package. CLR transformation was performed via `microbiome::transform()`.

### Statistical Testing of Community Structure

**Permutational Multivariate Analysis of Variance (PERMANOVA).** The effect of developmental stage on community composition was tested using `vegan::adonis2()` applied to the Bray-Curtis distance matrix with 999 permutations.

**Permutational Analysis of Multivariate Dispersions (PERMDISP).** Homogeneity of group dispersions was assessed using `vegan::betadisper()` followed by `vegan::permutest()` with 999 permutations. Distances to group centroids were visualized as boxplots overlaid with jittered points.

### Differential Abundance Analysis

**DESeq2.** Differentially abundant ASVs between Seeds and Seedlings were identified using the DESeq2 package (Love et al., 2014). Count data were modeled using a negative binomial generalized linear model with size factors estimated via the "poscounts" method to accommodate zero-inflated sparse count data. Statistical significance was assessed using the Wald test, and *P* values were adjusted for multiple comparisons using the Benjamini-Hochberg procedure. Results were visualized as volcano plots (log2 fold-change vs. −log10 adjusted *P* value), with significant taxa (|log2FC| > 1 and adjusted *P* < 0.05) colored by direction of change and the top 15 significant taxa labeled using ggrepel (Slowikowski, 2024).

**ANCOM-BC2.** As a complementary approach, Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) was performed using `microbiomeMarker::run_ancombc2()` (Lin and Peddada, 2020). Significant features were identified at an adjusted *P* value threshold of 0.05.

### Co-occurrence Network Analysis

Microbial co-occurrence patterns were assessed among the top 50 taxa by mean relative abundance. Pairwise Spearman rank correlations were computed, and *P* values were adjusted for multiple comparisons using the Benjamini-Hochberg method. Significant co-occurrence relationships were defined as those with |ρ| ≥ 0.6 and adjusted *P* < 0.05. The resulting network was constructed using the igraph package (Csardi and Nepusz, 2006), with edges colored by the sign of the correlation (red: positive; blue: negative) and edge width scaled by correlation strength. Node size was proportional to node degree. The network was visualized using the Fruchterman-Reingold force-directed layout algorithm, and the edge list was exported for downstream analysis.

### Visualization

All figures were generated using ggplot2 (Wickham, 2016) with custom theme functions adapted from the visual style of Simonin et al. (2022, *New Phytologist*). Three reusable themes were implemented: (i) a faceted theme (`theme_pipeline_faceted()`) with grey panel backgrounds (#eeeeee), white gridlines, and #d9dbdb strip backgrounds for faceted diversity and composition plots; (ii) a classic theme (`theme_pipeline_classic()`) with bold axis labels for scatter plots, ordination, and bar charts; and (iii) an occupancy theme (`theme_pipeline_occupancy()`) for abundance-occupancy scatter plots. Multi-panel figures were assembled using cowplot (Wilke, 2024). Color palettes were drawn from the wild-again project (Lobato et al., 2024, *Microbiome*). All figures were saved at 300 dpi resolution.

### Software and Packages

The following R packages were used: phyloseq (McMurdie and Holmes, 2013), vegan (Oksanen et al., 2024), DESeq2 (Love et al., 2014), microbiome (Lahti and Shetty, 2017), microbiomeMarker, microbiomeutilities, decontam (Davis et al., 2018), indicspecies (De Cáceres and Legendre, 2009), picante (Kembel et al., 2010), ape (Paradis and Schliep, 2019), GUniFrac, igraph (Csardi and Nepusz, 2006), umap, ggplot2 (Wickham, 2016), cowplot, ggpubr, patchwork, ggrepel, dplyr, tidyr, readr, stringr, purrr, and ranacapa. All analyses were performed using a config-driven workflow to ensure reproducibility across full and smoke-test runs.
Plan to implement                                                    │
│                                                                      │
│ Plan: Apply Microbiome Analysis Skills to Seedling Meta-Analysis     │
│                                                                      │
│ Context                                                              │
│                                                                      │
│ The seedling meta-analysis project has a working Python pipeline     │
│ (microbiome_analysis.py) that computes basic diversity metrics on    │
│ 2,456 samples across 129,070 ASVs. A full run has already completed  │
│ with outputs in analysis_results/. However, the pipeline lacks       │
│ several analyses critical for publication: Kruskal-Wallis/post-hoc   │
│ tests, differential abundance, cross-study meta-analysis, advanced   │
│ visualizations, and phylogenetic diversity metrics.                  │
│                                                                      │
│ The goal is to invoke the available microbiome analysis skills       │
│ (diversity-analysis, differential-abundance, taxonomic-analysis,     │
│ meta-analysis, microbiome-visualization) against the existing        │
│ datasets, then update microbiome_analysis.py with the enhanced       │
│ functionality those skills produce.                                  │
│                                                                      │
│ Existing Full-Run Data Available                                     │
│ Output: Alpha diversity (2,456 samples)                              │
│ Path: analysis_results/diversity/alpha_diversity.csv                 │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: Bray-Curtis distance matrix                                  │
│ Path: analysis_results/diversity/braycurtis_distance.csv             │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: Jaccard distance matrix                                      │
│ Path: analysis_results/diversity/jaccard_distance.csv                │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: PCoA coordinates                                             │
│ Path: analysis_results/diversity/braycurtis_pcoa.csv                 │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: PERMANOVA results                                            │
│ Path: analysis_results/diversity/permanova_results.csv               │
│ Status: Complete (99 perms)                                          │
│ ────────────────────────────────────────                             │
│ Output: PERMDISP results                                             │
│ Path: analysis_results/diversity/permdisp_results.csv                │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: Genus/Family/Phylum abundance                                │
│ Path: analysis_results/taxonomy/*.csv                                │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: QC summaries                                                 │
│ Path: analysis_results/qc/*.csv                                      │
│ Status: Complete                                                     │
│ ────────────────────────────────────────                             │
│ Output: Phylogenetic metrics                                         │
│ Path: Skipped (no scikit-bio)                                        │
│ Status: Missing                                                      │
│ Key data characteristics:                                            │
│ - 3 Developmental_Stage groups: Seeds, Seedlings, soil (48 samples) │
│ - 8 studies (study_alias), 17 species, 10 geographies                │
│ - 22 zero-library samples, 82% unassigned taxa                       │
│ - 2,000 features retained after filtering                            │
│                                                                      │
│ Step-by-Step Plan                                                    │
│                                                                      │
│ Step 1: Diversity Analysis Skill                                     │
│                                                                      │
│ Skill: diversity-analysis                                            │
│ Input: Pre-computed alpha diversity CSV, distance matrices,          │
│ metadata, tree                                                       │
│ Goal: Add Kruskal-Wallis tests, Dunn's post-hoc comparisons, NMDS    │
│ ordination, and increased permutation count for PERMANOVA/PERMDISP   │
│ Output: Enhanced statistical test results, NMDS coordinates          │
│ Update to script: Add kruskal_wallis(), dunns_posthoc(), nmds()      │
│ functions; increase default permutations to 999                      │
│                                                                      │
│ Step 2: Differential Abundance Skill                                 │
│                                                                      │
│ Skill: differential-abundance                                        │
│ Input: analysis_results/taxonomy/genus_relative_abundance.csv,       │
│ Datesets/metadata.csv                                                │
│ Goal: Identify taxa significantly enriched/depleted in Seeds vs      │
│ Seedlings, controlling for study effects                            │
│ Output: DA results table with log2 fold-changes, p-values, volcano   │
│ plot data                                                            │
│ Update to script: Add differential_abundance() function; add         │
│ differential_abundance/ output directory                             │
│                                                                      │
│ Step 3: Taxonomic Analysis Skill                                     │
│                                                                      │
│ Skill: taxonomic-analysis                                            │
│ Input: Datesets/taxa.csv, genus/family abundance tables, metadata    │
│ Goal: Core microbiome identification, indicator species analysis,    │
│ shared/unique taxa between groups                                    │
│ Output: Core microbiome tables, indicator species results            │
│ Update to script: Add core_microbiome(), indicator_species()         │
│ functions; extend taxonomy outputs                                   │
│                                                                      │
│ Step 4: Meta-Analysis Skill                                          │
│                                                                      │
│ Skill: meta-analysis                                                 │
│ Input: Alpha diversity with study labels, per-study PERMANOVA, DA    │
│ results                                                              │
│ Goal: Per-study effect sizes (Hedges' g) for Seeds vs Seedlings     │
│ diversity, random-effects pooling, heterogeneity (I-squared),        │
│ forest/funnel plots                                                  │
│ Output: Effect size tables, pooled estimates, heterogeneity stats    │
│ Update to script: Add meta_analysis/ output directory; add effect    │
│ size and pooling functions                                           │
│                                                                      │
│ Step 5: Microbiome Visualization Skill                               │
│                                                                      │
│ Skill: microbiome-visualization                                      │
│ Input: All computed outputs from Steps 1-4                           │
│ Goal: Publication-quality figures: alpha boxplots with significance  │
│ brackets, PCoA with ellipses, stacked bars, heatmaps, volcano plots, │
│  forest plots                                                        │
│ Output: High-DPI PNG/PDF figures in analysis_results/figures/        │
│ Update to script: Enhance existing plot_*() functions; add new       │
│ visualization functions                                              │
│                                                                      │
│ Step 6: Update microbiome_analysis.py                                │
│                                                                      │
│ Goal: Integrate the best patterns and functions from skill outputs   │
│ into the main pipeline                                               │
│ Changes:                                                             │
│ - Add new analysis functions (Kruskal-Wallis, Dunn's, DA, core       │
│ microbiome, meta-analysis)                                           │
│ - Add new output directories (differential_abundance/,               │
│ meta_analysis/, figures/)                                            │
│ - Update run_pipeline() to call new analysis steps                   │
│ - Update analysis_config.json with new parameters (permutations:     │
│ 999)                                                                 │
│ - Update requirements-microbiome.txt with new dependencies (scipy,   │
│ statsmodels)                                                         │
│                                                                      │
│ Step 7: Update Tests                                                 │
│                                                                      │
│ File: tests/test_microbiome_analysis.py                              │
│ - Add tests for new functions (Kruskal-Wallis, DA, core microbiome)  │
│ - Extend the tiny-dataset E2E test to cover new outputs              │
│ - Run full test suite to verify                                      │
│                                                                      │
│ Critical Files                                                       │
│ File: microbiome_analysis.py                                         │
│ Action: Extend with new analysis functions                           │
│ ────────────────────────────────────────                             │
│ File: analysis_config.json                                           │
│ Action: Add new parameters                                           │
│ ────────────────────────────────────────                             │
│ File: tests/test_microbiome_analysis.py                              │
│ Action: Add tests for new functions                                  │
│ ────────────────────────────────────────                             │
│ File: requirements-microbiome.txt                                    │
│ Action: Add scipy, statsmodels                                       │
│ ────────────────────────────────────────                             │
│ File: analysis_results/diversity/alpha_diversity.csv                 │
│ Action: Input to skills (read-only)                                  │
│ ────────────────────────────────────────                             │
│ File: analysis_results/taxonomy/genus_relative_abundance.csv         │
│ Action: Input to skills (read-only)                                  │
│ ────────────────────────────────────────                             │
│ File: Datesets/metadata.csv                                          │
│ Action: Input to skills (read-only)                                  │
│ Verification                                                         │
│                                                                      │
│ 1. Run updated tests: python -m pytest                               │
│ tests/test_microbiome_analysis.py -v                                 │
│ 2. Run smoke test: python microbiome_analysis.py --config            │
│ analysis_config.smoke.json                                           │
│ 3. Verify new output files exist in smoke results (kruskal_wallis,   │
│ DA, core microbiome, meta-analysis)                                  │
│ 4. Verify figures generate without errors                            │
│ 5. Spot-check statistical results for sanity (p-values in valid      │
│ range, effect sizes reasonable) 

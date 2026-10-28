# ExecPlan: follow the pasted methods text for generalists and specialists (16S)

> This ExecPlan is a living document maintained in accordance with ~/.claude/PLANS.md.
> Last updated: 2026-10-06 — IN PROGRESS. The user chose option A (everything in the text)
> for Paragraphs 1-3, then pasted three more: functional traits (4), gene content (5) and
> special classification groups (6), and on 2026-10-06 said "implement the plan", which I
> took as approval of the second decision too. See Progress for what is done.

## Purpose

Make the generalist/specialist analysis match the three methods paragraphs you pasted,
clause by clause, for 16S only (the paragraphs also mention ITS; there is none here).
Where the project currently does something different from the text, this plan says so
and says what it takes to close the gap.

## Context and Orientation

- The analysis lives in `seedling_analysis_master.qmd`: Part 1 chunk
  `generalist-specialist` computes the classes and tables; Part 2 chunk
  `fx-fig-tree-genus` draws the tree figures. `make_generalist_specialist.R` and
  `make_tree_figures.R` run those chunks alone.
- Classes come from the host "double screening" already implemented. At ASV level:
  **263 generalists, 1,120 specialists, 0 pseudo-generalists removed**; generalist
  threshold 6 host species.
- ASV sequences are in `rep-seqs.qza` (a zip holding `dna-sequences.fasta`). All 263
  generalist and all 1,120 specialist ASVs have a sequence there (checked).
- Homebrew is installed. None of MUSCLE, MAFFT, trimAl, IQ-TREE, FastTree or vsearch is
  installed, and neither is the R package `itol.toolkit`. conda is not installed.

## The pasted text, clause by clause

Status: **Done** = already matches the text. **Gap** = differs from the text now.

### Paragraph 1 — identification of generalists and specialists

| # | The text says | Status | What is done or needed |
| --- | --- | --- | --- |
| 1.1 | Two-stage screening on the host species-level presence/absence matrix | Done | ASV × host species matrix (16 hosts). |
| 1.2 | Stage 1: within-column permutation null preserving each host's richness T; 1,000 iterations; null distribution of cross-species prevalence | Done | 1,000 iterations (`screen_iterations`). |
| 1.3 | Threshold = smallest k at which the observed count first exceeds the null mean | Done, one reading to confirm | At k = 1 the observed count also exceeds the null (many single-host taxa), so the rule is applied above the expected breadth. Result here: **k = 6**. The "5" in the text is that study's value, not ours. |
| 1.4 | Stage 2: for each k ≥ 2, draw k host species 1,000 times; mean pairwise host phylogenetic distance (MPD); per-k null | Done | Host distances from your TimeTree file. |
| 1.5 | Standardised z-score, tested one-sided for dispersion and for clustering, each BH-FDR corrected | Done | Columns `host_mpd_z`, `q_clustered`, `q_dispersed`. |
| 1.6 | Generalist = n_hosts ≥ threshold and not significantly clustered (q_clus ≥ 0.05); specialist = n_hosts = 1; intermediate = the remainder | Done | Column `screen_class`. |
| 1.7 | Report pseudo-generalists removed and final counts | Done | 0 removed; 263 generalists; 1,120 specialists. The text's 32 / 1,103 / 5,830 must be replaced by these. |
| 1.8 | For the clustered taxa, tabulate the representative host genus (the genus most frequently hosting each taxon) | **Gap (small)** | The representative host genus is computed for every taxon, but no separate table of clustered taxa is written. Write `generalist_specialist_clustered_taxa.csv`. It will be empty at q < 0.05 (smallest q = 0.087); the file will say so. |

### Paragraph 2 — phylogenetic trees

| # | The text says | Status | What is needed |
| --- | --- | --- | --- |
| 2.1 | Separate maximum-likelihood trees for generalist and for specialist taxa | **Gap** | The two current figures prune the pipeline's existing tree. New trees must be inferred. |
| 2.2 | Generalists: MUSCLE (v5.3) → trimAl (v1.5, `-automated1`) → IQ-TREE under ModelFinder with 1,000 ultrafast bootstraps | **Gap** | Install MUSCLE, trimAl, IQ-TREE. Run on the 263 sequences. |
| 2.3 | Specialists: MAFFT (`--auto`) → trimAl → FastTree (`-nt -gtr -gamma`) | **Gap** | Install MAFFT, FastTree. Run on the 1,120 sequences. |
| 2.4 | Visualised in iTOL (v7) | **Gap** | iTOL is a website. I prepare the tree and annotation files; **you upload them**. I will not upload to an external service. |
| 2.5 | Annotation files generated with the `itol.toolkit` R package (v1.1.2) | **Gap** | Install `itol.toolkit` from CRAN. |
| 2.6 | Generalist tree: phylum/class colour strips | Gap in format only | Already in the R figures; to be written as iTOL colour-strip files. |
| 2.7 | BugBase-predicted trait rings: aerobic type, biofilm formation, Gram type | **Gap, uncertain** | See "BugBase" below. |
| 2.8 | Outer n_hosts bar | Gap in format only | Already in the R figure; to be written as an iTOL bar file. |
| 2.9 | Specialist tree annotated analogously | Gap in format only | Same files for the specialist tree. n_hosts is 1 for every specialist, so the bar is flat; I will add the host plant family as an extra strip and say that it is an addition. |

### Paragraph 3 — composition

| # | The text says | Status | What is done or needed |
| --- | --- | --- | --- |
| 3.1 | Fisher's exact tests of taxon representation across major ranks, BH-FDR | Done | `generalist_specialist_composition_tests.csv`, five ranks. |
| 3.2 | Taxa occurring exclusively in either group reported separately | **Gap (small)** | They are flagged in a column of the same table. Write them as their own table, `generalist_specialist_exclusive_taxa.csv`. |
| 3.3 | ASVs aggregated to genus; genera ranked by host breadth (n_hosts and n_hostg); the 25 broadest reported with mean relative abundance | Done | `generalist_specialist_top_genera.csv`. |
| 3.4 | Within the major ASV-rich phyla, host breadth against non-zero relative abundance by linear regression and one-sided Spearman | Done, one choice to confirm | The regression is on log10 abundance, which the text does not say. I will report both the raw and the log10 fit so the text is followed as written. "Major" = the five phyla with the most ASVs. |

### Paragraph 4 — functional traits and genes (added by the user after approval)

| # | The text says | Status | What is needed |
| --- | --- | --- | --- |
| 4.1 | All traits are in silico predictions from marker genes and are called "predicted traits" | To do | Naming convention in tables, figure legends and captions. |
| 4.2 | ASV × trait tables compiled from four resources | To do | One table per resource and one combined table, keyed by ASV. |
| 4.3 | BugBase (v0.1.0; 16S; Greengenes v13.8): aerobic type, biofilm formation, Gram type, stress tolerance, mobile elements, potential pathogenicity | **Gap, uncertain** | Same route as clause 2.7, now for all six trait groups. See "BugBase" below. |
| 4.4 | FAPROTAX (v1.2.12; 16S): ecological functions, including human association, plant pathogenicity, nitrate reduction, nitrogen respiration | **Gap, feasible** | FAPROTAX assigns functions from taxonomy names, which SILVA provides. Needs its database file downloaded. I will use the version the download gives and report it. |
| 4.5 | metaTraits (v0.1.14; 16S; GTDB-based): continuous traits (temperature, pH, salinity, genome size, gene count, GC content, coding density) and binary traits (carbon-source use, enzyme activity, motility, gas production, Gram staining, others) | **Gap, unknown to me** | I do not know this resource: not where it is obtained, nor how 16S ASVs are mapped to it. I will look for it; if I cannot find and verify it, I stop and ask you for the source. I will not substitute a different trait database. |
| 4.6 | FUNGuild (ITS) | Not applicable | ITS only; there are no fungal data here. |
| 4.7 | Trait columns classified as binary (TRUE/FALSE) or numeric | To do | Done when the combined table is built. |
| 4.8 | Binary traits analysed among annotated ASVs (NAs excluded) with binomial GLMs against log2(n_hosts + 1), Wald tests, BH-FDR | To do | `glm(trait ~ log2(n_hosts + 1), family = binomial)` per trait, on all ASVs with a value for it. |
| 4.9 | Traits with fewer than 100 annotated ASVs, or no variation, excluded | To do | Thresholds go in the config (`trait_min_annotated`: 100). |
| 4.10 | Effect sizes as rank-biserial correlations and odds ratios per doubling of host breadth, with 95% confidence intervals | To do | OR per doubling = exp(GLM coefficient), since the predictor is log2; Wald interval. Rank-biserial from the Mann-Whitney comparison of host breadth between ASVs with and without the trait. |
| 4.11 | Numeric traits tested by Spearman correlation with n_hosts, BH-FDR | To do | Per numeric trait, on ASVs with a value. |
| 4.12 | n_hosts capped at 20 to reduce long-tail influence | Done by construction | There are 16 host species here, so no ASV exceeds 16 and the cap changes nothing. The cap is still applied in code (`trait_max_hosts`: 20) and the methods should say it had no effect. |
| 4.13 | "115 binary traits and 15 continuous traits were screened" | To report | Those are the source study's counts. Ours are reported once the tables exist. |

### Paragraph 5 — gene content (KOs)

| # | The text says | Status | What is needed |
| --- | --- | --- | --- |
| 5.1 | KO presence/absence predicted for all 16S ASVs with PICRUSt2 v2.6.2 (`--max_nsti 0.15`) | **Gap, heavy** | PICRUSt2 is not installed and is normally installed through conda, which is not on this machine. Needs a conda distribution (Miniforge) and PICRUSt2 with its reference data, several GB. "All 16S ASVs" here = the 2,823 filtered ASVs. |
| 5.2 | KOs present in all or in no ASVs excluded | To do | Count reported (the text's 9,580 is that study's). |
| 5.3 | KO occurrence modelled against log2(n_hosts + 1) in the same logistic framework; likelihood-ratio tests; BH-FDR | To do | One GLM per KO, about ten thousand; minutes in R. |
| 5.4 | Wald odds ratios with 95% confidence intervals as effect sizes | To do | From the same fits. |
| 5.5 | Significant KOs (FDR < 0.05): generalist-enriched if log2 OR ≥ 0.5, specialist-enriched if log2 OR ≤ −0.5, weaker effects unassigned | To do | Thresholds in the config (`ko_fdr`, `ko_log2_or`). |
| 5.6 | KOs annotated to KEGG level-2 categories | **Gap** | Needs the KEGG KO hierarchy, downloaded from KEGG. KEGG's terms allow academic use; you would need to confirm that applies. |
| 5.7 | and to five-level plant-growth-promoting gene (PGPG) functions from the PGPg_finder database | **Gap, to verify** | I know the tool exists but have not seen the layout of its annotation table. First step is to fetch and check it. |

### Paragraph 6 — special classification groups

| # | The text says | Status | What is needed |
| --- | --- | --- | --- |
| 6.1 | (i) Mycorrhizal fungi from FUNGuild (ITS) | Not applicable | No fungal data. Three groups here, not four. |
| 6.2 | (ii) Anaerobic bacteria = BugBase "Anaerobic == TRUE" | Depends on 4.3 | Needs BugBase to work. |
| 6.3 | (iii) Nitrogen-fixing bacteria = FAPROTAX "nitrogen_fixation == TRUE" | Depends on 4.4 | Feasible once FAPROTAX is run. |
| 6.4 | (iv) Small-genome bacteria = the smallest 1% of non-NA metaTraits genome sizes | Depends on 4.5 | Needs metaTraits, which I have not located. |
| 6.5 | For each group: ASV counts, sample-level relative-abundance distributions, order-level taxonomic composition, generalist/specialist composition | To do | Tables per group. |
| 6.6 | Host loading: the number of group-member genera detected per host genus, and its correlation with group relative abundance | To do | One table per group. |
| 6.7 | Group ASVs annotated to the finest resolvable taxonomic level and aggregated into a label × host genus mean relative-abundance matrix, keeping labels detected in ≥ 3 host genera | To do | Finest level = genus where assigned, else family, order and so on. |
| 6.8 | Matrix log10-transformed, row-scaled (z-score), double-clustered (ward.D2) | To do | One heatmap per group, drawn to the layout rules below. A pseudocount is needed before log10; the text does not give one, so I will state the one used. |
| 6.9 | Per-genus × host-order Kruskal-Wallis tests, BH-FDR | To do | Host plant order is added to the host species table (Brassicales, Fabales, Apiales, Asterales, Solanales, Poales). |

These three paragraphs add three milestones, after M4 (the BugBase rings on the trees and
the BugBase trait columns are the same data):

- **M6 — predicted traits** (Paragraph 4): ASV × trait tables, binomial GLMs and Spearman
  tests. Outputs `predicted_traits_asv.csv`, `predicted_trait_breadth_tests.csv`.
- **M7 — gene content** (Paragraph 5): PICRUSt2, KO filtering, per-KO GLMs, enrichment
  classes, KEGG and PGPG annotation. Outputs `ko_breadth_tests.csv`, `ko_enrichment_summary.csv`.
- **M8 — special groups** (Paragraph 6): per-group tables and one clustered heatmap per
  group. Outputs `special_group_*.csv`, `figures/special_group_<name>_heatmap.tiff`.

What can be delivered depends on which resources can actually be obtained:

| Resource | Used by | Standing |
| --- | --- | --- |
| FAPROTAX | 4.4, 6.3 | Should work: assigns from taxonomy names. |
| BugBase + Greengenes 13.8 | 2.7, 4.3, 6.2 | Uncertain: file format not yet verified. |
| metaTraits | 4.5, 6.4 | Unknown to me: source not located. |
| PICRUSt2 | 5.1-5.5 | Heavy: needs conda and several GB. |
| KEGG hierarchy | 5.6 | Download; licence to confirm. |
| PGPg_finder annotations | 5.7 | To verify. |

### Paragraph 7 — statistical analysis statement (added 2026-10-06; to do AFTER M1-M8)

The user pasted the source text's "Statistical analysis" paragraph and asked to "update
stats ... after the previous task". To do once the trees, traits and gene content are
finished: write the statistical-analysis statement for THIS project, clause by clause
against the pasted paragraph, stating what was actually used. Known differences already:

| The text says | Here |
| --- | --- |
| R v4.5.2 | R 4.3.2 |
| MRM used ecodist | MRM is implemented in `host_phylosymbiosis.R` (ecodist is not installed); either install ecodist and switch, or say so |
| Hierarchical clustering and ARI; Cramer's V; chi-squared with Monte Carlo simulation; within-bin permutation | Not part of any paragraph pasted so far; no such analysis exists here yet |
| βNTI / RCbray with the iCAMP workflow on a server | Not done here; not requested so far |
| Heatmaps with ComplexHeatmap / pheatmap | Heatmaps are drawn with ggplot2 (both packages are installed if a switch is wanted) |
| 999 or 1,000 permutations; λ bootstrap B = 500 | Same |
| Sample sizes 4,107 (16S) and 2,367 (ITS) | 1,886 samples (16S) after quality control; no ITS |
| p < 0.001 reporting rule; α = 0.05; BH-FDR | Same |

## Second decision needed (Paragraphs 4-6)

Your approval of option A covered Paragraphs 1-3: the tree software, `itol.toolkit`,
BugBase and Greengenes. Paragraphs 4-6 need more, and I have not started any of it:

- **Install a conda distribution (Miniforge) and PICRUSt2**, several GB on disk.
- **Download FAPROTAX, the KEGG KO hierarchy and the PGPg_finder annotations.**
- **metaTraits:** tell me where it comes from (a link or a citation), or let me search
  for it and report back before using anything.

This is several sessions of work, not one. A sensible order is FAPROTAX first (cheap and
likely to work), then BugBase, then PICRUSt2, with metaTraits whenever its source is known.

## BugBase (clauses 2.7 and 4.3): what I know and do not know

BugBase predicts traits for OTUs that carry **Greengenes** identifiers. This project's
taxonomy is SILVA and its units are ASVs, so the ASVs must first be matched to Greengenes
OTUs (closed-reference, 97% identity, with vsearch), and the trait of the matched OTU
read from BugBase's precalculated trait table.

What I have **not** verified: the exact name, location and format of BugBase's
precalculated per-OTU trait file, and which Greengenes release it is keyed to. The first
step is to fetch the BugBase repository and check. If the file is not usable that way, I
will stop and tell you; I will not substitute another trait source without asking.

ASVs with no Greengenes match at 97% get "unassigned" in the trait rings, and the number
of them is reported.

## Plan of Work

Order matters: the gaps that need nothing installed first, then installs, then trees, then
annotation.

1. **Small gaps, no installs.** In chunk `generalist-specialist`: write
   `generalist_specialist_clustered_taxa.csv` (1.8) and
   `generalist_specialist_exclusive_taxa.csv` (3.2); add the raw-abundance regression
   columns beside the log10 ones (3.4).
2. **Install the tree software** with Homebrew: MUSCLE, MAFFT, trimAl, IQ-TREE, FastTree,
   vsearch. Record the exact version of each. If Homebrew offers a different version from
   the one in the text (for example FastTree 2.1 instead of 2.2.0), I use what installs and
   report it, because the manuscript must cite the versions actually run.
3. **New script `generalist_specialist_trees.sh`** (config-driven paths, no hardcoded
   thresholds): extract the two sets of sequences from `rep-seqs.qza`; run 2.2 and 2.3
   exactly as written; write alignments, trimmed alignments, trees and logs to
   `analysis_results/trees/`.
4. **Install `itol.toolkit`**; new script `generalist_specialist_itol.R` writes, for each
   tree, the iTOL annotation files (phylum strip, class strip, n_hosts bar, and the trait
   rings if BugBase works) into `analysis_results/itol/generalists/` and `.../specialists/`,
   with a short README saying which file to drag onto which tree in iTOL.
5. **BugBase** as described above, feeding step 4.
6. **R figures redrawn on the new trees.** The two existing figures
   (`tree_host_generalists.tiff`, `tree_host_specialists.tiff`) switch from the pruned
   pipeline tree to the new ML trees, with bootstrap support shown on the generalist tree,
   so there is a reproducible figure in the project as well as the iTOL files.

## Figure layout and legends (your instruction)

Applied already to the two current figures, and to be kept for every redraw:

- Sized to the journal page: 7.2 in (183 mm) wide; the tree fills the width.
- No dead bands: where a tree has no tip names (specialists), the subtitle and legends are
  pulled into the blank ring a fan layout leaves.
- Legends strong and complete: 8.5 pt text, 9.5 pt bold titles, outlined keys, and **a key
  for every ring and bar** (the outer bars and the host-family rings now have their own).
- Every legend inside the canvas; a class has the same colour in both trees.

## Milestones

M1-M5 cover Paragraphs 1-3; M6-M8 are described under Paragraph 6 above.

1. **M1 — small gaps.** Run `Rscript make_generalist_specialist.R analysis_config.json`.
   Verify: two new tables exist; existing columns of the genus and ASV tables unchanged;
   raw and log10 regression columns both present.
2. **M2 — trees.** Run `generalist_specialist_trees.sh`. Verify: each tree has exactly 263
   and 1,120 tips; IQ-TREE log names the selected model and reports 1,000 bootstraps;
   trimmed alignment length and gap fraction are reported (see Risks).
3. **M3 — iTOL package.** Verify: annotation files load without error in
   `itol.toolkit`; every tip in a file exists in its tree; README written.
4. **M4 — traits.** Verify: number of ASVs matched to Greengenes; trait rings written, or
   a clear report that BugBase could not be used.
5. **M5 — figures.** `Rscript make_tree_figures.R analysis_config.json` reports OK for all
   seven; both figures checked by eye against the layout rules above.
6. `python3 -m pytest tests/ -q` still passes (20 tests).

## Risks

- **The ASVs come from eight studies that sequenced different regions of the 16S gene.**
  A new alignment of fragments that do not overlap gives unreliable deep branches, and
  trimAl may leave very few columns. The text's method is followed as written, but if the
  trimmed alignment is shorter than 100 columns I stop and report before building trees.
  The pipeline's existing tree has the same limitation.
- **IQ-TREE with ModelFinder and 1,000 bootstraps on 263 short sequences** should take
  minutes, not hours; FastTree on 1,120 is fast.
- **Version mismatch with the text** (see step 2).
- **BugBase may not be usable** (see above).

## Interfaces and Dependencies

- New software (needs your approval): MUSCLE, MAFFT, trimAl, IQ-TREE, FastTree, vsearch
  via Homebrew; R package `itol.toolkit` from CRAN.
- External downloads (needs your approval): the BugBase repository and the Greengenes
  reference sequences it is keyed to.
- External service: iTOL, used by you, not by me.
- New config keys in `generalist_specialist`: `tree_bootstraps` (1000),
  `greengenes_identity` (0.97).
- New files: `generalist_specialist_trees.sh`, `generalist_specialist_itol.R`, outputs
  under `analysis_results/trees/` and `analysis_results/itol/`.

## Decision needed

How far to go:

- **A. Everything in the text**: installs, new ML trees, iTOL files, BugBase attempt.
- **B. Trees and iTOL files, no BugBase**: installs except vsearch; trait rings left out
  and stated as omitted.
- **C. No installs**: only the small gaps (M1); trees stay pruned from the pipeline tree
  and drawn in R, and the methods text is edited to say so.

## Progress

- [x] 2026-10-06 Decision: option A, everything in the text
- [x] 2026-10-06 M1 small gaps: `generalist_specialist_exclusive_taxa.csv` (471 taxa),
      `generalist_specialist_clustered_taxa.csv` (no rows: nothing clustered at q < 0.05),
      raw-abundance regression columns added
- [x] 2026-10-06 M2 trees: `generalist_specialist_trees.sh` (88 min, of which IQ-TREE 85).
      Generalists 263 tips, TIM3+I+R3, 291 columns, no gaps. Specialists 1,120 tips,
      466 columns, 26% gaps.
- [x] 2026-10-06 M3 iTOL package: `generalist_specialist_itol.R`; 6 annotation files for
      generalists and 7 for specialists in `analysis_results/itol/`; every row names a tip.
- [x] 2026-10-06 M4 traits: BugBase rings written into the iTOL files (see M6).
- [x] 2026-10-06 M5 figures: both ASV tree figures redrawn on the ML trees; seven tree
      figures render; 20 Python tests pass.
- [x] 2026-10-06 M6 predicted traits: `predicted_traits.py` + `predicted_trait_tests.R`.
      1,265 predicted traits (BugBase 9, FAPROTAX 92, metaTraits 1,164); 333 binary and
      29 numeric tested; 130 and 11 significant at FDR < 0.05.
- [x] 2026-10-06 M7 gene content: `ko_predict.sh` (PICRUSt2 2.6.2, 56 min) and
      `ko_breadth_tests.R`. 2,755 of 2,823 ASVs within NSTI 0.15; 7,883 of 13,793 KOs
      tested; 1,036 generalist-enriched, 1,451 specialist-enriched, 5,396 unassigned.
- [x] 2026-10-06 M9 statistical-analysis statement, and a methods paragraph for every
      new analysis, written into `analysis_results/methodology.md` at the user's request
      (ahead of M2 and M7, whose two paragraphs are marked "Run in progress" there)
- [x] 2026-10-06 M8 special groups, tables: `special_groups.R`. Anaerobic 154 ASVs,
      nitrogen-fixing 23, small-genome 22. Heatmap figures still to draw.

## Surprises and Discoveries

- 2026-10-06 — Homebrew core has MAFFT 7.526, IQ-TREE 3.1.4, trimAl 1.5.1 and vsearch
  2.32.0, but **no MUSCLE and no FastTree**. Those two have to come from their authors'
  own releases. Versions therefore differ from the text (IQ-TREE 3.1.4, not 3.0.1;
  trimAl 1.5.1, not 1.5); the manuscript must cite what was run.
- 2026-10-06 — **Homebrew no longer supports Intel Macs** and started building from
  source and upgrading cmake and xz, so that install was stopped. All six programs were
  taken from their own release pages into `~/.local/share/seedling_tools/` (outside the
  project; README there). This gave the exact versions in the text: MUSCLE 5.3,
  IQ-TREE 3.0.1, trimAl 1.5, FastTree 2.2.0, plus MAFFT 7.526 and vsearch 2.32.0.
  MUSCLE's binary needed its compiler runtime re-pointed to the installed gcc.
- 2026-10-06 — The 263 generalist ASVs are all 291-293 bp: one 16S region. Their
  alignment is 293 columns with 0.4% gaps, so the generalist tree does not have the
  mixed-region problem. The specialists span 252-471 bp and may.
- 2026-10-06 — **BugBase works**: its repository (v.0.1.0) holds a per-OTU trait table;
  2,754 of 2,823 ASVs match a Greengenes 13.8 97% OTU at >= 97% identity. Five of its
  nine columns are 0/1 (binary); the other four (biofilm, mobile elements, pathogenic,
  stress tolerant) are fractions and are treated as numeric, as the text's rule requires.
- 2026-10-06 — **metaTraits found**: metatraits.embl.de (Bork group, NAR 2026). ASVs are
  matched to its GTDB genus-level summary by genus name: 2,179 ASVs in 201 genera; 82
  genus names not found. The download carries no version number or licence statement.
- 2026-10-06 — itol.toolkit installed as 1.1.7 (the text cites 1.1.2).
- 2026-10-06 — No special group contains a single host generalist ASV.
- 2026-10-06 — **IQ-TREE's ultrafast bootstrap did not converge** in 1,000 iterations
  (correlation 0.981 against a target of 0.99); 54 of 260 branches have support >= 95%.
  The tree is used, with support described as approximate. A rerun with `-nm 2000`
  or more would settle it (about three hours).
- 2026-10-06 — itol.toolkit accepts a data frame for colour strips but fails on one for
  gradient and bar datasets, which it only takes as a tibble.
- 2026-10-06 — On the abundance itself (not log10), host breadth does not predict
  abundance in Proteobacteria (slope about zero, p = 0.96), while on log10 it does
  (p < 0.001). The two fits answer different questions; both are in the table.

## Decision Log

- 2026-10-06 — Written after the user asked for a plan in which "the pasted text is
  followed clearly". The earlier plan (`generalist-specialist-trees-composition.md`)
  adapted the text to what was installed; this one closes the gaps.

## Audit of the methodology against the code (2026-10-07)

The user asked whether `analysis_results/methodology.md` matches the code. Six mismatches
were found and the text corrected. One was in the analysis itself:

- **Stage two of the host screening could never call a taxon clustered.** It drew k hosts
  1,000 times, so the smallest possible p was 0.001, which cannot survive BH correction
  over 1,703 taxa. The pasted text reads as a z-score test; a normal approximation to z
  calls 195 ASVs clustered, but the null is skewed and discrete (120 host pairs at k = 2).
- **Decision (user: "address the mismatches"; my stated recommendation, not contradicted):**
  the null is now exact. Every set of k hosts is listed (at most 12,870) and p is the
  share of sets with an MPD no greater than observed. Result: 29 ASVs clustered, 22
  pseudo-generalists removed, **241 generalist ASVs** (was 263); genera unchanged (53).
  The z-normal q is kept in the tables as `q_clustered_z` for comparison.
- Consequences rerun: composition tables, special groups, the generalist tree (rebuilt on
  241 sequences), iTOL files, figures. Trait and KO tests use host breadth, not class,
  and are unaffected.
- Other text corrections: 16 host species, not 17 (the 17th was soil); `r.g` is the
  point-biserial coefficient, not phi; not all work is in R; the small-genome heatmap is
  not drawn; Faith's PD and UniFrac were computed on 5 October but are not used.

- 2026-10-07 — Generalist tree rebuilt on the 241 sequences (29 min; TIM3+I+R3; the
  bootstrap converged this time, after 800 iterations; 49 of 238 branches >= 95%). iTOL
  files rewritten, all seven tree figures redrawn, methodology and legends updated.

- 2026-10-07 — **One-tree main figure, at the user's request.** `tree_generalist_specialist.tiff`
  is now an ASV-level figure of both classes on one tree (1,354 tips), built on a third
  tree, `trees/combined/` (MAFFT, trimAl, FastTree; 251 columns, no gaps). Rings: class,
  phylum, BugBase oxygen use, six host-family rings; outer bars for host breadth. 600 dpi,
  ragg, no title in the artwork. The genus-level figure with the CLAM stage ring is kept
  as `tree_generalist_specialist_genus.tiff`. `generalist_specialist_trees.sh` takes a
  second argument naming which trees to build, so the combined tree does not need IQ-TREE
  rerun.

- 2026-10-07 — **Stage effect within strata, a fifth way of testing taxon-stage association.**
  Built first with study as the stratum (my recommendation), then changed at the user's
  instruction to plant family and cultivar ("make plant family or cultivar the variable not
  study"). Chunks `stage-effects-by-stratum` (Part 1) and `fx-fig-stage-forest` (Part 2);
  launcher `make_stage_forest.R`. Plant families pooled with REML + Knapp-Hartung;
  cultivars with a three-level model (study / cultivar), because 11 of 12 cultivars sit
  inside one study. No genus is significant after FDR either way; Bacillus is the most
  consistent. The study-stratified tables were removed as superseded.

- 2026-10-07 — Stage-effect analysis: cultivars taken within host species ("sativus" was
  both a carrot and a radish); a legend per forest panel; the methodology section cut to
  what was done, at the user's request; and a leave-one-stratum-out sensitivity analysis
  added (tables `stage_effect_leave_one_out_by_*_genus.csv`, summary columns `loo_*` in
  the pooled tables). The user's "implement how best it represents the analysis
  (sensitivity)" was read as a request for that sensitivity analysis.

- 2026-10-07 — At the user's request the forest plot became a heatmap
  (`stage_effect_heatmap.tiff`, chunk `fx-fig-stage-heatmap`, launcher renamed
  `make_stage_effects.R`), and `analysis_results/methodology.md` was cut to a concise
  version (about 1,300 words); the previous full text is kept as `methodology_detailed.md`.

- 2026-10-07 — **Taxon-stage associations reduced to two steps**, on the user's approval of
  my recommendation: screen (DESeq2 and IndVal.g) then the stratified meta-analysis. The
  top-five genus correlations were removed from the master (Part 1 chunk
  `top-5-genus-correlation-stats`, Part 2 chunk `fx-fig-correlation-heatmaps`, two
  captions); the removed code and its five output files are in
  `tasks/plans/pruned-2026-10-07/`. The co-occurrence network is kept, and described with
  composition instead of under stage associations (the user did not say whether to drop it).

- 2026-10-07 — "Update the code and the figures", read as bringing both into line with the
  two-step design: the screen is saved (`tables/stage_screen_asvs.csv`) and drawn
  (`stage_screen.tiff`, chunk `fx-fig-stage-screen`); the DESeq2 count printed by the
  pipeline now applies the fold-change threshold as well as the adjusted p; the indicator
  plot's axis and caption said point-biserial r but the values are IndVal.g, corrected.
  Found on the way: DESeq2 assigns an adjusted p to only 602 of 2,823 ASVs.

- 2026-10-07 — **Master brought up to date** at the user's request. Until now nine analyses
  ran only from their own scripts, so a render of the master did not produce the tables
  its Part 2 figures read. New Part 1 chunk `external-analyses` runs them in dependency
  order (host_phylosignal, host_phylosymbiosis, generalist_specialist_trees,
  predicted_traits, predicted_trait_tests, ko_predict, ko_breadth_tests, special_groups,
  generalist_specialist_itol). Skips, never fails, on a missing script, input, program
  or package; reruns a step only when the MD5 of its script or inputs changes
  (`cache/external_<step>.md5`, log beside it); adopts outputs found without a record.
  Config: `external_analyses.run` (false in the smoke config) and `.force`. Paths tested
  outside a render: adopted, up to date, rerun on missing output, rerun on changed script,
  not run. The master's overview list and Part 2 input table were updated to match.

## Outcomes and Retrospective

All nine milestones are done for 16S. Not done, and why:

- The R figures do not carry the BugBase trait rings; those are in the iTOL files only, as
  the text places them there.
- The two-panel generalist/specialist figure is not built; the user has not asked for it.
- iTOL itself: the files are written, the upload is the user's.
- metaTraits' version could not be confirmed (the download carries none).
- The full master render was not run.

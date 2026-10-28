> This ExecPlan is a living document maintained in accordance with ~/.claude/PLANS.md.
> Last updated: 2026-10-06 - EXECUTED. Approved by the user in plan mode on 2026-10-06.

# Plan: trees and composition of host generalists and specialists (16S)

## Context

You pasted two methods paragraphs to follow on from the host generalist/specialist
"double screening": (1) phylogenetic trees of the generalist and the specialist taxa with
annotation rings, and (2) a compositional comparison of the two groups. This plan adapts
both to this project. ITS is ignored throughout; the data here are 16S only.

**Where the double screening stands (already run, before plan mode):**

| Level | Generalist threshold | Generalists | Intermediate | Specialists (1 host) | Pseudo-generalists removed |
| --- | --- | --- | --- | --- | --- |
| ASV | 6 hosts | 263 | 1,440 | 1,120 | 0 |
| Genus | 7 hosts | 53 | 88 | 141 | 0 |

No taxon is significantly clustered on the host tree after FDR correction (smallest
q = 0.087), so the second screen removes nothing. An earlier figure I had on screen
(threshold 9) came from a bug that miscounted observed breadth; it is fixed and the table
above is the corrected result. These numbers are in
`analysis_results/tables/generalist_specialist_summary.csv`.

The two new paragraphs work on **ASVs** ("sequences", "ASVs were aggregated"), so this
plan uses the ASV-level classes: 263 generalists and 1,120 specialists.

## What cannot be done as the paragraphs describe, and what replaces it

| Paragraph says | Situation here | Plan |
| --- | --- | --- |
| New maximum-likelihood trees (MUSCLE, trimAl, IQ-TREE, MAFFT, FastTree) | None of these programs is installed | Prune the pipeline's existing ASV tree (`Datesets/tree.nexus`, built by QIIME2 from the same sequences) to the generalist ASVs and to the specialist ASVs. Same relationships, no new inference. This is what the other three tree figures already do. |
| Visualised in iTOL with `itol.toolkit` | `itol.toolkit` is not installed; iTOL is an external website | Draw both trees in R with `ggtree`/`ggtreeExtra`, as the project's other tree figures are, so they are reproducible from the master document. |
| BugBase trait rings (aerobic type, biofilm, Gram type) | BugBase needs Greengenes OTU IDs and an external tool; nothing in the project provides traits | **Omitted.** No trait rings. Stated in the caption and the report. |

If you want any of the three done as written, it needs software installed or an external
service, which I would ask you about separately. Say so when you review this plan.

## Part A — two tree figures (Part 2 of the master)

File: `seedling_analysis_master.qmd`, chunk `fx-fig-tree-genus`, two new `fig_section`s
after `tree_generalist_specialist`. `make_tree_figures.R` then draws them with no change
beyond its header comment.

Reads only persisted outputs: `setup/ps_clean.rds`, the tree (`TREE_FILE`), and
`tables/generalist_specialist_asv_host_screen.csv` (class, `n_hosts`,
`representative_host_genus` per ASV).

Reused, already in that chunk: `TREE_FILE`, `TREE_EXCLUDE_GENERA`, `TREE_PHYLUM_COLOURS`,
`TREE_TAXON_COLOURS`, the tip-label quote stripping, midpoint rooting, and the
ring-drawing pattern from `tree_generalist_specialist` (ring data must not share column
names with the tree data; `ggnewscale` between fills). `HOST_FAMILY_COLOURS` and
`tables/host_plant_species.csv` give host plant family.

1. **`tree_host_generalists.tiff`** — the 263 generalist ASVs.
   - Fan dendrogram, branch lengths not drawn.
   - Ring 1: phylum. Ring 2: class (the most frequent classes named, the rest "Other").
   - Outer bars: `n_hosts`, with a radial axis.
   - Tips of the most abundant ASVs labelled with their genus.
2. **`tree_host_specialists.tiff`** — the 1,120 specialist ASVs, annotated the same way.
   - Ring 1: phylum. Ring 2: class.
   - Outer ring: the one host each specialist was found on, coloured by host plant
     family (six colours), because `n_hosts` is 1 for every tip and a bar would say nothing.
   - Too many tips to label.

Standing rules from `tasks/lessons.md` applied from the start: sized to the journal page
(7.2 in wide), text at least 7 pt with bold legend titles, every legend inside the canvas
with outlined keys, solid and clearly distinct colours, no colour meaning two things,
navy and vermillion kept for seeds and seedlings only.

A caption for each goes in the caption registry at the top of the master.

## Part B — compositional comparison (Part 1 of the master)

File: `seedling_analysis_master.qmd`, chunk `generalist-specialist`, inside
`perform_generalist_specialist()` after the double screening, so
`make_generalist_specialist.R` produces everything in under a minute. All at ASV level.

1. **Fisher's exact tests by rank** → `tables/generalist_specialist_composition_tests.csv`
   - For each rank (phylum, class, order, family, genus) and each taxon in it: generalist
     vs specialist ASVs, inside vs outside the taxon; odds ratio, p, and BH-FDR within rank.
   - Column `exclusive_to` marks taxa found only among generalists or only among specialists.
2. **Genera ranked by host breadth** → `tables/generalist_specialist_top_genera.csv`
   - ASVs aggregated to genus; `n_hosts` (host species) and `n_hostg` (host genera);
     mean relative abundance across all samples. The 25 broadest genera.
3. **Host breadth against abundance** → `tables/generalist_specialist_breadth_abundance.csv`
   - Per ASV: `n_hosts` and non-zero relative abundance (mean over the samples where it occurs).
   - Within each of the phyla richest in ASVs: linear regression of log10 abundance on
     `n_hosts` (slope, R², p) and a one-sided Spearman correlation (does abundance rise
     with breadth), with BH-FDR across phyla.

Two new settings in the `generalist_specialist` section of `analysis_config.json` and
`analysis_config.smoke.json`, not hardcoded: `top_genera` (25) and `major_phyla` (5).

Summary rows for each are appended to `generalist_specialist_summary.csv`.

## Files touched

- `seedling_analysis_master.qmd` — Part 1 chunk `generalist-specialist`; Part 2 chunk
  `fx-fig-tree-genus`; two captions.
- `analysis_config.json`, `analysis_config.smoke.json` — two settings.
- `make_tree_figures.R` — header comment only.
- `tasks/plans/generalist-specialist-trees-composition.md` — a copy of this plan with
  progress and outcomes, per the project's ExecPlan rule.
- New outputs: three tables and two figures listed above.

Not touched: `analysis_results/manuscript/`, the validation path, the existing columns of
`generalist_specialist_genus.csv`, the existing tree figures and their tables, `CLAUDE.md`.

## Order of work

1. Part B tables (fast, and Part A's figure does not depend on them).
2. Part A generalist tree, then specialist tree, checking each image by eye.
3. Captions, launcher header, ExecPlan record.

## Verification

- `Rscript make_generalist_specialist.R analysis_config.json` completes; the existing
  columns of `generalist_specialist_genus.csv` are identical before and after.
- Fisher table: for every rank, generalist counts sum to 263 and specialist counts to 1,120.
- Top-genera table has 25 rows; `n_hostg` never exceeds `n_hosts`.
- `Rscript make_tree_figures.R analysis_config.json` reports OK for all seven figures; the
  three existing tree tables are byte-identical before and after.
- Tip counts: 263 and 1,120, less any ASV of the excluded genus (*Akkermansia*) or absent
  from the tree; the exact numbers go in the subtitles and are reported.
- Both figures opened and checked: every ring matches its legend, legends fully visible.
- `python3 -m pytest tests/ -q` still passes (20 tests).
- The full master render is not run (it takes hours); the new code is exercised through
  the two launcher scripts, which evaluate the master's own chunks.

## Things to know before approving

- **The trees are pruned, not re-inferred**, and carry **no BugBase traits**. This is the
  main departure from your paragraph.
- **Host is largely tied to study.** A specialist ASV may simply belong to a study that
  sequenced one host. This applies to every result here and will be stated in the captions.
- **"Major ASV-rich phyla" is set to the five with the most ASVs**; change `major_phyla`
  if the manuscript names a different set.
- **Nothing is committed**, as before.

## Progress

- [x] 2026-10-06 Part B tables (Fisher tests, top genera, breadth against abundance)
- [x] 2026-10-06 Part A generalist tree and specialist tree
- [x] 2026-10-06 Captions, launcher header, this record

## Surprises and Discoveries

- The first double-screening run reported a generalist threshold of 9 hosts. Inside
  `tibble()`, `observed = tabulate(n_hosts, ...)` read the new `n_hosts` column (1 to 16),
  not the per-taxon breadth, so every observed count was 1. Counting before the table is
  built gives the correct thresholds: 6 hosts (ASV) and 7 (genus).
- No taxon is significantly clustered on the host tree after FDR correction (smallest
  q = 0.087), so the second screen removes no pseudo-generalists at either level.

## Decision Log

- 2026-10-06 - Trees pruned from the pipeline's ASV tree, not re-inferred; drawn with
  ggtree, not iTOL; no BugBase rings. None of MUSCLE, MAFFT, trimAl, IQ-TREE, FastTree or
  itol.toolkit is installed, and nothing in the project provides traits.
- 2026-10-06 - **Departure from the approved plan:** the specialist tree shows host plant
  family as six thin rings marked by position, not as six colours. With seven phylum
  colours and eight class colours already in the figure, six more could not be kept
  clearly distinct, which `tasks/lessons.md` requires.
- 2026-10-06 - Class colours are assigned from the two trees together, so a class has the
  same colour in both figures.
- 2026-10-06 - On the generalist tree each genus is named once, on its most abundant ASV;
  naming the 30 most abundant ASVs repeated *Pantoea* and *Pseudomonas* many times.

## Outcomes and Retrospective

ASV level: 263 host generalists, 1,120 host specialists (1,113 on the tree after the
excluded genus and ASVs absent from the tree are dropped).
Composition: specialists include 15 phyla with no generalist at all; Bacteroidota (0
generalists, 98 specialists) and Acidobacteriota (0, 65) are the strongest contrasts, and
Proteobacteria and Firmicutes are over-represented among generalists.
Breadth against abundance: ASVs found on more hosts are more abundant where they occur in
Proteobacteria and Actinobacteriota (regression and Spearman both significant), weakly in
Firmicutes and Acidobacteriota (Spearman only), and not in Bacteroidota.

Not done: the full master render; BugBase traits; new ML trees; iTOL files.

# The five figures of the manuscript

Each figure is merged with patchwork from finished figures of the analysis, at the full
column width (183 mm) and 600 dpi. Files: `Figure_1.tiff` to `Figure_5.tiff` in this folder.
The arrangement follows the figure plan of `CELL-REPORTS-D-26-05893.pdf`, with a figure of
composition and beta diversity added as Figure 2 and the developmental-stage figure in the
place of that manuscript's clustering figure. The colours are that manuscript's too, sampled
from its figure pages and set in the `palette` section of `analysis_config.json`: seeds are
orange and seedlings teal throughout.

The panels are listed in `main_figures.csv` and set in the `main_figures` section of the
config. Rebuild with `Rscript make_main_figures.R`.

## Figure 1. Overview of the compiled dataset

*File:* `Figure_1.tiff` (study_overview)

**Overview of the compiled dataset.** The dataset comprises 1,886 samples from eight studies after filtering. **A**, Countries of origin, coloured as in B. **B**, Share of samples by country. **C**, Share of samples from seeds and from seedlings. **D**, Samples per study, by developmental stage (seeds / seedlings). **E, F**, Samples per host plant genus among seeds (E) and seedlings (F). Seeds are orange and seedlings teal throughout.

## Figure 2. Composition and beta diversity of the seed and seedling microbiota

*File:* `Figure_2.tiff` (composition_dendrogram + pcoa_css_by_stage)

**Composition and beta diversity of the seed and seedling microbiota.** **A**, Circular dendrogram of the host plant genera at each developmental stage (at least three samples), clustered on mean genus-level composition (Bray–Curtis, average linkage) separately for seeds and seedlings. Rings, from the centre outwards: stage, plant family and host genus. Bars: mean relative abundance (%) of bacterial classes (inner bars) and of bacterial genera (outer bars). The number beside each tip is its number of samples. **B–E**, Principal coordinates analysis of Bray–Curtis dissimilarities (cumulative sum scaling, genus level), computed separately for seeds (B, C) and seedlings (D, E) and coloured by host genus (B, D) and plant family (C, E), in the colours of A. R² and p of the PERMANOVA of each variable within each stage are given on the panels (999 permutations). PERMANOVA F was host genus, 11.7 (seeds) and 9.5 (seedlings); plant family, 15.6 (seeds) and 8.5 (seedlings).

## Figure 3. Alpha diversity: host phylogeny and developmental stage

*File:* `Figure_3.tiff` (host_phylosignal + alpha_stage_effect_heatmap)

**Alpha diversity: host phylogeny and developmental stage.** **A**, Time-calibrated phylogeny of the 16 host plant species (TimeTree 5), coloured by plant family, with the mean observed richness, Shannon index and Pielou's evenness of each host at each developmental stage after rarefaction to 500 reads. Point size gives the number of samples; hosts with fewer than three samples at a stage carry no point. Blomberg's *K* (999 permutations) and Pagel's λ are given above each panel. **B, C**, Effect of developmental stage on alpha diversity. Cells are Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. Observed richness, Shannon diversity and Pielou's evenness at 500 reads, at 1,000 reads and unrarefied, within plant families (B) and within host genera (C).

## Figure 4. Community structure and genera: host phylogeny and developmental stage

*File:* `Figure_4.tiff` (host_phylosymbiosis_scatter + stage_effect_heatmap)

**Community structure and genera: host phylogeny and developmental stage.** **A**, Divergence time between pairs of host species against the mean Bray–Curtis dissimilarity between their samples, at genus and at ASV level within each stage. Shade gives the share of sample pairs from different studies; the line is the least-squares fit. Mantel r and the R² of the multiple regression on distance matrices are given with permutation p-values (999 permutations). **B, C**, Effect of developmental stage on genera. Cells are Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. Centred log-ratio abundance of the 10 genera that hold an ASV significant in both DESeq2 and the indicator analysis, within plant families (B) and within host genera (C); the strip above the columns gives the plant family.

## Figure 5. Host generalists and specialists: phylogeny, host breadth and predicted functions

*File:* `Figure_5.tiff` (tree_host_generalists | tree_host_specialists; host_breadth_top_genera; trait_breadth_forest | ko_function_bars)

**Host generalists and specialists: phylogeny, host breadth and predicted functions.** **A**, Maximum-likelihood phylogeny of the 241 host generalist ASVs (at least six host species). Rings, from the centre outwards: phylum and class. Outer bars: number of host plant species on which the ASV occurs. Dots: branches with ultrafast bootstrap support of 95% or more. The 20 most abundant genera are named. **B**, Maximum-likelihood phylogeny of the 1,113 host specialist ASVs (a single host species). Rings, from the centre outwards: phylum, class, and one ring per plant family (Apiaceae, Asteraceae, Brassicaceae, Fabaceae, Poaceae and Solanaceae), dark where it is the family of the host. **C**, The 25 genera found on the most host species: number of host species and of host genera (bars) and mean relative abundance across all samples (diamonds, right axis). **D**, Predicted traits most strongly associated with host breadth: odds ratio per doubling of the number of host species, with its 95% confidence interval (binomial generalised linear models, q < 0.05; ten traits in each direction), coloured by source. **E**, Share of the KEGG orthologues of each PGPg_finder function that is enriched in host specialists (left) or in host generalists (right) (FDR < 0.05, absolute log2 odds ratio of at least 0.5), with the number of orthologues beside each bar; CE, competitive exclusion. Traits and orthologues are predicted from the 16S rRNA gene.

## Where the former figures went

The manuscript draft read on 2026-10-12 cites Figures 1–11 and S1–S7. With five main and five supplementary figures:

| Former number | Content | Now |
| --- | --- | --- |
| Figure 1 | Dataset overview | Figure 1 (A–F) |
| Figure 2 | Dominant phyla and classes (bar chart) | Figure 2A (class and genus composition of each host genus at each stage, dendrogram) |
| Figure 3 | Richness by host species and stage, 500 reads | Not in the figure sets; `02_alpha_diversity/alpha_richness_by_species_depth500.tiff` |
| Figure 4 | PCoA coloured by host species, all samples | Figure 2B–E (by stage; coloured by host genus and plant family) |
| Figure 5 | Screen for stage-associated ASVs | Not in the figure sets; `06_taxon_stage_associations/stage_screen.tiff` |
| Figure 6A, 6B | Stage effect on genera within plant family and cultivar | Figure 4B, 4C (within plant family and host genus) |
| Figure 7 | Host phylogeny and alpha diversity | Figure 3A |
| Figure 8 | Generalists and specialists on one phylogeny | Figure 5A, 5B (one tree per class) |
| Figure 9 | Predicted functions against host breadth (heatmaps) | Figure 5D, 5E (traits and plant-growth-promotion functions); Figure S4A (KEGG categories) |
| Figure 10 | Bacteria with distinctive predicted traits | Figure S4B, C |
| Figure 11 | Predicted functions between stages | Figure S5 |
| Figure S1 | Genus composition by plant family and stage (Sankey) | Figure 2A; the Sankey is `04_composition/fig_alluvial_genus_by_plant_family.tiff` |
| Figure S2 | Alpha diversity at 500 and 1,000 reads | Not in the figure sets; `02_alpha_diversity/alpha_indices_by_stage_rarefied.tiff` |
| Figure S3 | Stage contrast at genus and at ASV level | Not in the figure sets; `03_beta_diversity/pcoa_css_stage_unit_sensitivity.tiff` |
| Figure S4 | Co-occurrence network | Figure S1 |
| Figure S5 | Stage and host classes on the genus phylogeny | Not in the figure sets; `09_stage_generalists_specialists_clam/tree_generalist_specialist_genus.tiff` |
| Figure S6 | Phylogeny of host generalist ASVs | Figure 5A |
| Figure S7 | Phylogeny of host specialist ASVs | Figure 5B |
| (new) | Variance in community structure explained by each variable | Figure S2A, B |
| (new) | Seed and seedling communities by study | Figure S2C, D |
| (new) | Seed and seedling communities at ASV level, by study and host | Figure S2E–L |
| (new) | Null model behind the host-breadth threshold | Figure S3A |
| (new) | Taxonomic differences between host generalists and specialists | Figure S3B, C |
| (new) | Host divergence against community dissimilarity | Figure 4A |
| (new) | Stage effect on alpha diversity within plant families and host genera | Figure 3B, 3C |
| (new) | The 25 most widespread genera | Figure 5C |

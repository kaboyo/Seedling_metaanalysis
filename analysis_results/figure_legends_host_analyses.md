# Figure legends: dataset overview and host analyses

Legends of the component figures; the letters are placeholders. Each states what is drawn and no more. The files are in the analysis folders under `analysis_results/by_analysis/`.

## Figure A. Dataset overview

*File:* `study_overview.tiff`

**Overview of the compiled seed and seedling 16S rRNA gene amplicon dataset.** The dataset comprises 1,886 samples from eight studies after filtering. **A**, Countries of origin, coloured as in B. **B**, Share of samples by country. **C**, Share of samples from seeds and from seedlings. **D**, Samples per study, by developmental stage (seeds / seedlings). **E, F**, Samples per host plant genus among seeds (E) and seedlings (F). Seeds are orange and seedlings teal throughout.

*Source:* `setup/ps_clean.rds`, `tables/metadata_overview_publication.csv`.

## Figure B. Host plant phylogeny and alpha diversity

*File:* `host_phylosignal.tiff`

**Host plant phylogeny and the alpha diversity of the seed and seedling microbiota of each host.** Time-calibrated phylogeny of the 16 host plant species (TimeTree 5), coloured by plant family, with the mean observed richness, Shannon index and Pielou's evenness of each host at each developmental stage after rarefaction to 500 reads. Point size gives the number of samples; hosts with fewer than three samples at a stage carry no point. Blomberg's *K* (999 permutations) and Pagel's λ are given above each panel.

*Source:* `tables/host_phylosignal.csv`, `tables/host_phylosignal_species_means.csv`, `tables/host_phylosignal_tree.nwk`.

## Figure C. Host generalist and host specialist ASVs on one phylogeny

*File:* `tree_generalist_specialist.tiff`

**Host generalist and host specialist ASVs on one phylogeny.** Maximum-likelihood phylogeny of the 241 host generalist ASVs (purple; at least six host species) and the 1,113 host specialist ASVs (dark slate; a single host species). **Rings, from the centre:** host class; phylum; oxygen use predicted by BugBase; one ring per plant family, filled where the ASV was detected. **Outer bars**, number of host species. **Table**, linear regression of log10 relative abundance on the number of host species and one-sided Spearman correlation, for the five phyla richest in ASVs.

*Source:* `tables/tree_generalist_specialist.csv`, `tables/generalist_specialist_breadth_abundance.csv`, `trees/combined/tree.nwk`, `trees/tree_build_summary.csv`.

## Figure D. Phylogeny of host generalist ASVs

*File:* `tree_host_generalists.tiff`

**Phylogeny of the host generalist ASVs.** Maximum-likelihood phylogeny of the 241 host generalist ASVs (at least six host species). Rings, from the centre outwards: phylum and class. Outer bars: number of host plant species on which the ASV occurs. Dots: branches with ultrafast bootstrap support of 95% or more. The 20 most abundant genera are named.

*Source:* `tables/generalist_specialist_asv_host_screen.csv`, `trees/generalists/tree.nwk`, `trees/tree_build_summary.csv`.

## Figure E. Phylogeny of host specialist ASVs

*File:* `tree_host_specialists.tiff`

**Phylogeny of the host specialist ASVs.** Maximum-likelihood phylogeny of the 1,113 host specialist ASVs (a single host species). Rings, from the centre outwards: phylum, class, and one ring per plant family (Apiaceae, Asteraceae, Brassicaceae, Fabaceae, Poaceae and Solanaceae), dark where it is the family of the host.

*Source:* `tables/generalist_specialist_asv_host_screen.csv`, `trees/specialists/tree.nwk`, `trees/tree_build_summary.csv`.

## Figure F. Bacteria with distinctive predicted traits across host genera

*File:* `special_group_heatmaps.tiff`

**Distribution across host plant genera of bacteria with distinctive predicted traits.** **A**, Predicted anaerobic bacteria (BugBase). **B**, Predicted nitrogen-fixing bacteria (FAPROTAX). Colour is the mean relative abundance of each taxon on each host plant genus, log10-transformed and scaled by row; rows and columns are clustered by Ward's method. A star marks taxa that differ among host plant orders (Kruskal–Wallis, q < 0.05); f, family; o, order.

*Source:* `tables/special_group_heatmap_matrix.csv`, `tables/special_group_host_order_tests.csv`, `tables/special_group_summary.csv`.

## Figure H. Predicted functions against host breadth

*File:* `functional_host_breadth_heatmaps.tiff`

**Predicted functions against host breadth.** **A**, The 12 binary predicted traits most significantly associated with host breadth in each direction (binomial generalised linear models, q < 0.05): log2 odds ratio of carrying the trait, by number of host species, relative to all annotated ASVs; grey, fewer than five ASVs; OR, odds ratio per doubling of host number; B, BugBase; F, FAPROTAX; M, metaTraits. **B, C**, Share of the KEGG orthologues of each KEGG level-2 category (B) and PGPg_finder function (C) enriched in host generalists or specialists (FDR < 0.05, absolute log2 odds ratio of at least 0.5), with the number of orthologues in each cell; CE, competitive exclusion.

*Source:* `tables/predicted_trait_breadth_tests.csv`, `tables/predicted_traits_asv.csv`, `tables/functional_trait_heatmap_matrix.csv`, `tables/ko_enrichment_summary.csv`.

## Figure I. Predicted functions between developmental stages

*File:* `functional_stage_heatmaps.tiff`

**Predicted functions between developmental stages.** Cells are, within each plant family, Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. **A**, Predicted traits. **B, C**, KEGG level-2 and PGPg_finder categories. **D**, The three special groups.

*Source:* `tables/predicted_trait_stage_effect_{by_plant_family,pooled}.csv`, `tables/ko_category_stage_effect_{by_plant_family,pooled}.csv`, `tables/special_group_stage_effect_{by_plant_family,pooled}.csv`.

## Figure J. Effect of developmental stage on alpha diversity within plant families and host genera

*File:* `alpha_stage_effect_heatmap.tiff`

**Effect of developmental stage on alpha diversity within plant families and host genera.** Cells are Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. Observed richness, Shannon diversity and Pielou's evenness at 500 reads, at 1,000 reads and unrarefied, within plant families (B) and within host genera (C).

*Source:* `tables/alpha_stage_effect_by_stratum.csv`, `tables/alpha_stage_effect_pooled.csv`.

## Figure S. Screen for stage-associated ASVs

*File:* `stage_screen.tiff`

**Screen for stage-associated ASVs.** Each point is an ASV tested by DESeq2, placed by its log2 fold change (seedlings over seeds) and adjusted p-value; dashed lines mark adjusted p = 0.05 and a twofold change. ASVs also called by the indicator analysis are in the stage colours (seeds, n = 86; seedlings, n = 85); the others are grey.

*Source:* `tables/stage_screen_asvs.csv`.

## Figure G. Stage effect within each plant family and host genus

*File:* `stage_effect_heatmap.tiff`

**Stage effect on genera within each plant family (B) and each host genus (C).** Cells are Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. Centred log-ratio abundance of the 10 genera that hold an ASV significant in both DESeq2 and the indicator analysis, within plant families (B) and within host genera (C); the strip above the columns gives the plant family.

*Source:* `tables/stage_effect_by_plant_family_genus.csv`, `tables/stage_effect_by_host_genus_genus.csv`, and the matching `stage_effect_pooled_by_*_genus.csv` and `stage_effect_leave_one_out_by_*_genus.csv`.

## Supplementary figure. Generalist and specialist genera on the genus phylogeny

*File:* `tree_generalist_specialist_genus.tiff`

**Stage and host classes of bacterial genera on the genus phylogeny.** Dendrogram of 249 genera, pruned from the phylogeny of all ASVs; tip points are coloured by phylum and the 40 most abundant genera are named. **Inner bars**, stage class (CLAM): seed specialist (orange), seedling specialist (teal) or stage generalist (navy), full-length where the class recurs in at least two studies. **Second ring**, host class: generalist (purple), specialist (black) or intermediate (grey). **Outer heatmap**, mean relative abundance (%) on each host plant genus, on a logarithmic scale; grey, not detected.

*Source:* `tables/generalist_specialist_genus.csv`, `tables/tree_generalist_specialist_genus.csv`.

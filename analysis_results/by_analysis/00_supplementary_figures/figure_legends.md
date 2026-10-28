# Supplementary figures

Numbered as in the manuscript (Manuscript_revised.docx, 22 October 2026), in order of first citation. Each is
183 mm wide and 600 dpi. Files: `Figure_S1.tiff` to `Figure_S5.tiff` in this folder; the panels are set in the
`supplementary_figures` section of `analysis_config.json` and rebuilt with `Rscript make_main_figures.R`.

| Figure | Shows | Panels | Supplementary table |
| --- | --- | --- | --- |
| S1 | Co-occurrence networks of seeds and seedlings | cooccurrence_network | – |
| S2 | What structures the communities: explained variance, and study and host within each stage | beta_variance_explained + pcoa_css_by_stage_study + pcoa_css_asv_by_stage | – (values on the panels) |
| S3 | Host breadth: the null model and the taxa of host generalists and specialists | host_breadth_null + host_class_composition | S4, S5 |
| S4 | Predicted functions and host: orthologue categories and bacteria with distinctive traits | ko_kegg_category_heatmap + special_group_heatmaps | S7 |
| S5 | Predicted functions between developmental stages | functional_stage_heatmaps | S8 |

## Legends

**Figure S1. Co-occurrence networks of seeds and seedlings.** Networks of the 50 most abundant ASVs of each developmental stage, built separately for seeds (A; 1,282 samples) and seedlings (B; 602 samples) and drawn by connected module. Nodes are ASVs, coloured by genus and sized by degree; edges are Spearman correlations with |ρ| ≥ 0.6 and adjusted p < 0.05 within the stage. Each module is named after its genera.

**Figure S2. What structures the communities: explained variance, and study and host within each stage.** A, PERMANOVA R² of each variable alone (developmental stage, plant family, and host genus), at genus and at ASV level (999 permutations). B, Variance partitioning between host (as host genus or plant family) and developmental stage (adjusted R²), at genus and at ASV level. C, D, Principal coordinates analysis of Bray–Curtis dissimilarities (cumulative sum scaling, genus level), computed separately for seeds (C) and seedlings (D), as in Figure 2B–E, and coloured by study. E–J, The same at ASV level, for seeds (E–G) and seedlings (H–J), coloured by study (E, H), host genus (F, I) and plant family (G, J). R² and p of the PERMANOVA of the variable within the stage are given on each panel (999 permutations). Each stage has its own axes.

**Figure S3. Host breadth: the null model and the taxa of host generalists and specialists.** A, ASVs and genera found on each number of host plant species (bars) against the null expectation (points, mean; error bars, 95% range of 1,000 iterations). The dashed line marks the generalist threshold. B, Taxa whose share of ASVs differs between the two classes (Fisher's exact test, q < 0.05). C, Spearman correlation (ρ) between the number of host species and relative abundance, within the five phyla richest in ASVs, with the number of ASVs tested.

**Figure S4. Predicted functions and host: orthologue categories and bacteria with distinctive traits.** A, Share of the KEGG orthologues of each KEGG level-2 category that is enriched in host generalists or in host specialists (FDR < 0.05, absolute log2 odds ratio of at least 0.5), with the number of orthologues in each cell; the 15 categories with the most orthologues and, in the first row, all 7,883 orthologues. B, Predicted anaerobic bacteria (BugBase). C, Predicted nitrogen-fixing bacteria (FAPROTAX). Colour is the mean relative abundance of each taxon on each host plant genus, log10-transformed and scaled by row; rows and columns are clustered by Ward's method. A star marks taxa that differ among host plant orders (Kruskal–Wallis, q < 0.05); f, family; o, order. Orthologues, traits and special groups are predicted from the 16S rRNA gene.

**Figure S5. Predicted functions between developmental stages.** Cells are, within each plant family, Hedges' g, seedlings minus seeds (orange, higher in seeds; teal, higher in seedlings); the framed column is the random-effects pooled estimate, a dot marks a 95% confidence interval that excludes zero, and a cross a cell that could not be estimated. A, Predicted traits. B, C, KEGG level-2 and PGPg_finder categories. D, The three special groups. Orthologues, traits and special groups are predicted from the 16S rRNA gene.

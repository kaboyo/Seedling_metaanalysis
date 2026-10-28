# Manuscript Layout Plan

## Main Text

1. Dataset scope and metadata audit
   - Use the manuscript-ready overview table rather than a figure.
   - Report sample counts, studies, host species, cultivar coverage, geography, and developmental-stage coverage.

2. Community composition
   - Lead with `top_phyla_composition.tiff`.
   - Keep the phylum Sankey and rank-specific Sankey panels in the supplement unless they are needed to support a specific composition claim.

3. Alpha diversity
   - Use `alpha_richness_by_species_depth500.tiff` as the primary alpha-diversity figure.
   - Put `alpha_richness_by_species_depth1000.tiff` in the supplement as the rarefaction-depth sensitivity check; the stage contrast in richness depends on depth.
   - Keep developmental-stage faceting in both plots.

4. Beta diversity
   - Use the genus-level Bray-Curtis PCoA as the main ordination result.
   - Keep developmental-stage faceting for the assessment plots.
   - Report PERMANOVA and PERMDISP together so separation and dispersion are interpreted jointly.

5. Taxon-level associations
   - Keep the top-5 genus correlations stratified by family and cultivar.
   - Present the tables and compact plots together.

6. Core and interaction structure
   - Place core microbiome, indicator taxa, and co-occurrence network in the later Results section or Supplement, depending on length.

## Supplementary Material

- Sequencing depth and rarefaction diagnostics
- Full Sankey outputs across taxonomic ranks
- Correlation tables
- Network edge lists
- Any QC or removed-heatmap outputs

## Narrative Order

1. Dataset scope and filtering
2. Composition across developmental stage
3. Alpha diversity by species and cultivar
4. Beta diversity and genus-level clustering
5. Top genus associations with stage
6. Core microbiome and network structure

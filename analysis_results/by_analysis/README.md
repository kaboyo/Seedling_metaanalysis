# Results by analysis

One folder per analysis, with its tables and figures together.

- **Figures are kept here and nowhere else.** `figures/` is only a staging folder that
  the pipeline empties into these folders at the end of each run.
- **Tables are copies.** The pipeline reads its tables from `tables/`, and refreshes the
  copies here at the end of each run.

`index.csv` lists every file here and where it came from.

| Folder | Methodology paragraph | Analysis | Figures | Tables | Other |
| --- | --- | --- | --- | --- | --- |
| `01_data_and_filtering` | 1 | Data and filtering | 7 | 13 | 0 |
| `02_alpha_diversity` | 2 | Alpha diversity | 6 | 16 | 0 |
| `03_beta_diversity` | 2 | Beta diversity (CSS, Bray-Curtis, PERMANOVA) | 10 | 5 | 0 |
| `04_composition` | 2 | Taxonomic composition | 16 | 17 | 0 |
| `05_cooccurrence_network` | 2 | Co-occurrence network | 1 | 1 | 0 |
| `06_taxon_stage_associations` | 3 | Taxon-stage associations (two steps) | 4 | 10 | 0 |
| `07_host_phylogenetic_signal` | 4 | Host phylogenetic signal in alpha diversity | 2 | 2 | 2 |
| `08_phylosymbiosis` | 4 | Phylosymbiosis, marginal PERMANOVA and variance partitioning | 1 | 4 | 1 |
| `09_stage_generalists_specialists_clam` | 5 | Stage generalists and specialists (CLAM), genera | 1 | 2 | 0 |
| `10_host_generalists_specialists` | 6 | Host generalists and specialists (two-step screen) | 1 | 4 | 1 |
| `11_trees_and_class_comparison` | 7 | Trees of host generalists and specialists, and their comparison | 5 | 6 | 22 |
| `12_predicted_functions` | 8 | Predicted traits, gene content and special groups | 9 | 19 | 3 |

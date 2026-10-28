#!/usr/bin/env Rscript
# Redraw the taxonomic composition figures and tables without rendering the whole master
# document.
#
#   Rscript make_composition.R [analysis_config.json]
#
# This holds no analysis or plotting code of its own. It runs the Part 1 chunks `top-phyla`
# and `taxa-strips-by-study` of seedling_analysis_master.qmd on the filtered dataset the
# pipeline already saved, then the Part 2 chunks `fx-fig-top-phyla` and `fx-fig-alluvial-stage`
# (the Sankey diagrams, which use Part 1's `family_map` and palette) and the folder step, so
# what is written here is what a full render writes. The rank of the composition bars is
# the config's taxonomy_analysis.composition_rank (class). Edit the figures in the master.
#
# Reads  <outputs.directory>/setup/ps_clean.rds.
# Writes <outputs.directory>/tables/top_<rank>_summary.csv and taxa_strips_*.csv, and the
#        figures top_<rank>_composition.tiff and fig_taxa_strips_*.tiff, which the folder
#        step places in <outputs.directory>/by_analysis/04_composition/.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

MASTER <- "seedling_analysis_master.qmd"
lines <- readLines(MASTER, encoding = "UTF-8")
chunk_code <- function(label) {
    start <- grep(sprintf("^```\\{r %s[,}]", label), lines)
    if (length(start) != 1) stop("chunk not found exactly once in ", MASTER, ": ", label)
    end <- start + which(grepl("^```\\s*$", lines[(start + 1):length(lines)]))[1]
    lines[(start + 1):(end - 1)]
}
run_chunk <- function(label) eval(parse(text = chunk_code(label), keep.source = FALSE),
                                  envir = globalenv())

# Part 1's own setup: packages, palettes, captions, config, folders, themes.
for (label in c("load-packages", "color-palettes", "figure-captions", "figure-themes")) run_chunk(label)

# The objects Part 1 holds at this point, taken from the filtered dataset it saved.
ps_file <- file.path(setup_dir, "ps_clean.rds")
if (!file.exists(ps_file)) stop("no filtered dataset at ", ps_file, " - run the pipeline first")
ps <- readRDS(ps_file)
counts <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts <- t(counts)
taxa_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)[rownames(counts), ]
metadata_df <- as(sample_data(ps), "data.frame")
# Part 1 lists further metadata columns to carry along (chunk `metadata-summary`); the
# composition bars and tables do not use them.
plot_group_vars <- character()

for (label in c("top-phyla", "taxa-strips-by-study")) run_chunk(label)
for (label in c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core", "fx-fig-top-phyla",
                "fx-fig-alluvial-stage", "fx-fig-composition-dendrogram", "fx-analysis-folders")) run_chunk(label)

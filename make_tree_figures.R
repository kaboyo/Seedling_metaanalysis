#!/usr/bin/env Rscript
# Draw the eight phylogenetic tree figures without rendering the whole master document.
#
#   Rscript make_tree_figures.R [analysis_config.json]
#
# This holds no plotting code of its own. It runs the tree chunk of
# seedling_analysis_master.qmd (Part 2, `fx-fig-tree-genus`) together with the four Part 2
# setup chunks it depends on, so the figures written here are exactly the ones a full
# render writes. Edit the figures in the master, not here.
#
# Reads  <outputs.directory>/setup/ps_clean.rds and the tree named in the config; the
#        generalist/specialist tree also needs tables/generalist_specialist_genus.csv
#        (make_generalist_specialist.R) and is skipped without it.
# Writes <outputs.directory>/figures/tree_{genus,asv,asv_prevalent}_by_stage.tiff and the
#        matching tables in <outputs.directory>/tables/,
#        figures/tree_generalist_specialist.tiff (both ASV classes on one tree, the main
#        figure), tree_generalist_specialist_genus.tiff, tree_host_generalists.tiff and
#        tree_host_specialists.tiff (these four need make_generalist_specialist.R's tables;
#        the ASV trees also need generalist_specialist_trees.sh's trees, and the main figure
#        predicted_traits.py's table), and figures/host_phylosignal.tiff (the host plant
#        tree; needs host_phylosignal.R's tables). A figure whose inputs are missing is skipped.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

MASTER <- "seedling_analysis_master.qmd"
CHUNKS <- c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core", "fx-fig-tree-genus")

lines <- readLines(MASTER, encoding = "UTF-8")
chunk_code <- function(label) {
    start <- grep(sprintf("^```\\{r %s[,}]", label), lines)
    if (length(start) != 1) stop("chunk not found exactly once in ", MASTER, ": ", label)
    end <- start + which(grepl("^```\\s*$", lines[(start + 1):length(lines)]))[1]
    lines[(start + 1):(end - 1)]
}
for (label in CHUNKS) {
    eval(parse(text = chunk_code(label), keep.source = FALSE), envir = globalenv())
}

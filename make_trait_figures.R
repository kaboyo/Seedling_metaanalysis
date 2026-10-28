#!/usr/bin/env Rscript
# Draw the predicted-trait figures without rendering the whole master document.
#
#   Rscript make_trait_figures.R [analysis_config.json]
#
# This holds no plotting code of its own. It runs the chunks `fx-fig-special-groups`,
# `fx-fig-functional` and `fx-fig-functional-stage` of seedling_analysis_master.qmd (Part 2)
# with the Part 2 setup chunks
# they depend on, so the figures written here are exactly the ones a full render writes.
# Edit the figures in the master, not here.
#
# Reads  <outputs.directory>/tables/special_group_heatmap_matrix.csv and
#        special_group_host_order_tests.csv (written by special_groups.R).
#        For the functional heatmaps: tables/predicted_trait_breadth_tests.csv,
#        predicted_traits_asv.csv, generalist_specialist_asv_host_screen.csv and
#        ko_enrichment_summary.csv. For the stage heatmaps: the six
#        tables/*_stage_effect_*.csv files that functional_stage_effects.R writes.
# Writes special_group_heatmaps.tiff, functional_host_breadth_heatmaps.tiff and
#        functional_stage_heatmaps.tiff into their analysis folder under
#        <outputs.directory>/by_analysis/, and
#        tables/functional_trait_heatmap_matrix.csv.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

MASTER <- "seedling_analysis_master.qmd"
CHUNKS <- c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core", "fx-fig-special-groups",
            "fx-fig-functional", "fx-fig-functional-stage", "fx-fig-host-breadth-panels",
            "fx-fig-supplementary-panels")

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

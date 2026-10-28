#!/usr/bin/env Rscript
# Estimate the stage effect within each plant family and cultivar and draw its two figures
# (the screen and the heatmap), without rendering
# the whole master document.
#
#   Rscript make_stage_effects.R [analysis_config.json]
#
# This holds no analysis or plotting code of its own. It runs the Part 1 chunk
# `stage-effects-by-stratum` of seedling_analysis_master.qmd on the filtered dataset the
# pipeline already saved, then the Part 2 chunks `fx-fig-stage-screen` and `fx-fig-stage-heatmap`, so the tables and
# the figure written here are exactly the ones a full render writes. Edit them in the
# master, not here.
#
# Reads  <outputs.directory>/setup/ps_clean.rds,
#        <outputs.directory>/differential_abundance/deseq2_results.csv and
#        <outputs.directory>/tables/indicator_species_significant.csv.
# Writes <outputs.directory>/tables/stage_effect_{by,pooled_by}_{plant_family,cultivar}_genus.csv and
#        <outputs.directory>/tables/stage_screen_asvs.csv, and
#        <outputs.directory>/figures/stage_screen.tiff and stage_effect_heatmap.tiff.

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
for (label in c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core")) run_chunk(label)
if (is.null(ps)) stop("no filtered dataset at ", PS_CLEAN, " - run the pipeline first")
# The Part 1 names the analysis chunk uses, and Part 1's own config reader lifted from its
# setup chunk rather than copied.
group_col <- PRIMARY
tables_dir <- res("tables")
da_dir <- res("differential_abundance")
for (e in parse(text = chunk_code("figure-themes"), keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("get_config_value"))) eval(e, envir = globalenv())
}
# family_map (host species -> plant family) is defined in a Part 1 composition chunk
family_chunk <- sub("^```\\{r ([^,}]+).*", "\\1",
                    lines[max(grep("^```\\{r", lines)[grep("^```\\{r", lines) <
                                                          grep("^family_map <- ", lines)[1]])])
for (e in parse(text = chunk_code(family_chunk), keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("family_map"))) eval(e, envir = globalenv())
}
run_chunk("stage-effects-by-stratum")
if (!identical(stage_effects_result$status, "completed")) quit(status = 1)
for (label in c("fx-fig-stage-screen", "fx-fig-stage-heatmap")) run_chunk(label)

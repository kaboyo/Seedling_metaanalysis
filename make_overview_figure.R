#!/usr/bin/env Rscript
# Draw the dataset overview figure (panels A-F) without rendering the whole master document.
#
#   Rscript make_overview_figure.R [analysis_config.json]
#
# This holds no plotting code of its own. It runs the overview chunk of
# seedling_analysis_master.qmd (Part 2, `fx-fig-study-overview`) together with the four
# Part 2 setup chunks it depends on, so the figure written here is exactly the one a full
# render writes. Edit the figure in the master, not here.
#
# Reads  <outputs.directory>/setup/ps_clean.rds and
#        <outputs.directory>/tables/metadata_overview_publication.csv.
# Writes <outputs.directory>/figures/study_overview.tiff.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

MASTER <- "seedling_analysis_master.qmd"
CHUNKS <- c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core", "fx-fig-study-overview")

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

#!/usr/bin/env Rscript
# Merge the finished figures into the four figures of the manuscript, without rendering
# the whole master document.
#
#   Rscript make_main_figures.R [analysis_config.json]
#
# This holds no code of its own. It runs the Part 2 chunk `fx-fig-main-figures` of
# seedling_analysis_master.qmd after the Part 2 setup chunk. Which finished figure goes
# into which main figure is the `main_figures` section of the config.
#
# Reads  the figures in <outputs.directory>/by_analysis/<folder>/.
# Writes <outputs.directory>/by_analysis/00_main_figures/Figure_<n>.tiff and
#        main_figures.csv.

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
for (label in c("fx-setup", "fx-palettes", "fx-fig-main-figures")) {
    eval(parse(text = chunk_code(label), keep.source = FALSE), envir = globalenv())
}

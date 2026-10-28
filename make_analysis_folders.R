#!/usr/bin/env Rscript
# Deposit each analysis's tables and figures together, one folder per analysis, without
# rendering the whole master document.
#
#   Rscript make_analysis_folders.R [analysis_config.json]
#
# This holds no code of its own. It runs the `fx-analysis-folders` chunk of
# seedling_analysis_master.qmd (Part 2) after the Part 2 setup chunk, so the folders
# written here are exactly the ones a full render writes. Which file goes to which folder
# is the `analysis_folders` section of the config.
#
# Reads  <outputs.directory>/tables/, figures/ and the other files the config lists.
# Writes <outputs.directory>/<analysis_folders.directory>/<folder>/, with index.csv and
#        README.md beside the folders.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

MASTER <- "seedling_analysis_master.qmd"
CHUNKS <- c("fx-setup", "fx-analysis-folders")

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

#!/usr/bin/env Rscript
# Classify genera as generalists or specialists without rendering the whole master document.
#
#   Rscript make_generalist_specialist.R [analysis_config.json]
#
# This holds no analysis code of its own. It runs the Part 1 chunk `generalist-specialist`
# of seedling_analysis_master.qmd on the filtered dataset the pipeline already saved, so
# the tables written here are exactly the ones a full render writes. Edit the analysis in
# the master, not here. Draw the dendrogram afterwards with make_tree_figures.R.
#
# Reads  <outputs.directory>/setup/ps_clean.rds.
# Writes <outputs.directory>/tables/generalist_specialist_{genus,summary}.csv.

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
# Part 2's setup loads the config and ps_clean.rds, which is the `ps` Part 1 holds here.
for (label in c("fx-setup", "fx-load-core")) {
    eval(parse(text = chunk_code(label), keep.source = FALSE), envir = globalenv())
}
if (is.null(ps)) stop("no filtered dataset at ", PS_CLEAN, " - run the pipeline first")
# The two Part 1 names the chunk uses, and Part 1's own config reader lifted from its
# setup chunk rather than copied.
group_col <- PRIMARY
tables_dir <- res("tables")
for (e in parse(text = chunk_code("figure-themes"), keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("get_config_value"))) eval(e, envir = globalenv())
}
eval(parse(text = chunk_code("generalist-specialist"), keep.source = FALSE), envir = globalenv())
if (!identical(gs_result$status, "completed")) quit(status = 1)

#!/usr/bin/env Rscript
# Rebuild the alpha-diversity test tables without rendering the whole master document and
# without rarefying again.
#
#   Rscript make_alpha_tests.R [analysis_config.json]
#
# This holds no analysis code of its own. It runs the function `rarefied_alpha_tests()` of
# the Part 1 chunk `rarefied-alpha` and the whole Part 1 chunks `alpha-findings` and
# `alpha-stage-by-stratum` of
# seedling_analysis_master.qmd on the per-sample values the pipeline already saved, so the
# tables written here are the ones a full render writes. Rarefaction draws reads at
# random, so the per-sample tables are read back, never recomputed here.
#
# Reads  <outputs.directory>/setup/ps_clean.rds and tables/alpha_diversity.csv,
#        rarefied_alpha_diversity_depth<d>.csv.
# Writes <outputs.directory>/tables/rarefied_alpha_tests_depth<d>.csv,
#        alpha_diversity_findings.csv, alpha_stage_effect_by_stratum.csv and
#        alpha_stage_effect_pooled.csv, and the figure alpha_stage_effect_heatmap.tiff in
#        its analysis folder under <outputs.directory>/by_analysis/.

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
# One top-level assignment lifted from a chunk, rather than copied here.
lift <- function(label, name) {
    for (e in parse(text = chunk_code(label), keep.source = FALSE)) {
        if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(name))) {
            return(invisible(eval(e, envir = globalenv())))
        }
    }
    stop(name, " is not assigned at the top level of chunk ", label)
}
for (label in c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core")) run_chunk(label)
if (is.null(ps)) stop("no filtered dataset at ", PS_CLEAN, " - run the pipeline first")

# The Part 1 names the two chunks use.
group_col <- PRIMARY
tables_dir <- res("tables")
rare_depths <- sort(unique(as.integer(config$filters$rarefaction_depths)))
lift("top-phyla", "family_map")
lift("rarefied-alpha", "rarefied_alpha_tests")

for (depth in rare_depths) {
    file <- file.path(tables_dir, sprintf("rarefied_alpha_diversity_depth%d.csv", depth))
    if (!file.exists(file)) stop("run the pipeline first; missing ", file)
    tests <- rarefied_alpha_tests(read.csv(file, check.names = FALSE), depth)
    readr::write_csv(tests, file.path(tables_dir, sprintf("rarefied_alpha_tests_depth%d.csv", depth)))
    cat(sprintf("depth %d: %d Kruskal-Wallis tests\n", depth, nrow(tests)))
}

# The unrarefied indices, in the sample order of the filtered dataset. The saved table has
# no sample column, so the order is checked against the data before it is trusted.
unrarefied <- read.csv(file.path(tables_dir, "alpha_diversity.csv"), check.names = FALSE)
otu <- as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
if (nrow(unrarefied) != ncol(otu) || !all(unrarefied$Observed == colSums(otu > 0))) {
    stop("tables/alpha_diversity.csv does not follow the sample order of ps_clean.rds")
}
alpha_result <- list(alpha_df = unrarefied)
run_chunk("alpha-findings")
if (!identical(alpha_findings_result$status, "completed")) quit(status = 1)
# the stage effect within each plant family and cultivar, and across them
lift("figure-themes", "get_config_value")
run_chunk("alpha-stage-by-stratum")
if (!identical(alpha_stage_strata_result$status, "completed")) quit(status = 1)
# its figure, drawn with the heatmap panel of the functional stage figure
for (name in c("FUNCTIONAL_FAMILY_ORDER", "FUNCTIONAL_STAGE_LIMIT", "FUNCTIONAL_STAGE_MARKS",
               "stage_function_panel")) lift("fx-fig-functional-stage", name)
run_chunk("fx-fig-alpha-stage-strata")

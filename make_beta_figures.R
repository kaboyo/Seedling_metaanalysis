#!/usr/bin/env Rscript
# Ordinate seeds and seedlings apart and draw the figure of host structure within each
# stage, without rendering the whole master document.
#
#   Rscript make_beta_figures.R [analysis_config.json]
#
# This holds no analysis or plotting code of its own. It runs the Part 1 chunk
# `pcoa-by-stage` of seedling_analysis_master.qmd on the distances the pipeline cached,
# then the Part 2 chunk `fx-fig-pcoa-by-stage`, so the table and the figure written here
# are the ones a full render writes. Edit them in the master, not here.
#
# Reads  <outputs.directory>/setup/ps_clean.rds and cache/bray_css_genus.rds.
# Writes <outputs.directory>/diversity/pcoa_css_genus_by_stage.csv and the figure
#        pcoa_css_by_stage.tiff in its analysis folder under by_analysis/.

args <- commandArgs(trailingOnly = TRUE)
params <- list(config = if (length(args)) args[1] else "analysis_config.json")

# Every name of this launcher starts with .mk_, so that the master's code cannot overwrite it.
.mk_lines <- readLines("seedling_analysis_master.qmd", encoding = "UTF-8")
.mk_code <- function(label) {
    start <- grep(sprintf("^```\\{r %s[,}]", label), .mk_lines)
    if (length(start) != 1) stop("chunk not found exactly once: ", label)
    end <- start + which(grepl("^```\\s*$", .mk_lines[(start + 1):length(.mk_lines)]))[1]
    .mk_lines[(start + 1):(end - 1)]
}
.mk_run <- function(label) base::eval(base::parse(text = .mk_code(label), keep.source = FALSE),
                                      envir = globalenv())
# Part 1's own setup: packages, palettes, captions, config, folders, themes.
for (.mk_label in c("load-packages", "color-palettes", "figure-captions", "figure-themes")) .mk_run(.mk_label)
.mk_ps <- file.path(setup_dir, "ps_clean.rds")
if (!file.exists(.mk_ps)) stop("no filtered dataset at ", .mk_ps, " - run the pipeline first")
ps <- readRDS(.mk_ps)
# family_map (host species -> plant family) is defined in the composition chunk
for (.mk_e in base::parse(text = .mk_code("top-phyla"), keep.source = FALSE)) {
    if (is.call(.mk_e) && identical(.mk_e[[1]], as.name("<-")) && identical(.mk_e[[2]], as.name("family_map"))) {
        base::eval(.mk_e, envir = globalenv())
    }
}
.mk_run("pcoa-by-stage")
if (!identical(pcoa_by_stage_result$status, "completed")) quit(status = 1)
for (.mk_label in c("fx-setup", "fx-palettes", "fx-fig-registry", "fx-load-core", "fx-fig-pcoa-by-stage")) {
    .mk_run(.mk_label)
}

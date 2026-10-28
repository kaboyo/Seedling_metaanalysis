#!/usr/bin/env Rscript
# functional_stage_effects.R - do predicted functions differ between seeds and seedlings?
#
# Usage:
#   Rscript functional_stage_effects.R --config analysis_config.json
#
# Run predicted_trait_tests.R, ko_breadth_tests.R and special_groups.R first.
#
# The trait, gene-content and special-group analyses relate predicted functions to host
# breadth. This adds the study's primary contrast, developmental stage, to each of them,
# with the design of the taxon-stage analysis: the difference between seedlings and seeds
# is estimated as Hedges' g inside each plant family holding both stages, and the family
# estimates are pooled by random-effects meta-analysis (REML, Knapp-Hartung), BH-FDR
# within each scheme. Three per-sample measures:
#   traits          the share of a sample's annotated relative abundance held by ASVs
#                   carrying the trait (binary traits tested against host breadth)
#   KO categories   the share of a sample's predicted orthologue copies that falls in the
#                   category: ASV counts divided by the predicted 16S copy number, times
#                   the predicted copies of the category's orthologues (those tested by
#                   ko_breadth_tests.R, ASVs with NSTI <= max_nsti); KEGG level 2 and
#                   PGPg_finder level 3
#   special groups  the relative abundance of the group in the sample
# Host is largely tied to study here, and this design does not remove study: the
# per-family table records how many studies supplied each stage.
#
# Reads  <outputs.directory>/setup/ps_clean.rds, picrust2/, and tables/
#        predicted_traits_asv.csv, predicted_trait_breadth_tests.csv, ko_breadth_tests.csv,
#        special_group_sample_abundance.csv, host_plant_species.csv.
# Writes <outputs.directory>/tables/{predicted_trait,ko_category,special_group}_stage_effect_
#        {by_plant_family,pooled}.csv.

suppressPackageStartupMessages({
    library(jsonlite)
    library(phyloseq)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) default else args[i + 1]
}
config_path <- arg_value("--config", "analysis_config.json")
if (!file.exists(config_path)) stop("Configuration file not found: ", config_path)
config <- fromJSON(config_path)
base_dir <- dirname(normalizePath(config_path))
resolve_path <- function(path) if (startsWith(path, "/")) path else file.path(base_dir, path)
setting <- function(section, key, default) {
    value <- config[[section]][[key]]
    if (is.null(value)) default else value
}
STAGE_COL <- config$analysis$primary_group
STAGES <- config$differential_abundance$groups            # seeds first, seedlings second
HOST_COL <- setting("generalist_specialist", "host_group", "Scientific_Name")
STUDY_COL <- setting("generalist_specialist", "study_group", "study_alias")
MIN_SAMPLES <- setting("stage_effects", "min_samples_per_stage", 3)
MIN_STRATA <- setting("stage_effects", "min_strata", 2)
MAX_NSTI <- setting("gene_content", "max_nsti", 0.15)
FDR <- setting("predicted_traits", "fdr", 0.05)

out_dir <- resolve_path(config$outputs$directory)
tables_dir <- file.path(out_dir, "tables")

ps <- readRDS(file.path(out_dir, "setup", "ps_clean.rds"))
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
rel <- sweep(otu, 2, pmax(colSums(otu), 1), "/")
meta <- as(sample_data(ps), "data.frame")
hosts <- read_csv(file.path(tables_dir, "host_plant_species.csv"), show_col_types = FALSE)
family <- hosts$Family[match(meta[[HOST_COL]], hosts$Species)]
stage <- factor(meta[[STAGE_COL]], levels = STAGES)
study <- as.character(meta[[STUDY_COL]])

# ---- the three measures, each a feature x sample matrix -------------------------------
tested <- read_csv(file.path(tables_dir, "predicted_trait_breadth_tests.csv"),
                   show_col_types = FALSE) |>
    dplyr::filter(status == "tested", type == "binary")
traits <- read_csv(file.path(tables_dir, "predicted_traits_asv.csv"), show_col_types = FALSE,
                   na = "NA", guess_max = 5000)
carried <- as.matrix(traits[match(rownames(otu), traits$ASV), tested$trait])
storage.mode(carried) <- "logical"
annotated <- !is.na(carried)
carried[!annotated] <- FALSE
# sample x trait: abundance of carriers over abundance of annotated ASVs
trait_share <- t(crossprod(rel, carried * 1) / crossprod(rel, annotated * 1))

picrust <- file.path(out_dir, "picrust2")
ko_files <- list.files(picrust, pattern = "^(bac|arc)_KO_predicted\\.tsv\\.gz$", full.names = TRUE)
marker_file <- file.path(picrust, "combined_marker_predicted_and_nsti.tsv.gz")
if (!length(ko_files) || !file.exists(marker_file)) stop("Run ko_predict.sh first; nothing in ", picrust)
marker <- as_tibble(fread(marker_file)) |>
    dplyr::select(sequence = 1, copies_16s = dplyr::matches("16S"), nsti = dplyr::matches("NSTI"))
ko <- purrr::map(ko_files, function(f) fread(f)) |> rbindlist(fill = TRUE)
setnames(ko, 1, "sequence")
setnames(ko, sub("^ko:", "", names(ko)))
annotation <- read_csv(file.path(tables_dir, "ko_breadth_tests.csv"), show_col_types = FALSE) |>
    dplyr::select(KO, kegg_level2, pgpg_lv3)
sequence <- sub("^ASV [0-9]+ : ", "", rownames(otu))
usable <- sequence %in% ko$sequence & sequence %in% marker$sequence[marker$nsti <= MAX_NSTI]
copies <- as.matrix(ko[match(sequence[usable], ko$sequence), annotation$KO, with = FALSE])
copies[is.na(copies)] <- 0
# ASV counts per 16S copy, so that an ASV with many rRNA operons is not over-counted
per_genome <- otu[usable, , drop = FALSE] /
    pmax(marker$copies_16s[match(sequence[usable], marker$sequence)], 1)
all_copies <- as.vector(crossprod(per_genome, rowSums(copies)))          # per sample
category_share <- function(column, scheme) {
    member <- annotation |> dplyr::select(KO, category = dplyr::all_of(column)) |>
        dplyr::filter(!is.na(category)) |>
        tidyr::separate_longer_delim(category, "; ") |>
        dplyr::distinct()
    categories <- sort(unique(member$category))
    in_category <- sapply(categories, function(k) annotation$KO %in% member$KO[member$category == k])
    share <- t(crossprod(per_genome, copies %*% in_category) / all_copies)
    rownames(share) <- paste(scheme, categories, sep = "\t")
    share
}
ko_share <- rbind(category_share("kegg_level2", "KEGG level 2"),
                  category_share("pgpg_lv3", "PGPG level 3"))

group_share <- read_csv(file.path(tables_dir, "special_group_sample_abundance.csv"),
                        show_col_types = FALSE) |>
    dplyr::select(sample, group, relative_abundance) |>
    tidyr::pivot_wider(names_from = sample, values_from = relative_abundance) |>
    tibble::column_to_rownames("group") |> as.matrix()
group_share <- group_share[, colnames(otu), drop = FALSE]

# ---- Hedges' g within plant family, pooled across families ----------------------------
families <- names(which(tapply(stage, family, function(x) all(table(x) >= MIN_SAMPLES))))
cat(sprintf("%d plant families with at least %d samples of each stage: %s\n", length(families),
            MIN_SAMPLES, paste(families, collapse = ", ")))

stage_effects <- function(values) {
    by_family <- purrr::map(families, function(f) {
        inside <- !is.na(family) & family == f
        seeds <- values[, inside & stage == STAGES[1], drop = FALSE]
        seedlings <- values[, inside & stage == STAGES[2], drop = FALSE]
        n_seeds <- rowSums(!is.na(seeds))
        n_seedlings <- rowSums(!is.na(seedlings))
        # a trait may have no annotated ASV in any sample of a family: nothing to estimate
        ok <- n_seeds >= MIN_SAMPLES & n_seedlings >= MIN_SAMPLES
        if (!any(ok)) return(NULL)
        es <- metafor::escalc(
            measure = "SMD",
            m1i = rowMeans(seedlings[ok, , drop = FALSE], na.rm = TRUE),
            sd1i = apply(seedlings[ok, , drop = FALSE], 1, sd, na.rm = TRUE), n1i = n_seedlings[ok],
            m2i = rowMeans(seeds[ok, , drop = FALSE], na.rm = TRUE),
            sd2i = apply(seeds[ok, , drop = FALSE], 1, sd, na.rm = TRUE), n2i = n_seeds[ok])
        tibble(feature = rownames(values)[ok], plant_family = f,
               n_seeds = n_seeds[ok], n_seedlings = n_seedlings[ok],
               studies_with_seeds = dplyr::n_distinct(study[inside & stage == STAGES[1]]),
               studies_with_seedlings = dplyr::n_distinct(study[inside & stage == STAGES[2]]),
               studies_with_both = length(intersect(study[inside & stage == STAGES[1]],
                                                    study[inside & stage == STAGES[2]])),
               mean_seeds = rowMeans(seeds[ok, , drop = FALSE], na.rm = TRUE),
               mean_seedlings = rowMeans(seedlings[ok, , drop = FALSE], na.rm = TRUE),
               hedges_g = as.numeric(es$yi), variance = as.numeric(es$vi))
    }) |> purrr::list_rbind() |>
        # a feature absent from, or constant in, a family has no standard deviation there
        dplyr::filter(is.finite(hedges_g), is.finite(variance), variance > 0) |>
        dplyr::mutate(ci_low = hedges_g - qnorm(0.975) * sqrt(variance),
                      ci_high = hedges_g + qnorm(0.975) * sqrt(variance))
    pooled <- by_family |>
        dplyr::filter(dplyr::n() >= MIN_STRATA, .by = feature) |>
        dplyr::group_split(feature) |>
        purrr::map(function(d) {
            fit <- tryCatch(metafor::rma(yi = d$hedges_g, vi = d$variance, method = "REML",
                                         test = "knha"), error = function(e) NULL)
            if (is.null(fit)) return(NULL)
            tibble(feature = d$feature[1], n_families = nrow(d),
                   n_families_same_direction = sum(sign(d$hedges_g) == sign(fit$b[1])),
                   pooled_g = fit$b[1], ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval,
                   i2 = fit$I2)
        }) |> purrr::list_rbind()
    list(by_family = by_family, pooled = pooled)
}

#' Write one analysis's two tables. `scheme_of` gives the family of tests a feature belongs
#' to; BH-FDR is applied within it.
write_effects <- function(values, prefix, scheme_of, name_of = identity) {
    fit <- stage_effects(values)
    label <- function(d) d |>
        dplyr::mutate(scheme = scheme_of(feature), feature = name_of(feature), .before = 1)
    pooled <- label(fit$pooled) |>
        dplyr::mutate(q = p.adjust(p, method = "BH"), .by = scheme) |>
        dplyr::arrange(scheme, dplyr::desc(pooled_g))
    write_csv(label(fit$by_family),
              file.path(tables_dir, sprintf("%s_stage_effect_by_plant_family.csv", prefix)))
    write_csv(pooled, file.path(tables_dir, sprintf("%s_stage_effect_pooled.csv", prefix)))
    cat(sprintf("%-16s %4d features pooled; q < %s: %d (higher in seedlings %d, in seeds %d)\n",
                prefix, nrow(pooled), FDR, sum(pooled$q < FDR),
                sum(pooled$q < FDR & pooled$pooled_g > 0), sum(pooled$q < FDR & pooled$pooled_g < 0)))
    invisible(pooled)
}

resource_of <- setNames(tested$resource, tested$trait)
write_effects(trait_share, "predicted_trait", function(x) unname(resource_of[x]))
write_effects(ko_share, "ko_category", function(x) sub("\t.*$", "", x),
              function(x) sub("^.*\t", "", x))
write_effects(group_share, "special_group", function(x) rep("special group", length(x)))

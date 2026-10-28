#!/usr/bin/env Rscript
# host_phylosymbiosis.R - beta diversity and host relatedness, on the CSS-normalised data.
#
# Usage:
#   Rscript host_phylosymbiosis.R --config analysis_config.json
#
# Run host_phylosignal.R first: it writes the one-tip-per-species host tree used here.
#
# Two questions, each asked at genus level and at ASV level, on the pipeline's own
# CSS-normalised Bray-Curtis distances (cache/bray_css_genus.rds and
# cache/bray_css_all_features.rds; nothing is re-normalised here).
#
# 1. Phylosymbiosis: are the microbiota of closely related hosts more alike? Sample-level
#    Bray-Curtis is averaged over all sample pairs between two host species, within a
#    developmental stage, and set against the hosts' divergence time from the host tree:
#    Mantel tests (Pearson and Spearman) and multiple regression on distance matrices
#    (MRM). A second MRM adds the share of those sample pairs that come from two different
#    studies, because most host species were sequenced by one study only and two studies
#    differ for reasons that have nothing to do with the plant.
# 2. How much do host and stage each explain? Marginal PERMANOVA
#    (vegan::adonis2, by = "margin") and variance partitioning (vegan::varpart) on the
#    sample-level distances. Host family, host species and cultivar are nested, so the
#    marginal effect of one given another does not exist; they are tested in separate models,
#    the cultivar models on the samples with a cultivar recorded.
#
# Genus level is the one to read across studies: an ASV is an exact sequence and most
# ASVs occur in a single study, so ASV-level distances between studies are near their
# maximum whatever the host. ASV level is reported beside it as a sensitivity check.
#
# Reads  <outputs.directory>/cache/bray_css_{genus,all_features}.rds,
#        <outputs.directory>/setup/ps_clean.rds (sample metadata only),
#        <outputs.directory>/tables/host_phylosignal_tree.nwk,
#        <outputs.directory>/tables/host_plant_species.csv (plant family).
# Writes <outputs.directory>/tables/host_phylosymbiosis_tests.csv, _pairs.csv,
#        host_beta_permanova.csv and host_beta_variance_partition.csv.

suppressPackageStartupMessages({
    library(jsonlite)
    library(phyloseq)
    library(ape)
    library(vegan)
    library(dplyr)
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

setting <- function(key, default) {
    value <- config$host_phylosignal[[key]]
    if (is.null(value)) default else value
}
STAGE_COL   <- config$analysis$primary_group
STAGES      <- config$differential_abundance$groups
HOST_COL    <- if (is.null(config$generalist_specialist$host_group)) "Scientific_Name" else
    config$generalist_specialist$host_group
STUDY_COL   <- if (is.null(config$generalist_specialist$study_group)) "study_alias" else
    config$generalist_specialist$study_group
N_PERM      <- config$analysis$permutations
N_THREADS   <- if (is.null(config$analysis$threads)) 1L else as.integer(config$analysis$threads)
MIN_SAMPLES <- setting("min_samples_per_species", 3)
MIN_SPECIES <- setting("min_species", 4)

out_dir <- resolve_path(config$outputs$directory)
tables_dir <- file.path(out_dir, "tables")
DISTANCES <- c(genus = file.path(out_dir, "cache", "bray_css_genus.rds"),
               ASV = file.path(out_dir, "cache", "bray_css_all_features.rds"))
needed <- c(DISTANCES, file.path(out_dir, "setup", "ps_clean.rds"),
            file.path(tables_dir, "host_phylosignal_tree.nwk"),
            file.path(tables_dir, "host_plant_species.csv"))
if (any(!file.exists(needed))) {
    stop("missing input(s): ", paste(needed[!file.exists(needed)], collapse = ", "))
}

tip_name <- function(x) gsub(" ", "_", x)
tree <- read.tree(file.path(tables_dir, "host_phylosignal_tree.nwk"))
# divergence time of two hosts: half the path between their tips on a time tree
divergence <- cophenetic(tree) / 2
families <- read_csv(file.path(tables_dir, "host_plant_species.csv"), show_col_types = FALSE)

meta <- as(sample_data(readRDS(file.path(out_dir, "setup", "ps_clean.rds"))), "data.frame")
meta <- tibble(sample_id = rownames(meta),
               stage = meta[[STAGE_COL]],
               species = tip_name(meta[[HOST_COL]]),
               study = meta[[STUDY_COL]],
               genus = sub(" .*$", "", meta[[HOST_COL]]),
               # cultivar within host species ("sativus" is a cultivar of carrot and of radish);
               # missing where none was recorded
               cultivar = ifelse(is.na(meta$Cultivar) | meta$Cultivar %in% c("", "undefined"), NA,
                                 paste(meta$Cultivar, tip_name(meta[[HOST_COL]]), sep = " | "))) |>
    dplyr::mutate(family = families$Family[match(species, tip_name(families$Species))])
if (anyNA(meta$family)) stop("a host species has no plant family in host_plant_species.csv")

lower <- function(m) m[lower.tri(m)]

#' Multiple regression on distance matrices: OLS on the lower triangles, tested by
#' permuting the rows and columns of the response together (Lichstein 2007). `y` is a
#' square matrix, `xs` a named list of square matrices in the same order.
mrm <- function(y, xs, nperm) {
    X <- cbind(1, vapply(xs, lower, numeric(length(lower(y)))))
    fit <- function(v) {
        b <- qr.coef(qr(X), v)
        c(b[-1], r2 = 1 - sum((v - X %*% b)^2) / sum((v - mean(v))^2))
    }
    observed <- fit(lower(y))
    permuted <- replicate(nperm, {
        i <- sample(nrow(y))
        fit(lower(y[i, i]))
    })
    p <- (rowSums(abs(permuted) >= abs(observed)) + 1) / (nperm + 1)
    names(p) <- names(observed) <- c(names(xs), "r2")
    list(estimate = observed, p = p)
}

#' Mean of a sample-by-sample matrix over all sample pairs between two groups.
between_group_mean <- function(m, groups) {
    g <- model.matrix(~ 0 + factor(groups))
    g <- sweep(g, 2, colSums(g), "/")
    out <- t(g) %*% m %*% g
    dimnames(out) <- list(levels(factor(groups)), levels(factor(groups)))
    out
}

#' Phylosymbiosis tests for one distance matrix and one developmental stage.
phylosymbiosis_one <- function(d, unit, stage_level) {
    md <- meta |> dplyr::filter(sample_id %in% rownames(d), stage == stage_level) |>
        dplyr::filter(dplyr::n() >= MIN_SAMPLES, .by = species)
    hosts <- intersect(tree$tip.label, unique(md$species))
    if (length(hosts) < MIN_SPECIES) return(NULL)
    d <- d[md$sample_id, md$sample_id]
    bray <- between_group_mean(d, md$species)[hosts, hosts]
    # share of the sample pairs between two hosts that also lie between two studies
    between_study <- between_group_mean(outer(md$study, md$study, "!=") * 1,
                                        md$species)[hosts, hosts]
    host_time <- divergence[hosts, hosts] / 100                  # per 100 million years

    pearson <- mantel(as.dist(bray), as.dist(host_time), permutations = N_PERM)
    spearman <- mantel(as.dist(bray), as.dist(host_time), method = "spearman",
                       permutations = N_PERM)
    host_only <- mrm(bray, list(host = host_time), N_PERM)
    with_study <- mrm(bray, list(host = host_time, study = between_study), N_PERM)

    n_samples <- table(md$species)[hosts]
    pair <- which(lower.tri(bray), arr.ind = TRUE)
    list(
        tests = tibble(
            unit = unit, stage = stage_level, n_species = length(hosts), n_samples = nrow(md),
            mantel_r = pearson$statistic, mantel_p = pearson$signif,
            mantel_rho = spearman$statistic, mantel_rho_p = spearman$signif,
            mrm_slope_per_100my = host_only$estimate[["host"]],
            mrm_r2 = host_only$estimate[["r2"]], mrm_p = host_only$p[["host"]],
            mrm2_host_slope = with_study$estimate[["host"]],
            mrm2_host_p = with_study$p[["host"]],
            mrm2_study_slope = with_study$estimate[["study"]],
            mrm2_study_p = with_study$p[["study"]],
            mrm2_r2 = with_study$estimate[["r2"]]),
        pairs = tibble(
            unit = unit, stage = stage_level,
            species_a = hosts[pair[, "row"]], species_b = hosts[pair[, "col"]],
            divergence_my = 100 * host_time[pair], mean_bray = bray[pair],
            between_study_share = between_study[pair],
            n_a = as.integer(n_samples[pair[, "row"]]),
            n_b = as.integer(n_samples[pair[, "col"]])))
}

#' Sample metadata and distances for one unit, with host given as genus, family or
#' cultivar. Samples with no cultivar recorded are left out of the cultivar models only.
host_model_data <- function(d, host_level) {
    md <- meta |> dplyr::filter(sample_id %in% rownames(d)) |>
        dplyr::mutate(host = .data[[host_level]]) |> dplyr::filter(!is.na(host)) |> as.data.frame()
    list(md = md, dist = as.dist(d[md$sample_id, md$sample_id]))
}

#' Marginal PERMANOVA: host after stage, and stage after host. Country is not a term, at the
#' user's request (2026-10-19): it is kept as a description of the dataset only.
permanova_one <- function(d, unit, host_level) {
    x <- host_model_data(d, host_level)
    fit <- adonis2(x$dist ~ host + stage, data = x$md, by = "margin",
                   permutations = N_PERM, parallel = N_THREADS)
    as_tibble(as.data.frame(fit), rownames = "term") |>
        dplyr::transmute(unit = unit, host_level = host_level,
                         term = sub("^host$", host_level, term),
                         df = Df, R2, F, p = `Pr(>F)`, n = nrow(x$md))
}

#' Variance partitioning (distance-based RDA) between host and stage.
variance_partition_one <- function(d, unit, host_level) {
    x <- host_model_data(d, host_level)
    fractions <- varpart(x$dist, ~ host, ~ stage, data = x$md)$part$indfract
    adjusted <- grep("^Adj", names(fractions), value = TRUE)
    if (length(adjusted) != 1) stop("varpart: no adjusted R2 column among ",
                                    paste(names(fractions), collapse = ", "))
    labels <- c(paste(host_level, "alone"), "stage alone", paste(host_level, "and stage"), "unexplained")
    tibble(unit = unit, host_level = host_level, fraction = labels, df = fractions$Df,
           adj_R2 = fractions[[adjusted]], n = nrow(x$md))
}

show <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(where(is.double),
                                                                        ~ round(.x, 3)))),
                          row.names = FALSE)
set.seed(123)
distances <- purrr::map(DISTANCES, function(path) as.matrix(readRDS(path)))

# 1. Phylosymbiosis first: it takes seconds, and its tables are written before the slow part.
symbiosis <- purrr::imap(distances, function(d, unit) {
    purrr::map(STAGES, function(s) phylosymbiosis_one(d, unit, s))
}) |> purrr::list_flatten() |> purrr::compact()
tests <- purrr::map(symbiosis, "tests") |> purrr::list_rbind()
write_csv(tests, file.path(tables_dir, "host_phylosymbiosis_tests.csv"))
write_csv(purrr::map(symbiosis, "pairs") |> purrr::list_rbind(),
          file.path(tables_dir, "host_phylosymbiosis_pairs.csv"))
cat(sprintf("Phylosymbiosis (%d permutations)\n", N_PERM)); show(tests)

# 2. Marginal PERMANOVA and variance partitioning: minutes per model on ~1,800 samples.
# host as genus, plant family or cultivar: host species is not one of the host levels
# analysed (the user's choice, 2026-10-20); the phylosymbiosis tests above stay at species
# level, as the host phylogeny is a species tree
# host attributes are plant family and host genus only, at the user's request (2026-10-22)
models <- tidyr::expand_grid(unit = names(distances), host_level = c("genus", "family"))
timed <- function(what, f) function(unit, host_level) {
    started <- Sys.time()
    out <- f(distances[[unit]], unit, host_level)
    cat(sprintf("  %s, %s level, host as %s: %.1f min\n", what, unit, host_level,
                as.numeric(difftime(Sys.time(), started, units = "mins"))))
    out
}
partition <- purrr::pmap(models, timed("variance partition", variance_partition_one)) |>
    purrr::list_rbind()
write_csv(partition, file.path(tables_dir, "host_beta_variance_partition.csv"))
permanova <- purrr::pmap(models, timed("PERMANOVA", permanova_one)) |> purrr::list_rbind()
write_csv(permanova, file.path(tables_dir, "host_beta_permanova.csv"))
cat("\nMarginal PERMANOVA\n"); show(permanova)
cat("\nVariance partitioning (adjusted R2)\n"); show(partition)

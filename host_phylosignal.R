#!/usr/bin/env Rscript
# host_phylosignal.R - does host plant relatedness shape seed and seedling alpha diversity?
#
# Usage:
#   Rscript host_phylosignal.R --config analysis_config.json [--tree <newick file>]
#
# The host tree is NOT made here. It is a time-calibrated tree of the host species from
# TimeTree (timetree.org: upload the species list, download the Newick file), named in the
# config as inputs.host_tree; --tree overrides it. TimeTree returns subspecies and
# varieties as tips of their own (six for Brassica oleracea), so the tree is first reduced
# to one tip per species.
#
# For each rarefaction depth x developmental stage x alpha index, one mean per host
# species goes onto the tree and two measures of phylogenetic signal are computed with
# phytools::phylosig():
#   Blomberg's K  - tested by shuffling species on the tree (analysis.permutations times)
#   Pagel's lambda - fitted by maximum likelihood, tested against lambda = 0
# Species means rest on very unequal sample numbers, so lambda also gets a bootstrap 95%
# interval: samples are resampled with replacement within each species, lambda is
# re-estimated each time, and the 2.5th and 97.5th percentiles are reported.
#
# Reads  <outputs.directory>/tables/rarefied_alpha_diversity.csv (the pipeline's own
#        rarefied alpha diversity, so richness is compared at equal depth).
# Writes <outputs.directory>/tables/host_phylosignal.csv,
#        <outputs.directory>/tables/host_phylosignal_species_means.csv and
#        <outputs.directory>/tables/host_phylosignal_tree.nwk (the one-tip-per-species
#        tree the tests were run on, which the figure in the master document draws).
#
# Read the result with care: with 16 host species at most, K and lambda have little
# power, and host species is largely tied to study.

suppressPackageStartupMessages({
    library(jsonlite)
    library(ape)
    library(phytools)
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

tree_path <- arg_value("--tree")
if (is.null(tree_path) && !is.null(config$inputs$host_tree)) {
    tree_path <- resolve_path(config$inputs$host_tree)
}
if (is.null(tree_path) || !file.exists(tree_path)) {
    stop("No host tree. Pass --tree <newick file> or set inputs.host_tree in the config.")
}

setting <- function(key, default) {
    value <- config$host_phylosignal[[key]]
    if (is.null(value)) default else value
}
STAGE_COL   <- config$analysis$primary_group
HOST_COL    <- if (is.null(config$generalist_specialist$host_group)) "Scientific_Name" else
    config$generalist_specialist$host_group
N_PERM      <- config$analysis$permutations
N_BOOT      <- setting("bootstrap", 500)
MIN_SAMPLES <- setting("min_samples_per_species", 3)
MIN_SPECIES <- setting("min_species", 4)

tables_dir <- file.path(resolve_path(config$outputs$directory), "tables")
alpha_path <- file.path(tables_dir, "rarefied_alpha_diversity.csv")
if (!file.exists(alpha_path)) stop("Run the pipeline first; missing ", alpha_path)

# Tip labels and species names are matched with underscores for spaces.
tip_name <- function(x) gsub(" ", "_", gsub("'", "", x))

# Author names in the metadata columns are Latin-1, hence the encoding.
alpha <- read_csv(alpha_path, locale = locale(encoding = "latin1"), show_col_types = FALSE) |>
    dplyr::transmute(
        depth = rarefaction_depth,
        stage = .data[[STAGE_COL]],
        species = tip_name(.data[[HOST_COL]]),
        Observed = Observed_rarefied,
        Shannon = Shannon_rarefied,
        # Pielou's evenness; undefined for a sample with a single ASV
        Evenness = ifelse(Observed_rarefied > 1, Shannon_rarefied / log(Observed_rarefied),
                          NA_real_)) |>
    tidyr::pivot_longer(c(Observed, Shannon, Evenness), names_to = "index",
                        values_to = "value") |>
    dplyr::filter(!is.na(value))

tree <- if (grepl("\\.(nexus|nex)$", tree_path, ignore.case = TRUE)) {
    read.nexus(tree_path)
} else read.tree(tree_path)
tree$tip.label <- tip_name(tree$tip.label)
# One tip per species: the first two words of the label. Any one tip of a species stands
# for it at the same depth, provided the species' tips form a clade - checked here.
species_of <- function(x) sub("^([^_]+_[^_]+).*$", "\\1", x)
for (sp in unique(species_of(tree$tip.label))) {
    tips <- tree$tip.label[species_of(tree$tip.label) == sp]
    if (length(tips) > 1 && !is.monophyletic(tree, tips)) {
        stop("tips of ", sp, " do not form a clade in the host tree")
    }
}
tree <- keep.tip(tree, tree$tip.label[!duplicated(species_of(tree$tip.label))])
tree$tip.label <- species_of(tree$tip.label)
missing <- setdiff(unique(alpha$species), tree$tip.label)
if (length(missing)) stop("host species not in the tree: ", paste(missing, collapse = ", "))
tree <- keep.tip(tree, intersect(tree$tip.label, unique(alpha$species)))
write.tree(tree, file.path(tables_dir, "host_phylosignal_tree.nwk"))

#' Phylogenetic signal of one alpha index among the host species of one stage.
#' `d` holds one row per sample: `species` and `value`.
signal_one <- function(d) {
    d <- dplyr::filter(d, dplyr::n() >= MIN_SAMPLES, .by = species)
    if (dplyr::n_distinct(d$species) < MIN_SPECIES) return(NULL)
    tr <- keep.tip(tree, unique(d$species))
    by_species <- split(d$value, d$species)[tr$tip.label]
    means <- vapply(by_species, mean, numeric(1))

    fit_k <- phylosig(tr, means, method = "K", test = TRUE, nsim = N_PERM)
    fit_lambda <- phylosig(tr, means, method = "lambda", test = TRUE)
    boot <- replicate(N_BOOT, {
        resampled <- vapply(by_species, function(v) mean(sample(v, replace = TRUE)),
                            numeric(1))
        phylosig(tr, resampled, method = "lambda")$lambda
    })
    list(
        summary = tibble(n_species = length(means), n_samples = nrow(d),
                         K = fit_k$K, K_p = fit_k$P,
                         lambda = fit_lambda$lambda, lambda_p = fit_lambda$P,
                         lambda_ci_low = unname(quantile(boot, 0.025)),
                         lambda_ci_high = unname(quantile(boot, 0.975))),
        means = tibble(species = names(means), n = lengths(by_species), mean = means))
}

set.seed(123)
groups <- alpha |> dplyr::group_by(depth, stage, index) |> dplyr::group_split()
results <- purrr::map(groups, function(d) {
    out <- signal_one(d)
    if (is.null(out)) return(NULL)
    key <- dplyr::distinct(d, depth, stage, index)
    list(summary = dplyr::bind_cols(key, out$summary), means = dplyr::bind_cols(key, out$means))
}) |> purrr::compact()
if (!length(results)) stop("no stage has ", MIN_SPECIES, " host species with enough samples")

summary_tab <- purrr::map(results, "summary") |> purrr::list_rbind()
means_tab <- purrr::map(results, "means") |> purrr::list_rbind()
write_csv(summary_tab, file.path(tables_dir, "host_phylosignal.csv"))
write_csv(means_tab, file.path(tables_dir, "host_phylosignal_species_means.csv"))

cat(sprintf("Host phylogenetic signal: %d permutations for K, %d bootstraps for lambda\n",
            N_PERM, N_BOOT))
print(as.data.frame(dplyr::mutate(summary_tab, dplyr::across(where(is.double), ~ round(.x, 3)))),
      row.names = FALSE)

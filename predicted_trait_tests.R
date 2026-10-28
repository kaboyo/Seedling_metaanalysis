#!/usr/bin/env Rscript
# predicted_trait_tests.R - do predicted traits go with host breadth?
#
# Usage:
#   Rscript predicted_trait_tests.R --config analysis_config.json
#
# Run predicted_traits.py first: it writes the ASV x predicted-trait table read here.
#
# Binary traits: among the ASVs annotated for the trait (NAs excluded), a binomial GLM of
#   the trait on log2(n_hosts + 1), Wald test, BH-FDR over traits. Effect sizes: the odds
#   ratio per doubling of host breadth with its 95% interval (the predictor is log2, so
#   this is exp of the coefficient), and the rank-biserial correlation between the trait
#   and host breadth. Traits with fewer than `min_annotated` annotated ASVs, or with no
#   variation, are left out.
# Numeric traits: Spearman correlation with n_hosts, BH-FDR over traits.
# n_hosts is capped at `max_hosts` to limit the pull of a long tail.
#
# Reads  <outputs.directory>/tables/predicted_traits_asv.csv, predicted_traits_catalogue.csv
#        and generalist_specialist_asv_host_screen.csv (n_hosts per ASV).
# Writes <outputs.directory>/tables/predicted_trait_breadth_tests.csv.

suppressPackageStartupMessages({
    library(jsonlite)
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
    value <- config$predicted_traits[[key]]
    if (is.null(value)) default else value
}
MIN_ANNOTATED <- setting("min_annotated", 100)
MAX_HOSTS <- setting("max_hosts", 20)
FDR <- setting("fdr", 0.05)

tables_dir <- file.path(resolve_path(config$outputs$directory), "tables")
traits_path <- file.path(tables_dir, "predicted_traits_asv.csv")
if (!file.exists(traits_path)) stop("Run predicted_traits.py first; missing ", traits_path)

catalogue <- read_csv(file.path(tables_dir, "predicted_traits_catalogue.csv"),
                      show_col_types = FALSE)
traits <- read_csv(traits_path, show_col_types = FALSE, na = "NA", guess_max = 5000)
breadth <- read_csv(file.path(tables_dir, "generalist_specialist_asv_host_screen.csv"),
                    show_col_types = FALSE) |>
    dplyr::transmute(ASV, n_hosts = pmin(n_hosts, MAX_HOSTS))
data <- dplyr::inner_join(breadth, traits, dplyr::join_by(ASV))
cat(sprintf("%d ASVs; %d predicted traits; n_hosts capped at %d (largest observed: %d)\n",
            nrow(data), nrow(catalogue), MAX_HOSTS, max(breadth$n_hosts)))

#' One binary trait: GLM on log2 breadth, with odds ratio and rank-biserial correlation.
test_binary <- function(trait) {
    d <- tibble(n_hosts = data$n_hosts, value = as.logical(data[[trait]])) |>
        dplyr::filter(!is.na(value))
    base <- tibble(trait = trait, n_annotated = nrow(d), n_true = sum(d$value))
    if (nrow(d) < MIN_ANNOTATED) return(dplyr::mutate(base, status = "too few annotated"))
    if (base$n_true == 0 || base$n_true == nrow(d)) {
        return(dplyr::mutate(base, status = "no variation"))
    }
    fit <- suppressWarnings(glm(value ~ log2(n_hosts + 1), family = binomial, data = d))
    est <- summary(fit)$coefficients[2, ]
    # rank-biserial: +1 if every ASV with the trait has more hosts than every ASV without
    u <- suppressWarnings(wilcox.test(d$n_hosts[d$value], d$n_hosts[!d$value])$statistic)
    dplyr::mutate(base, status = "tested",
                  odds_ratio_per_doubling = exp(est[["Estimate"]]),
                  or_ci_low = exp(est[["Estimate"]] - qnorm(0.975) * est[["Std. Error"]]),
                  or_ci_high = exp(est[["Estimate"]] + qnorm(0.975) * est[["Std. Error"]]),
                  wald_z = est[["z value"]], p = est[["Pr(>|z|)"]],
                  rank_biserial = unname(2 * u / (base$n_true * (nrow(d) - base$n_true)) - 1))
}

#' One numeric trait: Spearman correlation with host breadth.
test_numeric <- function(trait) {
    d <- tibble(n_hosts = data$n_hosts, value = as.numeric(data[[trait]])) |>
        dplyr::filter(!is.na(value))
    base <- tibble(trait = trait, n_annotated = nrow(d))
    if (nrow(d) < MIN_ANNOTATED) return(dplyr::mutate(base, status = "too few annotated"))
    if (dplyr::n_distinct(d$value) < 2) return(dplyr::mutate(base, status = "no variation"))
    rho <- suppressWarnings(cor.test(d$n_hosts, d$value, method = "spearman"))
    dplyr::mutate(base, status = "tested", spearman_rho = unname(rho$estimate), p = rho$p.value)
}

by_type <- split(catalogue$trait, catalogue$type)
results <- dplyr::bind_rows(
    purrr::map(by_type$binary, test_binary) |> purrr::list_rbind() |>
        dplyr::mutate(type = "binary",
                      q = ifelse(status == "tested",
                                 p.adjust(ifelse(status == "tested", p, NA), method = "BH"), NA)),
    purrr::map(by_type$numeric, test_numeric) |> purrr::list_rbind() |>
        dplyr::mutate(type = "numeric",
                      q = ifelse(status == "tested",
                                 p.adjust(ifelse(status == "tested", p, NA), method = "BH"), NA))) |>
    dplyr::left_join(dplyr::select(catalogue, trait, resource), dplyr::join_by(trait)) |>
    dplyr::relocate(resource, type, .after = trait) |>
    dplyr::arrange(type, q, p)
write_csv(results, file.path(tables_dir, "predicted_trait_breadth_tests.csv"))

tested <- dplyr::filter(results, status == "tested")
cat(sprintf("Screened: %d binary and %d numeric predicted traits; tested: %d and %d\n",
            sum(results$type == "binary"), sum(results$type == "numeric"),
            sum(tested$type == "binary"), sum(tested$type == "numeric")))
cat(sprintf("Significant at FDR < %s: %d binary, %d numeric\n", FDR,
            sum(tested$type == "binary" & tested$q < FDR),
            sum(tested$type == "numeric" & tested$q < FDR)))
print(as.data.frame(dplyr::mutate(
    dplyr::select(dplyr::filter(tested, q < FDR), trait, type, n_annotated,
                  dplyr::any_of(c("n_true", "odds_ratio_per_doubling", "rank_biserial",
                                  "spearman_rho")), q),
    dplyr::across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

#!/usr/bin/env Rscript
# ko_breadth_tests.R - which predicted genes go with host breadth?
#
# Usage:
#   Rscript ko_breadth_tests.R --config analysis_config.json
#
# Run ko_predict.sh first: it runs PICRUSt2 and writes the per-ASV KO predictions.
#
# KO presence/absence (a predicted copy number above zero) is taken for the ASVs whose
# NSTI is at most `max_nsti`; KOs present in all or in none of them are dropped. Each KO
# is modelled against log2(n_hosts + 1) with a binomial GLM: likelihood-ratio test,
# BH-FDR over KOs, and the Wald odds ratio per doubling of host breadth with its 95%
# interval. A significant KO is generalist-enriched if log2 OR >= `log2_odds_ratio`,
# specialist-enriched if log2 OR <= -`log2_odds_ratio`, and unassigned otherwise.
# KOs are annotated with KEGG level-1 and level-2 categories (KEGG BRITE ko00001) and
# with the five-level plant-growth-promoting functions of the PGPg_finder database.
# All gene content here is predicted from the marker gene, not measured.
#
# Reads  <outputs.directory>/picrust2/ (PICRUSt2 output),
#        <outputs.directory>/tables/generalist_specialist_asv_host_screen.csv,
#        resources/kegg/ko00001.json and resources/pgpg/pathways_plabase.txt.
# Writes <outputs.directory>/tables/ko_breadth_tests.csv and ko_enrichment_summary.csv.

suppressPackageStartupMessages({
    library(jsonlite)
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
setting <- function(key, default) {
    value <- config$gene_content[[key]]
    if (is.null(value)) default else value
}
MAX_NSTI <- setting("max_nsti", 0.15)
FDR <- setting("fdr", 0.05)
LOG2_OR <- setting("log2_odds_ratio", 0.5)
MAX_HOSTS <- setting("max_hosts", 20)

out_dir <- resolve_path(config$outputs$directory)
tables_dir <- file.path(out_dir, "tables")
resources <- resolve_path(if (is.null(config$predicted_traits$resources_dir)) "resources" else
    config$predicted_traits$resources_dir)
picrust <- file.path(out_dir, "picrust2")
# PICRUSt2 places each ASV against a bacterial and an archaeal reference and keeps the
# better one: `combined_marker_predicted_and_nsti` holds that final NSTI, and there is one
# `<domain>_KO_predicted` table per domain that received ASVs.
ko_files <- list.files(picrust, pattern = "^(bac|arc)_KO_predicted\\.tsv\\.gz$", full.names = TRUE)
nsti_file <- file.path(picrust, "combined_marker_predicted_and_nsti.tsv.gz")
if (!length(ko_files) || !file.exists(nsti_file)) stop("Run ko_predict.sh first; nothing in ", picrust)

nsti <- as_tibble(fread(nsti_file)) |>
    dplyr::select(sequence = 1, nsti = dplyr::matches("NSTI"))
ko <- purrr::map(ko_files, function(f) fread(f)) |> rbindlist(fill = TRUE)
setnames(ko, 1, "sequence")
setnames(ko, sub("^ko:", "", names(ko)))
if (anyDuplicated(ko$sequence)) stop("an ASV has KO predictions in more than one domain")
kept <- nsti$sequence[nsti$nsti <= MAX_NSTI]
ko <- ko[sequence %in% kept]
cat(sprintf("%d ASVs predicted; %d with NSTI <= %s kept\n", nrow(nsti), nrow(ko), MAX_NSTI))

present <- as.matrix(ko[, -1]) > 0
present[is.na(present)] <- FALSE
rownames(present) <- ko$sequence
share <- colMeans(present)
present <- present[, share > 0 & share < 1, drop = FALSE]
cat(sprintf("%d KOs predicted; %d vary among ASVs and are tested\n", length(share),
            ncol(present)))

breadth <- read_csv(file.path(tables_dir, "generalist_specialist_asv_host_screen.csv"),
                    show_col_types = FALSE) |>
    dplyr::transmute(sequence = sub("^ASV [0-9]+ : ", "", ASV), n_hosts = pmin(n_hosts, MAX_HOSTS))
x <- log2(breadth$n_hosts[match(rownames(present), breadth$sequence)] + 1)
if (anyNA(x)) stop("an ASV with KO predictions has no host breadth")

design <- cbind(1, x)
test_ko <- function(y) {
    fit <- suppressWarnings(glm.fit(design, y, family = binomial()))
    se <- sqrt(diag(chol2inv(fit$qr$qr[seq_len(2), seq_len(2), drop = FALSE])))[2]
    beta <- fit$coefficients[[2]]
    c(beta = beta, se = se,
      lrt_p = pchisq(fit$null.deviance - fit$deviance, df = 1, lower.tail = FALSE))
}
fits <- t(apply(present, 2, test_ko))
results <- tibble(KO = colnames(present), n_asvs = nrow(present),
                  n_present = unname(colSums(present)),
                  log2_odds_ratio = fits[, "beta"] / log(2),
                  odds_ratio_per_doubling = exp(fits[, "beta"]),
                  or_ci_low = exp(fits[, "beta"] - qnorm(0.975) * fits[, "se"]),
                  or_ci_high = exp(fits[, "beta"] + qnorm(0.975) * fits[, "se"]),
                  lrt_p = fits[, "lrt_p"]) |>
    dplyr::mutate(q = p.adjust(lrt_p, method = "BH"),
                  enrichment = dplyr::case_when(
                      q < FDR & log2_odds_ratio >= LOG2_OR ~ "Generalist-enriched",
                      q < FDR & log2_odds_ratio <= -LOG2_OR ~ "Specialist-enriched",
                      TRUE ~ "Unassigned"))

# KEGG: a KO may sit in several level-2 categories; all are listed.
kegg <- fromJSON(file.path(resources, "kegg", "ko00001.json"), simplifyVector = FALSE)
kegg_rows <- purrr::map(kegg$children, function(a) {
    purrr::map(a$children, function(b) {
        kos <- unlist(purrr::map(b$children, function(pathway) {
            vapply(pathway$children %||% list(), function(k) sub(" .*$", "", k$name),
                   character(1))
        }))
        if (!length(kos)) return(NULL)
        tibble(KO = unique(kos), kegg_level1 = sub("^[0-9]+ ", "", a$name),
               kegg_level2 = sub("^[0-9]+ ", "", b$name))
    }) |> purrr::list_rbind()
}) |> purrr::list_rbind()
kegg_by_ko <- kegg_rows |>
    dplyr::summarise(kegg_level1 = paste(unique(kegg_level1), collapse = "; "),
                     kegg_level2 = paste(unique(kegg_level2), collapse = "; "), .by = KO)

# PGPg_finder: the KO is the last part of PGPT_ID ("PGPT0000001-nifA|anfA|vnfA-K02584").
pgpg_rows <- read_tsv(file.path(resources, "pgpg", "pathways_plabase.txt"),
                      show_col_types = FALSE) |>
    dplyr::mutate(KO = stringr::str_extract(PGPT_ID, "K[0-9]{5}")) |>
    dplyr::filter(!is.na(KO)) |>
    dplyr::distinct(KO, Lv1, Lv2, Lv3, Lv4, Lv5)
pgpg_by_ko <- pgpg_rows |>
    dplyr::summarise(dplyr::across(Lv1:Lv5, ~ paste(unique(.x), collapse = "; ")), .by = KO) |>
    dplyr::rename_with(~ paste0("pgpg_", tolower(.x)), Lv1:Lv5)

results <- results |>
    dplyr::left_join(kegg_by_ko, dplyr::join_by(KO)) |>
    dplyr::left_join(pgpg_by_ko, dplyr::join_by(KO)) |>
    dplyr::arrange(q, lrt_p)
write_csv(results, file.path(tables_dir, "ko_breadth_tests.csv"))

classes <- dplyr::select(results, KO, enrichment)
count_by <- function(annotation, scheme, level) {
    annotation |> dplyr::inner_join(classes, dplyr::join_by(KO)) |>
        dplyr::count(category = .data[[level]], enrichment, name = "n_kos") |>
        dplyr::mutate(scheme = scheme, .before = 1)
}
summary_tab <- dplyr::bind_rows(
    dplyr::count(classes, enrichment, name = "n_kos") |>
        dplyr::mutate(scheme = "all KOs", category = "all", .before = 1),
    count_by(kegg_rows, "KEGG level 2", "kegg_level2"),
    count_by(pgpg_rows, "PGPG level 2", "Lv2"),
    count_by(pgpg_rows, "PGPG level 3", "Lv3")) |>
    tidyr::pivot_wider(names_from = enrichment, values_from = n_kos, values_fill = 0)
write_csv(summary_tab, file.path(tables_dir, "ko_enrichment_summary.csv"))

cat(sprintf("\nKO enrichment at FDR < %s and |log2 OR| >= %s:\n", FDR, LOG2_OR))
print(table(results$enrichment))
cat(sprintf("KOs with a KEGG level-2 category: %d; with a PGPG function: %d\n",
            sum(!is.na(results$kegg_level2)), sum(!is.na(results$pgpg_lv1))))

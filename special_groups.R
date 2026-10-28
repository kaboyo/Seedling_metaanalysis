#!/usr/bin/env Rscript
# special_groups.R - where do bacteria with distinctive predicted traits occur?
#
# Usage:
#   Rscript special_groups.R --config analysis_config.json
#
# Run predicted_traits.py and make_generalist_specialist.R first.
#
# Three groups of ASVs are built from the predicted traits (the source text has a fourth,
# mycorrhizal fungi, which needs ITS data this project does not have):
#   anaerobic        BugBase  "Anaerobic == TRUE"
#   nitrogen_fixing  FAPROTAX "nitrogen_fixation == TRUE"
#   small_genome     the smallest 1% of non-NA metaTraits genome sizes
# For each group: ASV counts, the sample-level distribution of its relative abundance,
# its composition by bacterial order, its split into host generalists and specialists,
# and its host loading - the number of member genera detected on each host genus, and how
# that goes with the group's relative abundance there.
# For the heatmap, each ASV is named at the finest taxonomic level it resolves to; labels
# detected on at least `heatmap_min_host_genera` host genera form a label x host genus
# matrix of mean relative abundance, log10-transformed and scaled by row. Each label is
# also tested for a difference among host plant orders (Kruskal-Wallis, BH-FDR).
#
# Reads  <outputs.directory>/setup/ps_clean.rds and tables/predicted_traits_asv.csv,
#        generalist_specialist_asv_host_screen.csv, host_plant_species.csv.
# Writes <outputs.directory>/tables/special_group_*.csv.

suppressPackageStartupMessages({
    library(jsonlite)
    library(phyloseq)
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
    value <- config$special_groups[[key]]
    if (is.null(value)) default else value
}
SMALL_GENOME_SHARE <- setting("small_genome_share", 0.01)
MIN_HOST_GENERA <- setting("heatmap_min_host_genera", 3)
HOST_COL <- if (is.null(config$generalist_specialist$host_group)) "Scientific_Name" else
    config$generalist_specialist$host_group

out_dir <- resolve_path(config$outputs$directory)
tables_dir <- file.path(out_dir, "tables")

# Plant order of each host family (APG IV); the host table carries family only.
HOST_ORDER <- c(Brassicaceae = "Brassicales", Fabaceae = "Fabales", Apiaceae = "Apiales",
                Asteraceae = "Asterales", Solanaceae = "Solanales", Poaceae = "Poales")

ps <- readRDS(file.path(out_dir, "setup", "ps_clean.rds"))
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
rel <- sweep(otu, 2, pmax(colSums(otu), 1), "/")
tax <- as(tax_table(ps), "matrix")
meta <- as(sample_data(ps), "data.frame")
hosts <- read_csv(file.path(tables_dir, "host_plant_species.csv"), show_col_types = FALSE) |>
    dplyr::filter(!grepl("species$", Species))
if (!all(hosts$Family %in% names(HOST_ORDER))) stop("a host family has no plant order here")
sample_host <- tibble(sample = colnames(otu), species = meta[[HOST_COL]]) |>
    dplyr::mutate(host_genus = sub(" .*$", "", species),
                  host_order = unname(HOST_ORDER[hosts$Family[match(species, hosts$Species)]]))

traits <- read_csv(file.path(tables_dir, "predicted_traits_asv.csv"), show_col_types = FALSE,
                   na = "NA", guess_max = 5000)
screen <- read_csv(file.path(tables_dir, "generalist_specialist_asv_host_screen.csv"),
                   show_col_types = FALSE)
genome <- traits$metatraits_genome_size
groups <- list(
    anaerobic = list(definition = "BugBase Anaerobic == TRUE",
                     asvs = traits$ASV[dplyr::coalesce(traits$bugbase_Anaerobic, FALSE)]),
    nitrogen_fixing = list(definition = "FAPROTAX nitrogen_fixation == TRUE",
                           asvs = traits$ASV[dplyr::coalesce(traits$faprotax_nitrogen_fixation,
                                                             FALSE)]),
    small_genome = list(
        definition = sprintf("smallest %g%% of non-NA metaTraits genome sizes",
                             100 * SMALL_GENOME_SHARE),
        asvs = traits$ASV[!is.na(genome) &
                              genome <= quantile(genome, SMALL_GENOME_SHARE, na.rm = TRUE)]))

blank <- function(x) is.na(x) | x == "" |
    grepl("^(Unassigned|unidentified|uncultured|metagenome)", x, ignore.case = TRUE)
# The finest rank at which an ASV has a real name, e.g. "Pantoea" or "f: Erwiniaceae".
RANKS <- c(Genus = "", Family = "f: ", Order = "o: ", Class = "c: ", Phyla = "p: ")
finest_label <- function(asv) {
    vapply(asv, function(a) {
        for (rank in names(RANKS)) {
            if (!blank(tax[a, rank])) return(paste0(RANKS[[rank]], tax[a, rank]))
        }
        "Unclassified"
    }, character(1), USE.NAMES = FALSE)
}

one_group <- function(name) {
    asvs <- intersect(groups[[name]]$asvs, rownames(otu))
    classes <- table(factor(screen$screen_class[match(asvs, screen$ASV)],
                            levels = c("Generalist", "Intermediate", "Specialist")))
    sample_rel <- colSums(rel[asvs, , drop = FALSE])
    by_sample <- dplyr::mutate(sample_host, group = name, relative_abundance = sample_rel)

    summary <- tibble(group = name, definition = groups[[name]]$definition,
                      n_asvs = length(asvs), share_of_asvs = length(asvs) / nrow(otu),
                      generalists = classes[["Generalist"]],
                      intermediate = classes[["Intermediate"]],
                      specialists = classes[["Specialist"]],
                      samples_with_group = sum(sample_rel > 0),
                      median_relative_abundance = median(sample_rel),
                      mean_relative_abundance = mean(sample_rel),
                      max_relative_abundance = max(sample_rel))

    orders <- tibble(asv = asvs, bacterial_order = unname(tax[asvs, "Order"]),
                     mean_rel = rowMeans(rel[asvs, , drop = FALSE])) |>
        dplyr::mutate(bacterial_order = ifelse(blank(bacterial_order), "Unassigned",
                                               bacterial_order)) |>
        dplyr::summarise(n_asvs = dplyr::n(), mean_relative_abundance = sum(mean_rel),
                         .by = bacterial_order) |>
        dplyr::mutate(group = name, share_of_group_asvs = n_asvs / sum(n_asvs), .before = 1) |>
        dplyr::arrange(dplyr::desc(n_asvs))

    # host loading: member genera detected on each host genus, and the group's abundance there
    genus <- unname(tax[asvs, "Genus"])
    named <- !blank(genus)
    loading <- by_sample |>
        dplyr::summarise(n_samples = dplyr::n(),
                         group_relative_abundance = mean(relative_abundance),
                         .by = host_genus)
    loading$member_genera_detected <- vapply(loading$host_genus, function(h) {
        cols <- sample_host$sample[sample_host$host_genus == h]
        present <- rowSums(otu[asvs[named], cols, drop = FALSE]) > 0
        dplyr::n_distinct(genus[named][present])
    }, integer(1))
    loading <- dplyr::mutate(loading, group = name, .before = 1)
    rho <- suppressWarnings(cor.test(loading$member_genera_detected,
                                     loading$group_relative_abundance, method = "spearman"))
    loading_test <- tibble(group = name, n_host_genera = nrow(loading),
                           spearman_rho = unname(rho$estimate), p = rho$p.value)

    # heatmap matrix: label x host genus mean relative abundance
    label <- finest_label(asvs)
    label_rel <- rowsum(rel[asvs, , drop = FALSE], label)                # label x sample
    by_host <- sapply(split(seq_len(ncol(label_rel)), sample_host$host_genus),
                      function(i) rowMeans(label_rel[, i, drop = FALSE]))
    by_host <- matrix(by_host, nrow = nrow(label_rel),
                      dimnames = list(rownames(label_rel), sort(unique(sample_host$host_genus))))
    kept <- rownames(by_host)[rowSums(by_host > 0) >= MIN_HOST_GENERA]
    matrix_long <- NULL
    kruskal <- NULL
    if (length(kept) >= 2) {
        m <- by_host[kept, , drop = FALSE]
        # log10 needs a positive floor: half the smallest non-zero mean in the matrix
        pseudocount <- min(m[m > 0]) / 2
        z <- t(scale(t(log10(m + pseudocount))))
        matrix_long <- as_tibble(as.data.frame(as.table(m)), .name_repair = "minimal") |>
            setNames(c("label", "host_genus", "mean_relative_abundance")) |>
            dplyr::mutate(z = as.vector(z), group = name, pseudocount = pseudocount,
                          .before = 1)
        kruskal <- purrr::map(kept, function(l) {
            test <- kruskal.test(label_rel[l, ], factor(sample_host$host_order))
            tibble(group = name, label = l, statistic = unname(test$statistic),
                   df = unname(test$parameter), p = test$p.value)
        }) |> purrr::list_rbind() |> dplyr::mutate(q = p.adjust(p, method = "BH"))
    }
    cat(sprintf("%-16s %4d ASVs; %d labels on >= %d host genera\n", name, length(asvs),
                length(kept), MIN_HOST_GENERA))
    list(summary = summary, by_sample = by_sample, orders = orders, loading = loading,
         loading_test = loading_test, matrix = matrix_long, kruskal = kruskal)
}

results <- purrr::map(names(groups), one_group)
write_part <- function(part, file) {
    write_csv(purrr::map(results, part) |> purrr::compact() |> purrr::list_rbind(),
              file.path(tables_dir, file))
}
write_part("summary", "special_group_summary.csv")
write_part("by_sample", "special_group_sample_abundance.csv")
write_part("orders", "special_group_order_composition.csv")
write_part("loading", "special_group_host_loading.csv")
write_part("loading_test", "special_group_host_loading_tests.csv")
write_part("matrix", "special_group_heatmap_matrix.csv")
write_part("kruskal", "special_group_host_order_tests.csv")

print(as.data.frame(dplyr::mutate(
    purrr::map(results, "summary") |> purrr::list_rbind() |> dplyr::select(-definition),
    dplyr::across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

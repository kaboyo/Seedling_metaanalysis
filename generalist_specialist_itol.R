#!/usr/bin/env Rscript
# generalist_specialist_itol.R - iTOL annotation files for the generalist and specialist trees.
#
# Usage:
#   Rscript generalist_specialist_itol.R --config analysis_config.json
#
# Run generalist_specialist_trees.sh (the trees) and predicted_traits.py (the traits) first.
#
# iTOL (itol.embl.de) is a website: this script only writes the files. For each tree it
# writes the tree and one annotation file per ring, made with the itol.toolkit package:
#   1 phylum colour strip        2 class colour strip
#   3 predicted aerobic type     4 predicted Gram type       (BugBase, colour strips)
#   5 predicted biofilm formation                            (BugBase, gradient)
#   6 number of host species (bar)
#   7 host plant family (colour strip; specialist tree only, where every bar is 1)
# In iTOL: upload tree.nwk, then drag the numbered .txt files onto the tree in order.
#
# Reads  <outputs.directory>/trees/{generalists,specialists}/tree.nwk,
#        <outputs.directory>/tables/generalist_specialist_asv_host_screen.csv,
#        <outputs.directory>/tables/predicted_traits_asv.csv,
#        <outputs.directory>/tables/host_plant_species.csv, and setup/ps_clean.rds (class).
# Writes <outputs.directory>/itol/{generalists,specialists}/.

suppressPackageStartupMessages({
    library(jsonlite)
    library(phyloseq)
    library(ape)
    library(dplyr)
    library(readr)
    library(itol.toolkit)
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
out_dir <- resolve_path(config$outputs$directory)
tables_dir <- file.path(out_dir, "tables")

asv_hash <- function(x) sub("^ASV [0-9]+ : ", "", x)
blank <- function(x) is.na(x) | x == ""

tax_class <- as(tax_table(readRDS(file.path(out_dir, "setup", "ps_clean.rds"))), "matrix")[, "Class"]
host_family <- read_csv(file.path(tables_dir, "host_plant_species.csv"), show_col_types = FALSE) |>
    dplyr::filter(!grepl("species$", Species)) |>
    dplyr::transmute(host_genus = sub(" .*$", "", Species), host_family = Family) |>
    dplyr::distinct()
traits <- read_csv(file.path(tables_dir, "predicted_traits_asv.csv"), show_col_types = FALSE,
                   na = "NA", guess_max = 5000) |>
    dplyr::transmute(
        ASV,
        # BugBase gives three yes/no oxygen traits; one ring shows which applies
        aerobic_type = dplyr::case_when(bugbase_Anaerobic ~ "Anaerobic",
                                        bugbase_Facultatively_Anaerobic ~ "Facultatively anaerobic",
                                        bugbase_Aerobic ~ "Aerobic",
                                        !is.na(bugbase_Aerobic) ~ "None predicted",
                                        TRUE ~ "No Greengenes match"),
        gram_type = dplyr::case_when(bugbase_Gram_Positive ~ "Gram positive",
                                     bugbase_Gram_Negative ~ "Gram negative",
                                     !is.na(bugbase_Gram_Positive) ~ "None predicted",
                                     TRUE ~ "No Greengenes match"),
        biofilm = bugbase_Forms_Biofilms)
asvs <- read_csv(file.path(tables_dir, "generalist_specialist_asv_host_screen.csv"),
                 show_col_types = FALSE) |>
    dplyr::mutate(id = asv_hash(ASV), Class = unname(tax_class[ASV]),
                  Phylum = ifelse(blank(Phylum), "Unassigned", Phylum),
                  Class = ifelse(blank(Class), "Unassigned", Class)) |>
    dplyr::left_join(traits, dplyr::join_by(ASV)) |>
    dplyr::left_join(host_family, dplyr::join_by(representative_host_genus == host_genus))

#' Write the tree and its annotation files for one class of ASVs.
write_itol <- function(folder, with_host_family) {
    tree_path <- file.path(out_dir, "trees", folder, "tree.nwk")
    if (!file.exists(tree_path)) stop("Run generalist_specialist_trees.sh first; missing ", tree_path)
    target <- file.path(out_dir, "itol", folder)
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    file.copy(tree_path, file.path(target, "tree.nwk"), overwrite = TRUE)
    tips <- read.tree(tree_path)$tip.label
    d <- dplyr::filter(asvs, id %in% tips)
    if (nrow(d) != length(tips)) stop("tree tips without an ASV record in ", folder)

    rings <- list(
        list("1_phylum", "DATASET_COLORSTRIP", "Phylum"),
        list("2_class", "DATASET_COLORSTRIP", "Class"),
        list("3_predicted_aerobic_type", "DATASET_COLORSTRIP", "aerobic_type"),
        list("4_predicted_gram_type", "DATASET_COLORSTRIP", "gram_type"),
        list("5_predicted_biofilm_formation", "DATASET_GRADIENT", "biofilm"),
        list("6_n_hosts", "DATASET_SIMPLEBAR", "n_hosts"))
    if (with_host_family) {
        rings <- c(rings, list(list("7_host_plant_family", "DATASET_COLORSTRIP", "host_family")))
    }
    for (ring in rings) {
        values <- d |> dplyr::select(id, value = dplyr::all_of(ring[[3]])) |>
            dplyr::filter(!is.na(value))
        # itol.toolkit takes a data frame for colour strips but fails on one for numeric
        # datasets, which it only accepts as a tibble
        if (ring[[2]] == "DATASET_COLORSTRIP") values <- as.data.frame(values)
        unit <- create_unit(data = values, key = ring[[1]], type = ring[[2]], tree = tree_path)
        write_unit(unit, file.path(target, paste0(ring[[1]], ".txt")))
    }
    writeLines(c(
        sprintf("iTOL files for the host %s ASVs (%d tips).", folder, length(tips)),
        "", "1. Go to itol.embl.de and upload tree.nwk.",
        "2. Drag the numbered .txt files onto the tree, in order.", "",
        "Rings 3 to 5 are BugBase predictions from the marker gene, not measurements.",
        "Tips are named by the ASV's sequence hash, as in the project's tables."),
        file.path(target, "README.txt"))
    cat(sprintf("%s: %d tips, %d annotation files -> %s\n", folder, length(tips),
                length(rings), target))
}

write_itol("generalists", with_host_family = FALSE)
write_itol("specialists", with_host_family = TRUE)

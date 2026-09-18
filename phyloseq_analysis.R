#!/usr/bin/env Rscript
# phyloseq_analysis.R — Phylogenetic metrics & DESeq2 DA for seedling microbiome pipeline
#
# Usage:
#   Rscript phyloseq_analysis.R --config analysis_config.json \
#       [--filtered-features analysis_results/qc/filtered_feature_ids.csv]
#
# Reads the same analysis_config.json as the Python pipeline.
# Outputs are written into the config's output directory under diversity/ and
# differential_abundance/ subdirectories.

suppressPackageStartupMessages({
    library(jsonlite)
    library(phyloseq)
    library(ape)
    library(picante)
    library(vegan)
    library(ggplot2)
    library(RColorBrewer)
})

# ── Color Palettes (sourced from wild-again project) ─────────────────────────

if (file.exists("wildagain_joint_pellet.R")) {
    source("wildagain_joint_pellet.R")
} else {
    # Fallback definitions if source file is missing
    PRIMARY_PELLET <- c("#6159a4", "#feae01", "#0ea5e9", "#ff6682")
    CORE_PELLET <- c("Biomarker" = "#ff6682", "Core" = "#ce232a", "Flexible" = "#918bc0")
    HIGHLIGHT_PELLET <- c("#ff6682")
}

# Helper function to get appropriate color palette (wraps get_pellet)
get_palette <- function(n_groups, palette_type = "wildagain", context = NULL) {
    if (exists("get_pellet")) {
        pal <- get_pellet(n_groups, palette_type, context)
        # If palette is still too small, supplement with RColorBrewer
        if (length(pal) < n_groups) {
            additional <- RColorBrewer::brewer.pal(min(n_groups - length(pal), 12), "Set3")
            pal <- c(pal, additional[1:(n_groups - length(pal))])
        }
        return(pal)
    }
    # Fallback if get_pellet is missing
    return(RColorBrewer::brewer.pal(min(n_groups, 12), "Set3"))
}

theme_nature_publication <- function(base_size = 8) {
    theme_classic(base_size = base_size) %+replace%
        theme(
            text = element_text(color = "black"),
            axis.title = element_text(color = "black", size = base_size),
            axis.text = element_text(color = "black", size = base_size - 1),
            legend.text = element_text(color = "black", size = base_size - 1),
            legend.title = element_text(color = "black", size = base_size),
            strip.text = element_text(size = base_size - 1, face = "bold"),
            strip.background = element_rect(fill = "white", colour = "white"),
            panel.border = element_rect(fill = NA, colour = "#1f1f1f", linewidth = 0.35),
            panel.grid.major = element_line(linewidth = 0.2, colour = "#d9d9d9"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1),
            plot.subtitle = element_text(size = base_size - 1),
            plot.margin = margin(4, 4, 4, 4)
        )
}


# ── Helpers ──────────────────────────────────────────────────────────────────

load_config <- function(path) {
    config <- fromJSON(path)
    list(config = config, base_dir = dirname(normalizePath(path)))
}

resolve_path <- function(base_dir, configured_path) {
    if (startsWith(configured_path, "/")) {
        return(configured_path)
    }
    normalizePath(file.path(base_dir, configured_path), mustWork = FALSE)
}

# ── Data loading ─────────────────────────────────────────────────────────────

load_asv_counts <- function(path, max_samples = NULL, filtered_features = NULL) {
    read_limit <- if (is.null(max_samples)) -1 else max_samples
    if (requireNamespace("data.table", quietly = TRUE)) {
        counts_df <- data.table::fread(
            path,
            nrows = read_limit,
            data.table = FALSE,
            check.names = FALSE
        )
        rownames(counts_df) <- counts_df[[1]]
        counts_df[[1]] <- NULL
    } else {
        counts_df <- read.csv(path, row.names = 1, check.names = FALSE, nrows = read_limit)
    }

    counts <- as.matrix(counts_df)
    storage.mode(counts) <- "numeric"

    if (!is.null(filtered_features) && length(filtered_features) > 0) {
        common <- intersect(colnames(counts), filtered_features)
        if (length(common) > 0) {
            counts <- counts[, common, drop = FALSE]
        }
    }
    counts
}

build_phyloseq <- function(counts, taxa_path, meta_path, tree_path, filters) {
    # Taxonomy
    taxa_df <- read.csv(taxa_path, row.names = 1, check.names = FALSE)
    common_taxa <- intersect(colnames(counts), rownames(taxa_df))
    if (length(common_taxa) == 0) {
        stop("No overlapping taxa between counts and taxonomy table")
    }
    counts <- counts[, common_taxa, drop = FALSE]
    taxa_df <- taxa_df[common_taxa, , drop = FALSE]

    # Metadata
    meta_df <- read.csv(meta_path, check.names = FALSE)
    if ("Developmental_Stage" %in% colnames(meta_df)) {
        meta_df <- meta_df[meta_df$Developmental_Stage %in% c("Seeds", "Seedlings"), , drop = FALSE]
    }
    rownames(meta_df) <- meta_df$run_accession
    common_samples <- intersect(rownames(counts), rownames(meta_df))
    if (length(common_samples) == 0) {
        stop("No overlapping samples between counts and metadata table")
    }
    counts <- counts[common_samples, , drop = FALSE]
    meta_df <- meta_df[common_samples, , drop = FALSE]

    # OTU table: samples in rows
    otu <- otu_table(counts, taxa_are_rows = FALSE)
    tax <- tax_table(as.matrix(taxa_df))
    samp <- sample_data(meta_df)

    # Tree
    tree_raw <- tryCatch({
        if (grepl("\\.(nexus|nex)$", tree_path, ignore.case = TRUE)) {
            read.nexus(tree_path)
        } else {
            read.tree(tree_path)
        }
    }, error = function(e) {
        message("Warning: could not read tree: ", conditionMessage(e))
        return(NULL)
    })

    if (inherits(tree_raw, "multiPhylo")) {
        tree_raw <- tree_raw[[1]]
    }

    tree <- NULL
    if (!is.null(tree_raw)) {
        # Clean tip labels: remove single quotes that read.nexus might preserve
        tree_raw$tip.label <- gsub("^'|'$", "", tree_raw$tip.label)

        common_tips <- intersect(tree_raw$tip.label, colnames(counts))
        if (length(common_tips) >= 3) {
            tree <- keep.tip(tree_raw, common_tips)
        } else {
            message("Warning: fewer than 3 shared tips between tree and counts, skipping tree")
        }
    }

    # Build phyloseq
    if (!is.null(tree)) {
        ps <- merge_phyloseq(otu, tax, samp, phy_tree(tree))
    } else {
        ps <- merge_phyloseq(otu, tax, samp)
    }

    # Apply taxonomy filters using prune_taxa to avoid subset_taxa NSE issues
    tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
    keep_mask <- rep(TRUE, ntaxa(ps))

    keep_kingdom <- filters$keep_kingdom
    if (!is.null(keep_kingdom) && length(keep_kingdom) > 0 && "Kingdom" %in% colnames(tt)) {
        keep_mask <- keep_mask & (tt$Kingdom %in% keep_kingdom)
    }
    exclude_order <- filters$exclude_order
    if (!is.null(exclude_order) && length(exclude_order) > 0 && "Order" %in% colnames(tt)) {
        keep_mask <- keep_mask & !(tt$Order %in% exclude_order)
    }
    exclude_family <- filters$exclude_family
    if (!is.null(exclude_family) && length(exclude_family) > 0 && "Family" %in% colnames(tt)) {
        keep_mask <- keep_mask & !(tt$Family %in% exclude_family)
    }

    if (sum(keep_mask) < ntaxa(ps)) {
        ps <- prune_taxa(keep_mask, ps)
    }

    # Remove empty samples and taxa after alignment/filtering so downstream
    # distance-based methods and DESeq2 do not receive degenerate inputs.
    ps <- prune_samples(sample_sums(ps) > 0, ps)

    # Remove zero-sum taxa
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)

    ps
}

get_sample_columns <- function(ps) {
    colnames(as(sample_data(ps), "data.frame"))
}

get_available_tax_ranks <- function(ps) {
    colnames(as(tax_table(ps), "matrix"))
}

resolve_tax_rank <- function(ps, requested_rank) {
    available_ranks <- get_available_tax_ranks(ps)
    rank_aliases <- c(
        "Phylum" = "Phyla",
        "Phyla" = "Phylum"
    )

    if (requested_rank %in% available_ranks) {
        return(requested_rank)
    }

    alias <- unname(rank_aliases[[requested_rank]])
    if (!is.null(alias) && alias %in% available_ranks) {
        return(alias)
    }

    NULL
}

# ── Analysis functions ───────────────────────────────────────────────────────

compute_faiths_pd <- function(ps, output_path) {
    if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
        return(list(status = "skipped_no_tree"))
    }

    tryCatch({
        otu <- as(otu_table(ps), "matrix")
        if (taxa_are_rows(ps)) otu <- t(otu)
        tree <- phy_tree(ps)

        pd_result <- pd(otu, tree, include.root = FALSE)
        pd_df <- data.frame(
            sample_id = rownames(pd_result),
            faith_pd = pd_result$PD,
            phylogenetic_richness = pd_result$SR,
            stringsAsFactors = FALSE
        )

        write.csv(pd_df, output_path, row.names = FALSE)
        list(status = "completed", n_samples = nrow(pd_df))
    }, error = function(e) {
        message("Faith's PD error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

# Alpha diversity visualization (wild-again Figure1.Rmd approach)
plot_alpha_diversity <- function(ps, group_col, output_dir) {
    if (!(group_col %in% get_sample_columns(ps))) {
        return(list(status = "skipped_no_group"))
    }

    # Get unique groups and assign colors
    groups <- unique(sample_data(ps)[[group_col]])
    groups <- groups[!is.na(groups)]
    n_groups <- length(groups)
    palette <- get_palette(n_groups, "wildagain", context = "alpha_diversity")
    names(palette) <- groups

    # Estimate richness
    # Observed ASV richness plot (wild-again style)
    p_obs <- plot_richness(ps, x = group_col, measures = "Observed") +
        geom_boxplot(
            aes(fill = .data[[group_col]], color = .data[[group_col]]),
            outlier.shape = NA,
            alpha = 0.7,
            lwd = 0.5,
            show.legend = FALSE
        ) +
        geom_jitter(aes(color = .data[[group_col]]), size = 1.5, width = 0.2) +
        theme_nature_publication() +
        scale_fill_manual(values = palette, aesthetics = c("color", "fill")) +
        labs(x = group_col, y = "Observed ASV Richness") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "bottom"
        )

    # Shannon diversity plot
    p_sh <- plot_richness(ps, x = group_col, measures = "Shannon") +
        geom_boxplot(
            aes(fill = .data[[group_col]], color = .data[[group_col]]),
            outlier.shape = NA,
            alpha = 0.7,
            lwd = 0.5,
            show.legend = FALSE
        ) +
        geom_jitter(aes(color = .data[[group_col]]), size = 1.5, width = 0.2) +
        theme_nature_publication() +
        scale_fill_manual(values = palette, aesthetics = c("color", "fill")) +
        labs(x = group_col, y = "Shannon Diversity Index") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "bottom"
        )

    # Save plots
    ggsave(file.path(output_dir, "alpha_observed.png"), p_obs, width = 8, height = 6)
    ggsave(file.path(output_dir, "alpha_shannon.png"), p_sh, width = 8, height = 6)

    list(status = "completed", plots = c("alpha_observed.png", "alpha_shannon.png"))
}

compute_unifrac <- function(ps, output_dir) {
    if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
        return(list(status = "skipped_no_tree"))
    }

    tryCatch({
        uw <- UniFrac(ps, weighted = FALSE, normalized = TRUE, parallel = FALSE)
        w <- UniFrac(ps, weighted = TRUE, normalized = TRUE, parallel = FALSE)

        write.csv(as.matrix(uw), file.path(output_dir, "unifrac_unweighted_distance.csv"))
        write.csv(as.matrix(w), file.path(output_dir, "unifrac_weighted_distance.csv"))

        list(status = "completed", n_samples = nsamples(ps))
    }, error = function(e) {
        message("UniFrac error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

# Beta diversity visualization (wild-again Figure2.Rmd approach)
plot_beta_diversity <- function(ps, group_col, output_dir) {
    if (!(group_col %in% get_sample_columns(ps))) {
        return(list(status = "skipped_no_group"))
    }

    # Get unique groups and assign colors
    groups <- unique(sample_data(ps)[[group_col]])
    groups <- groups[!is.na(groups)]
    n_groups <- length(groups)
    palette <- get_palette(n_groups, "wildagain", context = "beta_diversity")
    names(palette) <- groups

    tryCatch({
        # PCoA with Bray-Curtis distance
        ord <- ordinate(ps, method = "PCoA", distance = "bray")

        # Extract variance explained
        var_explained <- round(ord$eig / sum(ord$eig) * 100, 2)

        # Plot PCoA (wild-again style)
        p_pcoa <- plot_ordination(ps, ord, type = "samples", color = group_col) +
            geom_point(size = 3, alpha = 0.7) +
            theme_nature_publication() +
            scale_color_manual(values = palette) +
            labs(
                x = paste0("PCoA1 (", var_explained[1], "%)"),
                y = paste0("PCoA2 (", var_explained[2], "%)"),
                color = group_col
            ) +
            theme(
                axis.text = element_text(size = 10),
                legend.position = "right"
            )

        ggsave(file.path(output_dir, "pcoa_bray.png"), p_pcoa, width = 8, height = 6)

        list(status = "completed", plots = c("pcoa_bray.png"))
    }, error = function(e) {
        message("Beta diversity plotting error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

# Taxonomic composition plot (wild-again Figure4.Rmd approach)
plot_taxonomy_composition <- function(ps, rank = "Phylum", group_col = NULL,
                                       n_taxa = 10, output_dir) {
    tryCatch({
        resolved_rank <- resolve_tax_rank(ps, rank)
        if (is.null(resolved_rank)) {
            return(list(
                status = "skipped_missing_rank",
                requested_rank = rank,
                available_ranks = get_available_tax_ranks(ps)
            ))
        }

        # Aggregate at specified taxonomic rank
        ps_agg <- tax_glom(ps, taxrank = resolved_rank)
        if (ntaxa(ps_agg) == 0 || nsamples(ps_agg) == 0) {
            return(list(status = "skipped_no_data"))
        }

        # Transform to relative abundance
        ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x) * 100)
        if (ntaxa(ps_rel) == 0 || nsamples(ps_rel) == 0) {
            return(list(status = "skipped_no_data"))
        }

        # Get top taxa
        taxa_totals <- taxa_sums(ps_rel)
        top_n <- min(n_taxa, length(taxa_totals))
        if (top_n == 0) {
            return(list(status = "skipped_no_taxa"))
        }
        top_taxa <- names(sort(taxa_totals, decreasing = TRUE)[seq_len(top_n)])

        # Filter to top taxa
        ps_top <- prune_taxa(top_taxa, ps_rel)
        if (ntaxa(ps_top) == 0 || nsamples(ps_top) == 0) {
            return(list(status = "skipped_no_taxa"))
        }

        # Create plot
        if (!is.null(group_col) && group_col %in% get_sample_columns(ps)) {
            p <- plot_bar(ps_top, fill = resolved_rank, x = group_col) +
                theme_nature_publication() +
                scale_fill_manual(values = get_palette(n_taxa, "wildagain", context = "composition")) +
                labs(x = group_col, y = "Relative Abundance [%]", fill = resolved_rank) +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    legend.position = "bottom"
                )
        } else {
            p <- plot_bar(ps_top, fill = resolved_rank, x = "Sample") +
                theme_nature_publication() +
                scale_fill_manual(values = get_palette(n_taxa, "wildagain", context = "composition")) +
                labs(x = "Sample", y = "Relative Abundance [%]", fill = resolved_rank) +
                theme(
                    axis.text.x = element_blank(),
                    legend.position = "bottom"
                )
        }

        ggsave(file.path(output_dir, paste0("composition_", tolower(rank), ".png")),
               p, width = 10, height = 6)

        list(status = "completed", plots = paste0("composition_", tolower(rank), ".png"))
    }, error = function(e) {
        message("Composition plotting error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

run_deseq2 <- function(ps, group_col, g1, g2, output_path) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        return(list(status = "skipped_missing_deseq2"))
    }

    tryCatch({
        sd <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
        if (!(group_col %in% colnames(sd))) {
            return(list(status = "skipped_no_group_column"))
        }

        keep <- sd[[group_col]] %in% c(g1, g2)
        if (sum(keep) < 4) {
            return(list(status = "skipped_too_few_samples"))
        }

        ps_sub <- prune_samples(keep, ps)
        ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
        ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
        if (nsamples(ps_sub) < 4 || ntaxa(ps_sub) == 0) {
            return(list(status = "skipped_insufficient_counts"))
        }
        retained_groups <- unique(as.character(sample_data(ps_sub)[[group_col]]))
        retained_groups <- retained_groups[!is.na(retained_groups)]
        if (length(retained_groups) < 2) {
            return(list(status = "skipped_single_group_after_filtering"))
        }

        sample_data(ps_sub)[[group_col]] <- factor(
            sample_data(ps_sub)[[group_col]],
            levels = c(g2, g1)
        )

        otu_mat <- as(otu_table(ps_sub), "matrix")
        if (taxa_are_rows(ps_sub)) {
            count_data <- otu_mat
        } else {
            count_data <- t(otu_mat)
        }
        storage.mode(count_data) <- "integer"

        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = count_data,
            colData = data.frame(sample_data(ps_sub)),
            design = as.formula(paste("~", group_col))
        )
        dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq2::DESeq(dds, quiet = TRUE)

        res <- DESeq2::results(dds, contrast = c(group_col, g1, g2))
        res_df <- as.data.frame(res)
        res_df$feature_id <- rownames(res_df)
        res_df <- res_df[order(res_df$padj), ]

        write.csv(res_df, output_path, row.names = FALSE)
        list(
            status = "completed",
            n_features = nrow(res_df),
            n_significant = sum(res_df$padj < 0.05, na.rm = TRUE)
        )
    }, error = function(e) {
        message("DESeq2 error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

# DESeq2 visualization (volcano plot with wild-again colors)
plot_deseq2_results <- function(res_df, output_path) {
    tryCatch({
        # Add significance labels
        res_df$significant <- ifelse(
            res_df$padj < 0.05 & !is.na(res_df$padj),
            "Significant",
            "Not Significant"
        )

        # Volcano plot (wild-again style with highlight palette)
        # Use pink/red for significant and gray for non-significant
        volcano_colors <- c(
            "Significant" = "#ff6682",
            "Not Significant" = "#a6a6a6"
        )

        p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
            geom_point(alpha = 0.6, size = 2) +
            scale_color_manual(values = volcano_colors) +
            theme_nature_publication() +
            labs(
                x = "Log2 Fold Change",
                y = "-Log10 Adjusted P-value",
                color = "Significance"
            ) +
            theme(
                legend.position = "right",
                axis.text = element_text(size = 10)
            )

        ggsave(output_path, p_volcano, width = 8, height = 6)

        list(status = "completed", plots = basename(output_path))
    }, error = function(e) {
        message("DESeq2 plotting error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

rarefaction_curves <- function(ps, output_path, n_depths = 10) {
    tryCatch({
        otu <- as(otu_table(ps), "matrix")
        if (taxa_are_rows(ps)) otu <- t(otu)

        max_depth <- max(rowSums(otu))
        depths <- unique(round(seq(100, max_depth, length.out = n_depths)))

        results <- list()
        for (d in depths) {
            for (i in seq_len(nrow(otu))) {
                samp_counts <- otu[i, ]
                total <- sum(samp_counts)
                if (total < d) next

                rarefied <- rarefy(samp_counts, sample = d)
                results[[length(results) + 1]] <- data.frame(
                    sample_id = rownames(otu)[i],
                    depth = d,
                    expected_richness = as.numeric(rarefied),
                    stringsAsFactors = FALSE
                )
            }
        }

        if (length(results) > 0) {
            rare_df <- do.call(rbind, results)
            write.csv(rare_df, output_path, row.names = FALSE)
            list(status = "completed", n_points = nrow(rare_df))
        } else {
            list(status = "skipped_no_data")
        }
    }, error = function(e) {
        message("Rarefaction error: ", conditionMessage(e))
        list(status = "failed", error = conditionMessage(e))
    })
}

# ── Main ─────────────────────────────────────────────────────────────────────

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    config_path <- "analysis_config.json"
    filtered_features_path <- NULL

    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--config" && i < length(args)) {
            config_path <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--filtered-features" && i < length(args)) {
            filtered_features_path <- args[i + 1]
            i <- i + 2
        } else {
            i <- i + 1
        }
    }

    cfg <- load_config(config_path)
    config <- cfg$config
    base_dir <- cfg$base_dir

    counts_path <- resolve_path(base_dir, config$inputs$counts)
    taxa_path <- resolve_path(base_dir, config$inputs$taxonomy)
    meta_path <- resolve_path(base_dir, config$inputs$metadata)
    tree_path <- resolve_path(base_dir, config$inputs$tree)
    output_dir <- resolve_path(base_dir, config$outputs$directory)
    max_samples <- config$analysis$max_samples
    group_col <- config$analysis$primary_group
    generate_plots <- !identical(config$outputs$generate_plots, FALSE)

    diversity_dir <- file.path(output_dir, "diversity")
    da_dir <- file.path(output_dir, "differential_abundance")
    figures_dir <- file.path(output_dir, "figures")
    dir.create(diversity_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(da_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

    # Load filtered feature IDs if provided
    filtered_features <- NULL
    if (!is.null(filtered_features_path) && file.exists(filtered_features_path)) {
        ff_df <- read.csv(filtered_features_path, stringsAsFactors = FALSE)
        filtered_features <- ff_df$OTUID
        message("Using ", length(filtered_features), " filtered feature IDs from Python pipeline")
    }

    # Load counts
    message("Loading ASV counts from: ", counts_path)
    counts <- load_asv_counts(counts_path, max_samples, filtered_features)
    message("Counts matrix: ", nrow(counts), " samples x ", ncol(counts), " features")

    # Build phyloseq
    message("Building phyloseq object...")
    ps <- build_phyloseq(counts, taxa_path, meta_path, tree_path, config$filters)
    message("Phyloseq: ", nsamples(ps), " samples, ", ntaxa(ps), " taxa")

    status <- list()

    # Faith's PD
    message("Computing Faith's PD...")
    status$faith_pd <- compute_faiths_pd(ps, file.path(diversity_dir, "faith_pd.csv"))
    message("  Faith's PD: ", status$faith_pd$status)

    # Alpha diversity plots (wild-again style)
    if (generate_plots && !is.null(group_col) && group_col %in% get_sample_columns(ps)) {
        message("Plotting alpha diversity...")
        status$alpha_plots <- plot_alpha_diversity(ps, group_col, figures_dir)
        message("  Alpha diversity plots: ", status$alpha_plots$status)
    }

    # UniFrac
    message("Computing UniFrac distances...")
    status$unifrac <- compute_unifrac(ps, diversity_dir)
    message("  UniFrac: ", status$unifrac$status)

    # Beta diversity plots (wild-again style)
    if (generate_plots && !is.null(group_col) && group_col %in% get_sample_columns(ps)) {
        message("Plotting beta diversity...")
        status$beta_plots <- plot_beta_diversity(ps, group_col, figures_dir)
        message("  Beta diversity plots: ", status$beta_plots$status)

        # Taxonomic composition
        message("Plotting taxonomic composition...")
        status$composition <- list(
            phylum = plot_taxonomy_composition_nature(ps, "Phylum", group_col, figures_dir, top_n = 12),
            family = plot_taxonomy_composition_nature(ps, "Family", group_col, figures_dir, top_n = 30),
            genus = plot_taxonomy_composition_nature(ps, "Genus", group_col, figures_dir, top_n = 30)
        )
        message("  Composition plots: ",
                paste(vapply(status$composition, `[[`, character(1), "status"), collapse = ", "))
    }

    # DESeq2
    da_groups <- config$differential_abundance$groups
    if (!is.null(da_groups) && length(da_groups) == 2) {
        message("Running DESeq2 (", da_groups[1], " vs ", da_groups[2], ")...")
        status$deseq2 <- run_deseq2(
            ps, group_col, da_groups[1], da_groups[2],
            file.path(da_dir, "deseq2_results.csv")
        )
        message("  DESeq2: ", status$deseq2$status)

        # Plot DESeq2 results if successful
        if (generate_plots && status$deseq2$status == "completed") {
            message("Plotting DESeq2 results...")
            res_df <- read.csv(file.path(da_dir, "deseq2_results.csv"))
            status$deseq2_plot <- plot_deseq2_results(
                res_df,
                file.path(figures_dir, "deseq2_volcano.png")
            )
            message("  DESeq2 volcano plot: ", status$deseq2_plot$status)
        }
    } else {
        status$deseq2 <- list(status = "skipped_no_groups")
    }

    # Rarefaction curves
    message("Computing rarefaction curves...")
    status$rarefaction <- rarefaction_curves(
        ps, file.path(diversity_dir, "rarefaction_data.csv")
    )
    message("  Rarefaction: ", status$rarefaction$status)

    # Print JSON status on last line for Python to parse
    status_json <- toJSON(status, auto_unbox = TRUE)
    cat(paste0("\n__R_STATUS__", status_json, "\n"))
}

## Nature-style taxa composition helper
## This is additive and does not change the legacy plotting path above.
plot_taxonomy_composition_nature <- function(ps, rank = "Phylum", group_col = NULL, output_dir,
                                             top_n = NULL) {
    resolved_rank <- resolve_tax_rank(ps, rank)
    if (is.null(resolved_rank)) {
        return(list(
            status = "skipped_missing_rank",
            requested_rank = rank,
            available_ranks = get_available_tax_ranks(ps)
        ))
    }

    ps_agg <- tax_glom(ps, taxrank = resolved_rank, NArm = FALSE)
    ps_rel <- transform_sample_counts(ps_agg, function(x) {
        total <- sum(x)
        if (total > 0) x / total * 100 else x
    })

    df <- psmelt(ps_rel)
    if (nrow(df) == 0) {
        return(list(status = "skipped_no_data"))
    }

    taxon_values <- as.character(df[[resolved_rank]])
    taxon_values[is.na(taxon_values) | taxon_values == ""] <- "Unassigned"
    df$taxon_label <- taxon_values

    if (is.null(top_n) || is.na(top_n) || top_n <= 0) {
        top_n <- if (identical(rank, "Phylum")) 12 else if (identical(rank, "Family")) 30 else 30
    }

    top_taxa <- df |>
        dplyr::group_by(taxon_label) |>
        dplyr::summarise(mean_abundance = mean(Abundance), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(mean_abundance)) |>
        dplyr::slice_head(n = min(top_n, dplyr::n())) |>
        dplyr::pull(taxon_label)

    df <- df |>
        dplyr::mutate(taxon_grouped = ifelse(taxon_label %in% top_taxa, taxon_label, "Other")) |>
        dplyr::group_by(Sample, taxon_grouped) |>
        dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

    order_tbl <- df |>
        dplyr::group_by(taxon_grouped) |>
        dplyr::summarise(mean_abundance = mean(Abundance), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(mean_abundance)) |>
        dplyr::pull(taxon_grouped)
    order_tbl <- c(setdiff(order_tbl, "Other"), "Other")
    df$taxon_grouped <- factor(df$taxon_grouped, levels = order_tbl)

    if (!is.null(group_col) && group_col %in% get_sample_columns(ps)) {
        meta_df <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
        meta_df$Sample <- rownames(meta_df)
        df <- dplyr::left_join(df, meta_df[, c("Sample", group_col)], by = "Sample")
        sample_order <- df |>
            dplyr::distinct(Sample, .data[[group_col]]) |>
            dplyr::arrange(.data[[group_col]], Sample) |>
            dplyr::pull(Sample)
        df$Sample <- factor(df$Sample, levels = sample_order)
    } else {
        df$Sample <- factor(df$Sample, levels = unique(df$Sample))
    }

    palette <- get_palette(length(order_tbl), "wildagain", context = "composition")
    names(palette) <- order_tbl

    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = taxon_grouped)) +
        geom_col(width = 0.85, color = "white", linewidth = 0.1) +
        scale_fill_manual(values = palette, drop = FALSE) +
        scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) +
        labs(
            x = if (!is.null(group_col) && group_col %in% get_sample_columns(ps)) group_col else "Sample",
            y = "Relative Abundance [%]",
            fill = resolved_rank
        ) +
        theme_nature_publication() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "bottom"
        )

    if (!is.null(group_col) && group_col %in% get_sample_columns(ps)) {
        p <- p + facet_grid(stats::as.formula(paste("~", group_col)), scales = "free_x", space = "free_x")
    }

    out_file <- file.path(output_dir, paste0("composition_", tolower(rank), "_nature.png"))
    ggsave(out_file, p, width = 12, height = 6, dpi = 300)

    list(status = "completed", plots = basename(out_file), rank = rank, top_taxa = top_taxa)
}

main()

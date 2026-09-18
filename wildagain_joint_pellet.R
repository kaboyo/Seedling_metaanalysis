# Wild-Again Joint Color Pellet Definition
# 
# Comprehensive color palette system from the wild-again project.
# All palettes organized by context and use case.
#
# Reference: Lobato, C. et al. Wild again: recovery of a beneficial Cannabis seed 
# endophyte from low domestication genotypes. Microbiome 12, 239 (2024).
# https://doi.org/10.1186/s40168-024-01951-5

# ══════════════════════════════════════════════════════════════════════════════
# JOINT PELLET SYSTEM
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. PRIMARY PELLET ────────────────────────────────────────────────────────
# Usage: Alpha diversity plots, general group comparisons (up to 4 groups)
# Context: domestication levels, genotype categories

PRIMARY_PELLET <- c(
    "#6159a4",  # Purple
    "#feae01",  # Yellow/Gold
    "#0ea5e9",  # Blue
    "#ff6682"   # Pink/Red
)


# ── 2. CORE PELLET ───────────────────────────────────────────────────────────
# Usage: Abundance-prevalence plots, biomarker visualization
# Context: Core/Flexible/Biomarker classification

CORE_PELLET <- c(
    "Biomarker" = "#ff6682",   # Pink/Red
    "Core" = "#ce232a",        # Deep Red
    "Flexible" = "#918bc0"     # Muted Purple
)


# ── 3. HEAT TREE PELLET ─────────────────────────────────────────────────────
# Usage: Network and heat tree visualizations

HEAT_TREE_PELLET <- c(
    "#918bc0",  # Muted Purple
    "#ff6682",  # Pink/Red
    "#ce232a"   # Deep Red
)


# ── 4. COMPOSITION PELLET ────────────────────────────────────────────────────
# Usage: Taxonomic composition bar plots (up to 3 taxa)

COMPOSITION_PELLET <- c(
    "#564787",  # Dark Purple
    "#9a93d9",  # Light Purple
    "#BEBAEC"   # Very Light Purple
)


# ── 5. EXTENDED COMPOSITION PELLET ───────────────────────────────────────────
# Usage: Taxonomic composition with multiple taxa (up to 7 taxa)

EXTENDED_COMPOSITION_PELLET <- c(
    "#564787",  # Dark Purple
    "#1BB367",  # Green
    "#0ea5e9",  # Blue
    "#ce232a",  # Deep Red
    "#E65E93",  # Pink
    "#FC8D62",  # Orange
    "#fed680"   # Light Yellow
)


# ── 6. BLUE GRADIENT PELLET ─────────────────────────────────────────────────
# Usage: Sequential data representation, neatmap/composition plots (26 colors)

BLUE_GRADIENT_PELLET <- c(
    # Blue gradient (dark to light)
    "#032230",  # Very Dark Blue
    "#053c55",  # Dark Blue
    "#07567a",  # Medium-Dark Blue
    "#0a719f",  # Medium Blue
    "#0c8bc4",  # Light-Medium Blue
    "#0ea5e9",  # Light Blue
    "#2cb5f2",  # Lighter Blue
    "#51c2f5",  # Very Light Blue
    "#76cff7",  # Pale Blue
    "#9bdcf9",  # Very Pale Blue
    # Red-Orange gradient
    "#db3f29",  # Red
    "#db4e29",  # Red-Orange
    "#db5d29",  # Orange-Red
    "#e1834b",  # Orange
    "#e69a6c",  # Light Orange
    "#ecb18e",  # Peach
    # Purple gradient
    "#444d7e",  # Dark Purple
    "#5964a4",  # Medium-Dark Purple
    "#6159a4",  # Medium Purple
    "#9a93d9",  # Light Purple
    "#b5b0e3",  # Very Light Purple
    # Yellow/Orange
    "#feae01",  # Yellow/Gold
    "#fec74f",  # Light Yellow
    # Pink/Red
    "#ff6682",  # Pink
    "#ff8da2",  # Light Pink
    # Gray
    "#a6a6a6"   # Gray
)


# ── 7. LANDRACE PELLET ───────────────────────────────────────────────────────
# Usage: Landrace-specific core/flexible features (up to 11 taxa)
# Associated taxa shown in comments

LANDRACE_PELLET <- c(
    "#fec74f",  # Sphingomonas
    "#2cb5f2",  # Pelomonas_c
    "#76cff7",  # Ralstonia_c
    "#9bdcf9",  # Ralstonia_c
    "#053c55",  # Burkholderia_c
    "#db3f29",  # Bacillus
    "#ecb18e",  # Staphylococcus
    "#51c2f5",  # Pseudomonas_c
    "#0ea5e9",  # Pantoea_c
    "#07567a",  # Enhydro_c
    "#9a93d9"   # Rhodococcus_c
)


# ── 8. SELECTED LINE PELLET ─────────────────────────────────────────────────
# Usage: Selected line core/flexible features (up to 12 taxa)

SELECTED_LINE_PELLET <- c(
    "#2cb5f2",  # Pelomonas_c
    "#76cff7",  # Ralstonia_c
    "#9bdcf9",  # Ralstonia_c
    "#053c55",  # Burkholderia_c
    "#69d1a2",  # Dyella
    "#db3f29",  # Bacillus
    "#ecb18e",  # Staphylococcus
    "#51c2f5",  # Pseudomonas_c
    "#0ea5e9",  # Pantoea_c
    "#07567a",  # Enhydro_c
    "#a6a6a6",  # BD7-11
    "#9a93d9"   # Rhodococcus_c
)


# ── 9. CROSS HYBRID PELLET ──────────────────────────────────────────────────
# Usage: Cross hybrid core/flexible features (up to 9 taxa)

CROSS_HYBRID_PELLET <- c(
    "#2cb5f2",  # Pelomonas_c
    "#76cff7",  # Ralstonia_c
    "#9bdcf9",  # Ralstonia_c
    "#053c55",  # Burkholderia_c
    "#db3f29",  # Bacillus
    "#51c2f5",  # Pseudomonas_c
    "#0ea5e9",  # Pantoea_c
    "#07567a",  # Enhydro_c
    "#9a93d9"   # Rhodococcus_c
)


# ── 10. INBRED LINE PELLET ──────────────────────────────────────────────────
# Usage: Inbred line core/flexible features (up to 11 taxa)

INBRED_LINE_PELLET <- c(
    "#2cb5f2",  # Pelomonas_c
    "#76cff7",  # Ralstonia_c
    "#9bdcf9",  # Ralstonia_c
    "#053c55",  # Burkholderia_c
    "#aadb84",  # Pseudomonas
    "#51c2f5",  # Pseudomonas_c
    "#0ea5e9",  # Pantoea_c
    "#0c8bc4",  # Kosakonia
    "#0a719f",  # Enterobacter
    "#07567a",  # Enhydro_c
    "#9a93d9"   # Rhodococcus_c
)


# ── 11. ALL DOMESTICATION PELLET ────────────────────────────────────────────
# Usage: Combined features across all domestication levels (up to 8 taxa)

ALL_DOMESTICATION_PELLET <- c(
    "#2cb5f2",  # Pelomonas_c
    "#76cff7",  # Ralstonia_c
    "#9bdcf9",  # Ralstonia_c
    "#053c55",  # Burkholderia_c
    "#51c2f5",  # Pseudomonas_c
    "#0ea5e9",  # Pantoea_c
    "#07567a",  # Enhydro_c
    "#9a93d9"   # Rhodococcus_c
)


# ── 12. FUNCTIONAL GENES PELLET ─────────────────────────────────────────────
# Usage: Gene class visualization (up to 8 categories)

FUNCTIONAL_GENES_PELLET <- c(
    "#db3f29",  # Red
    "#feae01",  # Yellow/Gold
    "#E65E93",  # Pink
    "#5964a4",  # Purple
    "#69d1a2",  # Green
    "#2cb5f2",  # Light Blue
    "#0a719f",  # Medium Blue
    "#cc0025"   # Deep Red
)


# ── 13. HIGHLIGHT PELLET ────────────────────────────────────────────────────
# Usage: Single biomarker highlighting

HIGHLIGHT_PELLET <- c(
    "#ff6682"   # Pink/Red - biomarker highlight
)


# ── 14. EXTENDED PELLET ─────────────────────────────────────────────────────
# Usage: General purpose extended palette (up to 8 groups)

EXTENDED_PELLET <- c(
    "#6159a4",  # Purple
    "#feae01",  # Yellow/Gold
    "#0ea5e9",  # Blue
    "#ff6682",  # Pink/Red
    "#918bc0",  # Muted Purple
    "#ce232a",  # Deep Red
    "#00a878",  # Green
    "#f97b22"   # Orange
)


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Get appropriate pellet based on context
# ══════════════════════════════════════════════════════════════════════════════

get_pellet <- function(n_groups, pellet_type = "wildagain", context = NULL) {
    if (pellet_type == "wildagain") {
        # Context-specific pellet selection
        if (!is.null(context)) {
            pal <- switch(context,
                "alpha_diversity" = PRIMARY_PELLET[1:min(n_groups, 4)],
                "composition" = if (n_groups <= 3) COMPOSITION_PELLET else EXTENDED_COMPOSITION_PELLET[1:min(n_groups, length(EXTENDED_COMPOSITION_PELLET))],
                "heatmap" = HEAT_TREE_PELLET,
                "landrace" = LANDRACE_PELLET[1:min(n_groups, length(LANDRACE_PELLET))],
                "selected" = SELECTED_LINE_PELLET[1:min(n_groups, length(SELECTED_LINE_PELLET))],
                "cross_hybrid" = CROSS_HYBRID_PELLET[1:min(n_groups, length(CROSS_HYBRID_PELLET))],
                "inbred" = INBRED_LINE_PELLET[1:min(n_groups, length(INBRED_LINE_PELLET))],
                "functional" = FUNCTIONAL_GENES_PELLET[1:min(n_groups, length(FUNCTIONAL_GENES_PELLET))],
                "highlight" = HIGHLIGHT_PELLET,
                PRIMARY_PELLET[1:min(n_groups, 4)]  # default fallback
            )
            return(pal)
        }
        
        # Default behavior based on number of groups
        if (n_groups <= 4) {
            return(PRIMARY_PELLET[1:n_groups])
        } else if (n_groups <= 8) {
            return(EXTENDED_PELLET[1:n_groups])
        } else {
            # For large numbers, use blue gradient
            return(BLUE_GRADIENT_PELLET[1:min(n_groups, length(BLUE_GRADIENT_PELLET))])
        }
    }
    
    # Fallback to RColorBrewer Set1
    return(brewer.pal(min(n_groups, 8), "Set1"))
}


# ══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ══════════════════════════════════════════════════════════════════════════════

# Example 1: Alpha diversity boxplot with primary pellet
# ggplot(alpha_data, aes(x = domestication, y = shannon, fill = domestication)) +
#   geom_boxplot() +
#   scale_fill_manual(values = get_pellet(4, context = "alpha_diversity")) +
#   theme_minimal()

# Example 2: Taxonomic composition bar plot
# ggplot(comp_data, aes(x = sample, y = abundance, fill = taxon)) +
#   geom_col(position = "fill") +
#   scale_fill_manual(values = get_pellet(5, context = "composition")) +
#   theme_minimal()

# Example 3: Volcano plot with biomarker highlighting
# ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
#   geom_point() +
#   scale_color_manual(values = get_pellet(2, context = "highlight")) +
#   theme_minimal()

# Example 4: Heat tree visualization
# ggraph(phylo_tree) +
#   geom_edge_link() +
#   geom_node_point(aes(color = abundance), size = 3) +
#   scale_color_manual(values = get_pellet(3, context = "heatmap")) +
#   theme_void()

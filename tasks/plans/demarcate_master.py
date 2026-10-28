"""Insert methodology-section banners into seedling_analysis_master.qmd (Parts 1-3).

Each executable chunk gets a banner as its first lines saying which paragraph of
analysis_results/methodology.md it implements. Comments only: no code changes.
Refuses to run twice (the banners are detected).
"""

import re
import sys
from pathlib import Path

MASTER = Path(sys.argv[1])
RULE = "# " + "=" * 90
MARK = "# METHODS "

TITLES = {
    1: "Data and filtering",
    2: "Diversity, composition and co-occurrence",
    3: "Taxon-stage associations",
    4: "Host phylogenetic signal and phylosymbiosis",
    5: "Stage generalists and specialists (CLAM)",
    6: "Host generalists and specialists (two-step screen)",
    7: "Trees and composition of host generalists and specialists",
    8: "Predicted traits, gene content and special groups",
    9: "Software, statistics and validation",
}


def methods(sections: int | tuple[int, ...], detail: str) -> list[str]:
    """Banner for a chunk that implements one or more methodology paragraphs."""
    nums = (sections,) if isinstance(sections, int) else sections
    head = ", ".join(str(n) for n in nums)
    title = TITLES[nums[0]] if len(nums) == 1 else "see the parts marked below"
    return [RULE, f"{MARK}{head} | {title}", f"#   {detail}", RULE]


def other(kind: str, detail: str) -> list[str]:
    """Banner for a chunk outside the methodology text (support code or unreported analysis)."""
    return [RULE, f"{MARK}- | {kind}", f"#   {detail}", RULE]


SUPPORT = "support code (no analysis)"
UNREPORTED = "NOT IN THE METHODOLOGY TEXT"

CHUNKS: dict[str, list[str]] = {
    # ---- Part 1 ----
    "setup": other(SUPPORT, "knitr options"),
    "load-packages": methods(9, "R packages"),
    "color-palettes": other(SUPPORT, "colour palettes"),
    "figure-captions": other(SUPPORT, "caption registry for every figure"),
    "figure-themes": other(SUPPORT, "figure themes and the TIFF writer"),
    "load-asv-counts": methods(1, "read the ASV count table"),
    "load-metadata-taxonomy": methods(1, "read sample metadata and taxonomy"),
    "build-phyloseq": methods(1, "compile the studies into one phyloseq object"),
    "meta-analysis-scope": methods(1, "studies, host species and countries in the compilation"),
    "meta-overview": methods(1, "overview table of the compiled studies"),
    "geography-map": methods(1, "where the samples came from"),
    "metadata-summary": methods(1, "sample counts per group"),
    "sequencing-depth": methods(1, "reads per sample (basis of the 500-read minimum)"),
    "taxonomy-filtering": methods(
        1, "remove non-bacterial, chloroplast, mitochondrial and cyanobacterial ASVs"),
    "decontamination": methods(
        1, "contaminant removal: not possible, no negative controls were deposited"),
    "data-cleaning": methods(
        1, "remove soil samples, samples < 500 reads, rare ASVs -> 1,886 samples, 2,823 ASVs"),
    "fig1-geographic-donut": methods(1, "geographic overview of the analysed samples"),
    "memory-aware-beta-objects": methods(
        2, "beta diversity: objects for the distances (100 most abundant ASVs for PERMDISP)"),
    "top-phyla": methods(2, "composition: mean of per-sample relative abundances, phyla"),
    "bubble-composition-helpers": methods(2, "composition: plotting helpers"),
    "metadata-overview-publication": methods(1, "per-study table of the filtered dataset"),
    "host-plant-species": methods(1, "table of the 16 host plant species"),
    "taxa-strips-by-study": methods(
        2, "composition: mean relative abundance by study, host species and plant family"),
    "top-genera": methods(2, "composition: mean of per-sample relative abundances, genera"),
    "core-microbiome": other(UNREPORTED, "core microbiome"),
    "indicator-species": methods(3, "step 1 of 2: indicator analysis (IndVal.g, 999 permutations)"),
    "generalist-specialist": methods(
        (5, 6, 7), "CLAM (5), two-step host screen (6), Fisher tests and breadth vs abundance (7)"),
    "alpha-diversity": methods(
        2, "alpha diversity: indices on unrarefied counts (reported values: chunk rarefied-alpha)"),
    "rarefied-alpha": methods(
        2, "alpha diversity: richness, Shannon, Pielou's evenness at 500 and 1,000 reads"),
    "alpha-findings": methods(
        2, "alpha diversity: Kruskal-Wallis tests and effect sizes, overall and within study"),
    "alpha-stage-by-stratum": methods(
        2, "alpha diversity: stage effect within each plant family and cultivar, and across them"),
    "faiths-pd": other(UNREPORTED, "Faith's phylogenetic diversity"),
    "unifrac": other(UNREPORTED, "UniFrac distances"),
    "permanova": other(
        UNREPORTED, "PERMANOVA on raw-count dissimilarities (reported: chunk beta-diversity-css)"),
    "beta-diversity-css": methods(
        2, "beta diversity: CSS, Bray-Curtis (ASV and genus), PCoA, PERMANOVA alone and marginal"),
    "pcoa-by-stage": methods(
        2, "beta diversity: seeds and seedlings ordinated apart, for the figure coloured by host"),
    "nmds": other(UNREPORTED, "NMDS ordination"),
    "pcoa": other(UNREPORTED, "PCoA of raw-count Bray-Curtis and Jaccard dissimilarities"),
    "variable-assessment": other(
        UNREPORTED, "per-variable screen on raw counts (reported: chunk beta-diversity-css)"),
    "umap": other(UNREPORTED, "UMAP ordination"),
    "deseq2": methods(3, "step 1 of 2: DESeq2 screen of ASVs (poscounts size factors)"),
    "volcano-plot": methods(3, "step 1 of 2: volcano plot of the DESeq2 screen"),
    "ancombc": other(UNREPORTED, "ANCOM-BC2"),
    "stage-effects-by-stratum": methods(
        3, "step 2 of 2: CLR, Hedges' g per plant family and cultivar, random-effects pooling"),
    "external-analyses": methods(
        (4, 7, 8), "runs the standalone scripts; each call below is marked with its paragraph"),
    "cooccurrence-network": methods(
        2, "co-occurrence: Spearman network of the 50 most abundant ASVs of each developmental stage"),
    "rarefaction": other(UNREPORTED, "rarefaction curves"),
    "nature-manuscript-package": other(SUPPORT, "figure and table plan for the manuscript"),
    "summary": other(SUPPORT, "run summary"),
    "rendered-output-figure-gallery": other(SUPPORT, "gallery of the rendered figures"),
    "restore-render-working-directory": other(SUPPORT, "restore the working directory"),
    "performance-benchmark": other(SUPPORT, "timing benchmark (not evaluated)"),
    # ---- Part 2 ----
    "fx-setup": other(SUPPORT, "Part 2 setup"),
    "fx-palettes": other(SUPPORT, "Part 2 palettes"),
    "fx-fig-registry": other(SUPPORT, "fig_section(): draw a figure, or record a skip"),
    "fx-load-core": other(SUPPORT, "load the filtered phyloseq object"),
    "fx-fig-library-sizes": methods(1, "figure: library_sizes"),
    "fx-fig-sequencing-depth": methods(1, "figure: sequencing_depth_histogram"),
    "fx-fig-rarefaction": other(UNREPORTED, "figure: rarefaction_curves"),
    "fx-fig-top-phyla": methods(2, "figure: top_phyla_composition"),
    "fx-fig-top-genera": methods(2, "figure: top_genus_composition"),
    "fx-fig-alluvial-stage": methods(2, "figures: composition alluvials and seedling-vs-seed"),
    "fx-fig-composition-dendrogram": methods(2, "figure: composition_dendrogram"),
    "fx-fig-alpha-indices": methods(2, "figure: alpha_indices_by_stage_rarefied"),
    "fx-fig-alpha-richness-faceted": methods(
        2, "figures: alpha richness and evenness by species, at 500 and 1,000 reads"),
    "fx-fig-tree-genus": methods(
        (4, 5, 6, 7, 8), "tree figures and host_phylosignal; each figure below is marked"),
    "fx-fig-umap": other(UNREPORTED, "figure: umap"),
    "fx-fig-pcoa-by-stage": methods(2, "figures: pcoa_css_by_stage, pcoa_css_by_stage_study, pcoa_css_asv_by_stage"),
    "fx-fig-indicator": methods(3, "figure: indicator_species_plot (step 1)"),
    "fx-fig-stage-screen": methods(3, "figure: stage_screen (step 1: DESeq2 and IndVal.g)"),
    "fx-fig-stage-heatmap": methods(3, "figure: stage_effect_heatmap (step 2: Hedges' g)"),
    "fx-fig-cooccurrence": methods(2, "figure: cooccurrence_network"),
    "fx-fig-special-groups": methods(8, "figure: special_group_heatmaps"),
    "fx-fig-functional": methods(8, "figure: functional_host_breadth_heatmaps"),
    "fx-fig-functional-stage": methods(8, "figure: functional_stage_heatmaps"),
    "fx-fig-alpha-stage-strata": methods(2, "figure: alpha_stage_effect_heatmap"),
    "fx-fig-supplementary-panels": methods(
        (2, 6, 7), "figures: beta_variance_explained, host_breadth_null, host_class_composition"),
    "fx-fig-host-breadth-panels": methods(
        (4, 7, 8), "figures: host_phylosymbiosis_scatter, host_breadth_top_genera, three volcano plots"),
    "fx-fig-geography": methods(1, "figures: sample provenance and geographic distribution"),
    "fx-fig-study-scope": methods(1, "figure: study_scope"),
    "fx-fig-study-overview": methods(1, "figure: study_overview (the overview figure)"),
    "fx-manuscript-figs": methods(9, "validation: index of the audit figures"),
    "fx-manifest": other(SUPPORT, "figure manifest: every figure written or skipped"),
    "fx-fig-main-figures": other(
        SUPPORT, "the four main figures, each merged from finished figures with patchwork"),
    "fx-analysis-folders": other(
        SUPPORT, "deposit each analysis's tables and figures together in its own folder"),
    # ---- Part 3 ----
    "ms-run-validation": methods(9, "validation: optionally rerun validation/verify_manuscript.py"),
    "ms-setup": methods(9, "validation: setup; sets MS_READY"),
    "ms-read-npz": methods(9, "validation: load the audit's raw data"),
    "ms-cascade": methods(9, "validation: filtering cascade recomputed"),
    "ms-aggregate": methods(9, "validation: aggregate to genus"),
    "ms-fig1-data": methods(9, "validation: composition recomputed"),
    "ms-fig1": methods(9, "validation: audit figure 1"),
    "ms-fig2-data": methods(9, "validation: alpha diversity recomputed"),
    "ms-fig2": methods(9, "validation: audit figure 2"),
    "ms-fig3-data": methods(9, "validation: beta diversity from the tables of record"),
    "ms-fig3": methods(9, "validation: audit figure 3"),
    "ms-fig4-data": methods(9, "validation: indicator taxa from the tables of record"),
    "ms-fig4": methods(9, "validation: audit figure 4"),
    "ms-validate": methods(9, "validation: reconcile recomputed values with the tables of record"),
    "master-session": methods(9, "session information"),
    "figure-caption-index": other(SUPPORT, "caption index"),
}

# Markers inside the three chunks that serve several paragraphs: (anchor, paragraph, detail).
# The anchor is the start of a line that occurs exactly once; the marker goes above it,
# with the anchor's indentation.
INNER: list[tuple[str, int | None, str]] = [
    ("#' The host tree with one tip per species", 6, "host tree, one tip per species"),
    ("#' The two screens.", 6, "breadth null (screen 1) and host-relatedness null (screen 2)"),
    ("#' How host generalist and host specialist ASVs differ", 7,
     "Fisher's exact tests by rank, top genera, breadth vs abundance"),
    ("perform_generalist_specialist <- function(", 5,
     "CLAM by stage; also calls the host screen (6) and the comparison (7)"),
    ('run_external("host_phylosignal"', 4, "Blomberg's K and Pagel's lambda"),
    ('run_external("host_phylosymbiosis"', 4,
     "Mantel, MRM, marginal PERMANOVA, variance partitioning"),
    ('run_external("generalist_specialist_trees"', 7, "alignment and maximum-likelihood trees"),
    ('run_external("predicted_traits"', 8, "BugBase, FAPROTAX, metaTraits"),
    ('run_external("predicted_trait_tests"', 8, "traits against host breadth"),
    ('run_external("ko_predict"', 8, "PICRUSt2 KEGG orthologue prediction"),
    ('run_external("ko_breadth_tests"', 8, "orthologues against host breadth; KEGG and PGPg"),
    ('run_external("special_groups"', 8, "anaerobic, nitrogen-fixing, small-genome groups"),
    ('run_external("generalist_specialist_itol"', None,
     "iTOL annotation files (trees are drawn in R)"),
    ('fig_section("tree_genus_by_stage"', None, "figure: tree_genus_by_stage"),
    ('fig_section("tree_asv_by_stage"', None, "figure: tree_asv_by_stage"),
    ('fig_section("tree_asv_prevalent_by_stage"', None, "figure: tree_asv_prevalent_by_stage"),
    ('fig_section("tree_generalist_specialist_genus"', 5,
     "figure: tree_generalist_specialist_genus (CLAM class and host screen, genus level)"),
    ("host_screen_tree <- function(", 7, "figures: tree_host_generalists, tree_host_specialists"),
    ('fig_section("tree_generalist_specialist",', 7,
     "figure: tree_generalist_specialist (both classes; rings from paragraphs 6 and 8)"),
    ('fig_section("host_phylosignal"', 4, "figure: host_phylosignal"),
]

lines = MASTER.read_text(encoding="utf-8").split("\n")
if any(line.lstrip().startswith(MARK) for line in lines):
    sys.exit("banners already present; nothing done")

archive_start = next(i for i, line in enumerate(lines) if line.startswith("# Part 4 "))
header = re.compile(r"^```\{r ([A-Za-z0-9_-]+)[,}]")

out: list[str] = []
labelled: list[str] = []
unlabelled: list[str] = []
inner_hits = {anchor: 0 for anchor, _, _ in INNER}
for i, line in enumerate(lines):
    if i < archive_start:
        for anchor, section, detail in INNER:
            if line.lstrip().startswith(anchor):
                indent = line[: len(line) - len(line.lstrip())]
                tag = f"{section} |" if section else f"- | {UNREPORTED}:"
                out.append(f"{indent}{MARK}{tag} {detail}")
                inner_hits[anchor] += 1
    out.append(line)
    match = header.match(line) if i < archive_start else None
    if match:
        label = match.group(1)
        if label in CHUNKS:
            out.extend(CHUNKS[label])
            labelled.append(label)
        else:
            unlabelled.append(label)

missing = sorted(set(CHUNKS) - set(labelled))
bad_inner = {a: n for a, n in inner_hits.items() if n != 1}
too_long = [text for text in out if text.lstrip().startswith(MARK) and len(text) > 100]
if missing or unlabelled or bad_inner or too_long:
    sys.exit(f"not written.\n missing chunks: {missing}\n unlabelled chunks: {unlabelled}\n"
             f" inner anchors not found exactly once: {bad_inner}\n over 100 chars: {too_long}")

MASTER.write_text("\n".join(out), encoding="utf-8")
print(f"{len(labelled)} chunk banners and {len(INNER)} inner markers written; "
      f"{len(out) - len(lines)} lines added")

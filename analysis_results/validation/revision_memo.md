# Revision memo — *Meta-analysis of the Seedling Microbiome*

Every headline statistic in `Manuscript.docx` was independently recomputed from `asv-counts (1).csv` (2,456 × 129,070) and the aligned metadata, using an independent reimplementation of the reported pipeline (`verify_manuscript.py`) rather than a re-run of the original scripts. This memo lists what changed, what survived, and what must be softened or dropped. The original `Manuscript.docx` was **not** modified; corrections are in `Manuscript_revised.docx`.

## Headline verdict

Of 22 audited claims: **6 verified exactly**, **4 verified only on a different analytical basis than the paper states**, **4 are outright errors**, **1 has the right number but the wrong stated mechanism**, **1 table is incomplete against its own inclusion rule**, **4 replicate but rest on strata or models too weak to interpret**, **1 caveat is stated then argued past**, and **1 design-level assumption fails**.

The paper's central ordering — host species and cultivar > developmental stage — **holds and strengthens** under correct analysis. What does not hold is the magnitude of the stage effect, the basis on which the beta-diversity tests were computed, one alpha-diversity test statistic, and the claim that ASV-level cross-study pooling is valid.

## 1. Errors requiring correction

**1.1 Sample count: 2,457 → 2,456.** The metadata lists 2,457 runs but `SRR18176465` (a *Glycine max* seedling from PRJNA811248) has no row in the count matrix. Every downstream count is unaffected; only the retrieval figure changes. *Fixed in Abstract and Methods.*

**1.2 Shannon Kruskal–Wallis H = 38.2 → 36.4 (n = 1,886 → 1,883).** The reported H = 38.234 is reproducible only by assigning Shannon = 0 to the three samples with no retained reads, where the index is undefined. Excluding them gives H = 36.357, p = 1.64 × 10⁻⁹ at n = 1,883. The conclusion is unchanged; the statistic is not. Note this is why the Shannon and Observed tests legitimately carry different n — a point the original text never explained. *Fixed in Abstract, Results and Methods.*

**1.3 Table 4 basis: genus-level → ASV-level raw counts.** This is the most consequential finding. The caption states genus-level Bray–Curtis, but all five PERMANOVA F values, all five R², and all five PERMDISP F values reproduce to four significant figures **only** on Bray–Curtis distances among the 100 most abundant **ASVs**, computed on **raw counts** without relative-abundance normalisation. Genus-level distances give materially different values (stage F = 37.70, R² = 0.020 rather than F = 19.26, R² = 0.012). The numbers in the table are arithmetically correct; the caption describing them is not. *Caption corrected; Methods now state which object each analysis used.*

**1.4 Table 4 sample sizes: 1,883/1,805 → 1,553/1,529.** The `n` column was populated from the genus-level object while the F and R² came from the ASV-level one. 330 samples carry no reads among the top 100 ASVs and are dropped from the distance matrix. *Fixed in the table and explained in the caption.*

**1.5 Filtering mechanism: 1,886 → 1,883 is not soil removal.** The manuscript states the 48 soil controls account for this drop. In the pipeline's own step order, soil is removed at step 5, *before* the 500-read depth floor at step 6 — so 1,886 is already soil-free. The 3-sample difference is samples with zero reads among the retained features. *Fixed in Methods and Results; a per-step accounting is provided as Supplementary Table S5.*

**1.6 Actinobacteriota prevalence: 937 → 940 samples.** All other Table 2 values verified exactly. *Fixed.*

**1.7 Table 5 is incomplete against its own rule.** The caption states all family/genus pairs with adjusted p < 0.05 are listed; 9 pairs qualify and 8 are listed. Asteraceae/*Sphingomonas* (ρ = 0.759, q = 0.031, n = 10) was omitted. *Row added.* All 8 originally listed rows reproduce exactly, including n and p_adj.

## 2. Claims that must be softened

**2.1 The stage effect on beta diversity (R² = 0.012) is an upper bound, not an estimate.** All five Table 4 tests are univariate, so each attributes to its variable everything the others would also have explained. In a marginal (type-III) model fitting study, stage, and host species simultaneously, stage falls to **R² = 0.003** — a 4.8-fold reduction — while host species holds at 0.077 and study at 0.029. The ordering survives; the number does not. *Table 4b added; Results and Abstract now report both.*

**2.2 The Observed-richness stage effect does not exist unadjusted.** Fitted alone, β = −0.15, p = 0.91. It becomes significant (β = −10.8, p = 3.0 × 10⁻¹⁴) only once log depth and study enter the model, with adjusted R² rising from ≈0 to 0.46. The unadjusted comparison cannot be presented as a stage effect. Shannon is more robust but still shrinks 39% under the same adjustment (β = −0.197 → −0.120). The 500-read depth floor is permissive and richness correlates with library size at ρ = 0.30. *Depth-controlled models added to Methods and Results.*

**2.3 The PERMDISP caveat is stated and then argued past.** The original text acknowledges dispersion significance, then proceeds to interpret stage-level location differences. With PERMDISP significant for every grouping (stage F = 21.49, p = 0.003), a significant PERMANOVA F does not license a centroid-shift reading. *Results and Discussion now decline that interpretation explicitly rather than noting the caveat and continuing.*

**2.4 Genus-level phrasing of indicator taxa conceals within-genus direction splits.** All six named genera replicate with the stated direction, but four contain significant ASVs in **both** directions: *Pantoea* 58 seedling-side vs 15 seed-side, *Pseudomonas* 103 vs 26, *Erwinia* 17 seed-side vs 4 seedling-side, *Exiguobacterium* 18 vs 1. Only *Enterobacter* (14) and *Kosakonia* (7) are unanimous, and both are small sets. A genus label describes a majority, not a property. Separately, the **twelve strongest indicators in the whole analysis (IndVal.g 0.521–0.594) are all seed-side** — the paper's seedling emphasis inverts the strength ordering. *Both points added to Results and Conclusion.*

**2.5 Two Table 5 rows are not interpretable as biology.** Both Asteraceae rows rest on 10 *Lactuca sativa* samples from one study (7 seeds, 3 seedlings) with the genus detected in only 2 samples; ρ = 0.759 is a two-point rank correlation. At cultivar level, the strongest reported value (*Raphanus sativus* cv. radish, *Erwinia* ρ = −0.995, q = 1.1 × 10⁻²⁸) has 32 samples of which 5 are seeds — near-perfect rank separation between very unequal groups, not a graded relationship. *Both flagged in caption and Results as hypothesis-generating.* Note the manuscript states ρ = −0.996; recomputation gives −0.9954.

**2.6 Study identity outranks developmental stage in every analysis and should be reported as such.** Alpha diversity ε² = 0.254 vs 0.019; univariate beta R² = 0.070 vs 0.012; marginal beta R² = 0.029 vs 0.003. Author and study are near-collinear (most studies have one first author) and should be read as one provenance effect. Study identity is also not separable from host species, since each study contributed a distinct host set — the marginal model apportions variance under a model assumption rather than resolving the confound. *Added to Results, Discussion and Limitations.*

## 3. Design-level finding: ASV-level cross-study pooling does not hold

This is not a numerical error but it conditions how the whole beta-diversity section should be read, and it is the single most important addition to the paper.

Of **119,135 ASVs** observed anywhere in the compilation, **117,696 (98.79%) occur in exactly one study**; 1,224 occur in two, 172 in three, 43 in four, and **none in more than four of the eight studies**. Median pairwise Jaccard similarity of ASV presence between studies is **0.000** (maximum 0.014). Per-study private-feature fractions run from 89.4% (PRJEB6729) to 100.0% (PRJNA735565, PRJNA784432).

Aggregated to genus, the same eight studies share substantially more (median Jaccard **0.168**, maximum 0.318).

Two consequences. First, some of this is real host specificity and biogeography, but it is equally what differing primer choices, amplicon regions, read lengths, and independent denoising runs produce — and those causes cannot be separated with the metadata the source studies deposited. Second, **the ASV-level distance matrix is exactly the basis Table 4 was computed on** (finding 1.3), so those tests are partly measuring study provenance. The defensible position is the one the revision now takes: report both levels, rely on the agreement in *ordering* between them, and treat genus-level aggregation as the appropriate resolution for cross-study inference in a compilation assembled this way. *Added as a framing result before the beta-diversity claims, and as Limitation (iv).*

## 4. Claims that verified cleanly

- 1,886 samples after filtering — reproduces exactly through all 10 pipeline steps.
- Observed-richness Kruskal–Wallis H = 13.1, p = 2.9 × 10⁻⁴, n = 1,886 — exact, and n = 1,886 is correct here because richness is defined at zero depth.
- Table 2 phylum composition — all values exact except the one prevalence noted in 1.6.
- Table 3 top-ten genera — all 10 rows exact (total abundance, mean, prevalence).
- Table 5 — all 8 listed rows exact to three decimals including n and p_adj.
- All 5 Table 4 F, R² and PERMDISP F values — arithmetically exact on the ASV basis, including PERMDISP computed as distances to the group spatial median (`betadisper(type="median")`), which is what reproduces 21.488 vs the reported 21.49.
- All six indicator genera replicate at p = 0.001 with the stated direction.

## 5. What was changed in `Manuscript_revised.docx`

| Location | Change |
|---|---|
| Abstract, methods sentence | 2,456 recoverable; ASV-level beta basis stated; depth covariate and marginal model added |
| Abstract, results | Shannon H corrected to 36.4/n=1,883; ε² added; marginal R² added; PERMDISP p corrected; within-genus split noted |
| Abstract, implication | Study-identity ordering and the 98.8% single-study diagnostic added |
| Methods, dataset assembly | SRR18176465 absence documented |
| Methods, analysis | Filter order and depth-floor position stated; ε² and depth-controlled models added; two distance objects distinguished; spatial-median PERMDISP specified; marginal model and comparability diagnostic added; indicator direction-reporting rule added |
| Results, dataset scope | Filtering mechanism corrected; comparability diagnostic added as a new paragraph |
| Results, composition | Mean-vs-median bimodality explained; 56.6% read share of top three genera added |
| Results, alpha | H corrected; ε² added; differing n explained; depth-controlled models added as a new paragraph |
| Results, alpha (host) | ε² values and the provenance/host non-separability added |
| Results, beta | Basis and n corrected; marginal model added; dispersion interpretation declined; genus-level sensitivity added (3 new paragraphs) |
| Results, PCoA | Basis mismatch between Figure 4 and Table 4 disclosed; genus-level values given |
| Results, indicators | 811/2,000 significant ASVs; within-genus direction counts; seed-side strength ordering (2 new paragraphs) |
| Results, correlations | Small-stratum rows flagged; omitted row noted; cv. radish fragility added (1 new paragraph) |
| Discussion, opening | Marginal-model R² and study-identity ordering added |
| Discussion, caveats | Expanded from two to three, with the comparability diagnostic as the second (2 new paragraphs) |
| Limitations | Expanded from five items to eight |
| Conclusion | Within-genus splits, seed-side strength ordering, and both beta R² values added |
| Table 2 | Actinobacteriota prevalence 937 → 940; caption notes bimodality |
| Table 3 | Caption states n and read share |
| Table 4 | Caption rewritten to ASV basis + spatial median; n corrected to 1,553/1,529; PERMDISP p 0.001 → 0.003 |
| Table 4b | **New** — marginal type-III PERMANOVA |
| Table 5 | Asteraceae *Sphingomonas* row added; caption flags the two-point correlation |

## 6. Not addressed

Figures in the manuscript file were not replaced; corrected versions are supplied separately as `fig1_composition_corrected.png` through `fig4_indicator_corrected.png`. Supplementary Tables S1–S5 are referenced by the revised text and supplied as the recomputed CSVs in `work/`. The PERMDISP p-values are reported as 0.003 (the 1/300 floor of the permutation scheme actually used) rather than 0.001; if 999 permutations are used for the final version this becomes 0.001 and the table should be updated to match.

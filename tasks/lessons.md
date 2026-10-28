# Lessons

## 2026-09-28 — ExecPlan skipped on a 3-file change

**What happened:** "Address the problems" from a code review led straight to editing
`seedling_analysis_master.qmd`, `analysis_config.json` and `analysis_config.smoke.json`
with no ExecPlan. The plan was written afterwards
(`tasks/plans/fix-richness-and-genus-figures.json`), so the one data-affecting choice —
restricting the genus figures to studies that sampled both stages — got no approval gate.

**Rule:** Any task touching 3+ files gets an ExecPlan in `tasks/plans/` BEFORE the first
edit, however the request is phrased. The user asked for this plan in JSON; confirm
whether JSON is the format for all future plans.

**Rule:** When a fix changes what a figure shows scientifically (which samples or studies
count), treat it as a decision to put to the user, not an implementation detail.

## 2026-10-02 — Methodology written as bullets and tables, not prose

**What happened:** Asked to write the methodology into the results folder, the first
version of `analysis_results/methodology.md` used numbered sections, bullet lists and
tables. The user asked for it to be rewritten as continuous prose.

**Rule:** Methods, results and other manuscript-type text for this project is written as
continuous paragraphs, with numbers worked into the sentences. Keep bullets and tables for
chat replies and for reference material, unless asked otherwise.

## 2026-10-02 — Guessed the user's writing voice from the wrong source

**What happened:** Asked to "make sure my voice is heard" in `analysis_results/methodology.md`,
the style was inferred from `analysis_results/manuscript/Manuscript.docx` (first-person "we",
American spelling). The user then supplied their own Introduction page as the reference.

**Rule:** When asked to write in the user's voice, ask for or use a sample they point to;
do not infer it from a file in the repo. Their voice: British spelling, impersonal register,
long semicolon-joined sentences, connectives such as "albeit", "amidst", "this, in part,".

## 2026-10-02 — Simpson and Fisher's alpha dropped from alpha diversity

**What happened:** The methodology and the pipeline reported five alpha indices. The user
asked for Simpson and Fisher's alpha to be removed from the calculations and the master
script; alpha diversity is now Observed richness, Shannon and Pielou's evenness only.

**Rule:** Do not reintroduce Simpson or Fisher's alpha in Parts 1–3 of
`seedling_analysis_master.qmd` or in manuscript text. Parts 4–9 are archival and untouched.

## 2026-10-02 — Tree rings could not be read against their legends

**What happened:** The tree figures drew phylum, genus and stage from one palette, so a
colour meant different things in different rings (purple = Proteobacteria in ring 1 and
Enterobacter in ring 2; navy = a phylum and also Seeds). The user reported that the rings
did not correspond to the legend, and asked that the graphic come from a code script.

**Rule:** In a figure with several colour legends, no colour may appear in two of them, and
the stage colours (navy, vermillion) are reserved for stage. Keep a taxon's colour fixed
across related figures.

**Rule:** Every figure must be reproducible from a committed script, not from scratch code.
`make_tree_figures.R` runs the master's tree chunk; do the same for other single figures.

## 2026-10-06 — Overview figure still had too much empty space after a first tidy-up

**What happened:** Asked to "optimize the spaces" in `study_overview.tiff`, the first pass
only trimmed the existing 11.5 in canvas (moved a legend, shortened one row). The user came
back with "there is too much space; optimize like Nature figures". The fix was to rebuild the
composite at double-column width (7.2 in / 183 mm) with 5-7 pt text and 8 pt panel letters,
in two rows instead of three.

**Rule:** When asked to tighten a composite figure, size it to the journal page first
(183 mm double column, 5-7 pt text) and fit the panels to that, rather than trimming gaps on
an oversized canvas. Check every in-bar and in-slice label still fits after shrinking.

## 2026-10-06 — Legends on the generalist/specialist tree were not clearly visible

**What happened:** While the generalist/specialist dendrogram was being built, the user
asked to "make sure legends are visible". One intermediate version had its bottom legends
running off both sides of the canvas, the legend text was 6 pt, and the palest key
(grey95, "Too rare to test") disappeared against the white page.

**Rule:** After every render, check the legends specifically: all of them inside the
canvas, text at least 7 pt with bold titles, and every key outlined so pale classes show.
With `ggnewscale`, `override.aes = list(fill = ...)` is silently ignored because `fill` is
renamed; outline keys with `colour` instead.

**Rule:** The user also asked for the fonts on that figure to be "slightly bigger". On tree
figures use tip labels of about 6 pt (`size = 2.1`), legend text 8 pt and legend titles
8.5 pt bold, not the 4.5-6 pt the earlier trees use.

## 2026-10-06 — Heatmap gradient on the generalist/specialist tree lacked contrast

**What happened:** The host heatmap drew mean relative abundance as a continuous
square-root gradient. Over half the detected cells are below 0.1% and the top one is 86%,
so nearly every cell came out the same pale shade. The user reported the gradient "is not
well contrasting".

**Rule:** For abundance heatmaps on this dataset, check the distribution of the values
first. Where they span several orders of magnitude, use a logarithmic scale and a
multi-hue palette (viridis, darkest = most abundant).

**Rule:** Keep the continuous scale bar. A first fix replaced it with four discrete classes
and the user asked for the scale bar back with a different colour contrast.

## 2026-10-06 — Class colours in the two inner rings of the tree could not be told apart

**What happened:** Stage classes without multi-study support were faded to 30% opacity.
Faded navy, red and tan looked alike and blended with the greys of the host-breadth ring
beside them. The user reported the colour distinction "is not clear" and asked for
"clearly distinct colours".

**Rule:** Never encode a second variable by fading a categorical colour. Keep every class
a solid colour that is distinct both within its ring and from the neighbouring ring (here:
navy, vermillion, yellow / magenta, black, grey, white), and carry the second variable by
shape or size (here: bar length).

## 2026-10-06 — Tree figures left blank bands and their legends did not explain everything

**What happened:** After the two host generalist/specialist trees were drawn, the user
asked again to arrange the figures to maximise space and to strengthen the legends. The
specialist tree sat in the middle of its canvas with blank bands above and below, and the
outer bars and host-family rings had no legend key, only a line in the subtitle.

**Rule:** A legend must carry a key for every ring and bar in a figure; a subtitle is not
a legend. Legend text 8.5 pt, titles 9.5 pt bold.

**Rule:** A fan tree keeps a blank ring for tip names. Where a tree has no tip names, pull
the subtitle and legends into that ring and let the tree take the full page width.

## 2026-10-07 — The methodology did not match the code in six places

**What happened:** The user asked whether `analysis_results/methodology.md` corresponds to
the code. A line-by-line check found six mismatches, some in paragraphs written the day
before: a test described as a z-test that the code ran as a rank test, "17 host species"
(the 17th was soil), a wrong statistic name, "all analyses in R" when some were Python and
external programs, a figure described that was not drawn, and analyses said not to have
run that had. The user then asked for the mismatches to be addressed.

**Rule:** Before handing over methods text, check every claim against the code, the config
and the output files: counts, test names, software, and which figures exist. Do it in the
same turn the text is written, not when asked.

**Rule:** A Monte Carlo test with n draws cannot give p below 1/(n+1). Before reporting
"nothing significant after FDR", compare that floor with the number of tests; where the
null can be enumerated, enumerate it.

## 2026-10-07 — Stratify the taxon-stage test by plant family and cultivar, not by study

**What happened:** For the within-stratum stage effect and its forest plot, the analysis
was first built with study as the stratum. The user asked to "make plant family or cultivar
the variable not study".

**Rule:** In this project the user's strata for taxon-stage association are plant family
and cultivar (as in the top-5 genus correlations). Use those; handle study as a random
effect or a recorded caveat, not as the stratifying variable, unless asked.

## 2026-10-07 — Methods text too long, and a forest plot where a heatmap was wanted

**What happened:** The methodology kept growing with reasons, caveats and descriptions of
figures. The user asked twice for less ("only write what was done", then "as concise as
possible"). The stage-effect result was drawn as a forest plot with 12 cultivar symbols;
the user asked for a heatmap instead.

**Rule:** Methods text for this user states what was done, with settings and software, and
nothing else: no rationale, no caveats, no figure walkthroughs. Keep those for captions,
legends and the discussion.

**Rule:** For an effect across many strata (genera by cultivars or families), draw a
heatmap with the strata as columns, not a forest plot with a symbol per stratum.

## 2026-10-12 — "effect size of the cultivars and families" meant the stage effect within them

- **What happened:** asked to "calculate the effect size of the cultivars and families (host)", I
  added plant family as one more grouping variable in the overall Kruskal-Wallis tests. The user
  corrected it: "the effect size was not tested across the different cultivars and plant families".
  They wanted the seed-versus-seedling effect estimated inside each cultivar and plant family and
  then tested across them (my reading of the correction).
- **Rule:** in this project, an "effect size of/across cultivars and plant families" is the
  developmental-stage effect stratified by them (per-stratum effect, pooled estimate, heterogeneity
  test), the design of `stage-effects-by-stratum` — not the variable's own overall effect. When the
  phrase is ambiguous, deliver the stratified stage effect and say which reading was taken.

## 2026-10-13 — figure legends were too long

- **What happened:** the legends I wrote for the main, supplementary and component figures carried
  results, method detail (alignments, models, thresholds) and repeated caveats. The user asked:
  "on all figures, reduce the excessive description of legends". The main legends were cut by half
  and the component legends by about 60%.
- **Rule:** a figure legend states what is drawn and how to read it (panels, marks, colours,
  n, the test named on the panel, abbreviations) and stops. Results go in the Results, method
  detail in the Methods, and a caveat is given once in the text, not in every legend. The caption
  registry in the master (`FIGURE_CAPTIONS`) follows the same rule.

## 2026-10-13 — no supplementary figure for data a main figure already shows

- **What happened:** when the composition dendrogram replaced the class Sankey in Figure 2, I
  moved the Sankey into the supplementary set, which already held the genus Sankeys and the
  trait and PGPg heatmaps (the same tests as Figure 5D, E). The user: "for data with figures
  already present in 00_main_figures, do not process it into supplementary".
- **Rule:** the supplementary set holds figures only for data that no main-figure panel shows.
  When a figure leaves a main figure, or a main panel is added, check the supplementary list in
  the config (`supplementary_figures`) for the same data and take it out. Component figures stay
  in their analysis folders. The same rule was given earlier for the supplementary tables.

## 2026-10-13 — keys inside a figure name the variable only

- **What happened:** the keys of the circular figures were titled "Stage (inner ring)",
  "Class (%; inner bars)", "Host family (outer rings)" and the like. The user: "instead of
  labelling: inner ring ... on the legend, add it in the captions". On the PCoA panels likewise:
  "add only P and R2 values, the rest into the legend caption".
- **Rule:** a key title is the variable's name and nothing else; a panel carries only the
  headline statistic (R2 and p). Where a ring or bar sits, what a short bar means, the other
  test statistics and similar reading instructions go in the written legend.

## 2026-10-13 — Nature figure style, one host colour scheme, one key

- **What happened:** the figures followed the Cell Reports draft (capital panel letters, a
  heading on each panel) and each figure coloured hosts on its own, so plant family was plum in
  one figure and navy in another, and Figure 2 carried a key under every PCoA column. The user:
  "all figures should be nature type alignment and style", "the colours on plant species and
  family on fig 2 should be the same as the pcoa", "use one legend".
- **Rule:** panel letters are lower case, bold, 8 pt, top left, with no heading where the legend
  can say it; legends cite panels as **a**, **b**. A host level has one colour in every figure
  (`host_colours()`), and a key is drawn once per figure for whatever is shared: a panel gets a
  key of its own only for a variable no other panel shows (cultivar in Figure 2).

## 2026-10-13 — panel letters are capitals, on every panel; related supplementary figures are joined

- **What happened:** asked for "nature type alignment and style", I set the panel letters in
  lower case and lettered only the top row of the PCoA grid. The user: "make capital all panel
  labels" and "make sure all pcoas have a panel label and update the caption". They also asked to
  join supplementary figures "that carry connecting messages".
- **Rule:** panel letters are capitals (config `figures.panel_letter_case`, applied by
  `panel_letter()`), and every panel of a grid has its own letter, cited in the legend. In the
  supplementary set, figures with one message are joined into one figure (config
  `supplementary_figures`), within one page.

## 2026-10-13 — alpha-diversity panels belong in one main figure

- **What happened:** alpha diversity was split over main Figures 3 (host phylogeny) and 4 (stage
  effects). The user: "make sure the figures in main with alpha diversity are kept together".
- **Rule:** group main-figure panels by what is measured (alpha diversity; community structure;
  taxa) before grouping them by the factor tested. Figure 3 is alpha diversity; Figure 4 is
  phylosymbiosis and the stage effect on genera.

## 2026-10-13 — merged figures: no loose space, every ordination with its test

- **What happened:** the merged figures kept 0.2 in of white between stacked panels and a white
  border, and the study PCoAs carried no statistics. The user: "reduce the distance between upper
  and lower panels", "crop the spaces around the figures, produce a ready to insert figure",
  "all pcoa should have stats".
- **Rule:** a merged figure has at most about 0.1 in of white between panels and no border
  (config `trim_margin_in`, `max_gap_in`, `crop_border_in`). Every ordination panel shows the
  PERMANOVA R2 and p of the variable it is coloured by, computed in Part 1 for that subset.
- **Cause worth remembering:** `patchwork::wrap_elements(full = ...)` keeps a margin of its own
  that `& theme(plot.margin = ...)` does not remove; nest the patchworks instead.

## 2026-10-20 — host levels are cultivar, host genus and plant family; no dispersion tests

- **What happened:** the grouping tests used host species beside cultivar and plant family. The
  user: "in analysis we don't have host species, we have cultivar, genus and plant family". Asked
  about PERMDISP F values, they then asked to "remove all dispersion associated results".
- **Rule:** in grouping tests (alpha-diversity tests, PERMANOVA, marginal PERMANOVA, variance
  partitioning) the host is taken as cultivar, host genus or plant family, never as species; host
  species remains only where the method needs it (the species time tree for phylogenetic signal
  and phylosymbiosis, and host breadth counted in species). Dispersion (PERMDISP) is not tested
  or reported. When adding a derived column such as host genus to a sample table, add it before
  any copy of that table is taken (a late addition was silently skipped once).

## 2026-10-22 — host attributes are plant family and host genus

- **What happened:** the user asked to "redo analyses, limiting host attributes to family and genus".
  Cultivar was tested and shown in many places (alpha tests, PERMANOVAs, stratified stage effects,
  PCoA panels, dendrogram tips). Agreed scope: remove cultivar everywhere; keep host species only
  where the method needs it (species time tree for phylogenetic signal and phylosymbiosis; host
  breadth counted in host species).
- **Rule:** stratify and group by plant family and host genus (config `stage_effects.strata`);
  do not reintroduce cultivar or host species as a tested or plotted host variable.

## 2026-10-22 — PCoA panels square and unstretched
- **Correction:** "make the pcoas perfect squares; please. no stretching.." PCoA panels had been
  drawn to fill a wide slot, so one axis was stretched.
- **Rule:** every ordination is a square panel with equal units on both axes (`coord_equal`, both
  axes spanning the wider of the stage's two ranges; chunk `fx-fig-pcoa-by-stage`). If squares do
  not fit the page, change the layout (e.g. one row of four), never the aspect.

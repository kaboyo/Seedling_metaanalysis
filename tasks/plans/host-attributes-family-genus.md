# ExecPlan: host attributes limited to plant family and host genus (2026-10-22)

## Purpose
The user asked to redo the analyses with host described only by plant family and host genus,
and to update all figures. Scope agreed: cultivar is removed everywhere it is tested or shown;
host species stays only where the method needs it (species time tree for phylogenetic signal
and phylosymbiosis; host breadth counted in host species).

## Changes
- Config (both files): `stage_effects.strata` = Plant_Family, Host_Genus; study random effect
  for Host_Genus; nothing taken within host.
- Master, Part 1: alpha-diversity tests without cultivar; alpha and genus stage effects within
  plant family and host genus; CSS PERMANOVA (alone and marginal) without cultivar; within-stage
  PERMANOVA without cultivar; raw-count screen without cultivar.
- `host_phylosymbiosis.R`: marginal PERMANOVA and variance partitioning with host as genus or family.
- Master, Part 2: dendrogram tips = host genus x stage; PCoA panels host genus and plant family
  (Figure 2 B-E, Figure S2 by study, genus, family); stage heatmaps (Figures 3, 4) within plant
  family and within host genus; Figure S2A/B without cultivar; captions.
- Legends, methods, workbook Tables S1 and S3, CLAUDE.md.
- Superseded cultivar outputs moved to `../Seedling_metaanalysis_removed_2026-10-19/`.

## Verification
Parse check, tests, workbook audit against sources, every figure looked at.

## Progress
- [x] edits  - [x] reruns  - [x] figures  - [x] text and workbook  - [x] checks (2026-10-22)

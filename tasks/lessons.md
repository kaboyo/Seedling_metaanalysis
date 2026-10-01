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

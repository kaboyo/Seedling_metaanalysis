# ExecPlan: Integrate `artifacts/` code and prune redundant files

Status: **EXECUTED — Phase 1 and Phase 2 complete; _quarantine/ awaits your sign-off**
Date: 2026-09-12

## Scope

1. Fold the code under `artifacts/` into the main config-driven pipeline.
2. Delete redundant files.

## Hard constraint

`analysis_results/manuscript/` is **off-limits**. Nothing in this plan reads, moves, edits, or
deletes anything under that path. Two consequences, both accepted:

- The stale Office lock files `~$nuscript_revised.docx` and
  `~$imi_Seedling_Meta_Analysis_Manuscript.docx` in that folder stay, even though they are junk.
- `analysis_results/manuscript/Manuscript_revised.docx` stays as the copy of record.

---

## BLOCKER 0 — the git repository is damaged

`git fsck` reports a broken object graph:

```
broken link from  commit a4ad7fd0c452548ff1efff0e875cd13826df8856
              to  commit e485a99512278840e6013e48753a35bdbca31401
broken link from    tree 282523d4f31850ce0408162127dafaed1ef78268
              to    blob 5008ddfcf53c02e82d7eee2e57c38e5672ef89f6      (+5 more blobs)
error: refs/heads/main: invalid reflog entry e485a995...
```

`git log` dies on the third commit. **Git cannot restore anything deleted in this repo.**

Therefore every deletion in Phase 2 moves files to `_quarantine/` rather than `rm`-ing them.
Nothing is permanently removed until you confirm after inspecting the quarantine folder.

Separately worth doing (not in this plan): re-clone from
`https://github.com/kaboyo/Seedling_metaanalysis.git` into a fresh directory and compare, or
`git push` what is recoverable before any further history operations.

---

## Phase 1 — Integrate `artifacts/`

`artifacts/` is **untracked** (`?? artifacts/`). It holds the output of a manuscript audit run in a
different environment (Claude Science project + a `~/Downloads` working copy), so neither of its two
code files runs from this repo as written.

### What is broken about it today

| File | Lines | Why it can't run here |
| --- | --- | --- |
| `artifacts/manuscript_audit_2026-08-04/verify_manuscript.py` | 558 | Lines 27–29 call `host.artifact_path("cade576b-…")` — a sandbox host API that does not exist in this repo. Checkpoints to a relative `work/` dir. |
| `artifacts/figures.qmd` | 625 | Line 49 hardcodes `ART <- "/Users/expeditoolimi/Downloads/artifacts"`. |

### Design decision — what "one code" means here, and what it must NOT mean

`verify_manuscript.py`'s own docstring states:

> *"Deliberately NOT a re-run of phyloseq_analysis.qmd: the statistics are reimplemented from the
> raw inputs in numpy/scipy so that agreement is evidence rather than tautology."*

Its value as an audit comes precisely from **not** sharing an implementation with the main pipeline.
If I merge its PERMANOVA/PERMDISP/Shannon code into `microbiome_analysis.py` so both call the same
helpers, the audit stops being independent verification and starts being a tautology — it would then
agree with the pipeline by construction, including where both are wrong.

So I propose integrating **the plumbing, not the statistics**:

- one repo, one `analysis_config.json`, one path convention, one output tree, one way to invoke — yes
- one shared implementation of the statistics — **no**, deliberately kept separate

If you actually want the statistics merged too, say so and I'll do it, but the audit's independence
guarantee is lost and `discrepancy_audit.csv` stops meaning what it currently means.

### Steps

1. **`verify_manuscript.py` → `validation/verify_manuscript.py`**
   - Replace `host.artifact_path(...)` with config-driven paths: read `inputs.counts`,
     `inputs.taxonomy`, `inputs.metadata` from `analysis_config.json`, resolved with the same
     base-dir rule the rest of the repo uses.
   - Add `--config analysis_config.json` argparse, matching `microbiome_analysis.py`.
   - Redirect the `work/` checkpoints to `<outputs.directory>/validation/` so smoke and full runs
     don't collide.
   - Keep every statistical function byte-for-byte as-is.

2. **`artifacts/figures.qmd` → `manuscript_figures.qmd`** (repo root, beside `phyloseq_analysis.qmd`)
   - Replace the hardcoded `ART` with Quarto `params: config: "analysis_config.json"`, matching how
     `phyloseq_analysis.qmd` already does it (`get_param()` / `load_config()` at lines 341–370).
   - Point `SUPP` and `OUT` at `<outputs.directory>/validation/`.

3. **Relocate the audit's data + results into the tracked tree**
   - `asv_sparse.npz`, `metadata_aligned.csv`, `taxonomy_aligned.csv` → `analysis_results/validation/`
     (derived checkpoints, regenerable by step 1)
   - `supplementary/` TableS1–S6 + the 3 PERMANOVA/PERMDISP CSVs → `analysis_results/validation/supplementary/`
   - `discrepancy_audit.csv`, `claim_evidence_register_verified.csv`, `revision_memo.md`,
     `README.txt` → `analysis_results/validation/`
   - `figures/` (4 PNGs, Python originals) and `figures_R/` (4 PNGs, R regeneration) →
     `analysis_results/validation/figures/` and `.../figures_R/`

4. **Update `CLAUDE.md`** to describe the validation path as a fourth execution path.

5. **Verify**: `python3 validation/verify_manuscript.py --config analysis_config.smoke.json`, then
   `quarto render manuscript_figures.qmd`. Confirm the `Validation` chunk at line 587 of the qmd
   still reconciles against the audit tables.

### Needs your decision before step 3

`artifacts/manuscript_audit_2026-08-04/Manuscript_revised.docx` (50,596 B) is **not** the same file
as `analysis_results/manuscript/Manuscript_revised.docx` (68,787 B, newer). Since the manuscript
folder is off-limits I will not touch the copy of record — but I also won't delete the audit's copy,
because it's a different document. **Default: leave it in place under quarantine, untouched.**

---

## Phase 2 — Prune redundant files

Everything below moves to `_quarantine/<original/path>` (see Blocker 0). Total ~1.87 GB.

### 2a. Byte-identical duplicates — md5-verified

| File | Size | Verified identical to |
| --- | --- | --- |
| `Datesets/asv-counts (1).csv` | 638,480,121 B | `Datesets/asv-counts.csv` — md5 `19e78715…` |
| `Datesets/asv-counts (1).pre_harmonize.csv` | 683,949,008 B | `…pre_harmonize.csv` — md5 `82648eb9…` |
| `Datesets/taxa (1).csv` | 14,830,215 B | `Datesets/taxa.csv` — md5 `bbed1aa0…` |
| `Datesets/taxonomy (1).tsv` | 8,948,523 B | `Datesets/taxonomy.tsv` — md5 `265bce21…` |
| `joined-scripts/project/strawberry-excalibur 2.qmd` | 33,517 B | `strawberry-excalibur_EXC.qmd` — md5 `9086f8f6…` |

Subtotal **~1.34 GB**. All five are untracked, so no git operation is involved.

### 2b. 38 duplicate output files — `* (1).*` under `analysis_results/` and `analysis_results_smoke/`

494.9 MB. Three pairs spot-checked byte-identical (`cooccurrence_network`,
`geographic_distribution_map`, `pcoa_jaccard_genus`). **I will `cmp` all 38 against their originals
first and quarantine only the ones that match**; any that differ get listed for you instead.
None are tracked in git.

### 2c. Regenerable caches and intermediates

| Path | Size | Note |
| --- | --- | --- |
| `.mypy_cache/` | 40 MB | regenerable |
| `tmp_debug/` | 89 MB | scratch |
| `tmp_debug_narrow/` | 8 MB | scratch |
| `__pycache__/`, `tests/__pycache__/` | 332 KB | `__pycache__/microbiome_analysis.cpython-314.pyc` **is tracked** — needs `git rm --cached` |
| `.pytest_cache/` | 16 KB | regenerable |
| `phyloseq_analysis 2.rmarkdown` | 182 KB | older Quarto knitr intermediate of `phyloseq_analysis.qmd` (184,054 B) |
| `phyloseq_analysis 3.rmarkdown` | 183 KB | same |
| `analysis_results/~$Figures.pptx` | 162 B | stale Office lock file (outside the manuscript folder) |
| `wild-again/`, `structures/` | 0 B | empty directories |

### 2d. Also worth adding

Extend `.gitignore` — it currently covers only `.env`, `node_modules/`, `dist/`, `.DS_Store`:
`__pycache__/`, `.mypy_cache/`, `.pytest_cache/`, `.venv/`, `.RData`, `.Rhistory`, `*_files/`,
`tmp_debug*/`, `~$*`, `* (1).*`.

---

## NOT deleting — flagged for your call

| File | Why I'm not touching it |
| --- | --- |
| `joined-scripts/diversity-analysis/diversity-analysis 2.qmd` (62,866 B) | Despite the `2` name it is **not** a copy — 900-line diff vs `_qwen.qmd`, and it's 20× larger than the 108-line canonical `diversity-analysis.qmd`. A distinct variant. |
| `seedling_meta_analysis.rmarkdown` | Looks like a regenerable intermediate but is **tracked in git** — deleting it is a history change, and history here is already damaged. |
| `.RData` (9.2 MB), `.Rhistory` (20 KB) | Live R session state; may hold unsaved work. |
| `Datesets/feature-table.tsv` (1.36 GB), `feature-table.biom` (8.2 MB) | Probably superseded by `asv-counts.csv`, but I have not proven redundancy. Biggest single win if you confirm. |
| `silva-138-99-seqs.qza` (97 MB), `rep-seqs.qza`, `table.qza`, `Datesets/tree.qza` | Upstream QIIME2 artifacts, tracked in git. Reference DB is re-downloadable but that's your call. |
| `presentation_files/` (2.4 MB), `phyloseq_analysis_files/` (876 KB) | Quarto render support dirs — deleting breaks the committed HTML. |

---

## Order of execution

1. Phase 2 quarantine (reversible, frees ~1.87 GB, unblocks everything else)
2. You inspect `_quarantine/`, confirm, then I `rm -rf` it
3. Phase 1 integration
4. Verify: `python3 -m pytest tests/ -q` (19 tests) + smoke run of both new entry points

## Decisions needed from you

1. Statistics kept independent (recommended) or genuinely merged into shared helpers?
2. Quarantine-then-confirm, or delete directly?
3. Anything in the "NOT deleting" table you want included — particularly `feature-table.tsv` (1.36 GB)?

---

# Progress

## 2026-09-12 — Phase 2 (prune) COMPLETE, quarantined not deleted

- Verified all 38 `* (1).*` output files with `cmp`: **38 identical, 0 differ, 0 missing original**
  (log: `_quarantine/VERIFY.txt`).
- Moved 3.1 GB to `_quarantine/files/`, paths preserved, 57 entries logged in
  `_quarantine/MANIFEST.txt`. Restore with `rsync -a _quarantine/files/ ./`.
- `analysis_results/manuscript/` confirmed untouched: 11 files before, 11 after.
- `structures/` was genuinely empty — `rmdir`'d.
- **Correction to the plan:** `wild-again/` was *not* an empty directory as stated. It held an
  aborted `git init` skeleton (`.git/objects`, `.git/refs`, zero regular files). Quarantined rather
  than deleted.
- `Datesets/feature-table.tsv` (1.36 GB) quarantined per your approval. Provenance verified first:
  `convert_seedling_tables.py` reads it to produce `asv-counts.csv`, and it is re-exportable from the
  retained `feature-table.biom` (8.2 MB).
- Tests after prune: **19 passed**.

## 2026-09-12 — Phase 1 (integrate) COMPLETE

- `validation/verify_manuscript.py` — `host.artifact_path()` replaced with config-driven resolution
  (`configure()` reading `inputs.*` and `outputs.directory`), `work/` redirected to
  `<outputs.directory>/validation/`, argparse `--config` added to match `microbiome_analysis.py`.
  Diff vs original is 104 lines confined to docstring, header, 8 path literals, a 2-line guard and
  the entry point. **Every statistical function is byte-identical** — verified by diff.
- `manuscript_figures.qmd` — created from `artifacts/figures.qmd`; the hardcoded
  `ART <- "/Users/expeditoolimi/Downloads/artifacts"` replaced with `params: config:` plumbing
  matching `phyloseq_analysis.qmd`. 44-line diff, confined to YAML + the 5 path lines.
- Audit data and results relocated to `analysis_results/validation/` (25 files).
- `artifacts/` no longer exists; its 4 superseded files are in quarantine.

### Bug found and fixed during verification

`run()`'s checkpoint branch tested only for `asv_sparse.npz` but then loaded `feature_ids.npy` and
`sample_ids.npy`. The artifact drop contained only the first, so the integrated script would have
crashed on resume. The guard now requires all three members and reparses otherwise, reporting which
were missing.

### Verification

- `quarto render manuscript_figures.qmd` → exit 0, `manuscript_figures.html` produced, all four PNGs
  regenerated into `analysis_results/validation/figures_R/`.
- The notebook's `validate` chunk reconciles against the audit tables:

  | quantity | recomputed | audit | agrees |
  | --- | --- | --- | --- |
  | samples retained | 1886 | 1886 | TRUE |
  | features retained | 2000 | 2000 | TRUE |
  | reads retained | 18,247,170 | 18,247,170 | TRUE |
  | Shannon n (stage test) | 1883 | 1883 | TRUE |
  | Proteobacteria prevalence (n) | 1770 | 1770 | TRUE |
  | Pantoea total reads | 4,468,686 | 4,468,686 | TRUE |

  Kruskal-Wallis H: 36.354 recomputed (R) vs 36.357 in TableS1 (Python) — rounding.
- `python3 -m pytest tests/ -q` → **19 passed**.
- `validation/verify_manuscript.py --config analysis_config.json`: config resolution, work-dir
  creation and the checkpoint guard all verified. **The full recomputation was not executed** — the
  638 MB reparse is ~5 min plus 999-permutation PERMANOVA, and `statsmodels` is absent so
  `TableS2_alpha_depth_models.csv` would be skipped. Existing audit tables were left as the
  tables of record.

## Still open

1. `_quarantine/` (3.1 GB) awaits your inspection, then `rm -rf`.
2. `statsmodels` not installed — needed only for `TableS2_alpha_depth_models.csv`. Not installed
   without approval (Ask First: new dependencies).
3. Git object database still damaged — unresolved, and unrelated to this work.
4. Untouched by choice, from the earlier flagged list: `diversity-analysis 2.qmd`,
   `seedling_meta_analysis.rmarkdown` (tracked), `.RData`/`.Rhistory`, the `.qza` artifacts.
5. `analysis_results/validation/Manuscript_revised.docx` (50,596 B) is the audit's own output and is
   **not** the same file as `analysis_results/manuscript/Manuscript_revised.docx` (68,787 B, newer).
   The manuscript folder's copy remains the copy of record; neither was edited.

## 2026-09-12 — Phase 3: single joint figures notebook

Added [figures_all.qmd](../../figures_all.qmd) — one document generating **all 42 analysis figures**.

**Design:** reads persisted pipeline outputs only, never the 609 MB count table. Inputs are
`setup/ps_clean.rds` (429 KB, 1,886 x 2,000), `cache/*_distance.rds`, `diversity/pcoa_*.csv`,
`diversity/umap_data.csv`, `diversity/rarefaction_data.csv`, `setup/sample_statistics.csv` and
`tables/*.csv`. Renders in about a minute. Computes no new statistics.

Each figure is wrapped in `fig_section(name, needs = <paths>)`: a missing input records a skip and the
render continues, then `tables/figure_manifest.csv` accounts for all 44 attempts.

**Result: 42 written, 0 failed, 2 skipped.**

Alpha diversity had to be recomputed from `ps_clean.rds` via `estimate_richness()` because
`tables/alpha_diversity.csv` carries no sample-id column and so cannot be joined to metadata.

### Three bugs fixed during the render loop

1. `as.data.frame(sample_data(ps))` keeps the phyloseq class; dplyr verbs reject it
   ("no applicable method for 'mutate'"). Now `as(sample_data(ps), "data.frame")` then `as_tibble()`.
   Cost 2 figures on the first render.
2. `kruskal.test(Observed ~ .data[[PRIMARY]])` — `.data` is not valid inside a formula. Now
   `reformulate(PRIMARY, "Observed")`.
3. The by-variable PCoA exports store the grouping column generically as `group`, not under the
   variable's own name, while the overall exports use real names. `pcoa_plot()` now accepts either.
   Cost 4 figures on the first render.

### Figure-folder redundancies solved

- `figures/umap_data.csv` -> `diversity/umap_data.csv` (data, not a figure; now the canonical path)
- `figures/pairwise_observed.rds`, `figures/pairwise_shannon.rds` quarantined — data objects misfiled
  in a figures directory; the `diversity/` copies are canonical
- 8 figures that `phyloseq_analysis.qmd` scattered into `diversity/` (7 PCoA) and `setup/`
  (`library_sizes.tiff`) quarantined; `figures_all.qmd` writes all of them to `figures/`
- Checked and **not** deleted: no two figures are byte-identical; all 7 `pcoa_*.csv` have distinct
  md5s; `alpha_richness_by_species` (faceted jitter) is a genuinely different figure from
  `alpha_observed_by_scientific_name` (boxplot), despite the overlapping names

`analysis_results/figures/` now holds all 42 `.tiff` and is the only figure directory. The 8 manuscript
audit `.png` files stay separate under `validation/` by design.

### The 2 legitimate skips

- `pcoa_bray_genus_by_timepoint` — `analysis_config.json` lists `Timepoint` in
  `analysis.secondary_groups` but `Datesets/metadata.csv` has no such column. **Config/data mismatch
  worth fixing in the config.**
- `pcoa_bray_genus_by_geography` — `Geography` is in the metadata but the pipeline never exports that
  ordination table. The fix belongs in `phyloseq_analysis.qmd`'s variable loop; a figures document
  should not invent an ordination.

## 2026-09-12 — Datesets/ restored from quarantine at user request

All 5 `Datesets/` files returned: `asv-counts (1).csv`, `asv-counts (1).pre_harmonize.csv`,
`taxa (1).csv`, `taxonomy (1).tsv`, `feature-table.tsv` (2.5 GB total). User states `Datesets/` holds
original files.

Verified after restore: `taxa (1).csv` and `taxonomy (1).tsv` are byte-identical to their non-`(1)`
counterparts (`cmp`), and the two 600 MB+ pairs were md5-verified before quarantine
(`19e78715`, `82648eb9`) — so the restored copies are intact, not truncated by the move.

**`Datesets/` is now off-limits for pruning**, alongside `analysis_results/manuscript/`.

Quarantine after restore: 623 files, 721 MB.

## 2026-09-12 — Quarantine deleted (721 MB, 623 files) after verifying script integration

User authorised deletion conditional on the `.qmd` scripts being integrated. Verified each of the
**5 code files** in quarantine before deleting — none was taken on trust:

| Quarantined | Integrated as | Evidence |
| --- | --- | --- |
| `artifacts/figures.qmd` | `manuscript_figures.qmd` | `diff` — the only lines present solely in the old copy are the 5 hardcoded `ART <- "~/Downloads/artifacts"` path lines, deliberately replaced |
| `artifacts/.../verify_manuscript.py` | `validation/verify_manuscript.py` | every `def` carried over (`comm` on sorted function lists returned empty) |
| `joined-scripts/project/strawberry-excalibur 2.qmd` | `strawberry-excalibur_EXC.qmd` (in repo) | `cmp` byte-identical |
| `phyloseq_analysis 2.rmarkdown` | `phyloseq_analysis.qmd` | every chunk label still present in the current `.qmd` |
| `phyloseq_analysis 3.rmarkdown` | `phyloseq_analysis.qmd` | every chunk label still present in the current `.qmd` |

Also checked before deleting, having previously recorded them as differing:
`figures/pairwise_{observed,shannon}.rds` vs the `diversity/` copies are the **same statistical
result** — identical p-values (0.00029 and 6.3e-10), same Wilcoxon method and Bonferroni adjustment,
differing only in the `data.name` label string. Nothing of substance lost.

`tmp_debug*/` held only derived subsets of `Datesets/` inputs (`asv-counts.small.csv`,
`taxa.csv`, `tree.nexus`, narrow configs); every source file is present in `Datesets/`.

**Audit trail preserved** at `tasks/plans/pruned-2026-09-12/`: `MANIFEST.txt` (every file + reason),
`VERIFY.txt` (the 38-file byte-comparison log), `deleted-file-list.txt` (all 623 paths),
`gitignore.before`.

Post-delete verification: 19 tests pass, 42 figures present in `analysis_results/figures/`,
`Datesets/` intact at 13 files, `analysis_results/manuscript/` intact at 11 files, all 7 integrated
scripts present.

### Caveat stated to the user

`manuscript_figures.qmd` is a **separate document** from `figures_all.qmd`, not folded into it. The
audit's Fig 1–4 code is preserved and runnable, but it is not part of the master figures notebook —
that separation is deliberate (see the independence note in `CLAUDE.md`).

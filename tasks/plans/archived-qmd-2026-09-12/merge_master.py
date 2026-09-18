"""Merge the three Quarto notebooks into one master analysis document.

Order matters: phyloseq_analysis (pipeline) -> figures_all (figures) ->
manuscript_figures (audit figures). R's later-definition-wins semantics mean the
figures part's helper definitions are the ones active while its own figures are
drawn, which is what each part was written against.
"""

import re
from pathlib import Path

ROOT = Path("/Users/expeditoolimi/Documents/Southampton/Southampton/Seedling_metaanalysis")
OUT = ROOT / "seedling_analysis_master.qmd"

SOURCES = [
    ("phyloseq_analysis.qmd", "pipeline", None),
    ("figures_all.qmd", "figures", "fx"),
    ("manuscript_figures.qmd", "manuscript", "ms"),
]


def split_yaml(text):
    """Return (yaml_block, body) for a Quarto file starting with ---."""
    if not text.startswith("---"):
        return "", text
    end = text.index("\n---", 3)
    return text[: end + 4], text[end + 4 :].lstrip("\n")


def relabel(body, prefix):
    """Prefix every chunk label so labels stay unique across the merged file."""
    if prefix is None:
        return body

    def sub(m):
        fence, lang, label, rest = m.group(1), m.group(2), m.group(3), m.group(4)
        return f"{fence}{{{lang} {prefix}-{label}{rest}}}"

    return re.sub(r"^(```)\{(r|python)[ ,]+([A-Za-z0-9_.-]+)([^}]*)\}",
                  sub, body, flags=re.M)


YAML = '''---
title: "Seedling Seed Microbiome Meta-analysis — Master Analysis"
subtitle: "Complete workflow: pipeline, all figures, and the manuscript audit figures"
date: today
format:
  html:
    code-fold: true
    toc: true
    toc-depth: 3
    theme: cosmo
    embed-resources: true
    fig-align: center
execute:
  echo: true
  warning: false
  message: false
knitr:
  opts_chunk:
    dev: png
    dpi: 150
params:
  config: "analysis_config.json"
  filtered_features: NULL
  run_validation: false
---

::: {.callout-note}
## One document, three parts

This is the single entry point for the whole analysis. It was assembled from three
notebooks that previously ran separately:

| Part | Source | What it does |
| --- | --- | --- |
| 1 | `phyloseq_analysis.qmd` | the pipeline — load, harmonise, filter, diversity, composition, associations |
| 2 | `figures_all.qmd` | every analysis figure, drawn from the outputs Part 1 persists |
| 3 | `manuscript_figures.qmd` | the 2026-08-04 manuscript audit's corrected Fig 1–4 |

Render it whole:

```
quarto render seedling_analysis_master.qmd -P config:analysis_config.json
```

**Part 1 is the expensive step** — it reads the 609 MB count table. For a fast pass,
point `-P config:analysis_config.smoke.json` at the 120-sample smoke config.

Part 3 needs the audit tables under `<outputs.directory>/validation/`. Set
`-P run_validation:true` to have this document regenerate them by invoking
`validation/verify_manuscript.py` first; otherwise it reads whatever is already
there and skips cleanly if nothing is.

Helper functions are intentionally defined more than once across the parts. R takes
the most recent definition, so each part runs against the helpers it was written
for — do not "tidy" these into a single definition without re-checking every figure.
:::

'''


def part_header(n, title, note):
    return f"\n\n{'#' * 1} Part {n} — {title}\n\n::: {{.callout-tip collapse=\"true\"}}\n## About this part\n{note}\n:::\n\n"


VALIDATION_CHUNK = '''

# Part 3 — Manuscript audit figures

::: {.callout-tip collapse="true"}
## About this part
Fig 1–4 of the manuscript, as corrected by the 2026-08-04 audit. Fig 1 and Fig 2a are
recomputed from the aligned raw inputs; Fig 2b/2c, 3 and 4 are read from the audit's
tables of record so the panels show published values rather than a fresh permutation
run. The `Validation` section at the end reconciles every recomputed quantity.

The recomputation behind those tables is `validation/verify_manuscript.py`, which
reimplements the statistics independently of Part 1 — that is deliberate, so agreement
between them is evidence rather than tautology. Do not refactor the two to share code.
:::

```{r ms-run-validation}
# Optionally regenerate the audit tables before drawing Fig 1-4. Guarded the same way
# the Python pipeline guards its R step: if the interpreter or script is missing, the
# existing tables are used and the render continues.
run_validation <- isTRUE(as.logical(if (exists("params")) params$run_validation else FALSE))
validation_script <- file.path(base_dir, "validation", "verify_manuscript.py")

if (run_validation) {
    python_bin <- Sys.which("python3")
    if (nzchar(python_bin) && file.exists(validation_script)) {
        message("running ", validation_script)
        status <- system2(python_bin, c(shQuote(validation_script),
                                        "--config", shQuote(config_path)))
        message("verify_manuscript.py exit status: ", status)
    } else {
        message("skipping validation run - python3 or verify_manuscript.py not found")
    }
} else {
    message("run_validation is FALSE - using the audit tables already on disk")
}
```

'''


def main():
    chunks_out = [YAML]

    for idx, (name, kind, prefix) in enumerate(SOURCES, start=1):
        raw = (ROOT / name).read_text()
        _, body = split_yaml(raw)
        body = relabel(body, prefix)

        if kind == "pipeline":
            chunks_out.append(part_header(
                1, "Analysis pipeline",
                "The full workflow from `phyloseq_analysis.qmd`: dataset audit and "
                "harmonisation, phyloseq construction, taxonomy filtering, decontamination, "
                "diversity, composition, core/flexible partitioning and association tests. "
                "This part writes the outputs Part 2 reads."))
            chunks_out.append(body)

        elif kind == "figures":
            chunks_out.append(part_header(
                2, "All analysis figures",
                "Every analysis figure, from `figures_all.qmd`. It reads the persisted "
                "outputs Part 1 just wrote rather than recomputing them, and wraps each "
                "figure in `fig_section()` so a missing input records a skip instead of "
                "aborting the render. Writes `tables/figure_manifest.csv`."))
            chunks_out.append(body)

        else:
            # Part 3 was written as a standalone document, so a missing validation
            # directory was a hard stop(). Inside the master that would abort the whole
            # render - fatal for a smoke run, where <outputs.directory>/validation does
            # not exist. Turn the stop into a flag and make every later ms- chunk
            # conditional on it, matching how Part 2 skips.
            body = body.replace(
                """for (needed in c(DATA, SUPP)) {
  if (!dir.exists(needed)) {
    stop("Missing input directory: ", needed,
         "\\nRun: python3 validation/verify_manuscript.py --config ", config_path)
  }
}""",
                """MS_READY <- all(dir.exists(c(DATA, SUPP)))
if (!MS_READY) {
  message("Part 3 skipped - missing ", paste(c(DATA, SUPP)[!dir.exists(c(DATA, SUPP))],
                                             collapse = ", "),
          "\\n  populate it with: python3 validation/verify_manuscript.py --config ",
          config_path, "\\n  or re-render with -P run_validation:true")
}""")
            assert "MS_READY <- all(" in body, "Part 3 guard not rewritten"

            # Gate every ms- chunk after the setup chunk on MS_READY.
            def gate(m):
                label, rest = m.group(1), m.group(2)
                if label in ("ms-setup", "ms-run-validation"):
                    return m.group(0)
                return f"```{{r {label}{rest}, eval = MS_READY}}"

            body = re.sub(r"^```\{r (ms-[A-Za-z0-9_.-]+)([^}]*)\}",
                          gate, body, flags=re.M)
            chunks_out.append(VALIDATION_CHUNK)
            chunks_out.append(body)

    text = "\n".join(chunks_out)

    # The three sources each ended with their own session-info chunk; keep only the last.
    text = re.sub(r"\n```\{r (?:fx-|ms-)?session[^}]*\}\nsessionInfo\(\)\n```\n",
                  "\n", text)
    text += """

# Session information

```{r master-session}
sessionInfo()
```
"""

    OUT.write_text(text)

    labels = re.findall(r"^```\{(?:r|python)[ ,]+([A-Za-z0-9_.-]+)", text, flags=re.M)
    dupes = {l for l in labels if labels.count(l) > 1}
    print(f"wrote {OUT.name}")
    print(f"  lines: {text.count(chr(10)) + 1}")
    print(f"  chunks: {len(labels)}")
    print(f"  duplicate labels: {sorted(dupes) if dupes else 'none'}")


if __name__ == "__main__":
    main()

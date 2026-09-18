"""Append every remaining .qmd into the master as non-executing reference parts.

The 23 archived notebooks cannot run in this repo: the joined-scripts set targets the
strawberry/EXC projects (`Data/` paths, `alpha_df`, `all_dunn_results`, `TreatmentCode.x`),
presentation.qmd is revealjs, and seedling_meta_analysis.qmd is an earlier duplicate of
Part 1. Their code is therefore carried with `eval = FALSE`: the content lives in the
master and renders as visible code, but it cannot break the working analysis.
"""

import re
from pathlib import Path

ROOT = Path("/Users/expeditoolimi/Documents/Southampton/Southampton/Seedling_metaanalysis")
SRC = ROOT / "tasks/plans/archived-qmd-2026-09-12/sources"
MASTER = ROOT / "seedling_analysis_master.qmd"
STAGED = ROOT / "seedling_analysis_master.staged.qmd"

# (relative path under SRC, short unique chunk prefix, part grouping)
FILES = [
    ("seedling_meta_analysis.qmd", "sma", "Part 4 — Earlier analysis version"),

    ("joined-scripts/project/strawberry-excalibur.qmd", "pj1", "Part 5 — Archived project notebooks"),
    ("joined-scripts/project/strawberry-excalibur_EXC.qmd", "pj2", "Part 5 — Archived project notebooks"),
    ("joined-scripts/project/strawberry-excalibur_phyloseq-legacy.qmd", "pj3", "Part 5 — Archived project notebooks"),

    ("joined-scripts/beta-diversity/Beta-diversity.qmd", "bd1", "Part 6 — Archived beta-diversity notebooks"),
    ("joined-scripts/beta-diversity/Beta-diversity_EXC.qmd", "bd2", "Part 6 — Archived beta-diversity notebooks"),
    ("joined-scripts/beta-diversity/Beta-diversity_functions.qmd", "bd3", "Part 6 — Archived beta-diversity notebooks"),
    ("joined-scripts/beta-diversity/Beta-diversity_functions_data_variant.qmd", "bd4", "Part 6 — Archived beta-diversity notebooks"),

    ("joined-scripts/diversity-analysis/diversity-analysis.qmd", "da1", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis 2.qmd", "da2", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis_core.qmd", "da3", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis_jaccard.qmd", "da4", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis_modified.qmd", "da5", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis_qwen.qmd", "da6", "Part 7 — Archived diversity-analysis notebooks"),
    ("joined-scripts/diversity-analysis/diversity-analysis_strawberry.qmd", "da7", "Part 7 — Archived diversity-analysis notebooks"),

    ("joined-scripts/reports/Jorn_class_2_2025.qmd", "rp1", "Part 8 — Archived report sections"),
    ("joined-scripts/reports/colour-pellet.qmd", "rp2", "Part 8 — Archived report sections"),
    ("joined-scripts/reports/manuscript.qmd", "rp3", "Part 8 — Archived report sections"),
    ("joined-scripts/reports/plots-section.qmd", "rp4", "Part 8 — Archived report sections"),
    ("joined-scripts/reports/statistics-section.qmd", "rp5", "Part 8 — Archived report sections"),

    ("rphylopic.qmd", "aux1", "Part 9 — Auxiliary documents"),
    ("shinyphaser.qmd", "aux2", "Part 9 — Auxiliary documents"),
    ("presentation.qmd", "aux3", "Part 9 — Auxiliary documents"),
]

PART_NOTES = {
    "Part 4 — Earlier analysis version":
        "An earlier, simpler implementation of the same workflow as Part 1 (metadata summary, "
        "sequencing depth, top phyla/genera, alpha, beta, PERMANOVA) against the same "
        "`Datesets/` inputs. Superseded by Part 1 and carried here for provenance.",
    "Part 5 — Archived project notebooks":
        "End-to-end notebooks from the strawberry/EXC projects.",
    "Part 6 — Archived beta-diversity notebooks":
        "Beta-diversity variants from the same archived projects.",
    "Part 7 — Archived diversity-analysis notebooks":
        "Seven diversity-analysis variants, largely alternatives to one another.",
    "Part 8 — Archived report sections":
        "Report and manuscript fragments; several contain prose only.",
    "Part 9 — Auxiliary documents":
        "Library tutorials and the reveal.js slide deck. Not analysis code.",
}

ARCHIVE_BANNER = """

::: {.callout-warning}
## Parts 4–9 are reference only — they do not execute

Every chunk below carries `eval = FALSE`. The code is here, in full and unmodified, so the
master document holds the whole project's content — but none of it runs, deliberately:

- The `joined-scripts` notebooks target the **strawberry / EXC projects**, not this dataset.
  They reference `Data/` paths, `alpha_df`, `all_dunn_results`, `TreatmentCode.x` and
  `attach()` on objects that do not exist here. Executing them fails immediately.
- Several are variants of one another (seven `diversity-analysis` files, four
  `Beta-diversity` files), so running them would also contradict each other.
- `seedling_meta_analysis.qmd` (Part 4) duplicates Part 1's analyses in an earlier form.
- Part 9 holds library tutorials and a reveal.js deck — a different output format entirely.

Parts 1–3 are the live analysis. To revive anything below, lift the chunk out and adapt it
to this repo's schema (`run_accession`, `Developmental_Stage`, `Datesets/`) first.
:::
"""


def split_yaml(text):
    if not text.startswith("---"):
        return "", text
    try:
        end = text.index("\n---", 3)
    except ValueError:
        return "", text
    return text[: end + 4], text[end + 4 :].lstrip("\n")


def prepare(path, prefix):
    """Relabel chunks uniquely and force eval = FALSE on every one."""
    text = path.read_text(errors="replace")
    yaml_block, body = split_yaml(text)

    # A label can repeat inside one source file (e.g. two chunks both called `echo`);
    # Quarto refuses to render on duplicates, so suffix repeats.
    seen = {}

    auto = [0]

    def sub(m):
        fence, lang, label, rest = m.group(1), m.group(2), m.group(3), m.group(4)
        rest = rest.rstrip()
        # `{r echo=FALSE}` has no label at all - the first token is an option name.
        # Prefixing it would emit `{r pj2-echo=FALSE}`, which knitr rejects.
        if rest.lstrip().startswith("="):
            auto[0] += 1
            return f"{fence}{{{lang} {prefix}-c{auto[0]}, {label}{rest}, eval = FALSE}}"
        seen[label] = seen.get(label, 0) + 1
        if seen[label] > 1:
            label = f"{label}{seen[label]}"
        return f"{fence}{{{lang} {prefix}-{label}{rest}, eval = FALSE}}"

    body = re.sub(r"^(```)\{(r|python)[ ,]+([A-Za-z0-9_.-]+)([^}]*)\}",
                  sub, body, flags=re.M)

    # Unlabelled chunks: give them a label too, still non-executing.
    counter = [0]

    def sub_unlabelled(m):
        counter[0] += 1
        return f"```{{{m.group(1)} {prefix}-u{counter[0]}, eval = FALSE}}"

    body = re.sub(r"^```\{(r|python)\s*\}", sub_unlabelled, body, flags=re.M)

    # Demote headings by two levels so they nest under the part heading - but only
    # outside fenced blocks, or every `# comment` inside an R chunk gets rewritten too.
    out_lines, in_fence = [], False
    for line in body.split("\n"):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            out_lines.append(line)
            continue
        if not in_fence:
            line = re.sub(r"^(#{1,4}) ",
                          lambda m: "#" * min(len(m.group(1)) + 2, 6) + " ", line)
        out_lines.append(line)
    body = "\n".join(out_lines)

    return yaml_block, body


def main():
    out = [MASTER.read_text(), ARCHIVE_BANNER]
    current_part = None
    n_files = n_chunks = 0

    for rel, prefix, part in FILES:
        path = SRC / rel
        if not path.exists():
            print(f"  MISSING, skipped: {rel}")
            continue

        if part != current_part:
            out.append(f"\n\n# {part}\n\n::: {{.callout-note collapse=\"true\"}}\n"
                       f"## About this part\n{PART_NOTES[part]}\n:::\n")
            current_part = part

        yaml_block, body = prepare(path, prefix)
        n_files += 1
        n_chunks += len(re.findall(r"^```\{(?:r|python) ", body, flags=re.M))

        out.append(f"\n## Source: `{rel}`\n")
        if yaml_block:
            out.append("\nOriginal front matter:\n\n``` yaml\n"
                       + yaml_block.strip("-\n") + "\n```\n")
        out.append("\n" + body.rstrip() + "\n")

    text = "\n".join(out)
    STAGED.write_text(text)

    labels = re.findall(r"^```\{(?:r|python)[ ,]+([A-Za-z0-9_.-]+)", text, flags=re.M)
    dupes = sorted({l for l in labels if labels.count(l) > 1})
    print(f"wrote {STAGED.name}")
    print(f"  files appended : {n_files}")
    print(f"  chunks appended: {n_chunks}")
    print(f"  total lines    : {text.count(chr(10)) + 1}")
    print(f"  total chunks   : {len(labels)}")
    print(f"  duplicate labels: {dupes if dupes else 'none'}")


if __name__ == "__main__":
    main()

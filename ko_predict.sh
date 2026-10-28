#!/usr/bin/env bash
# ko_predict.sh - predict KEGG orthologue (KO) content for the 16S ASVs with PICRUSt2.
#
# Usage:
#   bash ko_predict.sh [analysis_config.json]
#
# Run make_generalist_specialist.R first: it lists the ASVs. Gene content is predicted
# from the marker gene, not measured.
#
# PICRUSt2 is looked for in $PICRUST2_BIN (default: the picrust2 conda environment under
# ~/.local/share/seedling_tools/miniforge3) and then on PATH.
#
# Reads  <outputs.directory>/tables/generalist_specialist_asv_host_screen.csv and the ASV
#        sequences named in the config (inputs.rep_seqs).
# Writes <outputs.directory>/picrust2/ (PICRUSt2's own output; the per-ASV KO predictions
#        and NSTI values are read from there by ko_breadth_tests.R).
#
# PICRUSt2 wants an abundance table as well as sequences. Only the per-ASV predictions
# are used here, so each ASV is given a count of 1 in a single dummy sample.

set -euo pipefail

CONFIG="${1:-analysis_config.json}"
[ -f "$CONFIG" ] || { echo "Configuration file not found: $CONFIG" >&2; exit 1; }
PICRUST="${PICRUST2_BIN:-$HOME/.local/share/seedling_tools/miniforge3/envs/picrust2/bin}"
[ -d "$PICRUST" ] && PATH="$PICRUST:$PATH"
command -v picrust2_pipeline.py >/dev/null || { echo "PICRUSt2 not found" >&2; exit 1; }

read -r OUT_DIR REP_SEQS MAX_NSTI THREADS < <(python3 - "$CONFIG" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1]).resolve()
config = json.loads(path.read_text())
resolve = lambda p: p if p.startswith("/") else str(path.parent / p)
print(resolve(config["outputs"]["directory"]), resolve(config["inputs"]["rep_seqs"]),
      config.get("gene_content", {}).get("max_nsti", 0.15), config["analysis"].get("threads", 1))
PY
)
SCREEN="$OUT_DIR/tables/generalist_specialist_asv_host_screen.csv"
[ -f "$SCREEN" ] || { echo "Run make_generalist_specialist.R first; missing $SCREEN" >&2; exit 1; }

WORK="$OUT_DIR/picrust2_input"
mkdir -p "$WORK"
python3 - "$SCREEN" "$REP_SEQS" "$WORK" <<'PY'
import csv, re, sys, zipfile
screen, rep_seqs, work = sys.argv[1:4]
with open(screen, newline="") as handle:
    wanted = {re.sub(r"^ASV [0-9]+ : ", "", row["ASV"]): "" for row in csv.DictReader(handle)}
with zipfile.ZipFile(rep_seqs) as archive:
    fasta = next(n for n in archive.namelist() if n.endswith("data/dna-sequences.fasta"))
    name = None
    for line in archive.read(fasta).decode().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            name = name if name in wanted else None
        elif name is not None:
            wanted[name] += line.strip()
missing = [k for k, v in wanted.items() if not v]
if missing:
    sys.exit(f"{len(missing)} ASVs have no sequence, e.g. {missing[0]}")
with open(f"{work}/asvs.fasta", "w") as out:
    for name, sequence in wanted.items():
        out.write(f">{name}\n{sequence}\n")
with open(f"{work}/table.tsv", "w") as out:
    out.write("#OTU ID\tall_asvs\n")
    for name in wanted:
        out.write(f"{name}\t1\n")
print(f"{len(wanted)} ASVs written for PICRUSt2")
PY

# PICRUSt2 refuses to write into a folder that exists; this one holds only its own output.
PICRUST_OUT="$OUT_DIR/picrust2"
[ -d "$PICRUST_OUT" ] && rm -r "$PICRUST_OUT"
picrust2_pipeline.py -s "$WORK/asvs.fasta" -i "$WORK/table.tsv" -o "$PICRUST_OUT" \
    -p "$THREADS" --max_nsti "$MAX_NSTI" --in_traits KO --no_pathways --verbose
echo "PICRUSt2 $(picrust2_pipeline.py --version 2>&1 | tail -1) finished: $PICRUST_OUT"
ls "$PICRUST_OUT"

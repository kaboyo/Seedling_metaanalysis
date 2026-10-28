#!/usr/bin/env bash
# generalist_specialist_trees.sh - maximum-likelihood trees of the host generalist and the
# host specialist ASVs (16S), built as the methods text describes.
#
# Usage:
#   bash generalist_specialist_trees.sh [analysis_config.json] ["generalists specialists combined"]
#
# Run make_generalist_specialist.R first: it writes the ASV classes read here. The second
# argument names the trees to build; all three by default. IQ-TREE takes about half an hour,
# the others a minute.
#
#   Generalists: MUSCLE -> trimAl (-automated1) -> IQ-TREE (ModelFinder, ultrafast bootstrap)
#   Specialists: MAFFT (--auto) -> trimAl (-automated1) -> FastTree (-nt -gtr -gamma)
#   Combined:    both classes on one tree, built as the specialists are
#
# The programs are looked for in $SEEDLING_TOOLS_BIN (default
# ~/.local/share/seedling_tools/bin, see the README there) and then on PATH.
#
# Reads  <outputs.directory>/tables/generalist_specialist_asv_host_screen.csv and the ASV
#        sequences named in the config (inputs.rep_seqs, a QIIME2 .qza).
# Writes <outputs.directory>/trees/{generalists,specialists,combined}/: sequences, alignment,
#        trimmed alignment, tree and logs; and trees/tree_build_summary.csv.
#
# The ASVs come from studies that sequenced different regions of the 16S gene. An alignment
# of fragments that do not overlap gives unreliable deep branches, so the script reports
# how much of each alignment survives trimming and stops if too little does.

set -euo pipefail

CONFIG="${1:-analysis_config.json}"
TREE_SETS="${2:-generalists specialists combined}"
want() { case " $TREE_SETS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
[ -f "$CONFIG" ] || { echo "Configuration file not found: $CONFIG" >&2; exit 1; }
TOOLS="${SEEDLING_TOOLS_BIN:-$HOME/.local/share/seedling_tools/bin}"
[ -d "$TOOLS" ] && PATH="$TOOLS:$PATH"
for tool in muscle mafft trimal iqtree3 FastTree; do
    command -v "$tool" >/dev/null || { echo "missing program: $tool" >&2; exit 1; }
done

# Paths and settings come from the config; paths are relative to the config's directory.
read -r OUT_DIR REP_SEQS BOOTSTRAPS MIN_COLUMNS THREADS < <(python3 - "$CONFIG" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1]).resolve()
config = json.loads(path.read_text())
resolve = lambda p: p if p.startswith("/") else str(path.parent / p)
gs = config.get("generalist_specialist", {})
print(resolve(config["outputs"]["directory"]), resolve(config["inputs"]["rep_seqs"]),
      gs.get("tree_bootstraps", 1000), gs.get("tree_min_alignment_columns", 100),
      config["analysis"].get("threads", 1))
PY
)
SCREEN="$OUT_DIR/tables/generalist_specialist_asv_host_screen.csv"
[ -f "$SCREEN" ] || { echo "Run make_generalist_specialist.R first; missing $SCREEN" >&2; exit 1; }
[ -f "$REP_SEQS" ] || { echo "ASV sequences not found: $REP_SEQS" >&2; exit 1; }
TREES="$OUT_DIR/trees"
mkdir -p "$TREES/generalists" "$TREES/specialists" "$TREES/combined"

# One FASTA per class, named by the ASV's sequence hash (the tip labels of every tree here).
python3 - "$SCREEN" "$REP_SEQS" "$TREES" <<'PY'
import csv, re, sys, zipfile
screen, rep_seqs, trees = sys.argv[1:4]
wanted = {"Generalist": {}, "Specialist": {}}
with open(screen, newline="") as handle:
    for row in csv.DictReader(handle):
        if row["screen_class"] in wanted:
            wanted[row["screen_class"]][re.sub(r"^ASV [0-9]+ : ", "", row["ASV"])] = None
with zipfile.ZipFile(rep_seqs) as archive:
    fasta = next(n for n in archive.namelist() if n.endswith("data/dna-sequences.fasta"))
    name = None
    for line in archive.read(fasta).decode().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif name is not None:
            for group in wanted.values():
                if name in group:
                    group[name] = (group[name] or "") + line.strip()
for label, folder in (("Generalist", "generalists"), ("Specialist", "specialists")):
    missing = [k for k, v in wanted[label].items() if not v]
    if missing:
        sys.exit(f"{len(missing)} {label} ASVs have no sequence, e.g. {missing[0]}")
    with open(f"{trees}/{folder}/sequences.fasta", "w") as out:
        for name, sequence in wanted[label].items():
            out.write(f">{name}\n{sequence}\n")
    lengths = sorted(len(v) for v in wanted[label].values())
    print(f"{label}s: {len(lengths)} sequences, {lengths[0]}-{lengths[-1]} bp, "
          f"median {lengths[len(lengths) // 2]} bp")
with open(f"{trees}/combined/sequences.fasta", "w") as out:
    for group in wanted.values():
        for name, sequence in group.items():
            out.write(f">{name}\n{sequence}\n")
PY

# Columns and gap share of an aligned FASTA, and the number of sequences left with no bases.
alignment_stats() {
    python3 - "$1" <<'PY'
import sys
sequences, name = {}, None
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; sequences[name] = ""
    elif name:
        sequences[name] += line
columns = len(next(iter(sequences.values())))
gaps = sum(s.count("-") for s in sequences.values()) / (columns * len(sequences))
empty = sum(1 for s in sequences.values() if not s.replace("-", ""))
print(len(sequences), columns, f"{gaps:.3f}", empty)
PY
}
check_trimmed() {     # folder
    read -r n columns gaps empty < <(alignment_stats "$1/trimmed.fasta")
    echo "  trimmed alignment: $n sequences, $columns columns, gap share $gaps, $empty with no bases left"
    if [ "$columns" -lt "$MIN_COLUMNS" ]; then
        echo "  STOP: fewer than $MIN_COLUMNS columns survive trimming; no tree built" >&2; exit 2
    fi
}

if want generalists; then
    echo "== Generalists: MUSCLE -> trimAl -> IQ-TREE =="
    G="$TREES/generalists"
    muscle -align "$G/sequences.fasta" -output "$G/aligned.fasta" -threads "$THREADS" > "$G/muscle.log" 2>&1
    read -r n columns gaps empty < <(alignment_stats "$G/aligned.fasta")
    echo "  alignment: $n sequences, $columns columns, gap share $gaps"
    trimal -in "$G/aligned.fasta" -out "$G/trimmed.fasta" -automated1 > "$G/trimal.log" 2>&1
    check_trimmed "$G"
    iqtree3 -s "$G/trimmed.fasta" -m MFP -B "$BOOTSTRAPS" -T "$THREADS" -seed 123 \
        --prefix "$G/iqtree" -redo > "$G/iqtree.stdout.log" 2>&1
    cp "$G/iqtree.contree" "$G/tree.nwk"
    echo "  $BOOTSTRAPS ultrafast bootstraps; tree: $G/tree.nwk"
fi

# MAFFT -> trimAl -> FastTree for one folder of sequences.
fasttree_tree() {     # folder
    mafft --auto --thread "$THREADS" "$1/sequences.fasta" > "$1/aligned.fasta" 2> "$1/mafft.log"
    read -r n columns gaps empty < <(alignment_stats "$1/aligned.fasta")
    echo "  alignment: $n sequences, $columns columns, gap share $gaps"
    trimal -in "$1/aligned.fasta" -out "$1/trimmed.fasta" -automated1 > "$1/trimal.log" 2>&1
    check_trimmed "$1"
    FastTree -nt -gtr -gamma "$1/trimmed.fasta" > "$1/tree.nwk" 2> "$1/fasttree.log"
    echo "  tree: $1/tree.nwk"
}
if want specialists; then
    echo "== Specialists: MAFFT -> trimAl -> FastTree =="
    fasttree_tree "$TREES/specialists"
fi
if want combined; then
    echo "== Generalists and specialists combined: MAFFT -> trimAl -> FastTree =="
    fasttree_tree "$TREES/combined"
fi

{
    echo "group,sequences,aligned_columns,trimmed_columns,trimmed_gap_share,programs"
    for folder in generalists specialists combined; do
        [ -f "$TREES/$folder/tree.nwk" ] || continue
        read -r n a _ _ < <(alignment_stats "$TREES/$folder/aligned.fasta")
        read -r _ t gaps _ < <(alignment_stats "$TREES/$folder/trimmed.fasta")
        if [ "$folder" = generalists ]; then
            model=$(grep -m1 "Best-fit model" "$TREES/$folder/iqtree.log" | sed 's/.*: //; s/ chosen.*//')
            programs="MUSCLE $(muscle --version | awk '{print $2}' | head -1); trimAl -automated1; IQ-TREE $(iqtree3 --version | awk 'NR==1{print $3}') $model -B $BOOTSTRAPS"
        else
            programs="MAFFT $(mafft --version 2>&1 | awk '{print $1}') --auto; trimAl -automated1; FastTree 2.2.0 -nt -gtr -gamma"
        fi
        echo "$folder,$n,$a,$t,$gaps,\"$programs\""
    done
} > "$TREES/tree_build_summary.csv"
echo "Summary: $TREES/tree_build_summary.csv"

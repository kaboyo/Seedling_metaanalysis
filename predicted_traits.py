#!/usr/bin/env python3
"""Compile the ASV x predicted-trait table for the 16S ASVs.

Usage:
    python3 predicted_traits.py --config analysis_config.json

Every trait here is an in-silico prediction from the marker gene, not a measurement, and
is called a "predicted trait" throughout. Three resources are compiled:

BugBase (v0.1.0, Greengenes 13.8)
    BugBase's traits are keyed to Greengenes OTUs. Each ASV is matched to its closest
    Greengenes 97% reference sequence with vsearch (closed reference, identity from the
    config) and takes that OTU's precalculated trait values. Columns that only hold 0 and 1
    become binary traits; the fractional ones stay numeric.

FAPROTAX (v1.2.12)
    FAPROTAX assigns ecological functions from taxonomy names, so it runs on this
    project's SILVA taxonomy as it is. An ASV placed in at least one functional group is
    annotated: TRUE for its groups and FALSE for the rest. An ASV placed in none is NA.

metaTraits (metatraits.embl.de, GTDB genus-level summary)
    metaTraits summarises traits per taxon. An ASV takes the summary of the GTDB genus
    whose name equals its SILVA genus; an ASV with no genus, or whose genus name is not in
    metaTraits, is NA. Boolean traits use the consensus value ("No robust majority" is NA);
    traits with a numeric unit use the genus median. Categorical traits such as cell shape
    are neither binary nor numeric and are left out. Every ASV of a genus therefore
    carries the same values.

Reads   <outputs.directory>/tables/generalist_specialist_asv_host_screen.csv (the ASVs),
        inputs.rep_seqs (sequences), inputs.taxonomy, and the reference files under the
        resources directory (predicted_traits.resources_dir).
Writes  <outputs.directory>/tables/predicted_traits_asv.csv        one row per ASV
        <outputs.directory>/tables/predicted_traits_catalogue.csv  one row per trait
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

BINARY = "binary"
NUMERIC = "numeric"


def resolve(base_dir: Path, configured: str) -> Path:
    """Resolve a config path relative to the config file's directory."""
    path = Path(configured)
    return path if path.is_absolute() else base_dir / path


def find_program(name: str) -> str:
    """Locate an external program in the tools folder or on PATH."""
    tools = Path(os.environ.get("SEEDLING_TOOLS_BIN", Path.home() / ".local/share/seedling_tools/bin"))
    candidate = tools / name
    if candidate.exists():
        return str(candidate)
    found = shutil.which(name)
    if found is None:
        raise FileNotFoundError(f"program not found: {name}")
    return found


def asv_hash(asv: str) -> str:
    """'ASV 12 : 3e8c...' -> '3e8c...', the sequence hash used in the FASTA and the trees."""
    return re.sub(r"^ASV [0-9]+ : ", "", asv)


def read_asvs(screen_path: Path) -> dict[str, str]:
    """Sequence hash -> full ASV name, for every ASV in the host screen table."""
    with screen_path.open(newline="") as handle:
        return {asv_hash(row["ASV"]): row["ASV"] for row in csv.DictReader(handle)}


def read_sequences(rep_seqs: Path, wanted: set[str]) -> dict[str, str]:
    """Sequences of the wanted ASVs from the QIIME2 representative-sequence archive."""
    sequences: dict[str, str] = {}
    with zipfile.ZipFile(rep_seqs) as archive:
        fasta = next(n for n in archive.namelist() if n.endswith("data/dna-sequences.fasta"))
        name: str | None = None
        for line in archive.read(fasta).decode().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                name = name if name in wanted else None
            elif name is not None:
                sequences[name] = sequences.get(name, "") + line.strip()
    missing = wanted - set(sequences)
    if missing:
        raise ValueError(f"{len(missing)} ASVs have no sequence, e.g. {sorted(missing)[0]}")
    return sequences


def bugbase_traits(
    sequences: dict[str, str], resources: Path, identity: float, threads: int, work: Path
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    """BugBase trait values per ASV, through its closest Greengenes OTU.

    Returns (ASV hash -> {trait: value}, trait -> "binary" | "numeric").
    """
    reference = resources / "greengenes/gg_13_8_otus/rep_set/97_otus.fasta"
    table = resources / "bugbase/default_traits_precalculated.txt.gz"
    for path in (reference, table):
        if not path.exists():
            raise FileNotFoundError(f"BugBase resource missing: {path}")

    query = work / "asvs.fasta"
    query.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    hits = work / "greengenes_hits.tsv"
    subprocess.run(
        [find_program("vsearch"), "--usearch_global", str(query), "--db", str(reference),
         "--id", str(identity), "--strand", "both", "--top_hits_only", "--maxaccepts", "1",
         "--blast6out", str(hits), "--threads", str(threads), "--quiet"],
        check=True,
    )
    otu_of: dict[str, str] = {}
    with hits.open() as handle:
        for line in handle:
            fields = line.split("\t")
            otu_of.setdefault(fields[0], fields[1])

    wanted_otus = set(otu_of.values())
    values: dict[str, list[str]] = {}
    with gzip.open(table, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")[1:]
        column_values: list[set[str]] = [set() for _ in header]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            # the type of a column is judged on the whole BugBase table, not on our ASVs
            for seen, value in zip(column_values, fields[1:]):
                if len(seen) < 3:
                    seen.add(value)
            if fields[0] in wanted_otus:
                values[fields[0]] = fields[1:]
    kinds = {
        f"bugbase_{name}": BINARY if seen <= {"0", "1"} else NUMERIC
        for name, seen in zip(header, column_values)
    }
    traits: dict[str, dict[str, str]] = {}
    for asv, otu in otu_of.items():
        if otu not in values:
            continue
        row: dict[str, str] = {"greengenes_otu": otu}
        for name, value in zip(header, values[otu]):
            key = f"bugbase_{name}"
            row[key] = ("TRUE" if value == "1" else "FALSE") if kinds[key] == BINARY else value
        traits[asv] = row
    print(f"BugBase: {len(otu_of)} of {len(sequences)} ASVs matched a Greengenes OTU at "
          f">= {identity:.0%} identity; {len(traits)} have trait values")
    return traits, kinds


def faprotax_traits(
    asvs: dict[str, str], taxonomy_path: Path, resources: Path, work: Path
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    """FAPROTAX functional groups per ASV, from the project's taxonomy."""
    folder = resources / "faprotax/FAPROTAX_1.2.12"
    script, database = folder / "collapse_table.py", folder / "FAPROTAX.txt"
    for path in (script, database):
        if not path.exists():
            raise FileNotFoundError(f"FAPROTAX resource missing: {path}")

    ranks = ["Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species"]
    table = work / "faprotax_input.tsv"
    with taxonomy_path.open(newline="", encoding="utf-8", errors="replace") as source, \
            table.open("w") as out:
        out.write("#OTU ID\tcount\ttaxonomy\n")
        for row in csv.DictReader(source):
            if row["OTUID"] in asvs:
                out.write(f"{row['OTUID']}\t1\t{';'.join(row[r] for r in ranks)}\n")

    membership = work / "faprotax_groups2records.tsv"
    subprocess.run(
        [sys.executable, str(script), "-i", str(table), "-o", str(work / "faprotax_functions.tsv"),
         "-g", str(database), "-d", "taxonomy", "-c", "#", "--omit_columns", "0",
         "--column_names_are_in", "last_comment_line", "-r", str(work / "faprotax_report.txt"),
         "--out_groups2records_table", str(membership), "--force"],
        check=True, capture_output=True, text=True,
    )
    # collapse_table.py writes one row per input record, in input order, keyed by the
    # record's taxonomy string; the ASV is recovered from the row's position.
    with table.open() as handle:
        next(handle)
        records = [line.rstrip("\n").split("\t") for line in handle]
    with membership.open() as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]
    header, rows = rows[0], rows[1:]
    if len(rows) != len(records):
        raise ValueError(f"FAPROTAX returned {len(rows)} rows for {len(records)} ASVs")
    per_asv: dict[str, dict[str, str]] = {}
    for (asv, _, taxonomy), fields in zip(records, rows):
        if fields[0] != taxonomy:
            raise ValueError(f"FAPROTAX row order differs from the input at {asv}")
        member = [float(v) != 0 for v in fields[1:]]
        if not any(member):
            continue                      # in no functional group: not annotated
        per_asv[asv] = {f"faprotax_{group}": "TRUE" if flag else "FALSE"
                        for group, flag in zip(header[1:], member)}
    kinds = {name: BINARY for flags in per_asv.values() for name in flags}
    print(f"FAPROTAX: {len(per_asv)} of {len(asvs)} ASVs placed in at least one of "
          f"{len(kinds)} functional groups")
    return per_asv, kinds


def metatraits_traits(
    asvs: dict[str, str], taxonomy_path: Path, resources: Path
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    """metaTraits genus-level traits per ASV, matched on the genus name."""
    summary = resources / "metatraits/gtdb_genus_summary_all.tsv.gz"
    if not summary.exists():
        raise FileNotFoundError(f"metaTraits resource missing: {summary}")

    genus_of: dict[str, str] = {}
    with taxonomy_path.open(newline="", encoding="utf-8", errors="replace") as source:
        for row in csv.DictReader(source):
            if row["OTUID"] in asvs and row["Genus"] not in ("", "Unassigned"):
                genus_of[row["OTUID"]] = row["Genus"]
    wanted = set(genus_of.values())

    by_genus: dict[str, dict[str, str]] = {}
    kinds: dict[str, str] = {}
    with gzip.open(summary, "rt") as handle:
        column = {name: i for i, name in enumerate(handle.readline().rstrip("\n").split("\t"))}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genus = fields[column["taxon_name"]]
            if genus not in wanted:
                continue
            name = "metatraits_" + re.sub(r"[^a-z0-9]+", "_", fields[column["trait_name"]].lower()).strip("_")
            if fields[column["unit"]] == "boolean":
                consensus = fields[column["consensus_value"]]
                if consensus not in ("true", "false"):
                    continue
                kinds[name] = BINARY
                by_genus.setdefault(genus, {})[name] = consensus.upper()
            elif fields[column["median"]] != "NA":
                kinds[name] = NUMERIC
                by_genus.setdefault(genus, {})[name] = fields[column["median"]]
    per_asv = {asv: by_genus[genus] for asv, genus in genus_of.items() if genus in by_genus}
    print(f"metaTraits: {len(per_asv)} of {len(asvs)} ASVs belong to one of {len(by_genus)} "
          f"matched genera ({len(wanted) - len(by_genus)} genus names not found); "
          f"{len(kinds)} traits")
    return per_asv, kinds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", default="analysis_config.json")
    config_path = Path(parser.parse_args().config).resolve()
    config = json.loads(config_path.read_text())
    base_dir = config_path.parent
    settings = config.get("predicted_traits", {})

    tables = resolve(base_dir, config["outputs"]["directory"]) / "tables"
    resources = resolve(base_dir, settings.get("resources_dir", "resources"))
    asvs = read_asvs(tables / "generalist_specialist_asv_host_screen.csv")
    sequences = read_sequences(resolve(base_dir, config["inputs"]["rep_seqs"]), set(asvs))

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        bugbase, bugbase_kinds = bugbase_traits(
            sequences, resources, float(settings.get("greengenes_identity", 0.97)),
            int(config["analysis"].get("threads", 1)), work)
        taxonomy = resolve(base_dir, config["inputs"]["taxonomy"])
        faprotax, faprotax_kinds = faprotax_traits(asvs, taxonomy, resources, work)
        metatraits, metatraits_kinds = metatraits_traits(asvs, taxonomy, resources)

    sources = {"bugbase": (bugbase, "BugBase 0.1.0"), "faprotax": (faprotax, "FAPROTAX 1.2.12"),
               "metatraits": (metatraits, "metaTraits (GTDB genus summary)")}
    kinds = {**bugbase_kinds, **faprotax_kinds, **metatraits_kinds}
    columns = sorted(kinds, key=lambda k: (k.split("_")[0], k))
    with (tables / "predicted_traits_asv.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ASV", "greengenes_otu", *columns])
        for name, asv in asvs.items():
            row = {**bugbase.get(name, {}), **faprotax.get(name, {}), **metatraits.get(name, {})}
            writer.writerow([asv, row.get("greengenes_otu", "NA"),
                             *(row.get(c, "NA") for c in columns)])
    with (tables / "predicted_traits_catalogue.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["trait", "resource", "type", "n_annotated", "n_true"])
        for column in columns:
            source, resource = sources[column.split("_")[0]]
            present = [row[column] for row in source.values() if column in row]
            n_true = sum(v == "TRUE" for v in present) if kinds[column] == BINARY else "NA"
            writer.writerow([column, resource, kinds[column], len(present), n_true])
    print(f"{len(columns)} predicted traits for {len(asvs)} ASVs -> "
          f"{tables / 'predicted_traits_asv.csv'}")


if __name__ == "__main__":
    main()

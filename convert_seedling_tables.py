from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent
FEATURE_IN = BASE / "Datesets" / "feature-table.tsv"
TAX_IN = BASE / "Datesets" / "taxonomy.tsv"
COUNTS_OUT = BASE / "Datesets" / "asv-counts.csv"
TAXA_OUT = BASE / "Datesets" / "taxa.csv"
COUNTS_TMP = BASE / "Datesets" / "asv-counts.tmp.npy"

RANK_HEADERS = ["Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species"]
FALLBACKS = {
    "Kingdom": "Unassigned",
    "Phyla": "unassigned_phyla",
    "Class": "unassigned_class",
    "Order": "unassigned_order",
    "Family": "unassigned_family",
    "Genus": "unassigned_genus",
    "Species": "unassigned_species",
}
PREFIX_TO_HEADER = {
    "d": "Kingdom",
    "p": "Phyla",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}


def convert_feature_table() -> tuple[int, int]:
    with FEATURE_IN.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        header = next(reader)
        if not header or header[0] != "#OTU ID":
            raise ValueError("Unexpected feature-table.tsv header")
        sample_ids = header[1:]
        feature_ids = [row[0] for row in reader if row]

    sample_count = len(sample_ids)
    feature_count = len(feature_ids)
    matrix = np.memmap(COUNTS_TMP, dtype=np.int32, mode="w+", shape=(sample_count, feature_count))

    with FEATURE_IN.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader)
        if not first or not first[0].startswith("# Constructed from biom file"):
            raise ValueError("Unexpected feature-table.tsv first line")

        header = next(reader)
        if not header or header[0] != "#OTU ID":
            raise ValueError("Unexpected feature-table.tsv header")

        for column_idx, row in enumerate(reader):
            if not row:
                continue
            values = row[1:]
            if len(values) != sample_count:
                raise ValueError(f"Row length mismatch for feature {row[0]}")
            matrix[:, column_idx] = np.fromiter(
                (int(float(value)) for value in values),
                dtype=np.int32,
                count=sample_count,
            )

    matrix.flush()

    with COUNTS_OUT.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", *feature_ids])
        for sample_idx, sid in enumerate(sample_ids):
            writer.writerow([sid, *matrix[sample_idx, :].tolist()])

    del matrix
    COUNTS_TMP.unlink(missing_ok=True)
    return feature_count, sample_count


def convert_taxonomy() -> None:
    with TAX_IN.open(newline="") as f, TAXA_OUT.open("w", newline="") as out_f:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.writer(out_f)
        writer.writerow(["OTUID", *RANK_HEADERS])

        for row in reader:
            feature_id = row["Feature ID"]
            taxon = row["Taxon"].strip()
            values = {header: "" for header in RANK_HEADERS}

            if taxon == "Unassigned":
                values.update(FALLBACKS)
            else:
                for segment in [seg.strip() for seg in taxon.split(";") if seg.strip()]:
                    if "__" not in segment:
                        continue
                    prefix, label = segment.split("__", 1)
                    header_name = PREFIX_TO_HEADER.get(prefix)
                    if header_name is not None:
                        values[header_name] = label.strip()

                seen_missing = False
                for header_name in RANK_HEADERS:
                    if not values[header_name] or seen_missing:
                        seen_missing = True
                        values[header_name] = FALLBACKS[header_name]

            writer.writerow([feature_id, *(values[header] for header in RANK_HEADERS)])


def main() -> None:
    feature_rows, sample_rows = convert_feature_table()
    convert_taxonomy()
    print(f"feature_rows={feature_rows}")
    print(f"sample_rows={sample_rows}")
    print(f"counts_out={COUNTS_OUT}")
    print(f"taxa_out={TAXA_OUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate a QIIME2 paired-end manifest.tsv from files named:
  SAMPLE_1.fastq.gz
  SAMPLE_2.fastq.gz
"""

import os
import glob
import csv
import sys


def main():
    r1_files = sorted(glob.glob("*_1.fastq.gz"))

    if not r1_files:
        sys.exit("ERROR: No *_1.fastq.gz files found")

    rows = []

    for r1 in r1_files:
        r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")

        if not os.path.exists(r2):
            print("WARNING: Missing R2 for", r1)
            continue

        sample_id = os.path.basename(r1).replace("_1.fastq.gz", "")

        rows.append((
            sample_id,
            os.path.abspath(r1),
            os.path.abspath(r2)
        ))

    if not rows:
        sys.exit("ERROR: No complete paired-end samples found")

    with open("manifest.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample-id",
            "forward-absolute-filepath",
            "reverse-absolute-filepath"
        ])
        for row in rows:
            writer.writerow(row)

    print("OK: manifest.tsv written")
    print("Samples detected:", len(rows))


if __name__ == "__main__":
    main()

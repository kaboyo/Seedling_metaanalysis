#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate QIIME2 paired-end manifest files
Supports:
  *_1.clean.fastq.gz / *_2.clean.fastq.gz
  *_R1.clean.fastq.gz / *_R2.clean.fastq.gz
Python 2 and 3 compatible
"""

import argparse
import os
import csv
import sys
import glob


def get_sample_id(filename):
    """
    Remove read and cleaning suffixes to obtain sample ID
    """
    suffixes = [
        "_R1.clean.fastq.gz",
        "_R2.clean.fastq.gz",
        "_1.clean.fastq.gz",
        "_2.clean.fastq.gz",
        ".clean.fastq.gz",
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq"
    ]

    for s in suffixes:
        if filename.endswith(s):
            return filename[:-len(s)]

    return filename


def find_paired_end(input_dir):
    """
    Find paired-end FASTQ files with clean suffixes
    """
    r1_files = (
        glob.glob(os.path.join(input_dir, "*_R1.clean.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_1.clean.fastq.gz"))
    )

    pairs = []

    for r1 in sorted(r1_files):
        if r1.endswith("_R1.clean.fastq.gz"):
            r2 = r1.replace("_R1.clean.fastq.gz", "_R2.clean.fastq.gz")
        else:
            r2 = r1.replace("_1.clean.fastq.gz", "_2.clean.fastq.gz")

        if not os.path.exists(r2):
            print("WARNING: Missing R2 for", os.path.basename(r1))
            continue

        sample_id = get_sample_id(os.path.basename(r1))

        pairs.append((
            sample_id,
            os.path.abspath(r1),
            os.path.abspath(r2)
        ))

    return pairs


def write_manifest(output, pairs):
    with open(output, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "sample-id",
            "forward-absolute-filepath",
            "reverse-absolute-filepath"
        ])
        for row in pairs:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Generate QIIME2 paired-end manifest from cleaned FASTQ files"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output manifest TSV")

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        sys.exit("ERROR: Input directory does not exist")

    pairs = find_paired_end(args.input)

    if not pairs:
        sys.exit("ERROR: No paired-end FASTQ files found")

    write_manifest(args.output, pairs)

    print("OK: Manifest written to", args.output)
    print("Samples detected:", len(pairs))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Rename *_1.clean.fastq.gz and *_2.clean.fastq.gz to
*_1.fastq.gz and *_2.fastq.gz
"""

import os
import glob

def main():
    files = glob.glob("*_*.clean.fastq.gz")

    if not files:
        print("ERROR: No *.clean.fastq.gz files found")
        return

    for old in sorted(files):
        new = old.replace(".clean.fastq.gz", ".fastq.gz")

        if os.path.exists(new):
            print("SKIP (already exists):", new)
            continue

        os.rename(old, new)
        print("RENAMED:", old, "->", new)

    print("Done.")

if __name__ == "__main__":
    main()

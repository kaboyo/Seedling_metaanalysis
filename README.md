# Seedling Meta-analysis Project

This repository contains scripts and metadata for the meta-analysis of seedling data. The primary goal is to process raw sequencing data, prepare it for analysis with QIIME2, and perform subsequent analyses.

## Repository Contents

- `Metadata_seedlings.csv` / `Metadata_seedlings.xlsx`: The sample metadata for the project, detailing each sample's properties.
- `generate_qiime2_manifest.py`: A Python script to generate the manifest file required by QIIME2 to import data.
- `merge_fastp.sh`: A shell script for merging paired-end reads using `fastp`.
- `merge_with_fastp_bash`: Another script related to `fastp` merging.
- `merge_with_fastp_multisamples_bash`: A script for running `fastp` on multiple samples.
- `script_multiqc.sh`: A script to run MultiQC for quality control analysis on the processed data.
- `elisabeth_script.md`: Markdown file containing notes or scripts from Elisabeth.

## Workflow

1.  **Prepare Metadata**: Ensure `Metadata_seedlings.csv` is up-to-date.
2.  **Generate Manifest**: Run `generate_qiime2_manifest.py` to create the QIIME2 manifest file.
3.  **Process Raw Data**: Use the `merge_...` scripts to process raw sequencing reads.
4.  **Quality Control**: Run `script_multiqc.sh` to generate a MultiQC report.
5.  **QIIME2 Analysis**: Use the generated manifest to import data into QIIME2 and proceed with analysis.

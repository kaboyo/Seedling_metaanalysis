#!/bin/bash

# Description:
# This script runs fastp on all paired-end FASTQ files in the current directory
# and then aggregates the reports using MultiQC.

# Run MultiQC on the fastp reports
echo "Running MultiQC..."
multiqc /Users/expeditoolimi/Documents/Southampton/Meta-analysis/Seedling_metaanalysis/reports -o multiqc_report

echo "All done. MultiQC report saved in ./multiqc_report"
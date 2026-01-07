#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH --time=14:30:00
#SBATCH --ntasks-per-core=1

module load conda
conda init
conda activate fastp_env
set -euo pipefail

# =========================
# fastp paired-end merging
# =========================

# Check that fastp is available
if ! command -v fastp &> /dev/null; then
  echo "ERROR: fastp not found in PATH"
  exit 1
fi

# Create output directories
mkdir -p merged unmerged qc

# Loop over forward reads
for R1 in *_1.clean.fastq.gz
do
  # Skip if no files match
  [[ -e "$R1" ]] || { echo "No input files found"; exit 1; }

  SAMPLE=$(basename "$R1" _1.clean.fastq.gz)
  R2="${SAMPLE}_2.clean.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "WARNING: Missing R2 for ${SAMPLE}, skipping"
    continue
  fi

  echo "Processing sample: ${SAMPLE}"

  fastp \
    -i "$R1" \
    -I "$R2" \
    --merge \
    --merged_out "merged/${SAMPLE}.merged.fastq.gz" \
    --unpaired1 "unmerged/${SAMPLE}.unmerged_R1.fastq.gz" \
    --unpaired2 "unmerged/${SAMPLE}.unmerged_R2.fastq.gz" \
    --overlap_len_require 20 \
    --overlap_diff_limit 5 \
    --overlap_diff_percent_limit 20 \
    --disable_trim_poly_g \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --html "qc/${SAMPLE}.fastp.html" \
    --json "qc/${SAMPLE}.fastp.json"

done

echo "All samples processed successfully."

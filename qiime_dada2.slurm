#!/bin/bash
#SBATCH --job-name=qiime_dada2
#SBATCH --output=qiime_dada2_%j.out
#SBATCH --error=qiime_dada2_%j.err
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Initialise Conda properly
# -----------------------------
# replace this path with the actual conda base on your system
source /local/software/conda/miniconda-py3-new/etc/profile.d/conda.sh  

# Now activate your QIIME 2 environment
conda activate qiime2-2023.5

# -----------------------------
# Import into QIIME 2
# -----------------------------
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /scratch/eo2r24/SEEDLING_DATA/data/data-output/demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --o-table /scratch/eo2r24/SEEDLING_DATA/data/data-output/table.qza \
  --o-representative-sequences /scratch/eo2r24/SEEDLING_DATA/data/data-output/rep-seqs.qza \
  --o-denoising-stats /scratch/eo2r24/SEEDLING_DATA/data/data-output/dada2-stats.qza


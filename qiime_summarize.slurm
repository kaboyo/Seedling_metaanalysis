#!/bin/bash
#SBATCH --job-name=qiime_import
#SBATCH --output=qiime_import_%j.out
#SBATCH --error=qiime_import_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

echo "Job started on $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Working directory: $(pwd)"

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
qiime demux summarize \
  --i-data /scratch/eo2r24/SEEDLING_DATA/data/data-output/demux.qza \
  --o-visualization /scratch/eo2r24/SEEDLING_DATA/data/data-output/demux.qzv
echo "QIIME 2 sumarize completed"
echo "Job finished on $(date)"

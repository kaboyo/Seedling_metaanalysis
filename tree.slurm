#!/bin/bash
#SBATCH --job-name=qiime_tree
#SBATCH --output=qiime_tree_%j.out
#SBATCH --error=qiime_tree_%j.err
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


qiime phylogeny align-to-tree-mafft-fasttree --i-sequences /scratch/eo2r24/SEEDLING_DATA/data/data-output/rep-seqs.qza \
 --output-dir /scratch/eo2r24/SEEDLING_DATA/data/data-output/tree \
 --p-n-threads 50 
--verbose




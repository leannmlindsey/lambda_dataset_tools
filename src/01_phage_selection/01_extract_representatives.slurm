#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=01_phage_select
#SBATCH --output=01_phage_select_%j.out
#SBATCH --error=01_phage_select_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g

# Step 01: select one random representative per vclust cluster and split into
# train/dev/test by cluster (seeded, no cluster overlap).
#
# Outputs -> ${DATA_DIR}/inphared_dataset/  (accession lists + cluster tsvs),
# which the subsampling step (04) consumes for all of 2k/4k/8k.

set -euo pipefail
# SLURM copies this script to a spool dir, so $0/BASH_SOURCE cannot locate the
# repo. Source config by absolute path (override with LAMBDA_CONFIG if cloned
# elsewhere); config.sh then defines SRC_DIR and all other absolute paths.
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

OUTPUT_DIR="${DATA_DIR}/inphared_dataset"
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "Step 01: Phage representative selection"
echo "  vclust:  ${INPHARED_VCLUST}"
echo "  output:  ${OUTPUT_DIR}"
echo "  split:   ${SPLIT_RATIO}   seed: ${SEED}"
echo "  start:   $(date)"
echo "=========================================="

python "${SRC_DIR}/01_phage_selection/extract_representatives.py" \
    --clusters "${INPHARED_VCLUST}" \
    --output "${OUTPUT_DIR}" \
    --split-ratio "${SPLIT_RATIO}" \
    --seed "${SEED}"

echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=select_phage_only
#SBATCH --output=select_phage_only_%j.out
#SBATCH --error=select_phage_only_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g

# Select phage genomes for Phage-Only benchmark (no overlap with fine-tuning)
#
# This identifies clusters in vclust that are NOT used in the fine-tuning
# train/dev/test splits, and selects representatives from those clusters.
#
# Usage: sbatch select_phage_only.slurm

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Source data
VCLUST_FILE="/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv"

# Fine-tuning dataset directory (contains train/dev/test_accessions.txt)
FINETUNING_DIR="${SCRIPT_DIR}/inphared_dataset"

# Output directory
OUTPUT_DIR="${SCRIPT_DIR}/phage_only_dataset"

# Random seed
SEED=42

# ============================================
# END CONFIGURATION
# ============================================

echo "=========================================="
echo "Select Phage-Only Genomes"
echo "=========================================="
echo "Start time: $(date)"
echo "vclust file: ${VCLUST_FILE}"
echo "Fine-tuning dir: ${FINETUNING_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Seed: ${SEED}"
echo "=========================================="

# Check files exist
if [ ! -f "${VCLUST_FILE}" ]; then
    echo "Error: vclust file not found: ${VCLUST_FILE}"
    exit 1
fi

if [ ! -d "${FINETUNING_DIR}" ]; then
    echo "Error: Fine-tuning directory not found: ${FINETUNING_DIR}"
    exit 1
fi

# Load environment
module load anaconda3
conda activate phage_env

# Run selection script
python ${SCRIPT_DIR}/select_phage_only_genomes.py \
    --vclust "${VCLUST_FILE}" \
    --finetuning-dir "${FINETUNING_DIR}" \
    --output "${OUTPUT_DIR}" \
    --seed ${SEED}

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Selection completed successfully!"
    echo "End time: $(date)"
    echo "=========================================="

    echo ""
    echo "Output files:"
    ls -la ${OUTPUT_DIR}/
else
    echo "Error: Selection failed!"
    exit 1
fi

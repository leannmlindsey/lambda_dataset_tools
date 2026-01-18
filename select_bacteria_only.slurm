#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=select_bacteria_only
#SBATCH --output=select_bacteria_only_%j.out
#SBATCH --error=select_bacteria_only_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g

# Select bacterial genomes for Bacteria-Only benchmark (no overlap with fine-tuning)
#
# This identifies genera in GTDB that are NOT used in the fine-tuning
# train/dev/test splits, and selects representatives from those genera.
# Also excludes prophage-contaminated genomes.
#
# Usage: sbatch select_bacteria_only.slurm

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Source data
GTDB_METADATA="/projects/bfzj/llindsey1/black_and_white/data/gtdb/taxonomy/bac120_metadata.tsv"

# Fine-tuning dataset directory (contains train/dev/test_accessions.txt)
FINETUNING_DIR="${SCRIPT_DIR}/gtdb_dataset_filtered_v3"

# Contaminated accessions (prophage-containing genomes to exclude)
CONTAMINATED_FILE="${SCRIPT_DIR}/gtdb_dataset_filtered_v3/contaminated_accessions.txt"

# Output directory
OUTPUT_DIR="${SCRIPT_DIR}/bacteria_only_dataset"

# Quality filters (same as fine-tuning selection)
MIN_COMPLETENESS=95.0
MAX_CONTAMINATION=5.0

# Random seed
SEED=42

# ============================================
# END CONFIGURATION
# ============================================

echo "=========================================="
echo "Select Bacteria-Only Genomes"
echo "=========================================="
echo "Start time: $(date)"
echo "GTDB metadata: ${GTDB_METADATA}"
echo "Fine-tuning dir: ${FINETUNING_DIR}"
echo "Contaminated file: ${CONTAMINATED_FILE}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Min completeness: ${MIN_COMPLETENESS}%"
echo "Max contamination: ${MAX_CONTAMINATION}%"
echo "Seed: ${SEED}"
echo "=========================================="

# Check files exist
if [ ! -f "${GTDB_METADATA}" ]; then
    echo "Error: GTDB metadata not found: ${GTDB_METADATA}"
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
python ${SCRIPT_DIR}/select_bacteria_only_genomes.py \
    --metadata "${GTDB_METADATA}" \
    --finetuning-dir "${FINETUNING_DIR}" \
    --contaminated "${CONTAMINATED_FILE}" \
    --output "${OUTPUT_DIR}" \
    --min-completeness ${MIN_COMPLETENESS} \
    --max-contamination ${MAX_CONTAMINATION} \
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

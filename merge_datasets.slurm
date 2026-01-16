#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=merge_datasets
#SBATCH --output=merge_datasets_%j.out
#SBATCH --error=merge_datasets_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g

# Merge and shuffle INPHARED (phage) and GTDB (bacteria) datasets
#
# Usage: sbatch merge_datasets.slurm

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"
BASE_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"
OUTPUT_BASE="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/merged_datasets"

SEED=42

# ============================================

echo "=========================================="
echo "Merge and Shuffle Datasets"
echo "=========================================="
echo "Start time: $(date)"
echo "Output base: ${OUTPUT_BASE}"
echo "=========================================="

# Process each segment length
for SEG_LEN in 2k 4k 8k; do
    echo ""
    echo "=========================================="
    echo "Processing ${SEG_LEN} segment datasets"
    echo "=========================================="

    # Set input directories based on segment length
    if [ "${SEG_LEN}" == "2k" ]; then
        INPHARED_DIR="${BASE_DIR}/inphared_dataset"
        GTDB_DIR="${BASE_DIR}/gtdb_dataset"
    else
        INPHARED_DIR="${BASE_DIR}/inphared_dataset_${SEG_LEN}"
        GTDB_DIR="${BASE_DIR}/gtdb_dataset_${SEG_LEN}"
    fi

    OUTPUT_DIR="${OUTPUT_BASE}/${SEG_LEN}"
    mkdir -p ${OUTPUT_DIR}

    # Process each split
    for SPLIT in train dev test; do
        echo ""
        echo "Processing ${SEG_LEN} ${SPLIT}..."

        PHAGE_FASTA="${INPHARED_DIR}/${SPLIT}_segments.fasta"
        BACTERIA_FASTA="${GTDB_DIR}/${SPLIT}_segments.fasta"
        OUTPUT_FASTA="${OUTPUT_DIR}/${SPLIT}_merged.fasta"
        OUTPUT_LABELS="${OUTPUT_DIR}/${SPLIT}_labels.tsv"

        # Check input files exist
        if [ ! -f "${PHAGE_FASTA}" ]; then
            echo "Warning: Phage FASTA not found: ${PHAGE_FASTA}"
            echo "Skipping ${SEG_LEN} ${SPLIT}..."
            continue
        fi

        if [ ! -f "${BACTERIA_FASTA}" ]; then
            echo "Warning: Bacteria FASTA not found: ${BACTERIA_FASTA}"
            echo "Skipping ${SEG_LEN} ${SPLIT}..."
            continue
        fi

        python ${SCRIPT_DIR}/merge_and_shuffle.py \
            --phage-fasta "${PHAGE_FASTA}" \
            --bacteria-fasta "${BACTERIA_FASTA}" \
            --output-fasta "${OUTPUT_FASTA}" \
            --output-labels "${OUTPUT_LABELS}" \
            --seed ${SEED}

        if [ $? -eq 0 ]; then
            echo "${SEG_LEN} ${SPLIT} completed!"
        else
            echo "Error: ${SEG_LEN} ${SPLIT} failed!"
        fi
    done
done

echo ""
echo "=========================================="
echo "All merging completed!"
echo "End time: $(date)"
echo "=========================================="

# Show output summary
echo ""
echo "Output files:"
for SEG_LEN in 2k 4k 8k; do
    echo ""
    echo "${SEG_LEN} segments:"
    ls -lh ${OUTPUT_BASE}/${SEG_LEN}/*.fasta 2>/dev/null || echo "  No files found"
done

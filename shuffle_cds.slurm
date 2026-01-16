#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=shuffle_cds
#SBATCH --output=shuffle_cds_%j.out
#SBATCH --error=shuffle_cds_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g

# Shuffle nucleotides within CDS regions to create control sequences
#
# Usage: sbatch shuffle_cds.slurm

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Input: Pharokka annotation directory (contains GenBank files)
PHAROKKA_DIR="/u/llindsey1/llindsey/black_and_white/data/inphared/pharokka"

# Input: List of accessions to process
ACCESSION_LIST="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset/all_representatives.txt"

# Output: Directory for shuffled FASTA files
OUTPUT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/shuffled_genomes"

# GenBank file suffix (Pharokka naming pattern)
GBK_SUFFIX="_genome.gbk"

SEED=42

# ============================================

echo "=========================================="
echo "Shuffle CDS Nucleotides"
echo "=========================================="
echo "Start time: $(date)"
echo "Pharokka directory: ${PHAROKKA_DIR}"
echo "Accession list: ${ACCESSION_LIST}"
echo "Output directory: ${OUTPUT_DIR}"
echo "GenBank suffix: ${GBK_SUFFIX}"
echo "Seed: ${SEED}"
echo "=========================================="

# Check inputs exist
if [ ! -d "${PHAROKKA_DIR}" ]; then
    echo "Error: Pharokka directory not found: ${PHAROKKA_DIR}"
    exit 1
fi

if [ ! -f "${ACCESSION_LIST}" ]; then
    echo "Error: Accession list not found: ${ACCESSION_LIST}"
    exit 1
fi

# Show sample of what files we're looking for
echo ""
echo "Looking for GenBank files like: <accession>${GBK_SUFFIX}"
echo "Sample files in Pharokka directory:"
ls ${PHAROKKA_DIR}/*.gbk 2>/dev/null | head -3 || echo "  No .gbk files found in root"
ls ${PHAROKKA_DIR}/*/*.gbk 2>/dev/null | head -3 || echo "  No .gbk files found in subdirs"
echo ""

mkdir -p ${OUTPUT_DIR}

python ${SCRIPT_DIR}/shuffle_cds_nucleotides.py \
    --input-dir "${PHAROKKA_DIR}" \
    --accession-list "${ACCESSION_LIST}" \
    --output-dir "${OUTPUT_DIR}" \
    --shuffle-intergenic \
    --suffix "${GBK_SUFFIX}" \
    --seed ${SEED}

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Shuffling completed!"
    echo "End time: $(date)"
    echo "=========================================="
    echo ""
    echo "Output files:"
    ls -lh ${OUTPUT_DIR}/*.fasta 2>/dev/null | head -10
    echo "..."
    echo "Total shuffled files: $(ls ${OUTPUT_DIR}/*.fasta 2>/dev/null | wc -l)"
else
    echo ""
    echo "=========================================="
    echo "Shuffling failed!"
    echo "End time: $(date)"
    echo "=========================================="
    exit 1
fi

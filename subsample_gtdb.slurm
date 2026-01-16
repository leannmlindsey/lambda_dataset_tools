#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=gtdb_subsample
#SBATCH --output=gtdb_subsample_%j.out
#SBATCH --error=gtdb_subsample_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Subsample segments from GTDB bacterial genomes for train/dev/test splits
#
# Handles GTDB's nested directory structure and gzipped files.
# Uses threading to process multiple genomes in parallel.
#
# Usage: sbatch subsample_gtdb.slurm

# ============================================
# CONFIGURATION - Modify these paths as needed
# ============================================

# Directory containing the scripts
SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# GTDB base directory (contains database/GCA/... structure with .fna.gz files)
GTDB_BASE="/u/llindsey1/llindsey/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226"

# GTDB dataset directory (contains accession lists and will receive output)
GTDB_DATASET_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset"

# Number of parallel threads (adjust based on available CPUs)
THREADS=16

# Subsampling parameters
SEGMENT_LENGTH=2000
MIN_LENGTH=5000
SAMPLE_PER_BP=10000
SEED=42

# Optional: Set TOTAL_SEGMENTS to match phage segment count
# Leave empty to use SAMPLE_PER_BP instead
# TOTAL_SEGMENTS=20000
TOTAL_SEGMENTS=""

# ============================================
# END CONFIGURATION
# ============================================

echo "=========================================="
echo "GTDB Segment Subsampling"
echo "=========================================="
echo "Start time: $(date)"
echo "Script directory: ${SCRIPT_DIR}"
echo "GTDB base directory: ${GTDB_BASE}"
echo "GTDB dataset directory: ${GTDB_DATASET_DIR}"
echo "Threads: ${THREADS}"
echo "Segment length: ${SEGMENT_LENGTH}"
echo "Min genome length: ${MIN_LENGTH}"
if [ -n "${TOTAL_SEGMENTS}" ]; then
    echo "Total segments target: ${TOTAL_SEGMENTS}"
else
    echo "Sample per bp: ${SAMPLE_PER_BP}"
fi
echo "Seed: ${SEED}"
echo "=========================================="

# Check that GTDB directory exists
if [ ! -d "${GTDB_BASE}/database" ]; then
    echo "Error: GTDB database directory not found: ${GTDB_BASE}/database"
    exit 1
fi

# Load any required modules (uncomment/modify as needed)
# module load python/3.9
# conda activate lambda_tools

# Build the optional arguments
EXTRA_ARGS=""
if [ -n "${TOTAL_SEGMENTS}" ]; then
    EXTRA_ARGS="--total-segments ${TOTAL_SEGMENTS}"
fi

# Process each split
for SPLIT in train dev test; do
    echo ""
    echo "=========================================="
    echo "Processing ${SPLIT} split..."
    echo "=========================================="

    ACCESSION_FILE="${GTDB_DATASET_DIR}/${SPLIT}_accessions.txt"
    OUTPUT_FASTA="${GTDB_DATASET_DIR}/${SPLIT}_segments.fasta"
    OUTPUT_METADATA="${GTDB_DATASET_DIR}/${SPLIT}_segments.tsv"

    if [ ! -f "${ACCESSION_FILE}" ]; then
        echo "Warning: Accession file not found: ${ACCESSION_FILE}"
        echo "Skipping ${SPLIT} split..."
        continue
    fi

    N_ACCESSIONS=$(wc -l < "${ACCESSION_FILE}")
    echo "Accession file: ${ACCESSION_FILE} (${N_ACCESSIONS} accessions)"
    echo "Output FASTA: ${OUTPUT_FASTA}"
    echo "Output metadata: ${OUTPUT_METADATA}"

    python ${SCRIPT_DIR}/subsample_gtdb_segments.py \
        --gtdb-dir "${GTDB_BASE}" \
        --accession-list "${ACCESSION_FILE}" \
        --output "${OUTPUT_FASTA}" \
        --output-metadata "${OUTPUT_METADATA}" \
        --segment-length ${SEGMENT_LENGTH} \
        --min-length ${MIN_LENGTH} \
        --sample-per-bp ${SAMPLE_PER_BP} \
        --threads ${THREADS} \
        --seed ${SEED} \
        ${EXTRA_ARGS}

    if [ $? -eq 0 ]; then
        echo "${SPLIT} split completed successfully!"
    else
        echo "Error: ${SPLIT} split failed!"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "All splits completed!"
echo "End time: $(date)"
echo "=========================================="

# List output files
echo ""
echo "Output files:"
ls -lh ${GTDB_DATASET_DIR}/*_segments.fasta ${GTDB_DATASET_DIR}/*_segments.tsv 2>/dev/null

# Show segment counts
echo ""
echo "Segment counts:"
for SPLIT in train dev test; do
    if [ -f "${GTDB_DATASET_DIR}/${SPLIT}_segments.fasta" ]; then
        COUNT=$(grep -c "^>" "${GTDB_DATASET_DIR}/${SPLIT}_segments.fasta" 2>/dev/null || echo "0")
        echo "  ${SPLIT}: ${COUNT} segments"
    fi
done

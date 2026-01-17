#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=filter_subsample
#SBATCH --output=filter_subsample_%j.out
#SBATCH --error=filter_subsample_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Filter prophage-contaminated bacteria and subsample segments at 2k, 4k, 8k
#
# This script:
# 1. Combines v2 BLAST results
# 2. Filters out contaminated bacterial accessions
# 3. Subsamples segments at 2000, 4000, and 8000 bp lengths
#
# Usage: sbatch filter_and_subsample_gtdb.slurm

# ============================================
# CONFIGURATION
# ============================================

# Script directory
SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# BLAST results from v2 (with alignment length filter)
BLAST_RESULTS_DIR="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2"

# Original GTDB accession lists
GTDB_ACCESSIONS="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset"

# Base output directory for filtered datasets
OUTPUT_BASE="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# GTDB genomes directory
GTDB_BASE="/u/llindsey1/llindsey/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226"

# Filtering thresholds
MIN_IDENTITY=90
MIN_ALIGNMENT_LENGTH=200

# Subsampling parameters
THREADS=16
MIN_LENGTH=5000
SEED=42

# Target segment counts (should match INPHARED counts)
# Update these after checking your INPHARED dataset counts
TOTAL_SEGMENTS_TRAIN=27000
TOTAL_SEGMENTS_DEV=3400
TOTAL_SEGMENTS_TEST=3400

# ============================================

echo "=========================================="
echo "Filter and Subsample GTDB (Prophage-Free)"
echo "=========================================="
echo "Start time: $(date)"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# ============================================
# STEP 1: Combine BLAST results
# ============================================
echo ""
echo "Step 1: Combining BLAST results..."

COMBINED_BLAST="${BLAST_RESULTS_DIR}/blast_filtered_combined.tsv"

if [ ! -f "${COMBINED_BLAST}" ]; then
    echo "  Combining filtered results from all array tasks..."
    cat ${BLAST_RESULTS_DIR}/blast_filtered_*.tsv > ${COMBINED_BLAST}
    echo "  Combined: $(wc -l < ${COMBINED_BLAST}) total hits"
else
    echo "  Using existing combined file: ${COMBINED_BLAST}"
    echo "  Total hits: $(wc -l < ${COMBINED_BLAST})"
fi

# ============================================
# STEP 2: Filter contaminated accessions
# ============================================
echo ""
echo "Step 2: Filtering contaminated accessions..."

# Output directory for filtered accession lists
FILTERED_ACCESSIONS="${OUTPUT_BASE}/gtdb_dataset_filtered"
mkdir -p ${FILTERED_ACCESSIONS}

python ${SCRIPT_DIR}/filter_prophage_contaminated.py \
    --blast-results "${COMBINED_BLAST}" \
    --input-accessions-dir "${GTDB_ACCESSIONS}" \
    --output-accessions-dir "${FILTERED_ACCESSIONS}" \
    --min-identity ${MIN_IDENTITY} \
    --min-alignment-length ${MIN_ALIGNMENT_LENGTH} \
    --output-contaminated "${FILTERED_ACCESSIONS}/contaminated_accessions.txt" \
    --output-summary "${FILTERED_ACCESSIONS}/filtering_summary.txt"

if [ $? -ne 0 ]; then
    echo "Error: Filtering failed!"
    exit 1
fi

# ============================================
# STEP 3: Subsample at 2000 bp
# ============================================
echo ""
echo "=========================================="
echo "Step 3: Subsampling at 2000 bp..."
echo "=========================================="

SEGMENT_LENGTH=2000

for SPLIT in train dev test; do
    echo ""
    echo "Processing ${SPLIT} split (${SEGMENT_LENGTH} bp)..."

    ACCESSION_FILE="${FILTERED_ACCESSIONS}/${SPLIT}_accessions.txt"
    OUTPUT_FASTA="${FILTERED_ACCESSIONS}/${SPLIT}_segments.fasta"
    OUTPUT_METADATA="${FILTERED_ACCESSIONS}/${SPLIT}_segments.tsv"

    if [ "${SPLIT}" == "train" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TRAIN}
    elif [ "${SPLIT}" == "dev" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_DEV}
    else
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TEST}
    fi

    if [ ! -f "${ACCESSION_FILE}" ]; then
        echo "Warning: Accession file not found: ${ACCESSION_FILE}"
        continue
    fi

    N_ACCESSIONS=$(wc -l < "${ACCESSION_FILE}")
    echo "  Accessions: ${N_ACCESSIONS}"
    echo "  Target segments: ${TARGET_SEGMENTS}"

    python ${SCRIPT_DIR}/subsample_gtdb_segments.py \
        --gtdb-dir "${GTDB_BASE}" \
        --accession-list "${ACCESSION_FILE}" \
        --output "${OUTPUT_FASTA}" \
        --output-metadata "${OUTPUT_METADATA}" \
        --segment-length ${SEGMENT_LENGTH} \
        --min-length ${MIN_LENGTH} \
        --threads ${THREADS} \
        --seed ${SEED} \
        --total-segments ${TARGET_SEGMENTS}

    if [ $? -eq 0 ]; then
        ACTUAL=$(grep -c "^>" "${OUTPUT_FASTA}" 2>/dev/null || echo "0")
        echo "  ${SPLIT} complete: ${ACTUAL} segments"
    else
        echo "  Error: ${SPLIT} failed!"
    fi
done

# ============================================
# STEP 4: Subsample at 4000 bp
# ============================================
echo ""
echo "=========================================="
echo "Step 4: Subsampling at 4000 bp..."
echo "=========================================="

SEGMENT_LENGTH=4000
OUTPUT_4K="${OUTPUT_BASE}/gtdb_dataset_filtered_4k"
mkdir -p ${OUTPUT_4K}

for SPLIT in train dev test; do
    echo ""
    echo "Processing ${SPLIT} split (${SEGMENT_LENGTH} bp)..."

    ACCESSION_FILE="${FILTERED_ACCESSIONS}/${SPLIT}_accessions.txt"
    OUTPUT_FASTA="${OUTPUT_4K}/${SPLIT}_segments.fasta"
    OUTPUT_METADATA="${OUTPUT_4K}/${SPLIT}_segments.tsv"

    if [ "${SPLIT}" == "train" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TRAIN}
    elif [ "${SPLIT}" == "dev" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_DEV}
    else
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TEST}
    fi

    if [ ! -f "${ACCESSION_FILE}" ]; then
        echo "Warning: Accession file not found: ${ACCESSION_FILE}"
        continue
    fi

    python ${SCRIPT_DIR}/subsample_gtdb_segments.py \
        --gtdb-dir "${GTDB_BASE}" \
        --accession-list "${ACCESSION_FILE}" \
        --output "${OUTPUT_FASTA}" \
        --output-metadata "${OUTPUT_METADATA}" \
        --segment-length ${SEGMENT_LENGTH} \
        --min-length ${MIN_LENGTH} \
        --threads ${THREADS} \
        --seed ${SEED} \
        --total-segments ${TARGET_SEGMENTS}

    if [ $? -eq 0 ]; then
        ACTUAL=$(grep -c "^>" "${OUTPUT_FASTA}" 2>/dev/null || echo "0")
        echo "  ${SPLIT} complete: ${ACTUAL} segments"
    else
        echo "  Error: ${SPLIT} failed!"
    fi
done

# ============================================
# STEP 5: Subsample at 8000 bp
# ============================================
echo ""
echo "=========================================="
echo "Step 5: Subsampling at 8000 bp..."
echo "=========================================="

SEGMENT_LENGTH=8000
OUTPUT_8K="${OUTPUT_BASE}/gtdb_dataset_filtered_8k"
mkdir -p ${OUTPUT_8K}

for SPLIT in train dev test; do
    echo ""
    echo "Processing ${SPLIT} split (${SEGMENT_LENGTH} bp)..."

    ACCESSION_FILE="${FILTERED_ACCESSIONS}/${SPLIT}_accessions.txt"
    OUTPUT_FASTA="${OUTPUT_8K}/${SPLIT}_segments.fasta"
    OUTPUT_METADATA="${OUTPUT_8K}/${SPLIT}_segments.tsv"

    if [ "${SPLIT}" == "train" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TRAIN}
    elif [ "${SPLIT}" == "dev" ]; then
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_DEV}
    else
        TARGET_SEGMENTS=${TOTAL_SEGMENTS_TEST}
    fi

    if [ ! -f "${ACCESSION_FILE}" ]; then
        echo "Warning: Accession file not found: ${ACCESSION_FILE}"
        continue
    fi

    python ${SCRIPT_DIR}/subsample_gtdb_segments.py \
        --gtdb-dir "${GTDB_BASE}" \
        --accession-list "${ACCESSION_FILE}" \
        --output "${OUTPUT_FASTA}" \
        --output-metadata "${OUTPUT_METADATA}" \
        --segment-length ${SEGMENT_LENGTH} \
        --min-length ${MIN_LENGTH} \
        --threads ${THREADS} \
        --seed ${SEED} \
        --total-segments ${TARGET_SEGMENTS}

    if [ $? -eq 0 ]; then
        ACTUAL=$(grep -c "^>" "${OUTPUT_FASTA}" 2>/dev/null || echo "0")
        echo "  ${SPLIT} complete: ${ACTUAL} segments"
    else
        echo "  Error: ${SPLIT} failed!"
    fi
done

# ============================================
# SUMMARY
# ============================================
echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="

echo ""
echo "Filtered accession lists:"
ls -lh ${FILTERED_ACCESSIONS}/*_accessions.txt 2>/dev/null

echo ""
echo "2k segments (${FILTERED_ACCESSIONS}):"
for SPLIT in train dev test; do
    if [ -f "${FILTERED_ACCESSIONS}/${SPLIT}_segments.fasta" ]; then
        COUNT=$(grep -c "^>" "${FILTERED_ACCESSIONS}/${SPLIT}_segments.fasta" 2>/dev/null || echo "0")
        echo "  ${SPLIT}: ${COUNT} segments"
    fi
done

echo ""
echo "4k segments (${OUTPUT_4K}):"
for SPLIT in train dev test; do
    if [ -f "${OUTPUT_4K}/${SPLIT}_segments.fasta" ]; then
        COUNT=$(grep -c "^>" "${OUTPUT_4K}/${SPLIT}_segments.fasta" 2>/dev/null || echo "0")
        echo "  ${SPLIT}: ${COUNT} segments"
    fi
done

echo ""
echo "8k segments (${OUTPUT_8K}):"
for SPLIT in train dev test; do
    if [ -f "${OUTPUT_8K}/${SPLIT}_segments.fasta" ]; then
        COUNT=$(grep -c "^>" "${OUTPUT_8K}/${SPLIT}_segments.fasta" 2>/dev/null || echo "0")
        echo "  ${SPLIT}: ${COUNT} segments"
    fi
done

echo ""
echo "=========================================="
echo "Done!"
echo "End time: $(date)"
echo "=========================================="

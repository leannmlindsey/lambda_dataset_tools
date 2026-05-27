#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=derep_fragments
#SBATCH --output=derep_fragments_%j.out
#SBATCH --error=derep_fragments_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Fragment-level dereplication of dev/test against train, for BOTH classes.
#
# For each per-class dataset directory it:
#   1. builds a BLAST nucleotide DB from train_segments.fasta
#   2. BLASTs dev_segments.fasta and test_segments.fasta against it
#   3. removes any dev/test fragment near-identical to a train fragment
#      (>= MIN_IDENTITY% identity AND >= MIN_COVERAGE% query coverage)
#
# The originals are left untouched; dereplicated copies + reports are written
# to a sibling "<dataset>_derep" directory, so both the full and "hard"
# (dereplicated) test sets are available.
#
# After this runs, rebuild the merged datasets from the *_derep dev/test
# FASTAs (train is unchanged). See note at the end.
#
# Usage: sbatch dereplicate_fragments.slurm

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Per-class segment datasets to dereplicate (each must contain
# train_segments.fasta, dev_segments.fasta, test_segments.fasta).
DATASETS=(
    "${SCRIPT_DIR}/inphared_dataset"
    "${SCRIPT_DIR}/inphared_dataset_4k"
    "${SCRIPT_DIR}/inphared_dataset_8k"
    "${SCRIPT_DIR}/gtdb_dataset_filtered_v3"
    "${SCRIPT_DIR}/gtdb_dataset_filtered_v3_4k"
    "${SCRIPT_DIR}/gtdb_dataset_filtered_v3_8k"
)

# Leakage thresholds
MIN_IDENTITY=95
MIN_COVERAGE=50

# BLAST parameters
#   -perc_identity 90 pre-filters below the 95% cutoff (cheap; exact cut is
#   applied in the Python filter), -dust no avoids masking real matches.
BLAST_PERC_IDENTITY=90
NUM_THREADS=16
BLAST_OUTFMT="6 qseqid sseqid pident length qlen qcovhsp evalue bitscore"

# Scratch directory for BLAST databases
WORK_DIR="/work/hdd/bfzj/llindsey1/fragment_dereplication"

# ============================================

echo "=========================================="
echo "Fragment-level Dereplication"
echo "=========================================="
echo "Start time: $(date)"
echo "Leakage cutoff: >= ${MIN_IDENTITY}% identity AND >= ${MIN_COVERAGE}% query coverage"
echo "=========================================="

source ~/.bashrc
conda activate lambda_tools

mkdir -p "${WORK_DIR}"

for DATASET in "${DATASETS[@]}"; do
    NAME=$(basename "${DATASET}")
    DEREP_DIR="${DATASET}_derep"

    echo ""
    echo "##############################################"
    echo "# ${NAME}"
    echo "##############################################"

    TRAIN_FASTA="${DATASET}/train_segments.fasta"
    if [ ! -f "${TRAIN_FASTA}" ]; then
        echo "  Warning: train FASTA not found, skipping: ${TRAIN_FASTA}"
        continue
    fi

    mkdir -p "${DEREP_DIR}"

    # Train is unchanged; copy it through so the *_derep dir is self-contained.
    cp "${TRAIN_FASTA}" "${DEREP_DIR}/train_segments.fasta"

    # ---- Step 1: BLAST DB from train fragments ----
    DB="${WORK_DIR}/${NAME}_train_db"
    echo "  Building BLAST DB from $(grep -c '^>' ${TRAIN_FASTA}) train fragments..."
    makeblastdb -in "${TRAIN_FASTA}" -dbtype nucl -out "${DB}" >/dev/null

    # ---- Steps 2-3: search + filter each query split ----
    for SPLIT in dev test; do
        QUERY_FASTA="${DATASET}/${SPLIT}_segments.fasta"
        if [ ! -f "${QUERY_FASTA}" ]; then
            echo "  Warning: ${SPLIT} FASTA not found, skipping: ${QUERY_FASTA}"
            continue
        fi

        BLAST_OUT="${WORK_DIR}/${NAME}_${SPLIT}_vs_train.tsv"
        echo "  BLASTing ${SPLIT} ($(grep -c '^>' ${QUERY_FASTA}) fragments) vs train..."

        blastn \
            -db "${DB}" \
            -query "${QUERY_FASTA}" \
            -out "${BLAST_OUT}" \
            -outfmt "${BLAST_OUTFMT}" \
            -perc_identity ${BLAST_PERC_IDENTITY} \
            -dust no \
            -max_target_seqs 5 \
            -max_hsps 1 \
            -evalue 1e-5 \
            -num_threads ${NUM_THREADS}

        echo "    Raw BLAST hit lines: $(wc -l < ${BLAST_OUT})"

        python "${SCRIPT_DIR}/filter_leaked_fragments.py" \
            --query-fasta "${QUERY_FASTA}" \
            --blast-output "${BLAST_OUT}" \
            --output-fasta "${DEREP_DIR}/${SPLIT}_segments.fasta" \
            --output-leaked "${DEREP_DIR}/${SPLIT}_leaked.tsv" \
            --output-report "${DEREP_DIR}/${SPLIT}_dereplication_report.txt" \
            --min-identity ${MIN_IDENTITY} \
            --min-coverage ${MIN_COVERAGE}
    done

    echo "  Dereplicated output -> ${DEREP_DIR}"
done

# ============================================
# SUMMARY
# ============================================
echo ""
echo "=========================================="
echo "Summary (fragments kept / removed)"
echo "=========================================="
for DATASET in "${DATASETS[@]}"; do
    DEREP_DIR="${DATASET}_derep"
    NAME=$(basename "${DATASET}")
    [ -d "${DEREP_DIR}" ] || continue
    echo ""
    echo "${NAME}:"
    for SPLIT in dev test; do
        ORIG="${DATASET}/${SPLIT}_segments.fasta"
        NEW="${DEREP_DIR}/${SPLIT}_segments.fasta"
        if [ -f "${ORIG}" ] && [ -f "${NEW}" ]; then
            O=$(grep -c '^>' "${ORIG}")
            N=$(grep -c '^>' "${NEW}")
            echo "  ${SPLIT}: ${N}/${O} kept ($((O - N)) removed)"
        fi
    done
done

echo ""
echo "NOTE: train fragments are unchanged. To build the 'hard' merged datasets,"
echo "rerun merge_datasets.slurm pointing dev/test at the *_derep directories"
echo "(train can still point at the originals). Keep the originals to report the"
echo "full vs dereplicated metric gap."
echo ""
echo "=========================================="
echo "Done! End time: $(date)"
echo "=========================================="

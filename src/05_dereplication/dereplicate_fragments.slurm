#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=05_derep
#SBATCH --output=05_derep_%j.out
#SBATCH --error=05_derep_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Step 05: fragment-level dereplication of dev/test against train, both classes.
#
# For each per-class dataset it builds a BLAST db from train_segments.fasta,
# searches dev/test, and removes any dev/test fragment near-identical to a
# train fragment (>= DEREP_MIN_IDENTITY% id AND >= DEREP_MIN_COVERAGE% qcov).
# Originals are left intact; dereplicated copies + reports go to "<dataset>_derep".
#
# Genome/cluster splitting alone does not catch this: phages are mosaic and
# bacteria share conserved genes (rRNA, housekeeping), so fragments from
# genomes in different splits can still be near-identical. The model only sees
# fragments, so leakage must be removed at the fragment level.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

source ~/.bashrc
conda activate "${CONDA_ENV}"

DATASETS=("${PHAGE_DATASETS[@]}" "${BACTERIA_DATASETS[@]}")
BLAST_OUTFMT="6 qseqid sseqid pident length qlen qcovhsp evalue bitscore"
DB_WORK="${WORK_DIR}/fragment_dereplication"
mkdir -p "${DB_WORK}"

echo "=========================================="
echo "Step 05: Fragment-level dereplication"
echo "  cutoff: >= ${DEREP_MIN_IDENTITY}% id AND >= ${DEREP_MIN_COVERAGE}% qcov"
echo "  start:  $(date)"
echo "=========================================="

for DATASET in "${DATASETS[@]}"; do
    NAME="$(basename "${DATASET}")"
    DEREP_DIR="${DATASET}_derep"
    TRAIN_FASTA="${DATASET}/train_segments.fasta"

    echo ""
    echo "### ${NAME}"
    if [ ! -f "${TRAIN_FASTA}" ]; then
        echo "  Warning: train FASTA not found, skipping: ${TRAIN_FASTA}"
        continue
    fi

    mkdir -p "${DEREP_DIR}"
    cp "${TRAIN_FASTA}" "${DEREP_DIR}/train_segments.fasta"   # train unchanged

    DB="${DB_WORK}/${NAME}_train_db"
    echo "  Building train BLAST db ($(grep -c '^>' "${TRAIN_FASTA}") fragments)..."
    makeblastdb -in "${TRAIN_FASTA}" -dbtype nucl -out "${DB}" >/dev/null

    for SPLIT in dev test; do
        QUERY_FASTA="${DATASET}/${SPLIT}_segments.fasta"
        [ -f "${QUERY_FASTA}" ] || { echo "  Warning: missing ${QUERY_FASTA}"; continue; }

        BLAST_OUT="${DB_WORK}/${NAME}_${SPLIT}_vs_train.tsv"
        echo "  BLAST ${SPLIT} ($(grep -c '^>' "${QUERY_FASTA}") fragments) vs train..."
        blastn \
            -db "${DB}" \
            -query "${QUERY_FASTA}" \
            -out "${BLAST_OUT}" \
            -outfmt "${BLAST_OUTFMT}" \
            -perc_identity "${DEREP_BLAST_PERC_IDENTITY}" \
            -dust no \
            -max_target_seqs 5 \
            -max_hsps 1 \
            -evalue 1e-5 \
            -num_threads "${THREADS}"

        python "${HERE}/filter_leaked_fragments.py" \
            --query-fasta "${QUERY_FASTA}" \
            --blast-output "${BLAST_OUT}" \
            --output-fasta "${DEREP_DIR}/${SPLIT}_segments.fasta" \
            --output-leaked "${DEREP_DIR}/${SPLIT}_leaked.tsv" \
            --output-report "${DEREP_DIR}/${SPLIT}_dereplication_report.txt" \
            --min-identity "${DEREP_MIN_IDENTITY}" \
            --min-coverage "${DEREP_MIN_COVERAGE}"
    done
    echo "  -> ${DEREP_DIR}"
done

echo ""
echo "Summary (fragments kept / original):"
for DATASET in "${DATASETS[@]}"; do
    DEREP_DIR="${DATASET}_derep"; NAME="$(basename "${DATASET}")"
    [ -d "${DEREP_DIR}" ] || continue
    echo "  ${NAME}:"
    for SPLIT in dev test; do
        ORIG="${DATASET}/${SPLIT}_segments.fasta"; NEW="${DEREP_DIR}/${SPLIT}_segments.fasta"
        if [ -f "${ORIG}" ] && [ -f "${NEW}" ]; then
            O=$(grep -c '^>' "${ORIG}"); N=$(grep -c '^>' "${NEW}")
            echo "    ${SPLIT}: ${N}/${O} ($((O - N)) removed)"
        fi
    done
done

echo ""
echo "NOTE: rebuild merged 'hard' datasets (step 06) with dev/test pointed at"
echo "the *_derep dirs; keep originals to report the full-vs-dereplicated gap."
echo "Done: $(date)"

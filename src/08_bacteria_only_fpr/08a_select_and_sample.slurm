#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=08a_fpr_select
#SBATCH --output=08a_fpr_select_%j.out
#SBATCH --error=08a_fpr_select_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48g

# Step 08a: pick FPR_N_GENOMES fresh bacteria (HQ-filtered, excluding the step-02
# selection), then plain-sample FPR_TARGET_FRAGMENTS at each of 2k/4k/8k.
# Step 08b will do the prophage check + replacement, step 08c the PHROG-format CSV.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

mkdir -p "${FPR_SELECTION_DIR}" "${FPR_DATA_DIR}"
ACC_FILE="${FPR_SELECTION_DIR}/fpr_accessions.txt"

echo "Step 08a: FPR selection + sampling -- start $(date)"

# ---- 1. select fresh bacteria ----
echo ""
echo "--- selecting fresh bacterial pool ---"
python "${SRC_DIR}/08_bacteria_only_fpr/select_fresh_bacteria.py" \
    --metadata "${GTDB_METADATA}" \
    --exclude-accessions "${GTDB_SELECTION_DIR}/all_accessions.txt" \
    --output "${ACC_FILE}" \
    --n-genomes "${FPR_N_GENOMES}" \
    --min-completeness "${GTDB_MIN_COMPLETENESS}" \
    --max-contamination "${GTDB_MAX_CONTAMINATION}" \
    --seed "${SEED}"
[ -s "${ACC_FILE}" ] || { echo "ERROR: empty FPR accession list"; exit 1; }

# ---- 2. plain-sample fragments at 2k/4k/8k ----
SEG_LENS=(2000 4000 8000)
LEN_NAMES=(2k 4k 8k)
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"; LEN="${LEN_NAMES[$idx]}"
    OUTDIR="${FPR_DATA_DIR}/${LEN}"
    mkdir -p "${OUTDIR}"
    MIN_CONTIG=$(( 2 * SEG )); [ "${MIN_CONTIG}" -lt 10000 ] && MIN_CONTIG=10000
    echo ""
    echo "--- ${SEG} bp -> ${OUTDIR} (min contig ${MIN_CONTIG}, target ${FPR_TARGET_FRAGMENTS}) ---"
    python "${SRC_DIR}/04_sample_bacteria_fragments/subsample_gtdb_segments.py" \
        --gtdb-dir "${GTDB_GENOMES}" \
        --accession-list "${ACC_FILE}" \
        --output "${OUTDIR}/test_segments.fasta" \
        --output-metadata "${OUTDIR}/test_segments.tsv" \
        --segment-length "${SEG}" \
        --min-contig-length "${MIN_CONTIG}" \
        --total-segments "${FPR_TARGET_FRAGMENTS}" \
        --threads "${THREADS}" \
        --seed "${SEED}"
    N=$(grep -c '^>' "${OUTDIR}/test_segments.fasta")
    echo "  ${LEN}: ${N} fragments (target ${FPR_TARGET_FRAGMENTS})"
done
echo "Done: $(date)"

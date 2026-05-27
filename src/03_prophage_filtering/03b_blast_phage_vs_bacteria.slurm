#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03b_phage_vs_bacteria
#SBATCH --output=03b_phage_vs_bacteria_%A_%a.out
#SBATCH --error=03b_phage_vs_bacteria_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --array=0-9

# Step 03b: BLAST INPHARED phages (query) against the selected-bacteria DB.
# Parallelized 10 ways over the phage FASTA (task 0 splits it, others wait).
# Each task writes raw hits for its chunk; 03c combines/filters/masks/analyzes.
# INPHARED-only by design (curated; no metagenomic viral DBs); no learned model.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

NUM_SPLITS=10                 # must match --array above
SPLIT_DIR="${PROPHAGE_BLAST_DIR}/splits"
SPLIT_DONE="${PROPHAGE_BLAST_DIR}/.split_done"
# Full coordinates: subject (bacterial) coords for masking; query (phage) coords
# + qlen for the reverse-contamination analysis; slen for clamping.
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
mkdir -p "${SPLIT_DIR}"

ls "${BACTERIA_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: bacteria DB missing (run 03a first)"; exit 1; }

# Task 0 splits the phage FASTA into NUM_SPLITS chunks; others wait.
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    echo "Task 0: splitting $(grep -c '^>' "${INPHARED_FASTA}") phage sequences into ${NUM_SPLITS} chunks..."
    NUM_SEQS=$(grep -c "^>" "${INPHARED_FASTA}")
    PER=$(( (NUM_SEQS + NUM_SPLITS - 1) / NUM_SPLITS ))
    awk -v n="${PER}" -v dir="${SPLIT_DIR}" '
        /^>/ { if (c >= n) { f++; c = 0 } c++ }
        { print > (dir "/split_" f ".fasta") }
    ' "${INPHARED_FASTA}"
    touch "${SPLIT_DONE}"
else
    while [ ! -f "${SPLIT_DONE}" ]; do sleep 10; done
fi

SPLIT_FILE="${SPLIT_DIR}/split_${SLURM_ARRAY_TASK_ID}.fasta"
RAW_OUT="${PROPHAGE_BLAST_DIR}/blast_raw_${SLURM_ARRAY_TASK_ID}.tsv"
[ -f "${SPLIT_FILE}" ] || { echo "No split file for task ${SLURM_ARRAY_TASK_ID}; nothing to do."; exit 0; }

echo "BLASTing $(grep -c '^>' "${SPLIT_FILE}") phages vs ${BACTERIA_BLASTDB} (task ${SLURM_ARRAY_TASK_ID})..."
blastn \
    -db "${BACTERIA_BLASTDB}" \
    -query "${SPLIT_FILE}" \
    -out "${RAW_OUT}" \
    -outfmt "${OUTFMT}" \
    -max_target_seqs "${PROPHAGE_MAX_TARGET_SEQS}" \
    -evalue "${PROPHAGE_EVALUE}" \
    -num_threads "${THREADS}"

echo "  task ${SLURM_ARRAY_TASK_ID}: $(wc -l < "${RAW_OUT}") raw hits. Done: $(date)"

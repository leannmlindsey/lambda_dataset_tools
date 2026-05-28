#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03b_bacteria_vs_phage
#SBATCH --output=03b_bacteria_vs_phage_%A_%a.out
#SBATCH --error=03b_bacteria_vs_phage_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --array=0-9

# Step 03b: BLAST selected bacterial genomes (QUERY) against the INPHARED DB
# (subject). This direction is chosen so that -max_target_seqs caps per-query
# (per bacterial contig) rather than per-phage -- a widely-integrated temperate
# phage with thousands of lysogens never exceeds the cap, so masking is complete.
#
# Task 0 splits the combined bacterial FASTA into NUM_SPLITS chunks; the other
# tasks wait. Each task writes raw hits for its chunk; 03c combines/filters/
# masks/analyzes.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

NUM_SPLITS=10
SPLIT_DIR="${PROPHAGE_BLAST_DIR}/splits"
SPLIT_DONE="${PROPHAGE_BLAST_DIR}/.split_done"
# Full coords: bacterial coords on the query side (for masking); phage coords on
# the subject side + slen (for the reverse-contamination analysis).
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
mkdir -p "${SPLIT_DIR}"

ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: phage DB missing (run 03a first)"; exit 1; }
[ -s "${BACTERIA_FASTA_COMBINED}" ] || { echo "ERROR: combined bacterial FASTA missing (run 03a first)"; exit 1; }

if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    echo "Task 0: splitting $(grep -c '^>' "${BACTERIA_FASTA_COMBINED}") bacterial sequences into ${NUM_SPLITS} chunks..."
    NUM_SEQS=$(grep -c "^>" "${BACTERIA_FASTA_COMBINED}")
    PER=$(( (NUM_SEQS + NUM_SPLITS - 1) / NUM_SPLITS ))
    awk -v n="${PER}" -v dir="${SPLIT_DIR}" '
        /^>/ { if (c >= n) { f++; c = 0 } c++ }
        { print > (dir "/split_" f ".fasta") }
    ' "${BACTERIA_FASTA_COMBINED}"
    touch "${SPLIT_DONE}"
else
    while [ ! -f "${SPLIT_DONE}" ]; do sleep 10; done
fi

SPLIT_FILE="${SPLIT_DIR}/split_${SLURM_ARRAY_TASK_ID}.fasta"
RAW_OUT="${PROPHAGE_BLAST_DIR}/blast_raw_${SLURM_ARRAY_TASK_ID}.tsv"
[ -f "${SPLIT_FILE}" ] || { echo "No split file for task ${SLURM_ARRAY_TASK_ID}; nothing to do."; exit 0; }

echo "BLASTing $(grep -c '^>' "${SPLIT_FILE}") bacterial sequences vs ${PHAGE_BLASTDB} (task ${SLURM_ARRAY_TASK_ID})..."
blastn \
    -db "${PHAGE_BLASTDB}" \
    -query "${SPLIT_FILE}" \
    -out "${RAW_OUT}" \
    -outfmt "${OUTFMT}" \
    -max_target_seqs "${PROPHAGE_MAX_TARGET_SEQS}" \
    -evalue "${PROPHAGE_EVALUE}" \
    -num_threads "${THREADS}"

echo "  task ${SLURM_ARRAY_TASK_ID}: $(wc -l < "${RAW_OUT}") raw hits. Done: $(date)"

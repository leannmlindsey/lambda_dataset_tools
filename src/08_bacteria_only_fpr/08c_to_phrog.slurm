#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=08c_to_phrog
#SBATCH --output=08c_to_phrog_%j.out
#SBATCH --error=08c_to_phrog_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g

# Step 08c: convert the checked FPR bacterial fragments to PHROG-style CSVs
# (columns: seq_id, start, end, sequence, label) with label = 0 (bacteria).
# Output: ${FPR_DATA_DIR}/{2k,4k,8k}/bacteria_segments_{2k,4k,8k}.csv

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

LEN_NAMES=(2k 4k 8k)
echo "Step 08c: convert checked FPR bacterial fragments to PHROG-style CSV -- start $(date)"
for LEN in "${LEN_NAMES[@]}"; do
    FASTA="${FPR_DATA_DIR}/${LEN}/test_segments_checked.fasta"
    META="${FPR_DATA_DIR}/${LEN}/test_segments_checked.tsv"
    OUT="${FPR_DATA_DIR}/${LEN}/bacteria_segments_${LEN}.csv"
    [ -f "${FASTA}" ] || { echo "ERROR: missing ${FASTA} (run 08b first)"; exit 1; }
    [ -f "${META}" ]  || { echo "ERROR: missing ${META}"; exit 1; }
    python "${SRC_DIR}/08_bacteria_only_fpr/to_phrog_csv.py" \
        --fasta "${FASTA}" --metadata "${META}" --output "${OUT}" --label 0
    N=$(( $(wc -l < "${OUT}") - 1 ))
    [ "${N}" -gt 0 ] || { echo "ERROR: ${OUT} has no data rows"; exit 1; }
    echo "  ${LEN}: ${N} rows in ${OUT}"
done
echo "Done: $(date)"

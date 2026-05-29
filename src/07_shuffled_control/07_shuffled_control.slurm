#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=07_shuffled
#SBATCH --output=07_shuffled_%j.out
#SBATCH --error=07_shuffled_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g

# Step 07: per-row nucleotide shuffle of each merged TEST CSV (length, GC,
# labels preserved; k-mer / sequence structure destroyed). Used as a control
# in evaluation. Only the test split is shuffled (per spec).

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

LEN_NAMES=(2k 4k 8k)

echo "Step 07: shuffled-nucleotide control -- start $(date)"
for L in "${LEN_NAMES[@]}"; do
    IN="${MERGED_DIR}/${L}/test.csv"
    OUT="${SHUFFLED_DIR}/${L}/test_shuffled.csv"
    [ -f "${IN}" ] || { echo "ERROR: missing ${IN} (run step 06 first)"; exit 1; }
    mkdir -p "$(dirname "${OUT}")"
    echo "=== ${L}  ${IN} -> ${OUT} ==="
    python "${SRC_DIR}/07_shuffled_control/shuffle_segments.py" \
        --input-csv "${IN}" --output-csv "${OUT}" --seed "${SEED}"
    N=$(( $(wc -l < "${OUT}") - 1 ))
    [ "${N}" -gt 0 ] || { echo "ERROR: ${OUT} has no data rows"; exit 1; }
    echo "  ${N} rows"
done
echo "Done: $(date)"

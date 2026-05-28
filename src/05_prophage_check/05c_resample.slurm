#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=05c_resample
#SBATCH --output=05c_resample_%j.out
#SBATCH --error=05c_resample_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 05c: drop hit fragments, sample replacements from non-phage regions of
# the affected genomes, and write `${SPLIT}_segments_checked.fasta` (+ .tsv)
# alongside step 04's originals.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

MASKS="${PROPHAGE_CHECK_DIR}/affected_masks.tsv"
[ -f "${MASKS}" ] || { echo "ERROR: ${MASKS} missing (run 05b first)"; exit 1; }

BACT_DIRS_CSV="${BACTERIA_DATASETS[0]},${BACTERIA_DATASETS[1]},${BACTERIA_DATASETS[2]}"

echo "Step 05c: drop hit fragments + sample replacements -- start $(date)"
python "${SRC_DIR}/05_prophage_check/resample_replacements.py" \
    --gtdb-dir "${GTDB_GENOMES}" \
    --fragment-dirs "${BACT_DIRS_CSV}" \
    --seg-lengths "2000,4000,8000" \
    --splits "train,dev,test" \
    --hits-input-dir "${PROPHAGE_CHECK_DIR}" \
    --masks "${MASKS}" \
    --seed "${SEED}"

# Output validation: each *_checked.fasta should exist and be non-empty
for idx in 0 1 2; do
    OUTDIR="${BACTERIA_DATASETS[$idx]}"
    for SPLIT in train dev test; do
        F="${OUTDIR}/${SPLIT}_segments_checked.fasta"
        [ -s "${F}" ] || { echo "ERROR: missing/empty ${F}"; exit 1; }
        N=$(grep -c '^>' "${F}")
        echo "  ${idx} ${SPLIT}: ${N} fragments in ${F}"
    done
done

echo "Done: $(date)"

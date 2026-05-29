#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=06_merge
#SBATCH --output=06_merge_%j.out
#SBATCH --error=06_merge_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 06: merge phage + checked-bacteria fragments into 1:1 train/val/test CSVs.
# Reads:  phage from PHAGE_DATASETS, bacteria from BACTERIA_DATASETS
#         (the `*_segments_checked.fasta` files produced by step 05c).
# Writes: ${MERGED_DIR}/{2k,4k,8k}/{train,dev,test}.csv
#         columns: segment_id, sequence, label (1=phage,0=bacteria), source

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

SEG_LENS=(2000 4000 8000)
LEN_NAMES=(2k 4k 8k)
SPLITS=(train dev test)

echo "Step 06: merge fragments -- start $(date)"
for idx in 0 1 2; do
    LEN_NAME="${LEN_NAMES[$idx]}"
    OUTDIR="${MERGED_DIR}/${LEN_NAME}"
    mkdir -p "${OUTDIR}"
    echo "=== ${LEN_NAME} -> ${OUTDIR} ==="
    for SPLIT in "${SPLITS[@]}"; do
        PHAGE="${PHAGE_DATASETS[$idx]}/${SPLIT}_segments.fasta"
        BACT="${BACTERIA_DATASETS[$idx]}/${SPLIT}_segments_checked.fasta"
        OUT="${OUTDIR}/${SPLIT}.csv"
        [ -f "${PHAGE}" ] || { echo "ERROR: missing phage FASTA: ${PHAGE} (run step 03)"; exit 1; }
        [ -f "${BACT}" ]  || { echo "ERROR: missing checked-bacteria FASTA: ${BACT} (run step 05c)"; exit 1; }
        python "${SRC_DIR}/06_merge/merge_fragments.py" \
            --phage-fasta "${PHAGE}" \
            --bacteria-fasta "${BACT}" \
            --output "${OUT}" \
            --seed "${SEED}"
        N=$(( $(wc -l < "${OUT}") - 1 ))
        [ "${N}" -gt 0 ] || { echo "ERROR: ${OUT} has no data rows"; exit 1; }
        echo "  ${SPLIT}.csv: ${N} rows"
    done
done

echo "Done: $(date)"

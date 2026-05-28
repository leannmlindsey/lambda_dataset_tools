#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=04a_subsample_phage
#SBATCH --output=04a_subsample_phage_%j.out
#SBATCH --error=04a_subsample_phage_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 04a: subsample PHAGE segments from INPHARED genomes at 2k/4k/8k.
# Phage uses the length-proportional rule (~1 segment per SAMPLE_PER_BP bp) to
# represent the wide phage size range (the asymmetry vs bacteria is deliberate —
# 50:50 balance is achieved by matching the bacterial TOTAL to the phage total).
# No masking on the phage side (host-like content retained; see step 03 notes).

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

ACC_DIR="${PHAGE_DATASETS[0]}"          # step 01 wrote accession lists here
SEG_LENS=(2000 4000 8000)
MIN_LENS=("${MIN_LENGTH_2K}" "${MIN_LENGTH_4K}" "${MIN_LENGTH_8K}")

echo "Step 04a: phage subsampling — start $(date)"
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"; MINLEN="${MIN_LENS[$idx]}"; OUTDIR="${PHAGE_DATASETS[$idx]}"
    mkdir -p "${OUTDIR}"
    echo "=== ${SEG} bp -> ${OUTDIR} (min genome ${MINLEN}) ==="
    for SPLIT in train dev test; do
        ACC="${ACC_DIR}/${SPLIT}_accessions.txt"
        [ -f "${ACC}" ] || { echo "ERROR: missing ${ACC} (run step 01 first)"; exit 1; }
        python "${SRC_DIR}/04_subsampling/subsample_segments.py" \
            --input-fasta "${INPHARED_FASTA}" \
            --accession-list "${ACC}" \
            --output "${OUTDIR}/${SPLIT}_segments.fasta" \
            --output-metadata "${OUTDIR}/${SPLIT}_segments.tsv" \
            --segment-length "${SEG}" \
            --min-length "${MINLEN}" \
            --sample-per-bp "${SAMPLE_PER_BP}" \
            --seed "${SEED}"
        N=$(grep -c '^>' "${OUTDIR}/${SPLIT}_segments.fasta" 2>/dev/null || echo 0)
        [ "${N}" -gt 0 ] || { echo "ERROR: no phage segments for ${SEG} ${SPLIT}"; exit 1; }
        echo "  ${SEG} ${SPLIT}: ${N} segments"
    done
done
echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03_phage_fragments
#SBATCH --output=03_phage_fragments_%j.out
#SBATCH --error=03_phage_fragments_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 3: sample PHAGE fragments from INPHARED at 2k / 4k / 8k for train/dev/test.
# This step SETS THE COUNT — step 4 (bacterial fragment sampling) will be told
# to match the per-split fragment count produced here.
# No filtering / no masking on the phage side (deliberate: phage signal like
# integrases that resemble bacterial genes must be preserved).

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

ACC_DIR="${PHAGE_DATASETS[0]}"          # step 01 wrote train/dev/test_accessions.txt here
SEG_LENS=(2000 4000 8000)
MIN_LENS=("${MIN_LENGTH_2K}" "${MIN_LENGTH_4K}" "${MIN_LENGTH_8K}")

echo "Step 3: phage fragment sampling — start $(date)"
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"
    MINLEN="${MIN_LENS[$idx]}"
    OUTDIR="${PHAGE_DATASETS[$idx]}"
    mkdir -p "${OUTDIR}"
    echo "=== ${SEG} bp -> ${OUTDIR} (min genome ${MINLEN}; 1 fragment per ${SAMPLE_PER_BP} bp) ==="
    for SPLIT in train dev test; do
        ACC="${ACC_DIR}/${SPLIT}_accessions.txt"
        [ -f "${ACC}" ] || { echo "ERROR: missing ${ACC} (run step 01 first)"; exit 1; }
        python "${SRC_DIR}/03_sample_phage_fragments/subsample_segments.py" \
            --input-fasta "${INPHARED_FASTA}" \
            --accession-list "${ACC}" \
            --output "${OUTDIR}/${SPLIT}_segments.fasta" \
            --output-metadata "${OUTDIR}/${SPLIT}_segments.tsv" \
            --segment-length "${SEG}" \
            --min-length "${MINLEN}" \
            --sample-per-bp "${SAMPLE_PER_BP}" \
            --seed "${SEED}"
        N=$(grep -c '^>' "${OUTDIR}/${SPLIT}_segments.fasta" 2>/dev/null || echo 0)
        [ "${N}" -gt 0 ] || { echo "ERROR: no phage fragments produced for ${SEG} ${SPLIT}"; exit 1; }
        echo "  ${SEG} ${SPLIT}: ${N} fragments"
    done
done

# Report the per-split counts step 4 must match (at the 2k length; step 4 will
# read these to know how many bacterial fragments to sample per split).
echo ""
echo "Per-split phage fragment counts (step 4 must match these):"
for SPLIT in train dev test; do
    N=$(grep -c '^>' "${PHAGE_DATASETS[0]}/${SPLIT}_segments.fasta" 2>/dev/null || echo 0)
    echo "  ${SPLIT}: ${N}"
done

echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=04b_subsample_bacteria
#SBATCH --output=04b_subsample_bacteria_%j.out
#SBATCH --error=04b_subsample_bacteria_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Step 04b: subsample BACTERIAL segments from the step-02 selection at 2k/4k/8k,
# mask-aware. Samples from all contigs >= MIN_CONTIG (=max(10kb, 2*seg)), excludes
# the step-03 masked (phage-similar) regions, never crosses a masked region, and
# hits the per-split TOTAL (matched to phage for 50:50). Plasmids retained.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

[ -s "${PROPHAGE_MASK_INTERVALS}" ] || { echo "ERROR: masked intervals missing (run step 03 first): ${PROPHAGE_MASK_INTERVALS}"; exit 1; }

SEG_LENS=(2000 4000 8000)
SPLITS=(train dev test)
TARGETS=("${TARGET_SEGMENTS_TRAIN}" "${TARGET_SEGMENTS_DEV}" "${TARGET_SEGMENTS_TEST}")

echo "Step 04b: bacterial subsampling — start $(date)"
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"; OUTDIR="${BACTERIA_DATASETS[$idx]}"
    MIN_CONTIG=$(( 2 * SEG )); [ "${MIN_CONTIG}" -lt 10000 ] && MIN_CONTIG=10000   # 2k/4k->10k, 8k->16k
    mkdir -p "${OUTDIR}"
    echo "=== ${SEG} bp -> ${OUTDIR} (min contig ${MIN_CONTIG}) ==="
    for sidx in 0 1 2; do
        SPLIT="${SPLITS[$sidx]}"; TARGET="${TARGETS[$sidx]}"
        ACC="${GTDB_SELECTION_DIR}/${SPLIT}_accessions.txt"
        [ -f "${ACC}" ] || { echo "ERROR: missing ${ACC} (run step 02 first)"; exit 1; }
        python "${SRC_DIR}/04_subsampling/subsample_gtdb_segments.py" \
            --gtdb-dir "${GTDB_GENOMES}" \
            --accession-list "${ACC}" \
            --masked-intervals "${PROPHAGE_MASK_INTERVALS}" \
            --output "${OUTDIR}/${SPLIT}_segments.fasta" \
            --output-metadata "${OUTDIR}/${SPLIT}_segments.tsv" \
            --segment-length "${SEG}" \
            --min-contig-length "${MIN_CONTIG}" \
            --total-segments "${TARGET}" \
            --threads "${THREADS}" \
            --seed "${SEED}"
        N=$(grep -c '^>' "${OUTDIR}/${SPLIT}_segments.fasta" 2>/dev/null || echo 0)
        [ "${N}" -gt 0 ] || { echo "ERROR: no bacterial segments for ${SEG} ${SPLIT}"; exit 1; }
        echo "  ${SEG} ${SPLIT}: ${N} segments (target ${TARGET})"
    done
done
echo "Done: $(date)"

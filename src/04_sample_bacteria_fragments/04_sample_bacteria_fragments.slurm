#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=04_bacteria_fragments
#SBATCH --output=04_bacteria_fragments_%j.out
#SBATCH --error=04_bacteria_fragments_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Step 4: sample BACTERIAL fragments from the step-02 selection at 2k/4k/8k.
#
# Per the user's spec:
# - This step is PLAIN sampling (no mask, no per-fragment filter).
# - The per-split target count is taken LIVE from step 03's output so the
#   bacterial fragment count matches the phage fragment count for that split,
#   length-by-length and split-by-split. That's the 1:1 balance.
# - Fragments are sampled from ALL contigs >= max(10 kb, 2*seg) -- not just the
#   longest contig -- so plasmids and secondary contigs participate.
# - Reproducible (hashlib-seeded per-accession RNG; same SEED as everything else).
#
# Step 5 will BLAST these fragments against INPHARED at 70%/200 and replace any
# that hit. This step does NOT do that check.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

SEG_LENS=(2000 4000 8000)
SPLITS=(train dev test)

echo "Step 4: bacterial fragment sampling -- start $(date)"
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"
    OUTDIR="${BACTERIA_DATASETS[$idx]}"
    PHAGE_DIR="${PHAGE_DATASETS[$idx]}"
    # min contig length floor: max(10 kb, 2*seg) -> 2k/4k -> 10k, 8k -> 16k
    MIN_CONTIG=$(( 2 * SEG )); [ "${MIN_CONTIG}" -lt 10000 ] && MIN_CONTIG=10000
    mkdir -p "${OUTDIR}"
    echo "=== ${SEG} bp -> ${OUTDIR} (min contig ${MIN_CONTIG}) ==="
    for SPLIT in "${SPLITS[@]}"; do
        ACC="${GTDB_SELECTION_DIR}/${SPLIT}_accessions.txt"
        PHAGE_FASTA="${PHAGE_DIR}/${SPLIT}_segments.fasta"
        [ -f "${ACC}" ]         || { echo "ERROR: missing ${ACC} (run step 02 first)"; exit 1; }
        [ -f "${PHAGE_FASTA}" ] || { echo "ERROR: missing ${PHAGE_FASTA} (run step 03 first)"; exit 1; }

        # *** Target = phage fragment count for the same length+split ***
        TARGET=$(grep -c '^>' "${PHAGE_FASTA}")
        [ "${TARGET}" -gt 0 ]   || { echo "ERROR: phage FASTA is empty: ${PHAGE_FASTA}"; exit 1; }
        echo "  ${SEG} ${SPLIT}: matching phage count = ${TARGET}"

        python "${SRC_DIR}/04_sample_bacteria_fragments/subsample_gtdb_segments.py" \
            --gtdb-dir "${GTDB_GENOMES}" \
            --accession-list "${ACC}" \
            --output "${OUTDIR}/${SPLIT}_segments.fasta" \
            --output-metadata "${OUTDIR}/${SPLIT}_segments.tsv" \
            --segment-length "${SEG}" \
            --min-contig-length "${MIN_CONTIG}" \
            --total-segments "${TARGET}" \
            --threads "${THREADS}" \
            --seed "${SEED}"
        # NOTE: no --masked-intervals -> plain sampling.

        N=$(grep -c '^>' "${OUTDIR}/${SPLIT}_segments.fasta" 2>/dev/null || echo 0)
        [ "${N}" -gt 0 ] || { echo "ERROR: no bacterial fragments for ${SEG} ${SPLIT}"; exit 1; }
        echo "    -> produced ${N} fragments (target ${TARGET})"
    done
done
echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03c_postprocess
#SBATCH --output=03c_postprocess_%j.out
#SBATCH --error=03c_postprocess_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 03c (runs after the 03b array): combine raw BLAST hits, filter to the
# masking threshold, build per-contig masked intervals for step 04, and run the
# threshold-sensitivity + reverse-contamination analyses.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

RAW="${PROPHAGE_BLAST_DIR}/blast_raw_combined.tsv"

echo "=========================================="
echo "Step 03c: postprocess phage-vs-bacteria BLAST"
echo "  thresholds: >= ${PROPHAGE_MIN_IDENTITY}% id, >= ${PROPHAGE_MIN_ALIGN_LEN} bp; pad ${PROPHAGE_MASK_PAD} bp"
echo "  start:      $(date)"
echo "=========================================="

# 1. Combine array outputs
shopt -s nullglob
RAW_PARTS=("${PROPHAGE_BLAST_DIR}"/blast_raw_[0-9]*.tsv)
[ ${#RAW_PARTS[@]} -gt 0 ] || { echo "ERROR: no blast_raw_*.tsv found (did 03b run?)"; exit 1; }
cat "${RAW_PARTS[@]}" > "${RAW}"
echo "  combined raw hits: $(wc -l < "${RAW}")"

# 2. Filter to the masking threshold (col3 = pident, col4 = alignment length)
awk -v id="${PROPHAGE_MIN_IDENTITY}" -v ln="${PROPHAGE_MIN_ALIGN_LEN}" \
    '$3 >= id && $4 >= ln' "${RAW}" > "${PROPHAGE_HITS}"
echo "  filtered hits (>= ${PROPHAGE_MIN_IDENTITY}% / ${PROPHAGE_MIN_ALIGN_LEN} bp): $(wc -l < "${PROPHAGE_HITS}")"

# 3. Build per-contig masked intervals (padded, merged) for step 04
python "${SRC_DIR}/03_prophage_filtering/make_masking_intervals.py" \
    --blast-hits "${PROPHAGE_HITS}" \
    --output "${PROPHAGE_MASK_INTERVALS}" \
    --min-identity "${PROPHAGE_MIN_IDENTITY}" \
    --min-length "${PROPHAGE_MIN_ALIGN_LEN}" \
    --pad "${PROPHAGE_MASK_PAD}"

# 4. Descriptive analyses (no filtering of the data): threshold sweep + reverse contamination
python "${SRC_DIR}/03_prophage_filtering/analyze_blast_hits.py" \
    --blast-hits "${RAW}" \
    --rev-identity "${PROPHAGE_MIN_IDENTITY}" \
    --rev-length "${PROPHAGE_MIN_ALIGN_LEN}" \
    --output "${PROPHAGE_BLAST_DIR}/blast_analysis_report.txt"

# Output validation
[ -s "${PROPHAGE_MASK_INTERVALS}" ] || { echo "ERROR: masked intervals file is empty: ${PROPHAGE_MASK_INTERVALS}"; exit 1; }
echo "  masked-interval rows: $(( $(wc -l < "${PROPHAGE_MASK_INTERVALS}") - 1 ))"
echo "Done: $(date)"

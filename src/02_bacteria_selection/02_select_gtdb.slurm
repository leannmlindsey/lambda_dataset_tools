#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=02_bacteria_select
#SBATCH --output=02_bacteria_select_%j.out
#SBATCH --error=02_bacteria_select_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48g

# Step 02: select GTDB bacterial representatives.
#   - HQ filter: CheckM2 completeness >= GTDB_MIN_COMPLETENESS, contamination <= GTDB_MAX_CONTAMINATION
#   - one random genome per genus, ALL genera (~15,865) — no contiguity/completeness restriction
#   - split train/dev/test by genus (no genus spans splits)
# Plasmids are retained; prophage is handled later by alignment (steps 03/04).
# Outputs accession lists + metadata (incl. assembly provenance columns) to
# ${GTDB_SELECTION_DIR}.

set -euo pipefail
# SLURM copies this script to a spool dir, so $0/BASH_SOURCE cannot locate the
# repo. Source config by absolute path (override with LAMBDA_CONFIG if cloned
# elsewhere); config.sh then defines SRC_DIR and all other absolute paths.
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

mkdir -p "${GTDB_SELECTION_DIR}"

echo "=========================================="
echo "Step 02: Bacteria selection"
echo "  metadata: ${GTDB_METADATA}"
echo "  quality:  completeness >= ${GTDB_MIN_COMPLETENESS}, contamination <= ${GTDB_MAX_CONTAMINATION}"
echo "  output:   ${GTDB_SELECTION_DIR}"
echo "  split:    ${SPLIT_RATIO}   seed: ${SEED}"
echo "  start:    $(date)"
echo "=========================================="

python "${SRC_DIR}/02_bacteria_selection/select_gtdb_representatives.py" \
    --metadata "${GTDB_METADATA}" \
    --output "${GTDB_SELECTION_DIR}" \
    --min-completeness "${GTDB_MIN_COMPLETENESS}" \
    --max-contamination "${GTDB_MAX_CONTAMINATION}" \
    --split-ratio "${SPLIT_RATIO}" \
    --seed "${SEED}"

# Output validation — makes afterok meaningful (fail if any split is empty/missing).
for SPLIT in train dev test; do
    f="${GTDB_SELECTION_DIR}/${SPLIT}_accessions.txt"
    [ -s "${f}" ] || { echo "ERROR: empty or missing ${f}"; exit 1; }
    echo "  ${SPLIT}: $(wc -l < "${f}") genera/genomes"
done

echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03a_bacteria_blastdb
#SBATCH --output=03a_bacteria_blastdb_%j.out
#SBATCH --error=03a_bacteria_blastdb_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 03a: build a BLAST nucleotide DB from the step-02 SELECTED bacterial
# genomes (all_accessions.txt). The phages are then BLASTed against this DB in
# 03b. Skips rebuilding if the DB already exists.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

ACCESSION_LIST="${GTDB_SELECTION_DIR}/all_accessions.txt"
DB_DIR="$(dirname "${BACTERIA_BLASTDB}")"
COMBINED_FASTA="${DB_DIR}/gtdb_selected_genomes.fna"
mkdir -p "${DB_DIR}"

echo "=========================================="
echo "Step 03a: build selected-bacteria BLAST DB"
echo "  genomes:    ${GTDB_GENOMES}"
echo "  accessions: ${ACCESSION_LIST}"
echo "  db:         ${BACTERIA_BLASTDB}"
echo "  start:      $(date)"
echo "=========================================="

[ -f "${ACCESSION_LIST}" ] || { echo "ERROR: accession list not found (run step 02 first): ${ACCESSION_LIST}"; exit 1; }

if ls "${BACTERIA_BLASTDB}".n* >/dev/null 2>&1; then
    echo "BLAST DB already exists at ${BACTERIA_BLASTDB} — skipping rebuild."
    exit 0
fi

# GTDB accession -> gzipped genome path (nested GCA/123/456/789 layout)
accession_to_path() {
    local acc=$1
    acc=${acc#GB_}; acc=${acc#RS_}
    local prefix=${acc%%_*}
    local rest=${acc#*_}
    local number=${rest%%.*}
    number=$(printf "%09d" $((10#${number})))
    echo "${GTDB_GENOMES}/database/${prefix}/${number:0:3}/${number:3:3}/${number:6:3}/${acc}_genomic.fna.gz"
}

rm -f "${COMBINED_FASTA}"
N=$(wc -l < "${ACCESSION_LIST}"); COUNTER=0; FOUND=0; NOT_FOUND=0
while IFS= read -r accession; do
    [ -z "$accession" ] && continue
    COUNTER=$((COUNTER + 1))
    p=$(accession_to_path "$accession")
    if [ -f "$p" ]; then
        zcat "$p" >> "${COMBINED_FASTA}"; FOUND=$((FOUND + 1))
    else
        # Try the other RefSeq/GenBank prefix (GCA <-> GCF)
        alt=${accession#GB_}; alt=${alt#RS_}
        if [[ $alt == GCA_* ]]; then alt=${alt/GCA_/GCF_}; else alt=${alt/GCF_/GCA_}; fi
        ap=$(accession_to_path "$alt")
        if [ -f "$ap" ]; then zcat "$ap" >> "${COMBINED_FASTA}"; FOUND=$((FOUND + 1));
        else echo "  Not found: $accession"; NOT_FOUND=$((NOT_FOUND + 1)); fi
    fi
    [ $((COUNTER % 1000)) -eq 0 ] && echo "  processed ${COUNTER}/${N}..."
done < "${ACCESSION_LIST}"

echo "  found ${FOUND}, not found ${NOT_FOUND}; sequences: $(grep -c '^>' "${COMBINED_FASTA}")"

makeblastdb -in "${COMBINED_FASTA}" -dbtype nucl -out "${BACTERIA_BLASTDB}" \
    -title "GTDB selected bacterial genomes (step 02)" -parse_seqids

ls "${BACTERIA_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: makeblastdb produced no DB files"; exit 1; }
echo "Done: $(date)"

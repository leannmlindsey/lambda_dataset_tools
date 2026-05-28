#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=03a_prep_blast
#SBATCH --output=03a_prep_blast_%j.out
#SBATCH --error=03a_prep_blast_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

# Step 03a: prep for the bacteria-as-query BLAST.
# - Build a BLAST nucleotide DB from INPHARED (curated phages) if not present.
# - Gather the step-02 selected bacterial genomes into one combined FASTA, which
#   03b will split and BLAST against the phage DB.
# Both outputs are reused across reruns (skip if present), so iteration on
# downstream steps is cheap.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

ACCESSION_LIST="${GTDB_SELECTION_DIR}/all_accessions.txt"
PHAGE_DB_DIR="$(dirname "${PHAGE_BLASTDB}")"
mkdir -p "${PHAGE_DB_DIR}" "${PROPHAGE_BLAST_DIR}"

echo "=========================================="
echo "Step 03a: prep BLAST (build phage DB + gather bacteria FASTA)"
echo "  INPHARED FASTA: ${INPHARED_FASTA}"
echo "  phage DB:       ${PHAGE_BLASTDB}"
echo "  accessions:     ${ACCESSION_LIST}"
echo "  combined FASTA: ${BACTERIA_FASTA_COMBINED}"
echo "  start:          $(date)"
echo "=========================================="

[ -f "${ACCESSION_LIST}" ] || { echo "ERROR: accession list not found (run step 02 first): ${ACCESSION_LIST}"; exit 1; }
[ -f "${INPHARED_FASTA}" ] || { echo "ERROR: INPHARED FASTA not found: ${INPHARED_FASTA}"; exit 1; }

# ---- 1. INPHARED BLAST DB ----
if ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1; then
    echo "Phage BLAST DB already exists at ${PHAGE_BLASTDB} -- skipping rebuild."
else
    echo "Building phage BLAST DB from $(grep -c '^>' "${INPHARED_FASTA}") sequences..."
    makeblastdb -in "${INPHARED_FASTA}" -dbtype nucl -out "${PHAGE_BLASTDB}" \
        -title "INPHARED phages (curated)" -parse_seqids
    ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: makeblastdb produced no DB files"; exit 1; }
fi

# ---- 2. Combined bacterial FASTA (the BLAST query for 03b) ----
if [ -s "${BACTERIA_FASTA_COMBINED}" ]; then
    echo "Combined bacterial FASTA already present ($(grep -c '^>' "${BACTERIA_FASTA_COMBINED}") seqs) -- skipping gather."
else
    echo "Gathering selected bacterial genomes from ${GTDB_GENOMES} ..."

    # accession -> gzipped genome path (nested GCA/123/456/789 layout)
    accession_to_path() {
        local acc=$1
        acc=${acc#GB_}; acc=${acc#RS_}
        local prefix=${acc%%_*}
        local rest=${acc#*_}
        local number=${rest%%.*}
        number=$(printf "%09d" $((10#${number})))
        echo "${GTDB_GENOMES}/database/${prefix}/${number:0:3}/${number:3:3}/${number:6:3}/${acc}_genomic.fna.gz"
    }

    : > "${BACTERIA_FASTA_COMBINED}"
    N=$(wc -l < "${ACCESSION_LIST}"); COUNTER=0; FOUND=0; NOT_FOUND=0
    while IFS= read -r accession; do
        [ -z "$accession" ] && continue
        COUNTER=$((COUNTER + 1))
        p=$(accession_to_path "$accession")
        if [ -f "$p" ]; then
            zcat "$p" >> "${BACTERIA_FASTA_COMBINED}"; FOUND=$((FOUND + 1))
        else
            alt=${accession#GB_}; alt=${alt#RS_}
            if [[ $alt == GCA_* ]]; then alt=${alt/GCA_/GCF_}; else alt=${alt/GCF_/GCA_}; fi
            ap=$(accession_to_path "$alt")
            if [ -f "$ap" ]; then zcat "$ap" >> "${BACTERIA_FASTA_COMBINED}"; FOUND=$((FOUND + 1));
            else echo "  Not found: $accession"; NOT_FOUND=$((NOT_FOUND + 1)); fi
        fi
        [ $((COUNTER % 1000)) -eq 0 ] && echo "  processed ${COUNTER}/${N}..."
    done < "${ACCESSION_LIST}"

    echo "  found ${FOUND}, not found ${NOT_FOUND}; sequences in combined FASTA: $(grep -c '^>' "${BACTERIA_FASTA_COMBINED}")"
    [ -s "${BACTERIA_FASTA_COMBINED}" ] || { echo "ERROR: combined bacterial FASTA is empty"; exit 1; }
fi

echo "Done: $(date)"

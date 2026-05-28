#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=05b_blast_affected
#SBATCH --output=05b_blast_affected_%j.out
#SBATCH --error=05b_blast_affected_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Step 05b: for the genomes flagged by 05a (their bacterial fragments hit
# INPHARED), gather their full sequences and BLAST them against INPHARED at the
# same threshold. The merged padded hits become per-contig masked intervals that
# step 05c samples around for replacement fragments.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

AFFECTED="${PROPHAGE_CHECK_DIR}/affected_accessions.txt"
COMBINED="${PROPHAGE_CHECK_DIR}/affected_bacteria.fna"
RAW="${PROPHAGE_CHECK_DIR}/affected_blast_raw.tsv"
FILTERED="${PROPHAGE_CHECK_DIR}/affected_blast_filtered.tsv"
MASKS="${PROPHAGE_CHECK_DIR}/affected_masks.tsv"
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"

echo "Step 05b: affected-genome BLAST -- start $(date)"
[ -f "${AFFECTED}" ] || { echo "ERROR: ${AFFECTED} missing (run 05a first)"; exit 1; }
[ -s "${AFFECTED}" ] || { echo "No affected genomes; nothing to do."; exit 0; }
ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: phage DB missing (run 05a first)"; exit 1; }
N=$(wc -l < "${AFFECTED}"); echo "  affected genomes: ${N}"

# ---- 1. Gather affected genomes into one FASTA ----
if [ -s "${COMBINED}" ]; then
    echo "Combined affected FASTA already present ($(grep -c '^>' "${COMBINED}") seqs) -- skipping gather."
else
    echo "Gathering affected genomes from ${GTDB_GENOMES} ..."
    accession_to_path() {
        local acc=$1
        acc=${acc#GB_}; acc=${acc#RS_}
        local prefix=${acc%%_*}; local rest=${acc#*_}; local number=${rest%%.*}
        number=$(printf "%09d" $((10#${number})))
        echo "${GTDB_GENOMES}/database/${prefix}/${number:0:3}/${number:3:3}/${number:6:3}/${acc}_genomic.fna.gz"
    }
    : > "${COMBINED}"
    C=0; FOUND=0; NOT=0
    while IFS= read -r acc; do
        [ -z "$acc" ] && continue
        C=$((C + 1))
        p=$(accession_to_path "$acc")
        if [ -f "$p" ]; then zcat "$p" >> "${COMBINED}"; FOUND=$((FOUND + 1));
        else
            alt=${acc#GB_}; alt=${alt#RS_}
            if [[ $alt == GCA_* ]]; then alt=${alt/GCA_/GCF_}; else alt=${alt/GCF_/GCA_}; fi
            ap=$(accession_to_path "$alt")
            if [ -f "$ap" ]; then zcat "$ap" >> "${COMBINED}"; FOUND=$((FOUND + 1));
            else echo "  Not found: $acc"; NOT=$((NOT + 1)); fi
        fi
        [ $((C % 500)) -eq 0 ] && echo "  processed ${C}/${N}..."
    done < "${AFFECTED}"
    echo "  found ${FOUND}, not found ${NOT}; contigs in combined FASTA: $(grep -c '^>' "${COMBINED}")"
fi

# ---- 2. BLAST combined affected vs INPHARED ----
echo ""
echo "BLASTing affected bacterial sequence vs INPHARED..."
blastn \
    -db "${PHAGE_BLASTDB}" \
    -query "${COMBINED}" \
    -out "${RAW}" \
    -outfmt "${OUTFMT}" \
    -perc_identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
    -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" \
    -evalue "${BLAST_EVALUE}" \
    -num_threads "${THREADS}"
echo "  raw hits: $(wc -l < "${RAW}")"

awk -v id="${PROPHAGE_CHECK_MIN_IDENTITY}" -v ln="${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
    '$3 >= id && $4 >= ln' "${RAW}" > "${FILTERED}"
echo "  filtered hits (>= ${PROPHAGE_CHECK_MIN_IDENTITY}% / ${PROPHAGE_CHECK_MIN_ALIGN_LEN} bp): $(wc -l < "${FILTERED}")"

# ---- 3. Build per-contig masks ----
python "${SRC_DIR}/05_prophage_check/make_per_genome_masks.py" \
    --blast-hits "${FILTERED}" \
    --output "${MASKS}" \
    --min-identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
    --min-length "${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
    --pad "${PROPHAGE_CHECK_MASK_PAD}"

[ -s "${MASKS}" ] || { echo "ERROR: masks file empty: ${MASKS}"; exit 1; }
echo "  mask rows: $(( $(wc -l < "${MASKS}") - 1 ))"
echo "Done: $(date)"

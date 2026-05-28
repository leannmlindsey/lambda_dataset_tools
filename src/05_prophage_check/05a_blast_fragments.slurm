#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=05a_blast_fragments
#SBATCH --output=05a_blast_fragments_%j.out
#SBATCH --error=05a_blast_fragments_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

# Step 05a: BLAST each of the 9 bacterial fragment FASTAs (3 lengths x 3 splits)
# against the INPHARED phage DB, filter hits at >=PROPHAGE_CHECK_MIN_IDENTITY% /
# >=PROPHAGE_CHECK_MIN_ALIGN_LEN bp, and identify which fragments hit and
# which source genomes those fragments came from. The full-genome BLAST for the
# affected genomes happens in 05b; replacement sampling in 05c.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

mkdir -p "${PROPHAGE_CHECK_DIR}"
PHAGE_DB_DIR="$(dirname "${PHAGE_BLASTDB}")"
mkdir -p "${PHAGE_DB_DIR}"

OUTFMT="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
SEG_LENS=(2000 4000 8000)
SPLITS=(train dev test)

echo "Step 05a: fragment BLAST -- start $(date)"
echo "  filter: >= ${PROPHAGE_CHECK_MIN_IDENTITY}% identity AND >= ${PROPHAGE_CHECK_MIN_ALIGN_LEN} bp"

# ---- 1. INPHARED BLAST DB (build if not present) ----
if ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1; then
    echo "Phage BLAST DB already exists at ${PHAGE_BLASTDB} -- reusing."
else
    [ -f "${INPHARED_FASTA}" ] || { echo "ERROR: INPHARED FASTA not found: ${INPHARED_FASTA}"; exit 1; }
    echo "Building phage BLAST DB from $(grep -c '^>' "${INPHARED_FASTA}") sequences..."
    makeblastdb -in "${INPHARED_FASTA}" -dbtype nucl -out "${PHAGE_BLASTDB}" \
        -title "INPHARED phages" -parse_seqids
fi

# ---- 2. BLAST each fragment file ----
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"
    BACT_DIR="${BACTERIA_DATASETS[$idx]}"
    for SPLIT in "${SPLITS[@]}"; do
        QFASTA="${BACT_DIR}/${SPLIT}_segments.fasta"
        [ -f "${QFASTA}" ] || { echo "ERROR: missing ${QFASTA} (run step 04 first)"; exit 1; }
        RAW="${PROPHAGE_CHECK_DIR}/fragment_raw_${SEG}_${SPLIT}.tsv"
        FILTERED="${PROPHAGE_CHECK_DIR}/fragment_hits_${SEG}_${SPLIT}.tsv"

        echo ""
        echo "BLASTing ${SEG} ${SPLIT}: $(grep -c '^>' "${QFASTA}") fragments..."
        blastn \
            -db "${PHAGE_BLASTDB}" \
            -query "${QFASTA}" \
            -out "${RAW}" \
            -outfmt "${OUTFMT}" \
            -perc_identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
            -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" \
            -evalue "${BLAST_EVALUE}" \
            -num_threads "${THREADS}"

        # Filter at >= id AND >= length (perc_identity above already pre-filters at id; safe to filter again)
        awk -v id="${PROPHAGE_CHECK_MIN_IDENTITY}" -v ln="${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
            '$3 >= id && $4 >= ln' "${RAW}" > "${FILTERED}"
        echo "  raw hits:      $(wc -l < "${RAW}")"
        echo "  filtered hits: $(wc -l < "${FILTERED}")"
    done
done

# ---- 3. Identify hit fragments + affected genomes ----
echo ""
echo "Identifying affected fragments + genomes..."
BACT_DIRS_CSV="${BACTERIA_DATASETS[0]},${BACTERIA_DATASETS[1]},${BACTERIA_DATASETS[2]}"
python "${SRC_DIR}/05_prophage_check/find_affected.py" \
    --hits-dir "${PROPHAGE_CHECK_DIR}" \
    --fragment-dirs "${BACT_DIRS_CSV}" \
    --seg-lengths "2000,4000,8000" \
    --splits "train,dev,test" \
    --output-dir "${PROPHAGE_CHECK_DIR}"

AFFECTED="${PROPHAGE_CHECK_DIR}/affected_accessions.txt"
[ -s "${AFFECTED}" ] || echo "  (no affected genomes -- no replacement needed)"

echo "Done: $(date)"

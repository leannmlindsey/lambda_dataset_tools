#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=08b_fpr_check
#SBATCH --output=08b_fpr_check_%j.out
#SBATCH --error=08b_fpr_check_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Step 08b: run the same prophage-check pipeline as step 05 on the FPR
# fragments. BLAST fragments at 70%/200 -> identify hits + affected genomes
# -> gather affected -> BLAST affected genomes -> per-genome masks -> drop hit
# fragments + sample replacements from clean regions. Reuses the step 05
# helpers (find_affected.py, make_per_genome_masks.py, resample_replacements.py).

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

mkdir -p "${FPR_CHECK_DIR}"
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
ls "${PHAGE_BLASTDB}".n* >/dev/null 2>&1 || { echo "ERROR: phage DB missing -- run step 05a first"; exit 1; }

SEG_LENS=(2000 4000 8000)
LEN_NAMES=(2k 4k 8k)

# ---- 1. BLAST each FPR fragment FASTA vs INPHARED ----
echo "Step 08b: FPR prophage check -- start $(date)"
for idx in 0 1 2; do
    SEG="${SEG_LENS[$idx]}"; LEN="${LEN_NAMES[$idx]}"
    Q="${FPR_DATA_DIR}/${LEN}/test_segments.fasta"
    [ -f "${Q}" ] || { echo "ERROR: missing ${Q}; run 08a first"; exit 1; }
    RAW="${FPR_CHECK_DIR}/fragment_raw_${SEG}_test.tsv"
    FILT="${FPR_CHECK_DIR}/fragment_hits_${SEG}_test.tsv"
    echo "--- BLAST ${SEG} test: $(grep -c '^>' "${Q}") fragments ---"
    blastn -db "${PHAGE_BLASTDB}" -query "${Q}" -out "${RAW}" -outfmt "${OUTFMT}" \
        -perc_identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
        -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" \
        -evalue "${BLAST_EVALUE}" -num_threads "${THREADS}"
    awk -v id="${PROPHAGE_CHECK_MIN_IDENTITY}" -v ln="${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
        '$3 >= id && $4 >= ln' "${RAW}" > "${FILT}"
    echo "  raw: $(wc -l < "${RAW}")   filtered: $(wc -l < "${FILT}")"
done

# ---- 2. identify hit fragments + affected genomes ----
echo ""
echo "--- find affected ---"
FPR_DIRS="${FPR_DATA_DIR}/2k,${FPR_DATA_DIR}/4k,${FPR_DATA_DIR}/8k"
python "${SRC_DIR}/05_prophage_check/find_affected.py" \
    --hits-dir "${FPR_CHECK_DIR}" \
    --fragment-dirs "${FPR_DIRS}" \
    --seg-lengths "2000,4000,8000" \
    --splits "test" \
    --output-dir "${FPR_CHECK_DIR}"

AFFECTED="${FPR_CHECK_DIR}/affected_accessions.txt"
COMBINED="${FPR_CHECK_DIR}/affected_bacteria.fna"
RAW_A="${FPR_CHECK_DIR}/affected_blast_raw.tsv"
FILT_A="${FPR_CHECK_DIR}/affected_blast_filtered.tsv"
MASKS="${FPR_CHECK_DIR}/affected_masks.tsv"

if [ ! -s "${AFFECTED}" ]; then
    echo "No affected genomes -- copying step-08a outputs as-is to *_checked.*"
    for LEN in 2k 4k 8k; do
        cp "${FPR_DATA_DIR}/${LEN}/test_segments.fasta" "${FPR_DATA_DIR}/${LEN}/test_segments_checked.fasta"
        cp "${FPR_DATA_DIR}/${LEN}/test_segments.tsv"   "${FPR_DATA_DIR}/${LEN}/test_segments_checked.tsv"
    done
    echo "Done: $(date)"; exit 0
fi
echo "  affected genomes: $(wc -l < "${AFFECTED}")"

# ---- 3. gather affected genomes + BLAST against INPHARED ----
echo ""
echo "--- gather affected genomes ---"
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
    C=$((C+1))
    p=$(accession_to_path "$acc")
    if [ -f "$p" ]; then zcat "$p" >> "${COMBINED}"; FOUND=$((FOUND+1));
    else
        alt=${acc#GB_}; alt=${alt#RS_}
        if [[ $alt == GCA_* ]]; then alt=${alt/GCA_/GCF_}; else alt=${alt/GCF_/GCA_}; fi
        ap=$(accession_to_path "$alt")
        if [ -f "$ap" ]; then zcat "$ap" >> "${COMBINED}"; FOUND=$((FOUND+1));
        else echo "  Not found: $acc"; NOT=$((NOT+1)); fi
    fi
done < "${AFFECTED}"
echo "  found ${FOUND}, not found ${NOT}; contigs gathered: $(grep -c '^>' "${COMBINED}")"

echo ""
echo "--- BLAST affected genomes vs INPHARED ---"
blastn -db "${PHAGE_BLASTDB}" -query "${COMBINED}" -out "${RAW_A}" -outfmt "${OUTFMT}" \
    -perc_identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
    -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" \
    -evalue "${BLAST_EVALUE}" -num_threads "${THREADS}"
awk -v id="${PROPHAGE_CHECK_MIN_IDENTITY}" -v ln="${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
    '$3 >= id && $4 >= ln' "${RAW_A}" > "${FILT_A}"
echo "  raw: $(wc -l < "${RAW_A}")   filtered: $(wc -l < "${FILT_A}")"

# ---- 4. build per-contig masks ----
python "${SRC_DIR}/05_prophage_check/make_per_genome_masks.py" \
    --blast-hits "${FILT_A}" --output "${MASKS}" \
    --min-identity "${PROPHAGE_CHECK_MIN_IDENTITY}" \
    --min-length "${PROPHAGE_CHECK_MIN_ALIGN_LEN}" \
    --pad "${PROPHAGE_CHECK_MASK_PAD}"
[ -s "${MASKS}" ] || { echo "ERROR: masks file empty"; exit 1; }

# ---- 5. resample replacements (writes *_segments_checked.fasta/.tsv) ----
echo ""
echo "--- resample replacements ---"
python "${SRC_DIR}/05_prophage_check/resample_replacements.py" \
    --gtdb-dir "${GTDB_GENOMES}" \
    --fragment-dirs "${FPR_DIRS}" \
    --seg-lengths "2000,4000,8000" \
    --splits "test" \
    --hits-input-dir "${FPR_CHECK_DIR}" \
    --masks "${MASKS}" \
    --seed "${SEED}"

# Output validation
for LEN in 2k 4k 8k; do
    F="${FPR_DATA_DIR}/${LEN}/test_segments_checked.fasta"
    [ -s "${F}" ] || { echo "ERROR: ${F} missing/empty"; exit 1; }
    N=$(grep -c '^>' "${F}")
    echo "  ${LEN}: ${N} fragments in checked FASTA"
done
echo "Done: $(date)"

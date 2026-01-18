#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=validate_removal
#SBATCH --output=validate_removal_%j.out
#SBATCH --error=validate_removal_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Validate that prophage-contaminated genomes were actually removed
#
# This script:
# 1. Creates a BLAST database from the CLEAN (filtered) bacterial genomes
# 2. BLASTs INPHARED phages against this new database
# 3. Compares results - should have very few/no hits
#
# Usage: sbatch validate_prophage_removal.slurm

# ============================================
# CONFIGURATION
# ============================================

# GTDB base directory
GTDB_BASE="/u/llindsey1/llindsey/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226"

# Clean accession list (genomes WITHOUT prophage hits)
CLEAN_ACCESSIONS="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/clean_accessions.txt"

# INPHARED phage genomes (query)
PHAGE_FASTA="/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa"

# Output directory
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/prophage_validation"

# BLAST parameters (same as original)
MIN_ALIGNMENT_LENGTH=200
MIN_IDENTITY=90
NUM_THREADS=16

# ============================================

echo "=========================================="
echo "Validate Prophage Removal"
echo "=========================================="
echo "Start time: $(date)"
echo "Clean accessions: ${CLEAN_ACCESSIONS}"
echo "Phage query: ${PHAGE_FASTA}"
echo "Output directory: ${OUTPUT_DIR}"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# Create output directory
mkdir -p ${OUTPUT_DIR}

# ============================================
# STEP 1: Create combined FASTA from clean genomes
# ============================================
echo ""
echo "Step 1: Creating combined FASTA from clean genomes..."

CLEAN_FASTA="${OUTPUT_DIR}/clean_genomes.fna"
rm -f ${CLEAN_FASTA}

# Function to convert accession to file path
accession_to_path() {
    local acc=$1

    # Strip GB_ or RS_ prefix if present
    acc=${acc#GB_}
    acc=${acc#RS_}

    # Extract prefix (GCA or GCF) and number
    local prefix=${acc%%_*}
    local rest=${acc#*_}
    local number=${rest%%.*}

    # Pad to 9 digits and split into chunks
    number=$(printf "%09d" $((10#${number})))
    local chunk1=${number:0:3}
    local chunk2=${number:3:3}
    local chunk3=${number:6:3}

    echo "${GTDB_BASE}/database/${prefix}/${chunk1}/${chunk2}/${chunk3}/${acc}_genomic.fna.gz"
}

NUM_CLEAN=$(wc -l < ${CLEAN_ACCESSIONS})
echo "  Processing ${NUM_CLEAN} clean genomes..."

COUNTER=0
FOUND=0

while IFS= read -r accession; do
    [ -z "$accession" ] && continue

    COUNTER=$((COUNTER + 1))

    fna_path=$(accession_to_path "$accession")

    if [ -f "$fna_path" ]; then
        zcat "$fna_path" >> ${CLEAN_FASTA}
        FOUND=$((FOUND + 1))
    fi

    if [ $((COUNTER % 1000)) -eq 0 ]; then
        echo "    Processed ${COUNTER}/${NUM_CLEAN}..."
    fi

done < ${CLEAN_ACCESSIONS}

echo "  Genomes found: ${FOUND}"
echo "  Combined FASTA size: $(ls -lh ${CLEAN_FASTA} | awk '{print $5}')"
echo "  Total sequences: $(grep -c '^>' ${CLEAN_FASTA})"

# ============================================
# STEP 2: Create BLAST database
# ============================================
echo ""
echo "Step 2: Creating BLAST database from clean genomes..."

DB_NAME="${OUTPUT_DIR}/clean_bacteria_db"

makeblastdb \
    -in ${CLEAN_FASTA} \
    -dbtype nucl \
    -out ${DB_NAME} \
    -title "Clean Bacterial Genomes (no prophage)" \
    -parse_seqids

if [ $? -ne 0 ]; then
    echo "Error: makeblastdb failed!"
    exit 1
fi

echo "  Database created: ${DB_NAME}"

# ============================================
# STEP 3: BLAST phages against clean database
# ============================================
echo ""
echo "Step 3: BLASTing phages against clean bacterial database..."

BLAST_RAW="${OUTPUT_DIR}/validation_blast_raw.tsv"
BLAST_FILTERED="${OUTPUT_DIR}/validation_blast_filtered.tsv"

blastn \
    -db ${DB_NAME} \
    -query ${PHAGE_FASTA} \
    -out ${BLAST_RAW} \
    -outfmt "6 std" \
    -max_target_seqs 500 \
    -evalue 1e-5 \
    -num_threads ${NUM_THREADS}

echo "  Raw BLAST hits: $(wc -l < ${BLAST_RAW})"

# Filter by alignment length and identity (same thresholds as original)
awk -v min_len=${MIN_ALIGNMENT_LENGTH} -v min_id=${MIN_IDENTITY} \
    '$4 >= min_len && $3 >= min_id' ${BLAST_RAW} > ${BLAST_FILTERED}

echo "  Filtered BLAST hits: $(wc -l < ${BLAST_FILTERED})"

# ============================================
# STEP 4: Analyze results
# ============================================
echo ""
echo "Step 4: Analyzing validation results..."

RAW_HITS=$(wc -l < ${BLAST_RAW})
FILTERED_HITS=$(wc -l < ${BLAST_FILTERED})

echo ""
echo "  Raw BLAST hits (all):        ${RAW_HITS}"
echo "  Filtered hits (≥${MIN_ALIGNMENT_LENGTH}bp, ≥${MIN_IDENTITY}%): ${FILTERED_HITS}"

# Analyze alignment length distribution of RAW hits
echo ""
echo "  Alignment length distribution (raw hits):"
echo "    <100bp:    $(awk '$4 < 100' ${BLAST_RAW} | wc -l)"
echo "    100-199bp: $(awk '$4 >= 100 && $4 < 200' ${BLAST_RAW} | wc -l)"
echo "    200-499bp: $(awk '$4 >= 200 && $4 < 500' ${BLAST_RAW} | wc -l)"
echo "    500-999bp: $(awk '$4 >= 500 && $4 < 1000' ${BLAST_RAW} | wc -l)"
echo "    ≥1000bp:   $(awk '$4 >= 1000' ${BLAST_RAW} | wc -l)"

# Check max alignment length in raw hits
if [ ${RAW_HITS} -gt 0 ]; then
    MAX_LEN=$(cut -f4 ${BLAST_RAW} | sort -n | tail -1)
    echo ""
    echo "  Maximum alignment length in raw hits: ${MAX_LEN} bp"
fi

# Count unique phages and genomes with filtered hits
PHAGES_WITH_HITS=$(cut -f1 ${BLAST_FILTERED} | sort -u | wc -l)
CONTIGS_WITH_HITS=$(cut -f2 ${BLAST_FILTERED} | sort -u | wc -l)

echo ""
echo "  Unique phages with hits (≥${MIN_ALIGNMENT_LENGTH}bp): ${PHAGES_WITH_HITS}"
echo "  Unique contigs with hits (≥${MIN_ALIGNMENT_LENGTH}bp): ${CONTIGS_WITH_HITS}"

# ============================================
# STEP 5: Compare with original results
# ============================================
echo ""
echo "Step 5: Comparison with original BLAST..."

ORIGINAL_BLAST="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/blast_filtered_combined.tsv"
ORIGINAL_HITS=$(wc -l < ${ORIGINAL_BLAST})
ORIGINAL_PHAGES=$(cut -f1 ${ORIGINAL_BLAST} | sort -u | wc -l)
ORIGINAL_CONTIGS=$(cut -f2 ${ORIGINAL_BLAST} | sort -u | wc -l)

NEW_HITS=$(wc -l < ${BLAST_FILTERED})

echo ""
echo "=========================================="
echo "VALIDATION SUMMARY"
echo "=========================================="
echo ""
echo "Original BLAST (all genomes):"
echo "  Total filtered hits:    ${ORIGINAL_HITS}"
echo "  Unique phages with hits: ${ORIGINAL_PHAGES}"
echo "  Unique contigs with hits: ${ORIGINAL_CONTIGS}"
echo ""
echo "Validation BLAST (clean genomes only):"
echo "  Total filtered hits:    ${NEW_HITS}"
echo "  Unique phages with hits: ${PHAGES_WITH_HITS}"
echo "  Unique contigs with hits: ${CONTIGS_WITH_HITS}"
echo ""

echo "Interpretation:"
echo ""
if [ ${NEW_HITS} -eq 0 ]; then
    echo "  SUCCESS: No hits ≥${MIN_ALIGNMENT_LENGTH}bp found in clean genomes!"
    echo "  The filtering worked correctly - all long prophage matches were removed."
else
    echo "  WARNING: ${NEW_HITS} hits ≥${MIN_ALIGNMENT_LENGTH}bp still found."
    echo "  This may indicate issues with the filtering process."
fi
echo ""
echo "  Note: Short hits (<${MIN_ALIGNMENT_LENGTH}bp) are expected and acceptable."
echo "  These represent partial/weak similarities that don't indicate prophages."
echo ""
if [ ${RAW_HITS} -gt 0 ]; then
    echo "  Raw hits found: ${RAW_HITS} (includes short alignments)"
    echo "  Max alignment in clean genomes: ${MAX_LEN} bp"
    if [ ${MAX_LEN} -lt ${MIN_ALIGNMENT_LENGTH} ]; then
        echo "  VERIFIED: Maximum alignment (${MAX_LEN}bp) < threshold (${MIN_ALIGNMENT_LENGTH}bp)"
    fi
fi

# Save summary
SUMMARY="${OUTPUT_DIR}/validation_summary.txt"
cat > ${SUMMARY} << EOF
Prophage Removal Validation Summary
===================================
Date: $(date)

Clean genomes used: ${FOUND}
Original BLAST hits: ${ORIGINAL_HITS}
Validation BLAST hits: ${NEW_HITS}

Original unique phages with hits: ${ORIGINAL_PHAGES}
Validation unique phages with hits: ${PHAGES_WITH_HITS}

Original unique contigs with hits: ${ORIGINAL_CONTIGS}
Validation unique contigs with hits: ${CONTIGS_WITH_HITS}

Filtering thresholds:
  Min alignment length: ${MIN_ALIGNMENT_LENGTH} bp
  Min identity: ${MIN_IDENTITY}%

Output files:
  ${CLEAN_FASTA}
  ${DB_NAME}.*
  ${BLAST_RAW}
  ${BLAST_FILTERED}
EOF

echo ""
echo "Summary saved to: ${SUMMARY}"

# Show top hits if any remain
if [ ${NEW_HITS} -gt 0 ]; then
    echo ""
    echo "Top 10 remaining hits:"
    head -10 ${BLAST_FILTERED}
fi

echo ""
echo "=========================================="
echo "Validation complete!"
echo "End time: $(date)"
echo "=========================================="

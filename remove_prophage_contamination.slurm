#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=remove_prophage
#SBATCH --output=remove_prophage_%j.out
#SBATCH --error=remove_prophage_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# Remove prophage-contaminated bacterial genomes from GTDB dataset
#
# Strategy:
# 1. BLAST INPHARED phages against GTDB selected bacteria
# 2. Filter hits by coverage (80%) and identity (90%)
# 3. Remove bacterial accessions with prophage hits
# 4. Create filtered accession lists for resampling
#
# Usage: sbatch remove_prophage_contamination.slurm

# ============================================
# CONFIGURATION
# ============================================

# INPHARED phage genomes (query)
PHAGE_FASTA="/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa"

# GTDB selected bacteria database (subject)
BACTERIA_DB="/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected"

# Original GTDB accession lists
GTDB_DATASET_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset"

# Output directory
OUTPUT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered"

# Filtering thresholds
MIN_COVERAGE=80      # Minimum query coverage (%)
MIN_IDENTITY=90      # Minimum sequence identity (%)
EVALUE="1e-10"       # E-value threshold

# BLAST parameters
NUM_THREADS=16

# ============================================

echo "=========================================="
echo "Remove Prophage Contamination from GTDB"
echo "=========================================="
echo "Start time: $(date)"
echo "Phage query: ${PHAGE_FASTA}"
echo "Bacteria DB: ${BACTERIA_DB}"
echo "Thresholds: ${MIN_COVERAGE}% coverage, ${MIN_IDENTITY}% identity"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check inputs
if [ ! -f "${PHAGE_FASTA}" ]; then
    echo "Error: Phage FASTA not found: ${PHAGE_FASTA}"
    exit 1
fi

if [ ! -f "${BACTERIA_DB}.nin" ] && [ ! -f "${BACTERIA_DB}.nal" ]; then
    echo "Error: Bacteria database not found: ${BACTERIA_DB}"
    exit 1
fi

# Count input phages
NUM_PHAGES=$(grep -c "^>" ${PHAGE_FASTA})
echo "Input phages: ${NUM_PHAGES}"

# Step 1: BLAST phages against bacteria
echo ""
echo "Step 1: BLASTing phages against bacterial genomes..."
BLAST_RAW="${OUTPUT_DIR}/prophage_blast_raw.tsv"

# Use -qcov_hsp_perc to filter at BLAST level (more efficient)
# Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
blastn \
    -db ${BACTERIA_DB} \
    -query ${PHAGE_FASTA} \
    -out ${BLAST_RAW} \
    -outfmt "6 std" \
    -qcov_hsp_perc ${MIN_COVERAGE} \
    -max_target_seqs 500 \
    -num_threads ${NUM_THREADS}

echo "  Raw BLAST hits: $(wc -l < ${BLAST_RAW})"

# Step 2: Filter by identity (coverage already filtered by -qcov_hsp_perc)
echo ""
echo "Step 2: Filtering hits (≥${MIN_IDENTITY}% identity)..."
BLAST_FILTERED="${OUTPUT_DIR}/prophage_blast_filtered.tsv"

# Filter: column 3 = pident (coverage already handled by BLAST)
awk -v min_id=${MIN_IDENTITY} '$3 >= min_id' ${BLAST_RAW} > ${BLAST_FILTERED}

echo "  Filtered hits: $(wc -l < ${BLAST_FILTERED})"

# Step 3: Extract bacterial accessions with prophage hits
echo ""
echo "Step 3: Extracting contaminated bacterial accessions..."
CONTAMINATED="${OUTPUT_DIR}/contaminated_accessions.txt"

# Column 2 = sseqid (bacterial sequence ID)
# Need to extract just the accession part (e.g., GCA_000005845.2 from longer ID)
cut -f2 ${BLAST_FILTERED} | \
    sed 's/_genomic.*//' | \
    sort -u > ${CONTAMINATED}

NUM_CONTAMINATED=$(wc -l < ${CONTAMINATED})
echo "  Bacterial genomes with prophage hits: ${NUM_CONTAMINATED}"

# Step 4: Get unique phages that hit bacteria
PHAGES_WITH_HITS="${OUTPUT_DIR}/phages_with_bacterial_hits.txt"
cut -f1 ${BLAST_FILTERED} | sort -u > ${PHAGES_WITH_HITS}
NUM_PHAGES_HITTING=$(wc -l < ${PHAGES_WITH_HITS})
echo "  Phages with bacterial hits: ${NUM_PHAGES_HITTING}"

# Step 5: Create filtered accession lists
echo ""
echo "Step 4: Creating filtered accession lists..."

for SPLIT in train dev test all; do
    if [ "${SPLIT}" == "all" ]; then
        ORIG_FILE="${GTDB_DATASET_DIR}/all_accessions.txt"
    else
        ORIG_FILE="${GTDB_DATASET_DIR}/${SPLIT}_accessions.txt"
    fi

    if [ ! -f "${ORIG_FILE}" ]; then
        echo "  Skipping ${SPLIT} - file not found"
        continue
    fi

    FILTERED_FILE="${OUTPUT_DIR}/${SPLIT}_accessions.txt"

    # Remove contaminated accessions
    # Need to handle the GB_/RS_ prefix in GTDB accessions
    grep -v -F -f ${CONTAMINATED} ${ORIG_FILE} > ${FILTERED_FILE}

    ORIG_COUNT=$(wc -l < ${ORIG_FILE})
    FILTERED_COUNT=$(wc -l < ${FILTERED_FILE})
    REMOVED=$((ORIG_COUNT - FILTERED_COUNT))

    echo "  ${SPLIT}: ${ORIG_COUNT} -> ${FILTERED_COUNT} (removed ${REMOVED})"
done

# Step 6: Generate summary
echo ""
echo "Step 5: Generating summary..."
SUMMARY="${OUTPUT_DIR}/prophage_removal_summary.txt"

cat > ${SUMMARY} << EOF
Prophage Contamination Removal Summary
======================================
Date: $(date)

Input:
  Phage sequences: ${NUM_PHAGES}
  Bacterial database: ${BACTERIA_DB}

Filtering thresholds:
  Minimum query coverage: ${MIN_COVERAGE}%
  Minimum identity: ${MIN_IDENTITY}%
  E-value: ${EVALUE}

Results:
  Raw BLAST hits: $(wc -l < ${BLAST_RAW})
  Filtered hits: $(wc -l < ${BLAST_FILTERED})
  Phages with bacterial hits: ${NUM_PHAGES_HITTING}
  Contaminated bacterial genomes: ${NUM_CONTAMINATED}

Accession list changes:
EOF

for SPLIT in train dev test all; do
    ORIG_FILE="${GTDB_DATASET_DIR}/${SPLIT}_accessions.txt"
    FILTERED_FILE="${OUTPUT_DIR}/${SPLIT}_accessions.txt"
    if [ -f "${ORIG_FILE}" ] && [ -f "${FILTERED_FILE}" ]; then
        ORIG_COUNT=$(wc -l < ${ORIG_FILE})
        FILTERED_COUNT=$(wc -l < ${FILTERED_FILE})
        echo "  ${SPLIT}: ${ORIG_COUNT} -> ${FILTERED_COUNT}" >> ${SUMMARY}
    fi
done

cat >> ${SUMMARY} << EOF

Output files:
  ${BLAST_RAW}
  ${BLAST_FILTERED}
  ${CONTAMINATED}
  ${PHAGES_WITH_HITS}
  ${OUTPUT_DIR}/*_accessions.txt

Next steps:
  1. Review the contaminated accessions
  2. Re-run subsampling with filtered accession lists
  3. Update TOTAL_SEGMENTS in GTDB subsampling scripts if needed
EOF

cat ${SUMMARY}

# Show top contaminated bacteria (most hits)
echo ""
echo "Top 20 most contaminated bacterial genomes (by hit count):"
cut -f2 ${BLAST_FILTERED} | sed 's/_genomic.*//' | sort | uniq -c | sort -rn | head -20

echo ""
echo "=========================================="
echo "Prophage removal analysis complete!"
echo "End time: $(date)"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="

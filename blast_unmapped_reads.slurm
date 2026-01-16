#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=blast_unmapped
#SBATCH --output=blast_unmapped_%j.out
#SBATCH --error=blast_unmapped_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# BLAST unmapped reads against phage and bacteria databases
# Subsamples reads from each FASTQ file for quick overview
#
# Usage: sbatch blast_unmapped_reads.slurm

# ============================================
# CONFIGURATION
# ============================================

# Input directory with FASTQ.gz files
INPUT_DIR="/work/hdd/bfzj/llindsey1/Noahs_unmapped_data"

# Output directory for results
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/unmapped_blast_results"

# BLAST databases
PHAGE_DB="/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages"
BACTERIA_DB="/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected"

# Subsampling: number of reads to take from EACH file
READS_PER_FILE=10000

# BLAST parameters
EVALUE="1e-10"
NUM_THREADS=16
MAX_TARGET_SEQS=1

# ============================================

echo "=========================================="
echo "BLAST Unmapped Reads Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Reads per file: ${READS_PER_FILE}"
echo "E-value threshold: ${EVALUE}"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check databases exist
if [ ! -f "${PHAGE_DB}.nin" ] && [ ! -f "${PHAGE_DB}.nal" ]; then
    echo "Warning: Phage database not found: ${PHAGE_DB}"
    echo "Run create_phage_blastdb.slurm first"
fi

if [ ! -f "${BACTERIA_DB}.nin" ] && [ ! -f "${BACTERIA_DB}.nal" ]; then
    echo "Warning: Bacteria database not found: ${BACTERIA_DB}"
    echo "Run create_gtdb_selected_blastdb.slurm first"
fi

# Step 1: Subsample and convert FASTQ to FASTA
echo ""
echo "Step 1: Subsampling reads from each FASTQ file..."
COMBINED_FASTA="${OUTPUT_DIR}/subsampled_reads.fasta"
rm -f ${COMBINED_FASTA}

TOTAL_READS=0
FILE_COUNT=0

for FASTQ_GZ in ${INPUT_DIR}/*.fastq.gz ${INPUT_DIR}/*.fq.gz; do
    [ -f "$FASTQ_GZ" ] || continue

    FILE_COUNT=$((FILE_COUNT + 1))
    BASENAME=$(basename ${FASTQ_GZ} .gz)
    BASENAME=${BASENAME%.fastq}
    BASENAME=${BASENAME%.fq}

    echo "  Processing: $(basename ${FASTQ_GZ})"

    # Subsample and convert to FASTA
    # seqtk would be ideal but using zcat + awk for compatibility
    zcat ${FASTQ_GZ} | awk -v n=${READS_PER_FILE} -v prefix="${BASENAME}" '
        BEGIN { srand(42); count=0 }
        NR%4==1 {
            header=$0
            getline seq
            getline plus
            getline qual
            if (count < n) {
                # Convert FASTQ header to FASTA
                gsub(/^@/, ">", header)
                print header
                print seq
                count++
            }
        }
    ' >> ${COMBINED_FASTA}

    TOTAL_READS=$((TOTAL_READS + READS_PER_FILE))
done

echo "  Files processed: ${FILE_COUNT}"
echo "  Total subsampled reads: ${TOTAL_READS}"
echo "  Combined FASTA: ${COMBINED_FASTA}"
ACTUAL_READS=$(grep -c "^>" ${COMBINED_FASTA})
echo "  Actual reads in FASTA: ${ACTUAL_READS}"

# Step 2: BLAST against phage database
echo ""
echo "Step 2: BLASTing against phage database..."
PHAGE_RESULTS="${OUTPUT_DIR}/phage_blast_results.tsv"

if [ -f "${PHAGE_DB}.nin" ] || [ -f "${PHAGE_DB}.nal" ]; then
    blastn \
        -task megablast \
        -db ${PHAGE_DB} \
        -query ${COMBINED_FASTA} \
        -out ${PHAGE_RESULTS} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -evalue ${EVALUE} \
        -max_target_seqs ${MAX_TARGET_SEQS} \
        -num_threads ${NUM_THREADS}

    PHAGE_HITS=$(cut -f1 ${PHAGE_RESULTS} | sort -u | wc -l)
    echo "  Reads with phage hits: ${PHAGE_HITS}"
else
    echo "  Skipped - phage database not found"
    PHAGE_HITS=0
fi

# Step 3: BLAST against bacteria database
echo ""
echo "Step 3: BLASTing against bacteria database..."
BACTERIA_RESULTS="${OUTPUT_DIR}/bacteria_blast_results.tsv"

if [ -f "${BACTERIA_DB}.nin" ] || [ -f "${BACTERIA_DB}.nal" ]; then
    blastn \
        -task megablast \
        -db ${BACTERIA_DB} \
        -query ${COMBINED_FASTA} \
        -out ${BACTERIA_RESULTS} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -evalue ${EVALUE} \
        -max_target_seqs ${MAX_TARGET_SEQS} \
        -num_threads ${NUM_THREADS}

    BACTERIA_HITS=$(cut -f1 ${BACTERIA_RESULTS} | sort -u | wc -l)
    echo "  Reads with bacteria hits: ${BACTERIA_HITS}"
else
    echo "  Skipped - bacteria database not found"
    BACTERIA_HITS=0
fi

# Step 4: Summarize results
echo ""
echo "Step 4: Generating summary..."
SUMMARY_FILE="${OUTPUT_DIR}/summary.txt"

# Get unique read IDs from each result
cut -f1 ${PHAGE_RESULTS} 2>/dev/null | sort -u > ${OUTPUT_DIR}/phage_read_ids.txt
cut -f1 ${BACTERIA_RESULTS} 2>/dev/null | sort -u > ${OUTPUT_DIR}/bacteria_read_ids.txt

# Find reads hitting both
comm -12 ${OUTPUT_DIR}/phage_read_ids.txt ${OUTPUT_DIR}/bacteria_read_ids.txt > ${OUTPUT_DIR}/both_read_ids.txt
BOTH_HITS=$(wc -l < ${OUTPUT_DIR}/both_read_ids.txt)

# Find reads hitting only phage
comm -23 ${OUTPUT_DIR}/phage_read_ids.txt ${OUTPUT_DIR}/bacteria_read_ids.txt > ${OUTPUT_DIR}/phage_only_read_ids.txt
PHAGE_ONLY=$(wc -l < ${OUTPUT_DIR}/phage_only_read_ids.txt)

# Find reads hitting only bacteria
comm -13 ${OUTPUT_DIR}/phage_read_ids.txt ${OUTPUT_DIR}/bacteria_read_ids.txt > ${OUTPUT_DIR}/bacteria_only_read_ids.txt
BACTERIA_ONLY=$(wc -l < ${OUTPUT_DIR}/bacteria_only_read_ids.txt)

# Reads with no hits
NO_HITS=$((ACTUAL_READS - PHAGE_HITS - BACTERIA_ONLY))
if [ ${NO_HITS} -lt 0 ]; then
    NO_HITS=$((ACTUAL_READS - PHAGE_ONLY - BACTERIA_ONLY - BOTH_HITS))
fi

# Write summary
cat > ${SUMMARY_FILE} << EOF
Unmapped Reads BLAST Analysis Summary
=====================================
Date: $(date)
Input directory: ${INPUT_DIR}
Files processed: ${FILE_COUNT}
Reads sampled per file: ${READS_PER_FILE}
Total reads analyzed: ${ACTUAL_READS}
E-value threshold: ${EVALUE}

Results:
--------
Phage hits (total):      ${PHAGE_HITS} ($(echo "scale=1; ${PHAGE_HITS}*100/${ACTUAL_READS}" | bc)%)
Bacteria hits (total):   ${BACTERIA_HITS} ($(echo "scale=1; ${BACTERIA_HITS}*100/${ACTUAL_READS}" | bc)%)

Breakdown:
  Phage only:            ${PHAGE_ONLY} ($(echo "scale=1; ${PHAGE_ONLY}*100/${ACTUAL_READS}" | bc)%)
  Bacteria only:         ${BACTERIA_ONLY} ($(echo "scale=1; ${BACTERIA_ONLY}*100/${ACTUAL_READS}" | bc)%)
  Both phage & bacteria: ${BOTH_HITS} ($(echo "scale=1; ${BOTH_HITS}*100/${ACTUAL_READS}" | bc)%)
  No hits:               ${NO_HITS} ($(echo "scale=1; ${NO_HITS}*100/${ACTUAL_READS}" | bc)%)

Output files:
  ${PHAGE_RESULTS}
  ${BACTERIA_RESULTS}
  ${SUMMARY_FILE}
EOF

# Print summary to stdout
cat ${SUMMARY_FILE}

# Step 5: Top hits breakdown
echo ""
echo "Top 10 phage hits:"
if [ -f "${PHAGE_RESULTS}" ]; then
    cut -f2 ${PHAGE_RESULTS} | sort | uniq -c | sort -rn | head -10
fi

echo ""
echo "Top 10 bacteria hits:"
if [ -f "${BACTERIA_RESULTS}" ]; then
    cut -f2 ${BACTERIA_RESULTS} | sort | uniq -c | sort -rn | head -10
fi

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "End time: $(date)"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="

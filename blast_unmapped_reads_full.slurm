#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=blast_unmapped_full
#SBATCH --output=blast_unmapped_full_%j.out
#SBATCH --error=blast_unmapped_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g

# BLAST ALL unmapped reads against phage and bacteria databases
# No subsampling - processes all reads
#
# Usage: sbatch blast_unmapped_reads_full.slurm

# ============================================
# CONFIGURATION
# ============================================

# Input directory with FASTQ.gz files
INPUT_DIR="/work/hdd/bfzj/llindsey1/Noahs_unmapped_data"

# Output directory for results
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/unmapped_blast_results_full"

# BLAST databases
PHAGE_DB="/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages"
BACTERIA_DB="/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected"

# BLAST parameters
EVALUE="1e-10"
NUM_THREADS=16
MAX_TARGET_SEQS=1

# ============================================

echo "=========================================="
echo "BLAST All Unmapped Reads Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "E-value threshold: ${EVALUE}"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Convert ALL FASTQ files to FASTA (no subsampling)
echo ""
echo "Step 1: Converting all FASTQ files to FASTA..."
COMBINED_FASTA="${OUTPUT_DIR}/all_unmapped_reads.fasta"
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

    # Convert FASTQ to FASTA (all reads)
    zcat ${FASTQ_GZ} | awk -v prefix="${BASENAME}" '
        NR%4==1 {
            # Convert FASTQ header to FASTA, add file prefix for tracking
            gsub(/^@/, ">" prefix "_", $0)
            print
        }
        NR%4==2 { print }
    ' >> ${COMBINED_FASTA}

done

TOTAL_READS=$(grep -c "^>" ${COMBINED_FASTA})
echo ""
echo "  Files processed: ${FILE_COUNT}"
echo "  Total reads: ${TOTAL_READS}"
echo "  Combined FASTA: ${COMBINED_FASTA}"
echo "  File size: $(ls -lh ${COMBINED_FASTA} | awk '{print $5}')"

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
NO_HITS=$((TOTAL_READS - PHAGE_ONLY - BACTERIA_ONLY - BOTH_HITS))

# Write summary
cat > ${SUMMARY_FILE} << EOF
Unmapped Reads BLAST Analysis Summary (FULL DATA)
==================================================
Date: $(date)
Input directory: ${INPUT_DIR}
Files processed: ${FILE_COUNT}
Total reads analyzed: ${TOTAL_READS}
E-value threshold: ${EVALUE}

Results:
--------
Phage hits (total):      ${PHAGE_HITS} ($(echo "scale=1; ${PHAGE_HITS}*100/${TOTAL_READS}" | bc)%)
Bacteria hits (total):   ${BACTERIA_HITS} ($(echo "scale=1; ${BACTERIA_HITS}*100/${TOTAL_READS}" | bc)%)

Breakdown:
  Phage only:            ${PHAGE_ONLY} ($(echo "scale=1; ${PHAGE_ONLY}*100/${TOTAL_READS}" | bc)%)
  Bacteria only:         ${BACTERIA_ONLY} ($(echo "scale=1; ${BACTERIA_ONLY}*100/${TOTAL_READS}" | bc)%)
  Both phage & bacteria: ${BOTH_HITS} ($(echo "scale=1; ${BOTH_HITS}*100/${TOTAL_READS}" | bc)%)
  No hits:               ${NO_HITS} ($(echo "scale=1; ${NO_HITS}*100/${TOTAL_READS}" | bc)%)

Output files:
  ${PHAGE_RESULTS}
  ${BACTERIA_RESULTS}
  ${OUTPUT_DIR}/phage_read_ids.txt
  ${OUTPUT_DIR}/bacteria_read_ids.txt
  ${OUTPUT_DIR}/phage_only_read_ids.txt
  ${OUTPUT_DIR}/bacteria_only_read_ids.txt
  ${OUTPUT_DIR}/both_read_ids.txt
  ${SUMMARY_FILE}
EOF

# Print summary to stdout
cat ${SUMMARY_FILE}

# Step 5: Top hits breakdown
echo ""
echo "Top 20 phage hits:"
if [ -f "${PHAGE_RESULTS}" ]; then
    cut -f2 ${PHAGE_RESULTS} | sort | uniq -c | sort -rn | head -20
fi

echo ""
echo "Top 20 bacteria hits:"
if [ -f "${BACTERIA_RESULTS}" ]; then
    cut -f2 ${BACTERIA_RESULTS} | sort | uniq -c | sort -rn | head -20
fi

# Step 6: Extract reads with no hits (for further analysis)
echo ""
echo "Step 5: Extracting reads with no hits..."
NO_HITS_FASTA="${OUTPUT_DIR}/no_hits_reads.fasta"

# Get all read IDs
grep "^>" ${COMBINED_FASTA} | sed 's/^>//' | sort > ${OUTPUT_DIR}/all_read_ids.txt

# Get IDs with any hit
cat ${OUTPUT_DIR}/phage_read_ids.txt ${OUTPUT_DIR}/bacteria_read_ids.txt | sort -u > ${OUTPUT_DIR}/any_hit_ids.txt

# Get IDs with no hits
comm -23 ${OUTPUT_DIR}/all_read_ids.txt ${OUTPUT_DIR}/any_hit_ids.txt > ${OUTPUT_DIR}/no_hit_ids.txt

echo "  Reads with no hits: $(wc -l < ${OUTPUT_DIR}/no_hit_ids.txt)"

# Extract no-hit sequences to FASTA (for eukaryotic virus BLAST)
echo "  Extracting no-hit sequences to FASTA..."
awk 'NR==FNR {ids[$1]; next} /^>/ {p=0; id=substr($1,2); if(id in ids) p=1} p' \
    ${OUTPUT_DIR}/no_hit_ids.txt ${COMBINED_FASTA} > ${NO_HITS_FASTA}

echo "  No-hits FASTA: ${NO_HITS_FASTA}"
echo "  Sequences: $(grep -c "^>" ${NO_HITS_FASTA})"
echo "  File size: $(ls -lh ${NO_HITS_FASTA} | awk '{print $5}')"

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "End time: $(date)"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
echo ""
echo "Next step: BLAST ${NO_HITS_FASTA} against eukaryotic virus database on Biowulf"

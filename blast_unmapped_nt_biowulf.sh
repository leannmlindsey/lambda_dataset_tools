#!/bin/bash
#SBATCH --job-name=blast_nt
#SBATCH --output=blast_nt_%j.out
#SBATCH --error=blast_nt_%j.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --partition=norm
#SBATCH --gres=lscratch:200

# BLAST unmapped reads against NCBI nt database (all sequences)
# For NIH Biowulf HPC
#
# Before running:
#   1. Copy your no_hits_reads.fasta to Biowulf
#   2. Update INPUT_FASTA path below
#   3. Submit with: sbatch blast_unmapped_nt_biowulf.sh
#
# Database info:
#   nt = complete NCBI nucleotide collection (~100GB+ database)
#   This will find matches to bacteria, viruses, fungi, plants, animals, etc.

# ============================================
# CONFIGURATION - UPDATE THESE PATHS
# ============================================

# Input FASTA (reads with no hits to phage/bacteria databases)
# Copy from: /work/hdd/bfzj/llindsey1/unmapped_blast_results_full/no_hits_reads.fasta
INPUT_FASTA="/data/$USER/unmapped_reads/no_hits_reads.fasta"

# Output directory
OUTPUT_DIR="/data/$USER/unmapped_reads/nt_blast_results"

# BLAST database (nt = all NCBI nucleotides)
# Check available databases: ls /fdb/blastdb/
BLAST_DB="/fdb/blastdb/nt"

# BLAST parameters
EVALUE="1e-10"
MAX_TARGET_SEQS=5
NUM_THREADS=32

# Use local scratch for temp files (faster I/O)
export TMPDIR=/lscratch/$SLURM_JOB_ID

# ============================================

echo "=========================================="
echo "BLAST Unmapped Reads vs NCBI nt Database"
echo "=========================================="
echo "Start time: $(date)"
echo "Input: ${INPUT_FASTA}"
echo "Database: ${BLAST_DB}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${NUM_THREADS}"
echo "=========================================="

# Load BLAST module
module load blast/2.14.0

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check input exists
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "Error: Input FASTA not found: ${INPUT_FASTA}"
    echo ""
    echo "Please copy your reads file to Biowulf first:"
    echo "  scp /work/hdd/bfzj/llindsey1/unmapped_blast_results_full/no_hits_reads.fasta \\
        ${USER}@biowulf.nih.gov:/data/${USER}/unmapped_reads/"
    exit 1
fi

NUM_READS=$(grep -c "^>" ${INPUT_FASTA})
echo "Input reads: ${NUM_READS}"

# ============================================
# Run BLAST against nt
# ============================================
echo ""
echo "Running BLASTn against nt database..."
echo "This may take several hours..."

BLAST_OUTPUT="${OUTPUT_DIR}/nt_blast_results.tsv"

# Output format includes taxonomy info
# Columns: qseqid sseqid stitle staxids sscinames sskingdoms pident length mismatch gapopen qstart qend sstart send evalue bitscore
blastn \
    -task megablast \
    -db ${BLAST_DB} \
    -query ${INPUT_FASTA} \
    -out ${BLAST_OUTPUT} \
    -outfmt "6 qseqid sseqid stitle staxids sscinames sskingdoms pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue ${EVALUE} \
    -max_target_seqs ${MAX_TARGET_SEQS} \
    -num_threads ${NUM_THREADS}

echo "BLAST complete!"
echo "Total hits: $(wc -l < ${BLAST_OUTPUT})"

# ============================================
# Analyze results by kingdom/domain
# ============================================
echo ""
echo "=========================================="
echo "Analyzing results by taxonomy..."
echo "=========================================="

# Unique reads with hits
READS_WITH_HITS=$(cut -f1 ${BLAST_OUTPUT} | sort -u | wc -l)
READS_NO_HITS=$((NUM_READS - READS_WITH_HITS))
echo "Reads with hits: ${READS_WITH_HITS}"
echo "Reads still unidentified: ${READS_NO_HITS}"

# Breakdown by kingdom (column 6 = sskingdoms)
echo ""
echo "Hits by Kingdom/Domain:"
echo "----------------------------------------"

# Count unique reads per kingdom
for KINGDOM in Bacteria Viruses Eukaryota Archaea; do
    COUNT=$(grep -i "${KINGDOM}" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l)
    echo "  ${KINGDOM}: ${COUNT} reads"
done

# Unclassified/other
OTHER=$(grep -v -i "Bacteria\|Viruses\|Eukaryota\|Archaea" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l)
echo "  Other/Unclassified: ${OTHER} reads"

# ============================================
# Detailed virus breakdown
# ============================================
echo ""
echo "=========================================="
echo "Virus hits breakdown:"
echo "=========================================="

VIRUS_HITS="${OUTPUT_DIR}/virus_hits.tsv"
grep -i "Viruses" ${BLAST_OUTPUT} > ${VIRUS_HITS}
VIRUS_COUNT=$(cut -f1 ${VIRUS_HITS} | sort -u | wc -l)
echo "Total reads with virus hits: ${VIRUS_COUNT}"

echo ""
echo "Top 20 virus species (by read count):"
cut -f5 ${VIRUS_HITS} | sort | uniq -c | sort -rn | head -20

# Categorize viruses
echo ""
echo "Virus categories:"
echo "  Bacteriophage: $(grep -i "phage" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Human viruses: $(grep -i "Human\|Homo sapiens" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Herpesvirus: $(grep -i "herpes" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Papillomavirus: $(grep -i "papilloma" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Adenovirus: $(grep -i "adeno" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Retrovirus/HIV: $(grep -i "retro\|HIV\|HTLV" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Coronavirus: $(grep -i "corona\|SARS\|COVID" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"
echo "  Influenza: $(grep -i "influenza\|Orthomyxo" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)"

# ============================================
# Detailed eukaryote breakdown
# ============================================
echo ""
echo "=========================================="
echo "Eukaryote hits breakdown:"
echo "=========================================="

EUK_HITS="${OUTPUT_DIR}/eukaryote_hits.tsv"
grep -i "Eukaryota" ${BLAST_OUTPUT} > ${EUK_HITS}
EUK_COUNT=$(cut -f1 ${EUK_HITS} | sort -u | wc -l)
echo "Total reads with eukaryote hits: ${EUK_COUNT}"

echo ""
echo "Top 20 eukaryote species:"
cut -f5 ${EUK_HITS} | sort | uniq -c | sort -rn | head -20

# Check for human sequences (potential host contamination)
HUMAN_HITS=$(grep -i "Homo sapiens" ${EUK_HITS} | cut -f1 | sort -u | wc -l)
echo ""
echo "Human sequences (potential host contamination): ${HUMAN_HITS}"

# ============================================
# Bacteria hits (unexpected - should have been caught earlier)
# ============================================
echo ""
echo "=========================================="
echo "Bacteria hits (check these):"
echo "=========================================="

BACT_HITS="${OUTPUT_DIR}/bacteria_hits.tsv"
grep -i "Bacteria" ${BLAST_OUTPUT} > ${BACT_HITS}
BACT_COUNT=$(cut -f1 ${BACT_HITS} | sort -u | wc -l)
echo "Reads with bacteria hits: ${BACT_COUNT}"
echo "(These weren't caught by your GTDB database - may be species not in GTDB)"

echo ""
echo "Top 10 bacteria species:"
cut -f5 ${BACT_HITS} | sort | uniq -c | sort -rn | head -10

# ============================================
# Generate comprehensive summary
# ============================================
echo ""
echo "=========================================="
echo "Generating summary files..."
echo "=========================================="

SUMMARY="${OUTPUT_DIR}/nt_blast_summary.txt"
cat > ${SUMMARY} << EOF
NCBI nt BLAST Analysis Summary
==============================
Date: $(date)
Input: ${INPUT_FASTA}
Total input reads: ${NUM_READS}
Database: ${BLAST_DB}
E-value threshold: ${EVALUE}

Overall Results:
----------------
Reads with hits: ${READS_WITH_HITS} ($(echo "scale=1; ${READS_WITH_HITS}*100/${NUM_READS}" | bc)%)
Reads unidentified: ${READS_NO_HITS} ($(echo "scale=1; ${READS_NO_HITS}*100/${NUM_READS}" | bc)%)

By Kingdom/Domain:
------------------
$(for K in Bacteria Viruses Eukaryota Archaea; do
    C=$(grep -i "${K}" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l)
    echo "  ${K}: ${C} reads"
done)

Virus Breakdown:
----------------
Total virus reads: ${VIRUS_COUNT}
  Bacteriophage: $(grep -i "phage" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)
  Human viruses: $(grep -i "Human\|Homo sapiens" ${VIRUS_HITS} | cut -f1 | sort -u | wc -l)

Eukaryote Breakdown:
--------------------
Total eukaryote reads: ${EUK_COUNT}
  Human (host contamination?): ${HUMAN_HITS}

Output Files:
-------------
${BLAST_OUTPUT} - All BLAST hits
${VIRUS_HITS} - Virus hits only
${EUK_HITS} - Eukaryote hits only
${BACT_HITS} - Bacteria hits only
${SUMMARY} - This summary

EOF

cat ${SUMMARY}

# Save top hits for each category
echo ""
echo "Saving top hits files..."

cut -f2,3,5,6 ${BLAST_OUTPUT} | sort | uniq -c | sort -rn | head -100 > ${OUTPUT_DIR}/top_100_hits_overall.txt
cut -f2,3,5 ${VIRUS_HITS} | sort | uniq -c | sort -rn | head -50 > ${OUTPUT_DIR}/top_50_virus_hits.txt
cut -f2,3,5 ${EUK_HITS} | sort | uniq -c | sort -rn | head -50 > ${OUTPUT_DIR}/top_50_eukaryote_hits.txt

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "End time: $(date)"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="

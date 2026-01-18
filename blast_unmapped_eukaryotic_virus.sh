#!/bin/bash
#SBATCH --job-name=blast_eukvirus
#SBATCH --output=blast_eukvirus_%j.out
#SBATCH --error=blast_eukvirus_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --partition=norm

# BLAST unmapped reads against eukaryotic virus database
# Designed for NIH Biowulf - adjust paths/modules as needed
#
# Usage: sbatch blast_unmapped_eukaryotic_virus.sh

# ============================================
# CONFIGURATION - ADJUST THESE FOR BIOWULF
# ============================================

# Input: the reads that had no hits to phage or bacteria
# You may need to copy this file to Biowulf first
INPUT_FASTA="/path/to/unmapped_reads/no_hits_reads.fasta"

# Or use the full subsampled reads file
# INPUT_FASTA="/path/to/unmapped_reads/subsampled_reads.fasta"

# Output directory
OUTPUT_DIR="/path/to/output/eukaryotic_virus_blast"

# BLAST database options on Biowulf (check with: ls /fdb/blastdb/)
# Option 1: RefSeq viral representative genomes (includes eukaryotic viruses)
VIRUS_DB="/fdb/blastdb/ref_viruses_rep_genomes"

# Option 2: If you want to search nt with viral taxonomy filter
# VIRUS_DB="/fdb/blastdb/nt"
# Add: -taxids 10239  (Viruses taxid)

# BLAST parameters
EVALUE="1e-10"
MAX_TARGET_SEQS=5
NUM_THREADS=16

# ============================================

echo "=========================================="
echo "BLAST Unmapped Reads vs Eukaryotic Viruses"
echo "=========================================="
echo "Start time: $(date)"
echo "Input: ${INPUT_FASTA}"
echo "Database: ${VIRUS_DB}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================="

# Load modules (Biowulf-specific)
module load blast

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check input exists
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "Error: Input FASTA not found: ${INPUT_FASTA}"
    echo "Please copy your reads file to Biowulf first."
    exit 1
fi

NUM_READS=$(grep -c "^>" ${INPUT_FASTA})
echo "Input reads: ${NUM_READS}"

# Run BLAST
echo ""
echo "Running BLASTn against viral database..."

BLAST_OUTPUT="${OUTPUT_DIR}/eukaryotic_virus_blast.tsv"

blastn \
    -db ${VIRUS_DB} \
    -query ${INPUT_FASTA} \
    -out ${BLAST_OUTPUT} \
    -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue ${EVALUE} \
    -max_target_seqs ${MAX_TARGET_SEQS} \
    -num_threads ${NUM_THREADS}

echo "BLAST complete!"
echo "Total hits: $(wc -l < ${BLAST_OUTPUT})"

# Analyze results
echo ""
echo "=========================================="
echo "Analyzing results..."
echo "=========================================="

# Unique reads with hits
READS_WITH_HITS=$(cut -f1 ${BLAST_OUTPUT} | sort -u | wc -l)
echo "Reads with viral hits: ${READS_WITH_HITS}"

# Top virus hits (by subject title - column 3)
echo ""
echo "Top 20 viral hits:"
TOP_HITS="${OUTPUT_DIR}/top_virus_hits.txt"
cut -f2,3 ${BLAST_OUTPUT} | sort | uniq -c | sort -rn | head -20 | tee ${TOP_HITS}

# Categorize by virus type (look for keywords)
echo ""
echo "=========================================="
echo "Categorization by keywords:"
echo "=========================================="

echo "Bacteriophage-related:"
grep -i "phage\|bacteriophage" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Human virus-related:"
grep -i "human\|homo sapiens" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Herpesvirus:"
grep -i "herpes" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Papillomavirus:"
grep -i "papilloma" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Adenovirus:"
grep -i "adeno" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Retrovirus:"
grep -i "retro\|HIV\|HTLV" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Parvovirus:"
grep -i "parvo" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

echo "Polyomavirus:"
grep -i "polyoma" ${BLAST_OUTPUT} | cut -f1 | sort -u | wc -l

# Generate summary
SUMMARY="${OUTPUT_DIR}/eukaryotic_virus_summary.txt"
cat > ${SUMMARY} << EOF
Eukaryotic Virus BLAST Summary
==============================
Date: $(date)
Input: ${INPUT_FASTA}
Total reads: ${NUM_READS}
Database: ${VIRUS_DB}
E-value threshold: ${EVALUE}

Results:
  Total BLAST hits: $(wc -l < ${BLAST_OUTPUT})
  Reads with hits: ${READS_WITH_HITS}
  Percent with hits: $(echo "scale=1; ${READS_WITH_HITS} * 100 / ${NUM_READS}" | bc)%

Top 20 virus hits:
$(cat ${TOP_HITS})

Output files:
  ${BLAST_OUTPUT}
  ${TOP_HITS}
  ${SUMMARY}
EOF

echo ""
echo "Summary saved to: ${SUMMARY}"
cat ${SUMMARY}

echo ""
echo "=========================================="
echo "Done!"
echo "End time: $(date)"
echo "=========================================="

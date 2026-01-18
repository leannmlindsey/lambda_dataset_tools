#!/bin/bash

# Debug: Investigate why some hits remain after filtering
# Check if the contigs with remaining hits were properly mapped

CONTIG_MAP="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contig_to_genome_map.tsv"
ORIGINAL_BLAST="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/blast_filtered_combined.tsv"
CONTAMINATED="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contaminated_accessions.txt"
CLEAN="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/clean_accessions.txt"

echo "=========================================="
echo "Debug: Investigating Missed Hits"
echo "=========================================="

# The contigs with remaining hits
CONTIGS=("NZ_CP045798.1" "NZ_JPST01000043.1")

for CONTIG in "${CONTIGS[@]}"; do
    echo ""
    echo "Investigating contig: ${CONTIG}"
    echo "----------------------------------------"

    # Check if contig is in the mapping
    echo "1. Checking contig-to-genome mapping:"
    GENOME=$(grep "^${CONTIG}" ${CONTIG_MAP} | cut -f2)
    if [ -n "$GENOME" ]; then
        echo "   Contig maps to genome: ${GENOME}"
    else
        echo "   WARNING: Contig NOT FOUND in mapping!"
        # Try partial match
        echo "   Trying partial match..."
        grep "${CONTIG}" ${CONTIG_MAP} | head -3
        continue
    fi

    # Check if this genome was in original BLAST hits
    echo ""
    echo "2. Checking original BLAST for this contig:"
    ORIGINAL_HITS=$(grep "${CONTIG}" ${ORIGINAL_BLAST} | wc -l)
    echo "   Original BLAST hits for this contig: ${ORIGINAL_HITS}"
    if [ ${ORIGINAL_HITS} -gt 0 ]; then
        echo "   Sample hits:"
        grep "${CONTIG}" ${ORIGINAL_BLAST} | head -3
    fi

    # Check if genome is in contaminated list
    echo ""
    echo "3. Checking if genome ${GENOME} is in contaminated list:"
    if grep -q "^${GENOME}$" ${CONTAMINATED}; then
        echo "   YES - genome IS in contaminated list"
        echo "   ERROR: This genome should have been removed!"
    else
        echo "   NO - genome is NOT in contaminated list"
        echo "   This explains why it wasn't filtered out."
    fi

    # Check if genome is in clean list
    echo ""
    echo "4. Checking if genome ${GENOME} is in clean list:"
    if grep -q "^${GENOME}$" ${CLEAN}; then
        echo "   YES - genome IS in clean list"
    else
        echo "   NO - genome is NOT in clean list"
    fi
done

# Check the phage MZ092003 in original BLAST
echo ""
echo "=========================================="
echo "Checking phage MZ092003 in original BLAST:"
echo "=========================================="
echo "Total hits for MZ092003 in original BLAST:"
grep "^MZ092003" ${ORIGINAL_BLAST} | wc -l

echo ""
echo "Unique contigs hit by MZ092003:"
grep "^MZ092003" ${ORIGINAL_BLAST} | cut -f2 | sort -u | head -10

echo ""
echo "Sample hits for MZ092003:"
grep "^MZ092003" ${ORIGINAL_BLAST} | head -5

echo ""
echo "=========================================="
echo "Debug complete"
echo "=========================================="

#!/bin/bash

# Identify what organisms the top bacterial contigs belong to
# Outputs deduplicated results to file

# ============================================
# CONFIGURATION
# ============================================

BACTERIA_RESULTS="/work/hdd/bfzj/llindsey1/unmapped_blast_results/bacteria_blast_results.tsv"
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/unmapped_blast_results"
CONTIG_MAP="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contig_to_genome_map.tsv"

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/top_bacteria_analysis.txt"

# ============================================

echo "=========================================="
echo "Identify Top Bacterial Hits"
echo "=========================================="
echo "Output file: ${OUTPUT_FILE}"
echo ""

# Start output file
cat > ${OUTPUT_FILE} << EOF
Top Bacterial Hits Analysis
===========================
Date: $(date)
Input: ${BACTERIA_RESULTS}

EOF

# Extract top 50 bacterial contigs by hit count (deduplicated)
echo "Extracting top 50 bacterial contigs by hit count..."

TOP_CONTIGS="${OUTPUT_DIR}/top_bacteria_contigs.txt"

# Clean up contig IDs (remove prefixes like emb|, ref|, gb|, tpe|, etc.)
cut -f2 ${BACTERIA_RESULTS} | \
    sed 's/^[a-z]*|//; s/|$//' | \
    sort | uniq -c | sort -rn | head -50 > ${TOP_CONTIGS}

echo "" >> ${OUTPUT_FILE}
echo "Top 50 Bacterial Contigs by Hit Count:" >> ${OUTPUT_FILE}
echo "=======================================" >> ${OUTPUT_FILE}
echo "Hits    Contig_ID" >> ${OUTPUT_FILE}
cat ${TOP_CONTIGS} >> ${OUTPUT_FILE}

# Extract unique accession prefixes (genome project IDs)
echo "" >> ${OUTPUT_FILE}
echo "" >> ${OUTPUT_FILE}
echo "Unique Genome Projects (by accession prefix):" >> ${OUTPUT_FILE}
echo "==============================================" >> ${OUTPUT_FILE}

# Group by genome project (first part of accession before numbers)
cut -f2 ${BACTERIA_RESULTS} | \
    sed 's/^[a-z]*|//; s/|$//' | \
    sed 's/[0-9]*\.[0-9]*$//' | \
    sort | uniq -c | sort -rn | head -30 >> ${OUTPUT_FILE}

# Map contigs to GTDB genomes if mapping exists
if [ -f "${CONTIG_MAP}" ]; then
    echo "" >> ${OUTPUT_FILE}
    echo "" >> ${OUTPUT_FILE}
    echo "Mapping to GTDB Genome Accessions:" >> ${OUTPUT_FILE}
    echo "===================================" >> ${OUTPUT_FILE}
    echo "(Contigs found in our GTDB selected database)" >> ${OUTPUT_FILE}
    echo "" >> ${OUTPUT_FILE}
    echo "Hits    Contig_ID -> Genome_Accession" >> ${OUTPUT_FILE}

    # Create temp file for genome hit counts
    GENOME_COUNTS="${OUTPUT_DIR}/genome_hit_counts.tmp"
    > ${GENOME_COUNTS}

    while read line; do
        count=$(echo $line | awk '{print $1}')
        contig=$(echo $line | awk '{print $2}')

        # Look up genome in mapping
        genome=$(grep "^${contig}	" ${CONTIG_MAP} | cut -f2)

        if [ -n "$genome" ]; then
            echo "${count}    ${contig} -> ${genome}" >> ${OUTPUT_FILE}
            echo "${count} ${genome}" >> ${GENOME_COUNTS}
        fi
    done < ${TOP_CONTIGS}

    # Aggregate by genome (sum hits across contigs from same genome)
    echo "" >> ${OUTPUT_FILE}
    echo "" >> ${OUTPUT_FILE}
    echo "Aggregated Hits by Genome (deduplicated):" >> ${OUTPUT_FILE}
    echo "==========================================" >> ${OUTPUT_FILE}
    echo "Total_Hits    Genome_Accession" >> ${OUTPUT_FILE}

    awk '{sum[$2]+=$1} END {for (g in sum) print sum[g], g}' ${GENOME_COUNTS} | \
        sort -rn | head -30 >> ${OUTPUT_FILE}

    rm -f ${GENOME_COUNTS}
fi

# Count total unique contigs and genomes
echo "" >> ${OUTPUT_FILE}
echo "" >> ${OUTPUT_FILE}
echo "Summary Statistics:" >> ${OUTPUT_FILE}
echo "===================" >> ${OUTPUT_FILE}

TOTAL_HITS=$(wc -l < ${BACTERIA_RESULTS})
UNIQUE_CONTIGS=$(cut -f2 ${BACTERIA_RESULTS} | sed 's/^[a-z]*|//; s/|$//' | sort -u | wc -l)
UNIQUE_READS=$(cut -f1 ${BACTERIA_RESULTS} | sort -u | wc -l)

echo "Total BLAST hits: ${TOTAL_HITS}" >> ${OUTPUT_FILE}
echo "Unique contigs hit: ${UNIQUE_CONTIGS}" >> ${OUTPUT_FILE}
echo "Unique reads with bacteria hits: ${UNIQUE_READS}" >> ${OUTPUT_FILE}

# Print summary to screen
echo ""
echo "Analysis complete!"
echo ""
echo "Summary:"
echo "  Total BLAST hits: ${TOTAL_HITS}"
echo "  Unique contigs hit: ${UNIQUE_CONTIGS}"
echo "  Unique reads with bacteria hits: ${UNIQUE_READS}"
echo ""
echo "Top 10 contigs:"
head -10 ${TOP_CONTIGS}
echo ""
echo "Full results saved to: ${OUTPUT_FILE}"
echo ""
cat ${OUTPUT_FILE}

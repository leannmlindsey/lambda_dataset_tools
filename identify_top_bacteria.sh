#!/bin/bash

# Identify what organisms the top bacterial contigs belong to
# Uses NCBI's Entrez Direct (EDirect) tools to look up accession info

# ============================================
# CONFIGURATION
# ============================================

BACTERIA_RESULTS="/work/hdd/bfzj/llindsey1/unmapped_blast_results/bacteria_blast_results.tsv"
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/unmapped_blast_results"

# ============================================

echo "=========================================="
echo "Identify Top Bacterial Hits"
echo "=========================================="

# Extract top 20 bacterial contigs by hit count
echo "Top 20 bacterial contigs by hit count:"
TOP_CONTIGS="${OUTPUT_DIR}/top_bacteria_contigs.txt"

cut -f2 ${BACTERIA_RESULTS} | \
    sed 's/.*|//' | sed 's/|$//' | \
    sort | uniq -c | sort -rn | head -20 > ${TOP_CONTIGS}

cat ${TOP_CONTIGS}

echo ""
echo "=========================================="
echo "Looking up accession information..."
echo "=========================================="

# Extract just the accession IDs (clean format)
ACCESSIONS="${OUTPUT_DIR}/top_bacteria_accessions.txt"
awk '{print $2}' ${TOP_CONTIGS} | sed 's/\.[0-9]*$//' > ${ACCESSIONS}

echo "Unique accession prefixes (first 4 chars indicate genome project):"
cut -c1-6 ${ACCESSIONS} | sort -u

echo ""
echo "Full accessions to look up:"
cat ${ACCESSIONS}

# Try to get organism info using efetch if available
if command -v efetch &> /dev/null; then
    echo ""
    echo "=========================================="
    echo "Fetching organism information from NCBI..."
    echo "=========================================="

    OUTPUT_INFO="${OUTPUT_DIR}/bacteria_organism_info.txt"
    > ${OUTPUT_INFO}

    while read acc; do
        echo "Looking up: ${acc}"
        # Get nucleotide record and extract organism
        efetch -db nuccore -id "${acc}" -format docsum 2>/dev/null | \
            xtract -pattern DocumentSummary -element Caption,Title,Organism >> ${OUTPUT_INFO}
        sleep 1  # Be nice to NCBI
    done < ${ACCESSIONS}

    echo ""
    echo "Results saved to: ${OUTPUT_INFO}"
    cat ${OUTPUT_INFO}
else
    echo ""
    echo "EDirect tools (efetch) not found."
    echo "To install: conda install -c bioconda entrez-direct"
    echo ""
    echo "Alternative: Look up these accessions manually at NCBI:"
    echo "https://www.ncbi.nlm.nih.gov/nuccore/"
    echo ""
    cat ${ACCESSIONS}
fi

# Also check if any of these are in our GTDB contig mapping
CONTIG_MAP="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contig_to_genome_map.tsv"

if [ -f "${CONTIG_MAP}" ]; then
    echo ""
    echo "=========================================="
    echo "Checking GTDB contig mapping..."
    echo "=========================================="

    while read line; do
        count=$(echo $line | awk '{print $1}')
        contig=$(echo $line | awk '{print $2}')

        # Clean up contig ID (remove prefixes like emb|, ref|, etc.)
        clean_contig=$(echo $contig | sed 's/.*|//' | sed 's/|$//')

        genome=$(grep "^${clean_contig}" ${CONTIG_MAP} | cut -f2)

        if [ -n "$genome" ]; then
            echo "${count} hits: ${clean_contig} -> ${genome}"
        fi
    done < ${TOP_CONTIGS}
fi

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

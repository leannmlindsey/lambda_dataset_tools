#!/bin/bash

# Check the format of BLAST subject IDs

BLAST_RESULTS="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/blast_filtered_combined.tsv"

echo "First 10 lines of BLAST results:"
head -10 ${BLAST_RESULTS}

echo ""
echo "Unique subject ID formats (column 2, first 20):"
cut -f2 ${BLAST_RESULTS} | head -100 | sort -u | head -20

echo ""
echo "Sample subject IDs:"
cut -f2 ${BLAST_RESULTS} | shuf | head -10

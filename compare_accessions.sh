#!/bin/bash

# Compare original accessions with BLAST hit accessions
# to see how many clean genomes remain

# ============================================
# CONFIGURATION
# ============================================

# Original accession list (the selected genomes)
ORIGINAL_ACCESSIONS="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset/all_accessions.txt"

# Combined BLAST results from v2
BLAST_RESULTS="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/blast_filtered_combined.tsv"

# Output directory
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/prophage_blast_results_v2"

# ============================================

echo "=========================================="
echo "Compare Original vs BLAST Hit Accessions"
echo "=========================================="

# Step 1: Combine BLAST results if not already done
if [ ! -f "${BLAST_RESULTS}" ]; then
    echo "Combining BLAST filtered results..."
    cat ${OUTPUT_DIR}/blast_filtered_*.tsv > ${BLAST_RESULTS}
fi

echo "BLAST results file: ${BLAST_RESULTS}"
echo "Total BLAST hits: $(wc -l < ${BLAST_RESULTS})"
echo ""

# Step 2: Extract unique bacterial accessions from BLAST hits
# Column 2 = sseqid (bacterial sequence ID)
# Extract GCA_XXXXXXXXX.Y or GCF_XXXXXXXXX.Y pattern
echo "Extracting unique bacterial accessions from BLAST hits..."
BLAST_ACCESSIONS="${OUTPUT_DIR}/blast_hit_accessions.txt"

cut -f2 ${BLAST_RESULTS} | \
    grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | \
    sort -u > ${BLAST_ACCESSIONS}

NUM_BLAST_ACCESSIONS=$(wc -l < ${BLAST_ACCESSIONS})
echo "Unique accessions with phage hits: ${NUM_BLAST_ACCESSIONS}"
echo ""

# Step 3: Normalize original accessions (strip GB_/RS_ prefix)
echo "Normalizing original accession list..."
ORIGINAL_NORMALIZED="${OUTPUT_DIR}/original_accessions_normalized.txt"

sed 's/^GB_//; s/^RS_//' ${ORIGINAL_ACCESSIONS} | sort -u > ${ORIGINAL_NORMALIZED}

NUM_ORIGINAL=$(wc -l < ${ORIGINAL_NORMALIZED})
echo "Original accessions (normalized, unique): ${NUM_ORIGINAL}"
echo ""

# Step 4: Find accessions that are in original but NOT in BLAST hits (clean genomes)
echo "Finding clean genomes (no phage hits)..."
CLEAN_ACCESSIONS="${OUTPUT_DIR}/clean_accessions.txt"

comm -23 ${ORIGINAL_NORMALIZED} ${BLAST_ACCESSIONS} > ${CLEAN_ACCESSIONS}

NUM_CLEAN=$(wc -l < ${CLEAN_ACCESSIONS})
echo "Clean accessions (no phage hits): ${NUM_CLEAN}"
echo ""

# Step 5: Find accessions that are in BLAST hits but NOT in original (shouldn't happen)
echo "Checking for BLAST hits not in original list..."
UNEXPECTED="${OUTPUT_DIR}/unexpected_accessions.txt"

comm -13 ${ORIGINAL_NORMALIZED} ${BLAST_ACCESSIONS} > ${UNEXPECTED}

NUM_UNEXPECTED=$(wc -l < ${UNEXPECTED})
if [ ${NUM_UNEXPECTED} -gt 0 ]; then
    echo "WARNING: ${NUM_UNEXPECTED} BLAST hit accessions not in original list!"
    echo "First 10:"
    head -10 ${UNEXPECTED}
else
    echo "All BLAST hit accessions are in original list (as expected)"
fi
echo ""

# Step 6: Summary
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Original accessions:     ${NUM_ORIGINAL}"
echo "Contaminated (with hits): ${NUM_BLAST_ACCESSIONS}"
echo "Clean (no hits):         ${NUM_CLEAN}"
echo ""

PERCENT_CLEAN=$(echo "scale=1; ${NUM_CLEAN} * 100 / ${NUM_ORIGINAL}" | bc)
PERCENT_CONTAMINATED=$(echo "scale=1; ${NUM_BLAST_ACCESSIONS} * 100 / ${NUM_ORIGINAL}" | bc)

echo "Percent clean:           ${PERCENT_CLEAN}%"
echo "Percent contaminated:    ${PERCENT_CONTAMINATED}%"
echo ""

if [ ${NUM_CLEAN} -eq 0 ]; then
    echo "WARNING: All genomes have at least one phage hit!"
    echo "You may need to adjust filtering thresholds or accept some contamination."
else
    echo "You have ${NUM_CLEAN} clean genomes to work with."
fi

echo ""
echo "Output files:"
echo "  ${BLAST_ACCESSIONS}"
echo "  ${ORIGINAL_NORMALIZED}"
echo "  ${CLEAN_ACCESSIONS}"
echo "=========================================="

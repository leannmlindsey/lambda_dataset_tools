#!/bin/bash
# Copy Pharokka annotations for test set phages to a separate directory
# This prepares the Phage-Only benchmark dataset for Zenodo upload
#
# Usage: bash copy_pharokka_test_annotations.sh

SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Input
TEST_ACCESSIONS="${SCRIPT_DIR}/inphared_dataset/test_accessions.txt"
PHAROKKA_DIR="/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka"

# Output
OUTPUT_DIR="${SCRIPT_DIR}/phage_only_annotations"

echo "=========================================="
echo "Copy Pharokka Annotations for Test Set"
echo "=========================================="
echo ""
echo "Test accessions: ${TEST_ACCESSIONS}"
echo "Pharokka source: ${PHAROKKA_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check inputs exist
if [ ! -f "${TEST_ACCESSIONS}" ]; then
    echo "Error: Test accessions file not found: ${TEST_ACCESSIONS}"
    exit 1
fi

if [ ! -d "${PHAROKKA_DIR}" ]; then
    echo "Error: Pharokka directory not found: ${PHAROKKA_DIR}"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Count accessions
TOTAL=$(wc -l < "${TEST_ACCESSIONS}")
echo "Total test accessions: ${TOTAL}"
echo ""

# Copy annotations for each accession
COPIED=0
MISSING=0

while read -r ACCESSION; do
    # Skip empty lines
    [ -z "${ACCESSION}" ] && continue

    # Check for pharokka output directory for this accession
    # Pharokka typically creates a directory per genome
    PHAROKKA_ACCESSION_DIR="${PHAROKKA_DIR}/${ACCESSION}"

    if [ -d "${PHAROKKA_ACCESSION_DIR}" ]; then
        cp -r "${PHAROKKA_ACCESSION_DIR}" "${OUTPUT_DIR}/"
        ((COPIED++))
    else
        # Try without version number (e.g., NC_001416 instead of NC_001416.1)
        ACCESSION_BASE=$(echo "${ACCESSION}" | sed 's/\.[0-9]*$//')
        PHAROKKA_ACCESSION_DIR="${PHAROKKA_DIR}/${ACCESSION_BASE}"

        if [ -d "${PHAROKKA_ACCESSION_DIR}" ]; then
            cp -r "${PHAROKKA_ACCESSION_DIR}" "${OUTPUT_DIR}/"
            ((COPIED++))
        else
            echo "  Missing: ${ACCESSION}"
            ((MISSING++))
        fi
    fi
done < "${TEST_ACCESSIONS}"

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Total test accessions: ${TOTAL}"
echo "Copied: ${COPIED}"
echo "Missing: ${MISSING}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# List what was copied
if [ ${COPIED} -gt 0 ]; then
    echo "Contents:"
    ls "${OUTPUT_DIR}" | head -20
    if [ $(ls "${OUTPUT_DIR}" | wc -l) -gt 20 ]; then
        echo "... and $(($(ls "${OUTPUT_DIR}" | wc -l) - 20)) more"
    fi
fi

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

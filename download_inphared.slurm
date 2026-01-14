#!/bin/bash
#SBATCH --job-name=inphared_download
#SBATCH --output=inphared_download_%j.out
#SBATCH --error=inphared_download_%j.err
#SBATCH --time=24:00:00
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --gpus-per-node=1

# Download INPHARED phage genome database
# Includes all files, PHROGs HMM, with decompression

# Set output directory - modify as needed
OUTPUT_DIR="/data/lindseylm/INPHARED_DATA"

# Path to the download script - modify as needed
SCRIPT_DIR="/data/lindseylm/PROPHAGE_DETECTION_GENOME_BIOLOGY/lambda/05_metrics"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "INPHARED Database Download"
echo "=========================================="
echo "Start time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo "=========================================="

# Load any required modules (uncomment/modify as needed)
# module load python/3.9

# Run the download script
python ${SCRIPT_DIR}/download_inphared.py \
    --output ${OUTPUT_DIR} \
    --category all \
    --include-phrogs \
    --use-wget \
    --decompress

# Check exit status
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Download completed successfully!"
    echo "End time: $(date)"
    echo "=========================================="
else
    echo "=========================================="
    echo "Download failed with errors"
    echo "End time: $(date)"
    echo "=========================================="
    exit 1
fi

# List downloaded files
echo ""
echo "Downloaded files:"
ls -lh ${OUTPUT_DIR}

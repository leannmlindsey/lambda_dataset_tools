#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=phage_vs_gtdb
#SBATCH --output=phage_vs_gtdb_%A_%a.out
#SBATCH --error=phage_vs_gtdb_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --array=0-9

# BLAST INPHARED phages against GTDB selected bacteria
# Uses SLURM array to split into 10 parallel jobs
#
# Usage: sbatch blast_phage_vs_gtdb.slurm

# ============================================
# CONFIGURATION
# ============================================

# INPHARED phage genomes (query)
PHAGE_FASTA="/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa"

# GTDB selected bacteria database (subject)
BACTERIA_DB="/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected"

# Output directory
OUTPUT_DIR="/work/hdd/bfzj/llindsey1/prophage_blast_results"

# Number of splits (must match --array above)
NUM_SPLITS=10

# BLAST parameters
MIN_COVERAGE=80
NUM_THREADS=16

# ============================================

echo "=========================================="
echo "BLAST Phages vs GTDB Bacteria"
echo "=========================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "=========================================="

# Activate conda environment
source ~/.bashrc
conda activate lambda_tools

# Create output and temp directories
mkdir -p ${OUTPUT_DIR}
SPLIT_DIR="${OUTPUT_DIR}/splits"
mkdir -p ${SPLIT_DIR}

# Step 1: Split the input FASTA (only task 0 does this, others wait)
SPLIT_DONE="${OUTPUT_DIR}/.split_done"

if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    echo "Task 0: Splitting input FASTA into ${NUM_SPLITS} chunks..."

    # Count sequences
    NUM_SEQS=$(grep -c "^>" ${PHAGE_FASTA})
    echo "  Total sequences: ${NUM_SEQS}"

    # Calculate sequences per split
    SEQS_PER_SPLIT=$(( (NUM_SEQS + NUM_SPLITS - 1) / NUM_SPLITS ))
    echo "  Sequences per split: ~${SEQS_PER_SPLIT}"

    # Split using awk
    awk -v n=${SEQS_PER_SPLIT} -v dir="${SPLIT_DIR}" '
        BEGIN { file_num = 0; seq_count = 0 }
        /^>/ {
            if (seq_count >= n) {
                file_num++
                seq_count = 0
            }
            seq_count++
        }
        { print > (dir "/split_" file_num ".fasta") }
    ' ${PHAGE_FASTA}

    echo "  Split complete!"
    ls -lh ${SPLIT_DIR}/

    # Signal that splitting is done
    touch ${SPLIT_DONE}
else
    # Wait for task 0 to finish splitting
    echo "Waiting for split to complete..."
    while [ ! -f ${SPLIT_DONE} ]; do
        sleep 10
    done
    echo "Split complete, proceeding..."
fi

# Step 2: Run BLAST on this task's chunk
SPLIT_FILE="${SPLIT_DIR}/split_${SLURM_ARRAY_TASK_ID}.fasta"
OUTPUT_FILE="${OUTPUT_DIR}/blast_results_${SLURM_ARRAY_TASK_ID}.tsv"

if [ ! -f "${SPLIT_FILE}" ]; then
    echo "No split file for task ${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

NUM_QUERIES=$(grep -c "^>" ${SPLIT_FILE})
echo ""
echo "Running BLAST on ${NUM_QUERIES} sequences..."
echo "  Input: ${SPLIT_FILE}"
echo "  Output: ${OUTPUT_FILE}"

blastn \
    -db ${BACTERIA_DB} \
    -query ${SPLIT_FILE} \
    -out ${OUTPUT_FILE} \
    -outfmt "6 std" \
    -qcov_hsp_perc ${MIN_COVERAGE} \
    -max_target_seqs 500 \
    -num_threads ${NUM_THREADS}

echo ""
echo "BLAST complete!"
echo "  Hits: $(wc -l < ${OUTPUT_FILE})"

# Step 3: If this is the last task to finish, combine results
# (Simple approach: each task checks if all outputs exist)
sleep 5  # Small delay to avoid race condition

ALL_DONE=true
for i in $(seq 0 $((NUM_SPLITS - 1))); do
    if [ ! -f "${OUTPUT_DIR}/blast_results_${i}.tsv" ]; then
        ALL_DONE=false
        break
    fi
done

if [ "${ALL_DONE}" = true ] && [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    echo ""
    echo "All tasks complete. Combining results..."
    COMBINED="${OUTPUT_DIR}/blast_results_combined.tsv"
    cat ${OUTPUT_DIR}/blast_results_*.tsv > ${COMBINED}

    echo "Combined results: $(wc -l < ${COMBINED}) hits"
    echo "Output: ${COMBINED}"

    # Cleanup split files (optional)
    # rm -rf ${SPLIT_DIR}
    # rm ${SPLIT_DONE}
fi

echo ""
echo "=========================================="
echo "Task ${SLURM_ARRAY_TASK_ID} complete!"
echo "End time: $(date)"
echo "=========================================="

#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=02_contiguity_report
#SBATCH --output=02_contiguity_report_%j.out
#SBATCH --error=02_contiguity_report_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g

# Decision-support (read-only): report how many genera survive each bacterial
# isolate/assembly-contiguity policy, so the step-02 selection cutoff can be set
# from real numbers. Reads GTDB metadata via config.sh; writes a report here.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

source ~/.bashrc
conda activate "${CONDA_ENV}"

echo "Metadata: ${GTDB_METADATA}"
echo "Quality:  completeness >= ${GTDB_MIN_COMPLETENESS}, contamination <= ${GTDB_MAX_CONTAMINATION}"

python "${HERE}/analyze_bacteria_contiguity.py" \
    --metadata "${GTDB_METADATA}" \
    --min-completeness "${GTDB_MIN_COMPLETENESS}" \
    --max-contamination "${GTDB_MAX_CONTAMINATION}" \
    --output "${HERE}/contiguity_retention_report.txt"

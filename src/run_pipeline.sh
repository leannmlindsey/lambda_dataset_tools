#!/bin/bash
# =============================================================================
# LAMBDA master orchestrator.
#
# Submits each pipeline step to SLURM in dependency order (sbatch
# --dependency=afterok), so the whole dataset is reconstructed with one
# command and fixed seeds:
#
#   bash run_pipeline.sh
#
# Every step sources config.sh, so seeds/paths/thresholds live in one place.
# Job IDs are recorded under logs/ for provenance.
#
# Steps are activated as they are ported into src/ (see STEP_NOTES.md in each
# step directory). Ported + verified steps below are live; the rest are shown
# as the intended dependency wiring and will be uncommented when ported.
# =============================================================================

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/config.sh"

mkdir -p "${LOG_DIR}" "${MANIFEST_DIR}"
RUN_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).jobids"
echo "# LAMBDA pipeline submission $(date)" > "${RUN_LOG}"
echo "# seed=${SEED} data_dir=${DATA_DIR}" >> "${RUN_LOG}"

# submit <step_name> <slurm_script> [dependency_jobid]
# Echoes the new job id; records it in the run log.
submit() {
    local name="$1" script="$2" dep="${3:-}"
    local opts=(--parsable)
    [ -n "${dep}" ] && opts+=(--dependency=afterok:"${dep}")
    local jid
    jid="$(sbatch "${opts[@]}" "${script}" | cut -d';' -f1)"
    echo "  submitted ${name}: job ${jid}${dep:+ (after ${dep})}"
    printf '%-28s %s\n' "${name}" "${jid}${dep:+ after=${dep}}" >> "${RUN_LOG}"
    echo "${jid}"
}

echo "=========================================="
echo "LAMBDA pipeline — submitting jobs"
echo "  seed:     ${SEED}"
echo "  data dir: ${DATA_DIR}"
echo "  log:      ${RUN_LOG}"
echo "=========================================="

# ---- 01  Phage selection (vclust reps -> train/dev/test) -------------------
JID_PHAGE_SELECT="$(submit 01_phage_selection \
    "${SRC_DIR}/01_phage_selection/01_extract_representatives.slurm")"

# ---- 02  Bacteria selection (GTDB HQ reps, 1/genus, all ~15,865 genera) -----
JID_BACT_SELECT="$(submit 02_bacteria_selection \
    "${SRC_DIR}/02_bacteria_selection/02_select_gtdb.slurm")"

# ---- 03  Prophage BLAST (build DB -> BLAST array -> postprocess) ------------
# Alignment-based, INPHARED-only; produces per-contig masked intervals for step 04.
JID_BACTERIA_DB="$(submit 03a_bacteria_blastdb \
    "${SRC_DIR}/03_prophage_filtering/03a_create_bacteria_blastdb.slurm" "${JID_BACT_SELECT}")"
JID_PHAGE_BLAST="$(submit 03b_phage_vs_bacteria \
    "${SRC_DIR}/03_prophage_filtering/03b_blast_phage_vs_bacteria.slurm" "${JID_BACTERIA_DB}")"
JID_PROPHAGE="$(submit 03c_postprocess \
    "${SRC_DIR}/03_prophage_filtering/03c_postprocess_blast.slurm" "${JID_PHAGE_BLAST}")"

# ---- 04  Subsampling (phage + bacteria segments, 2k/4k/8k) -----------------
# TODO(port): fixes = all-contig sampling, reproducible per-accession RNG,
#             count valid genomes before distributing target segments.
# Depends on phage selection AND prophage filtering.
# JID_SUBSAMPLE="$(submit 04_subsampling \
#     "${SRC_DIR}/04_subsampling/04_subsample.slurm" "${JID_PROPHAGE}")"

# ---- 05  Fragment dereplication (remove dev/test leakage vs train) ---------
# PORTED + self-tested. Activates once 04 is ported (needs the segment FASTAs).
# JID_DEREP="$(submit 05_dereplication \
#     "${SRC_DIR}/05_dereplication/dereplicate_fragments.slurm" "${JID_SUBSAMPLE}")"

# ---- 06  Merge + training CSVs (full and dereplicated "hard" sets) ---------
# JID_MERGE="$(submit 06_merge "${SRC_DIR}/06_merge/06_merge.slurm" "${JID_DEREP}")"

# ---- 07  Controls (GC shuffle; planned dinucleotide-preserving shuffle) ----
# JID_CONTROLS="$(submit 07_controls "${SRC_DIR}/07_controls/07_controls.slurm" "${JID_DEREP}")"

# ---- 08  Benchmarks (phage-only, bacteria-only prophage-free, bacterial CDS)
# JID_BENCH="$(submit 08_benchmarks "${SRC_DIR}/08_benchmarks/08_benchmarks.slurm" "${JID_PROPHAGE}")"

# ---- 09  Package for release -----------------------------------------------
# JID_PACKAGE="$(submit 09_package "${SRC_DIR}/09_package/09_package.slurm" "${JID_MERGE}")"

echo ""
echo "Submitted. Job IDs recorded in: ${RUN_LOG}"
echo "Monitor with: squeue -u \$USER"

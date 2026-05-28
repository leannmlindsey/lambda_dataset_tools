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

# ---- 03  Prophage BLAST (bacteria-as-query vs INPHARED): prep -> BLAST array -> postprocess
# Alignment-based, INPHARED-only. Bacteria-as-query keeps -max_target_seqs from
# capping out for widespread temperate phages, so masking is complete.
JID_PREP="$(submit 03a_prep_blast \
    "${SRC_DIR}/03_prophage_filtering/03a_prep_blast.slurm" "${JID_BACT_SELECT}")"
JID_BLAST="$(submit 03b_bacteria_vs_phage \
    "${SRC_DIR}/03_prophage_filtering/03b_blast_bacteria_vs_phage.slurm" "${JID_PREP}")"
JID_PROPHAGE="$(submit 03c_postprocess \
    "${SRC_DIR}/03_prophage_filtering/03c_postprocess_blast.slurm" "${JID_BLAST}")"

# ---- 04  Subsampling (phage from step 01; bacteria mask-aware from steps 02+03)
JID_SUB_PHAGE="$(submit 04a_subsample_phage \
    "${SRC_DIR}/04_subsampling/04a_subsample_phage.slurm" "${JID_PHAGE_SELECT}")"
JID_SUB_BACT="$(submit 04b_subsample_bacteria \
    "${SRC_DIR}/04_subsampling/04b_subsample_bacteria.slurm" "${JID_PROPHAGE}")"

# ---- 05  Fragment dereplication (remove dev/test leakage vs train) ---------
# Depends on BOTH subsampling jobs (afterok on a colon-list waits for all).
JID_DEREP="$(submit 05_dereplication \
    "${SRC_DIR}/05_dereplication/dereplicate_fragments.slurm" "${JID_SUB_PHAGE}:${JID_SUB_BACT}")"

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

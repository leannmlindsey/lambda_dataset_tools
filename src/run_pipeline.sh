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

# ---- 03  Sample PHAGE fragments (sets the per-split count) ------------------
JID_PHAGE_FRAG="$(submit 03_sample_phage_fragments \
    "${SRC_DIR}/03_sample_phage_fragments/03_sample_phage_fragments.slurm" "${JID_PHAGE_SELECT}")"

# ---- 04  Sample BACTERIAL fragments (matches phage count) -------------------
JID_BACT_FRAG="$(submit 04_sample_bacteria_fragments \
    "${SRC_DIR}/04_sample_bacteria_fragments/04_sample_bacteria_fragments.slurm" \
    "${JID_PHAGE_FRAG}:${JID_BACT_SELECT}")"

# ---- 05  Prophage check on bacterial fragments + replacement sampling ------
# Per-fragment BLAST at 70%/200 -> for hits, full-genome BLAST -> replace from
# clean regions of the same genome. No phage-side masking (deliberate).
JID_05A="$(submit 05a_blast_fragments \
    "${SRC_DIR}/05_prophage_check/05a_blast_fragments.slurm" "${JID_BACT_FRAG}")"
JID_05B="$(submit 05b_blast_affected_genomes \
    "${SRC_DIR}/05_prophage_check/05b_blast_affected_genomes.slurm" "${JID_05A}")"
JID_PROPHAGE_CHECK="$(submit 05c_resample \
    "${SRC_DIR}/05_prophage_check/05c_resample.slurm" "${JID_05B}")"

# ---- 06  Merge into 1:1 train/val/test CSVs --------------------------------
JID_MERGE="$(submit 06_merge \
    "${SRC_DIR}/06_merge/06_merge.slurm" "${JID_PROPHAGE_CHECK}")"

# ---- 07  Shuffled-nucleotide control (test set only) -----------------------
JID_CONTROLS="$(submit 07_shuffled_control \
    "${SRC_DIR}/07_shuffled_control/07_shuffled_control.slurm" "${JID_MERGE}")"

# ---- 08  Bacteria-only FPR dataset (PHROG-sized) ---------------------------
JID_08A="$(submit 08a_fpr_select \
    "${SRC_DIR}/08_bacteria_only_fpr/08a_select_and_sample.slurm" "${JID_BACT_SELECT}")"
JID_08B="$(submit 08b_fpr_check \
    "${SRC_DIR}/08_bacteria_only_fpr/08b_prophage_check.slurm" "${JID_08A}:${JID_PROPHAGE_CHECK}")"
JID_BENCH="$(submit 08c_fpr_to_phrog \
    "${SRC_DIR}/08_bacteria_only_fpr/08c_to_phrog.slurm" "${JID_08B}")"

# ---- 09  Package for release -----------------------------------------------
# JID_PACKAGE="$(submit 09_package "${SRC_DIR}/09_package/09_package.slurm" "${JID_MERGE}")"

echo ""
echo "Submitted. Job IDs recorded in: ${RUN_LOG}"
echo "Monitor with: squeue -u \$USER"

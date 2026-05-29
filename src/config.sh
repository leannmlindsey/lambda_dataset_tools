#!/bin/bash
# =============================================================================
# LAMBDA pipeline — single source of truth for paths, seeds, and thresholds.
#
# Every step script sources this file, so changing a value here changes it
# everywhere. Nothing downstream should hardcode a path, seed, or threshold.
#
#   source /path/to/src/config.sh
#
# Override any value without editing this file by exporting it first, e.g.:
#   SEED=7 bash run_pipeline.sh
# (each assignment below uses :=, so a pre-exported value wins).
# =============================================================================

# ----------------------------------------------------------------------------
# Reproducibility
# ----------------------------------------------------------------------------
# Single global seed used by EVERY randomized step (selection, splitting,
# subsampling, shuffling). Do not introduce per-script seeds.
: "${SEED:=42}"

# Pin Python hash randomization so any hash()-based logic is reproducible
# across runs/processes (see review: GTDB subsampling previously used the
# salted built-in hash()). Exported so it reaches every child process.
export PYTHONHASHSEED=0

# Conda environments
: "${CONDA_ENV:=lambda_tools}"        # pandas/biopython/blast
: "${CONDA_ENV_BAKTA:=bakta_env}"     # bacterial annotation
: "${CONDA_ENV_GENOMAD:=genomad_env}" # prophage detection (planned, review #1)

# ----------------------------------------------------------------------------
# Pipeline locations
# ----------------------------------------------------------------------------
# Root of this repository on the cluster. Override if cloned elsewhere.
: "${PIPELINE_ROOT:=/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools}"
: "${SRC_DIR:=${PIPELINE_ROOT}/src}"

# All generated datasets are written under DATA_DIR (kept separate from the
# scripts so a clean run never clobbers source code).
: "${DATA_DIR:=${PIPELINE_ROOT}/lambda_data}"

# High-capacity scratch for large intermediates (BLAST dbs, raw hits).
: "${WORK_DIR:=/work/hdd/bfzj/llindsey1/lambda_work}"

# Per-step logs and output checksums (reproducibility evidence).
: "${LOG_DIR:=${PIPELINE_ROOT}/src/logs}"
: "${MANIFEST_DIR:=${PIPELINE_ROOT}/src/manifests}"

# ----------------------------------------------------------------------------
# Source data — INPHARED (phage)
# ----------------------------------------------------------------------------
: "${INPHARED_FASTA:=/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa}"
: "${INPHARED_VCLUST:=/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv}"
: "${INPHARED_METADATA:=/projects/bfzj/llindsey1/black_and_white/data/inphared/14Apr2025_data.tsv}"
: "${INPHARED_PHAROKKA:=/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka}"

# ----------------------------------------------------------------------------
# Source data — GTDB (bacteria)
# ----------------------------------------------------------------------------
# NOTE: the running scripts use the /u/... copy of the genomes; PATHS.md lists a
# /projects/... copy. Confirm which is canonical and keep only one.
: "${GTDB_GENOMES:=/u/llindsey1/llindsey/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226}"
: "${GTDB_METADATA:=/projects/bfzj/llindsey1/black_and_white/data/gtdb/metadata/bac120_metadata.tsv}"
: "${GTDB_RELEASE:=r226}"   # verify against source data (docs disagree r220 vs r226)

# External tool databases
: "${BAKTA_DB:=/u/llindsey1/scratch/bakta_db/db}"

# ----------------------------------------------------------------------------
# BLAST databases
# ----------------------------------------------------------------------------
: "${PHAGE_BLASTDB:=${WORK_DIR}/inphared_blastdb/inphared_phages}"   # built by step 03a if not present

# ----------------------------------------------------------------------------
# Split / selection parameters
# ----------------------------------------------------------------------------
: "${SPLIT_RATIO:=80:10:10}"          # train:dev:test, applied at cluster/genus level
: "${GTDB_MIN_COMPLETENESS:=95.0}"    # CheckM2
: "${GTDB_MAX_CONTAMINATION:=5.0}"    # CheckM2
# Step 02 output: selected bacterial accessions + metadata (one per genus, all genera)
: "${GTDB_SELECTION_DIR:=${DATA_DIR}/gtdb_selection}"

# ----------------------------------------------------------------------------
# Segment subsampling
# ----------------------------------------------------------------------------
: "${SEGMENT_LENGTHS:=2000 4000 8000}"   # the three datasets
: "${SAMPLE_PER_BP:=10000}"              # 1 segment per N bp of genome
: "${MIN_LENGTH_2K:=5000}"
: "${MIN_LENGTH_4K:=5000}"
: "${MIN_LENGTH_8K:=10000}"              # 8k excludes 5-10kb genomes (see methods)
# Target bacterial segment counts (matched to phage for ~50:50 balance)
: "${TARGET_SEGMENTS_TRAIN:=27000}"
: "${TARGET_SEGMENTS_DEV:=3400}"
: "${TARGET_SEGMENTS_TEST:=3400}"

# ----------------------------------------------------------------------------
# Prophage filtering (BLAST phages vs bacteria)
# ----------------------------------------------------------------------------
# Step 06 -- merged 1:1 train/val/test CSVs (phage + checked-bacteria, labeled, shuffled)
: "${MERGED_DIR:=${DATA_DIR}/merged_datasets}"

# Step 07 -- shuffled-nucleotide control derived from the merged test CSVs
: "${SHUFFLED_DIR:=${DATA_DIR}/shuffled_controls}"

# Step 05 -- prophage check on bacterial fragments (BLAST fragments vs INPHARED;
# replace any fragment that hits with a fresh sample from the same genome's
# non-phage regions).
: "${PROPHAGE_CHECK_MIN_IDENTITY:=70}"        # filter threshold for fragment hits AND for full-genome hits
: "${PROPHAGE_CHECK_MIN_ALIGN_LEN:=200}"
: "${PROPHAGE_CHECK_MASK_PAD:=500}"           # bp pad around masked phage regions in affected genomes
: "${PROPHAGE_CHECK_DIR:=${WORK_DIR}/prophage_check}"   # step 05 working dir (raw hits, masks, gathered FASTA)
# BLAST execution parameters (shared)
: "${BLAST_MAX_TARGET_SEQS:=2000}"
: "${BLAST_EVALUE:=1e-5}"

# ----------------------------------------------------------------------------
# Fragment-level dereplication (step 05) — removes dev/test leakage vs train
# ----------------------------------------------------------------------------
: "${DEREP_MIN_IDENTITY:=95}"            # % identity to count as leakage
: "${DEREP_MIN_COVERAGE:=50}"            # % query coverage per HSP
: "${DEREP_BLAST_PERC_IDENTITY:=90}"     # BLAST pre-filter (below the cut)

# Per-class segment dataset directories (one train/dev/test FASTA set each).
# Used by the dereplication step and the merge step.
PHAGE_DATASETS=(
    "${DATA_DIR}/inphared_dataset"
    "${DATA_DIR}/inphared_dataset_4k"
    "${DATA_DIR}/inphared_dataset_8k"
)
BACTERIA_DATASETS=(
    "${DATA_DIR}/gtdb_dataset"
    "${DATA_DIR}/gtdb_dataset_4k"
    "${DATA_DIR}/gtdb_dataset_8k"
)

# ----------------------------------------------------------------------------
# Compute resources (defaults for SBATCH; steps may override)
# ----------------------------------------------------------------------------
: "${SLURM_ACCOUNT:=bfzj-dtai-gh}"
: "${SLURM_PARTITION:=ghx4}"
: "${THREADS:=16}"

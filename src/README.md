# LAMBDA reproducible pipeline (`src/`)

A single, ordered, seed-controlled pipeline that reconstructs the LAMBDA
datasets from source data. Everything is driven by one config file and one
master script so the build is reproducible and auditable end-to-end.

## How it works

- **`config.sh`** — the only place paths, the random seed, and thresholds are
  defined. Every step sources it. `PYTHONHASHSEED=0` is exported here so even
  hash-based logic is reproducible.
- **`run_pipeline.sh`** — submits each step to SLURM in dependency order
  (`sbatch --dependency=afterok`). One command runs the whole chain.
- **`NN_step/`** — each numbered directory is one pipeline stage: its
  script(s), a `*.slurm` wrapper that sources `config.sh`, and `STEP_NOTES.md`
  documenting purpose, inputs/outputs, and any pending methodological fixes.
- **`logs/`**, **`manifests/`** — submitted job IDs per run, and (planned)
  md5 checksums of step outputs so a rerun can be verified bit-identical.

## Reproduce

```bash
# 1. Build the environment (versions pinned in environment.yml)
conda env create -f environment.yml      # creates lambda_tools
# (bakta_env / genomad_env are created separately; see environment.yml)

# 2. Edit config.sh if your paths differ, then launch the whole chain:
bash run_pipeline.sh

# ...or run a single step:
sbatch 01_phage_selection/01_extract_representatives.slurm
```

All randomness derives from `SEED` (default 42) in `config.sh`; override per
run with `SEED=<n> bash run_pipeline.sh`.

## Steps

| # | Stage | Status |
|---|-------|--------|
| 01 | Phage selection (vclust reps → train/dev/test, **random representative**) | ✅ ported |
| 02 | Bacteria selection (GTDB HQ reps, 1/genus, all ~15,865 genera) | ✅ ported (run-verified) |
| 03 | Prophage BLAST (**bacteria-as-query** vs INPHARED) → masked intervals + threshold/reverse-contam analyses | ✅ ported |
| 04 | Subsampling (phage proportional; bacteria mask-aware, all-contig ≥10 kb) | ✅ ported (self-tested) |
| 05 | Fragment dereplication (remove dev/test leakage vs train) | ✅ ported (self-tested) |
| 06 | Merge + training CSVs (full **and** dereplicated "hard" sets) | ⏳ to port |
| 07 | Controls (GC shuffle; + planned dinucleotide-preserving) | ⏳ to port |
| 08 | Benchmarks (phage-only, bacteria-only prophage-free, bacterial CDS) | ⏳ to port |
| 09 | Package for release | ⏳ to port |

Status reflects migration into this clean pipeline; the original top-level
scripts remain in place until each stage is ported and verified here. See
each step's `STEP_NOTES.md` and `../METHODS_DATA_SELECTION.md` for rationale.

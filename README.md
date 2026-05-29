# LAMBDA — LAnguage Model Bacteriophage Detection Assessment

Reproducible pipeline that builds the LAMBDA benchmark datasets — paired phage
vs bacterial sequence fragments at 2 kb / 4 kb / 8 kb, balanced 1:1, with
matching diagnostic controls (shuffled-nucleotide), a held-out bacteria-only
false-positive test set sized to PHROG, and PHROG-format CSV exports.

## Where the pipeline lives

Everything is under [`src/`](src/). Start there.

- **`src/README.md`** — overview of the pipeline + per-step status table.
- **`src/config.sh`** — single source of truth for paths, seeds, thresholds.
  Override any value by pre-exporting before running.
- **`src/run_pipeline.sh`** — master orchestrator: submits every step to SLURM
  in dependency order (`sbatch --dependency=afterok`). One command from a fresh
  clone reproduces the whole release.
- **`src/NN_step/`** — each numbered directory is one pipeline stage with its
  scripts, a `*.slurm` wrapper, and a `STEP_NOTES.md` describing inputs /
  outputs / design choices.
- **`src/environment.yml`** — pinned conda env (`lambda_tools`).
- **`src/logs/`**, **`src/manifests/`** — submitted job IDs per run, and (when
  enabled) MD5 checksums of step outputs.

## Steps (high-level)

| # | Stage |
|---|---|
| 01 | Phage selection (INPHARED vclust 95% ANI reps, 80/10/10 by cluster) |
| 02 | Bacteria selection (GTDB HQ reps, 1 per genus, 80/10/10 by genus) |
| 03 | Phage fragment sampling (2 k / 4 k / 8 k, sets the per-split count) |
| 04 | Bacterial fragment sampling (matches the per-split phage count → 1:1) |
| 05 | Prophage check + replacement (BLAST fragments vs INPHARED at 70%/200, replace hits from clean regions of the same genome) |
| 06 | Merge into 1:1 train / val / test CSVs (`segment_id, sequence, label, source`) |
| 07 | Shuffled-nucleotide control from the test set (composition preserved, sequence destroyed) |
| 08 | Bacteria-only FPR dataset, fresh pool, PHROG-sized + PHROG CSV format |
| 09 | Package release (Zenodo / HuggingFace-ready directory + checksums) |

See [`src/README.md`](src/README.md) for the full status table and each
[`src/NN_step/STEP_NOTES.md`](src/) for design + reproducibility details.

## Quick start

```bash
# 1. build the env (versions pinned)
conda env create -f src/environment.yml

# 2. edit src/config.sh if your paths differ, then launch the whole chain:
bash src/run_pipeline.sh

# ...or submit any one step on its own:
sbatch src/06_merge/06_merge.slurm
```

All randomness is driven by `SEED` in `src/config.sh` (default 42), with
`PYTHONHASHSEED=0` exported so hash-based logic is reproducible.

## Release output

The packaged release directory (step 09) looks like this:

```
LAMBDA_v1/
├── README.md                              # dataset card (Zenodo + HF)
├── train_val_test/{2k,4k,8k}/{train,val,test}.csv
├── shuffled_controls/{2k,4k,8k}/test_shuffled.csv
├── fpr_test/{2k,4k,8k}/bacteria_segments_*.csv
├── metadata/
│   ├── phage_accessions/{train,val,test}.txt
│   ├── bacteria_accessions/{train,val,test}.txt
│   ├── fpr_bacteria_accessions.txt
│   └── pipeline_version.txt
└── checksums.md5
```

CSV schemas:

- **train_val_test / shuffled_controls:** `segment_id, sequence, label, source`
  - `label`: 1 = phage, 0 = bacteria
  - `source`: `inphared` or `gtdb`
- **fpr_test:** `seq_id, start, end, sequence, label` (PHROG-compatible format)
  - `seq_id`: NCBI contig accession (e.g. `NC_001234.1`)
  - `start` / `end`: 1-based, inclusive
  - `label`: 0 throughout (bacteria-only)

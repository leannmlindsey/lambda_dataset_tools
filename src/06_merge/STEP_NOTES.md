# Step 06 — Merge + training CSVs

**Status:** ⏳ to port

**Purpose:** Combine phage (label 1) and bacteria (label 0) segments, shuffle,
and write `train/dev/test.csv` (`segment_id, sequence, label, source`).

**Implements (current, top-level):** `merge_and_shuffle.py`,
`create_training_csv.py`, `merge_datasets.slurm`

**Inputs (from config.sh):** `PHAGE_DATASETS`, `BACTERIA_DATASETS`, `SEED`

**Pending fixes (review):**
- Build **two** merged variants: `full` (from original dev/test) and `hard`
  (from the `*_derep` dev/test). Report the metric gap between them.
- Class balance is ~2% off 50:50 (phage > bacteria). Either trim to exact
  balance or standardize on balanced-accuracy / per-class metrics and state it.

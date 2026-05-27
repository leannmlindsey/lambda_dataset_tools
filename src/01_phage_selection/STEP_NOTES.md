# Step 01 — Phage selection

**Status:** ✅ ported & live in `run_pipeline.sh`

**Purpose:** Select one representative per vclust (~95% ANI) phage cluster and
split into train/dev/test by cluster (no cluster spans splits).

**Implements:** `extract_representatives.py` → `01_extract_representatives.slurm`

**Inputs (from config.sh):** `INPHARED_VCLUST`, `SPLIT_RATIO`, `SEED`
**Outputs:** `${DATA_DIR}/inphared_dataset/{train,dev,test}_accessions.txt`,
`*_with_clusters.tsv`, `summary.txt`

**Fixes applied vs original:**
- Representative is now chosen **at random per cluster** (seeded), not "first in
  cluster", to avoid bias from vclust's within-cluster ordering.

**Still open (analysis, not code):**
- Quantify max inter-split ANI/Mash similarity between representatives (review
  #5) — different 95% ANI clusters can still be highly similar. Fragment-level
  leakage from this is handled in step 05.

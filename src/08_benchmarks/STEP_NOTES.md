# Step 08 — Benchmark datasets (full genomes + annotations)

**Status:** ⏳ to port

**Purpose:** Annotated, full-genome benchmarks for error analysis.

**Implements (current, top-level):**
- Phage-only: `copy_pharokka_test_annotations.sh` (869 test genomes + Pharokka)
- Bacteria-only: `select_bacteria_only_balanced.py` **and**
  `find_prophage_free_bacteria.py` / `select_prophage_free_bacteria.sh`,
  `run_bakta_bacteria.slurm`
- Bacterial CDS: `extract_bacterial_cds.py` / `extract_bacterial_cds.sh`

**Inputs (from config.sh):** `INPHARED_PHAROKKA`, `GTDB_GENOMES`,
`GTDB_METADATA`, `BAKTA_DB` (`CONDA_ENV_BAKTA`), `SEED`

**Pending fixes (review #7):**
- **Two competing bacteria-only selectors.** Standardize on the **prophage-free**
  selector (`find_prophage_free_bacteria.py`); deprecate/remove the plain
  balanced one so the wrong set can't be used.
- Same independent-RNG fix as step 02 (`sample(random_state=SEED)` reuse).
- `extract_bacterial_cds.py` reseeds inside the per-genome loop; and excludes
  `integrase/transposase/recombinase` (also bacterial MGE genes) — conservative;
  document.
- "Prophage-free" here still means "no INPHARED hit"; re-check with the
  orthogonal detector from step 03.

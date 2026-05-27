# Step 05 — Fragment-level dereplication

**Status:** ✅ ported & self-tested (filter logic validated on synthetic input;
BLAST runs on the cluster)

**Purpose:** Remove dev/test fragments that are near-identical to any train
fragment — the leakage that genome/cluster splitting cannot catch (mosaic
phages; conserved bacterial genes like rRNA/housekeeping). The model only sees
fragments, so leakage is judged at the fragment level.

**Implements:** `filter_leaked_fragments.py` → `dereplicate_fragments.slurm`

**Inputs (from config.sh):** `PHAGE_DATASETS`, `BACTERIA_DATASETS`,
`DEREP_MIN_IDENTITY` (95), `DEREP_MIN_COVERAGE` (50), `DEREP_BLAST_PERC_IDENTITY`
**Outputs:** `<dataset>_derep/{train,dev,test}_segments.fasta`,
`*_leaked.tsv`, `*_dereplication_report.txt` — originals untouched (keep both
the full and the "hard" sets).

**Method:** build BLAST db of train fragments; BLAST dev/test; drop any query
fragment with a hit ≥95% identity AND ≥50% query coverage. Applied to **both**
classes.

**Depends on:** step 04 outputs (segment FASTAs). Activate in `run_pipeline.sh`
once 04 is ported.

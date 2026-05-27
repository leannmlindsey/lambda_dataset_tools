# Step 02 — Bacteria selection

**Status:** ✅ ported (all genera, RNG fix, provenance metadata)

**Purpose:** Select GTDB bacterial representatives passing CheckM2 quality,
one genome per genus, split train/dev/test by genus.

**Implements (current, top-level):** `select_gtdb_representatives.py`

**Inputs (from config.sh):** `GTDB_METADATA`, `GTDB_MIN_COMPLETENESS`,
`GTDB_MAX_CONTAMINATION`, `SPLIT_RATIO`, `SEED`
**Outputs:** `${DATA_DIR}/gtdb_dataset/{train,dev,test}_accessions.txt` + metadata

**Fixes applied (in top-level `select_gtdb_representatives.py`, to be carried into the port):**
- ✅ Per-genus selection now uses a single seeded `np.random.RandomState(SEED)`
  threaded across draws (was `random_state=SEED` reused per genus → not
  independent).

**Decided (retention analysis ran — see `contiguity_retention_report.txt`):**
- Keep the **FULL HQ set**: completeness ≥95%, contamination ≤5%, one per genus →
  **all ~15,865 genera**. NO complete/isolate/contiguity filter at selection.
  Rationale: the HQ set is 86% draft MAGs (median longest contig ~360 kb), so
  requiring complete isolates would drop genera 15,865 → ~2,251 (−86%). Because we
  sample only short fragments, assembly *contiguity* is irrelevant to *fragment*
  quality — quality is controlled downstream (steps 03/04), not by completeness.
- Add `ncbi_assembly_level`, `ncbi_genome_category`, `contig_count`, `n50_contigs`,
  `longest_contig` to the output metadata for provenance.
- Plasmids are **retained** (valid bacterial negatives); no plasmid-removal step
  anywhere. Phage-like content is removed by the step-03 phage BLAST.

**Pending before porting:**
- Document the one-per-genus design (flattens abundance; tilts toward rare/
  uncultured genera) as an explicit, justified choice.

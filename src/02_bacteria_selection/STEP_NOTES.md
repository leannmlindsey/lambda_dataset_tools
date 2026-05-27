# Step 02 — Bacteria selection

**Status:** ⏳ to port

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

**Decided:**
- Negatives are **chromosome-only** (plasmids/secondary contigs excluded; see
  step 04). This is only sound for closed/isolate genomes, so selection must
  control assembly contiguity — CheckM2 completeness ≠ single-contig assembly.

**Pending before porting:**
- Run `analyze_bacteria_contiguity.py` on the real `bac120_metadata.tsv` to get
  genus-retention numbers, then add an isolate/assembly-level filter
  (`ncbi_genome_category` excludes MAG/SAG; `ncbi_assembly_level`) at the chosen
  cutoff. Add `ncbi_assembly_level`, `ncbi_genome_category`, `contig_count` to
  the output metadata for provenance.
- Document the one-per-genus design (flattens abundance; tilts toward rare/
  uncultured genera) as an explicit, justified choice.

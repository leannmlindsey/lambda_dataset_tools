# Step 03 — Prophage filtering

**Status:** ⏳ to port (largest methodological changes pending)

**Purpose:** Ensure bacterial negatives do not contain phage sequence.

**Implements (current, top-level):**
- `blast_phage_vs_gtdb_selected_v3.slurm` — BLAST phages vs selected bacteria
- `filter_prophage_contaminated_v2.py` — **⚠️ referenced but NOT present in repo; locate or rewrite**
- `validate_prophage_removal_v3.slurm` — re-BLAST clean genomes
- `create_contig_to_genome_map.py` — map BLAST contig IDs → genome accessions

**Inputs (from config.sh):** `INPHARED_FASTA`, `BACTERIA_BLASTDB`,
`PROPHAGE_MIN_IDENTITY`, `PROPHAGE_MIN_ALIGN_LEN`, `PROPHAGE_MAX_TARGET_SEQS`

**Pending fixes (review — highest priority):**
- **#1 Residual prophages:** BLAST vs INPHARED only catches prophages ≥90% id to
  *known* phages. Add an orthogonal, reference-free detector (geNomad/CheckV/
  VIBRANT) to find & **mask** divergent prophages. (`CONDA_ENV_GENOMAD`)
- **Circular validation:** current validation re-uses the same BLAST threshold.
  Validate with the orthogonal tool instead.
- **#4 Taxonomic skew:** whole-genome removal of hit genomes may deplete
  prophage-rich taxa. Prefer masking regions over dropping genomes; report
  phylum/genus distribution of removed vs retained.
- **Threshold sensitivity:** sweep identity/length (e.g. 80/90/95 × 200/500/1000).

# Step 03 — Prophage filtering

**Status:** ✅ ported (03a DB build, 03b BLAST array, 03c filter + masking intervals + analyses)

**Purpose:** Ensure bacterial negatives do not contain phage sequence.

**Implements (current, top-level):**
- `blast_phage_vs_gtdb_selected_v3.slurm` — BLAST phages vs selected bacteria
- `filter_prophage_contaminated_v2.py` — **⚠️ referenced but NOT present in repo; locate or rewrite**
- `validate_prophage_removal_v3.slurm` — re-BLAST clean genomes
- `create_contig_to_genome_map.py` — map BLAST contig IDs → genome accessions

**Inputs (from config.sh):** `INPHARED_FASTA`, `BACTERIA_BLASTDB`,
`PROPHAGE_MIN_IDENTITY`, `PROPHAGE_MIN_ALIGN_LEN`, `PROPHAGE_MAX_TARGET_SEQS`

**Decided (alignment-only — NO learned classifier in data curation):**
- Prophage detection stays **BLAST vs INPHARED only**. No geNomad/CheckV/VIBRANT:
  using a learned classifier to curate training data bakes its decision boundary
  into the labels (circularity) and is attackable. Keep INPHARED (advisor: RefSeq
  viral / IMG/VR carry too much uncurated metagenomic sequence); do NOT broaden.
- **Filter at the SEGMENT level** (in step 04), not whole-genome: drop sampled
  segments overlapping a phage hit (≥90% id, ≥200 bp), reusing the existing BLAST
  results. This preserves the genome and its genus (fixes the taxonomic-skew
  concern from whole-genome removal) and stays fully alignment-based.
- **Plasmids are retained** (valid bacterial negatives); phage-like plasmid content
  is removed by the same BLAST filter. No plasmid-removal step.

**Accepted limitation (document, do not "fix" with a model):**
- BLAST-vs-INPHARED catches only prophages similar to known phages; divergent/novel
  prophages may remain in the negatives. This is the deliberate cost of model-free,
  reproducible curation. Optionally quantify the residual rate descriptively.

**Ported (03a / 03b / 03c):**
- 03a: build the bacteria BLAST DB from the step-02 selection (skips if it exists).
- 03b: BLAST INPHARED → bacteria (10-way array), emitting full coords (subject for
  masking, query+qlen for the reverse-contamination analysis).
- 03c: filter to 90% / 200 bp, build padded per-contig masked intervals
  (`make_masking_intervals.py`, step-04 input), and run the threshold-sensitivity
  sweep + reverse-contamination stat (`analyze_blast_hits.py`).
- `filter_prophage_contaminated_v2.py` is no longer needed — we mask regions at the
  segment level in step 04 rather than dropping whole genomes.
- The old circular "prophage-free validation" is dropped; step 04 instead QCs that no
  sampled segment overlaps a masked interval.

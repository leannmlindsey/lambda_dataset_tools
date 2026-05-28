# Step 04 — Bacterial fragment sampling

**Role in the workflow:** plain (unfiltered, unmasked) sampling of bacterial
fragments from the step-02 selection, sized to match the per-split phage
fragment counts produced by step 03 (1:1 balance, length-by-length and
split-by-split).

**No prophage check is performed here.** Step 05 BLASTs each sampled bacterial
fragment against INPHARED at ≥70% / ≥200 bp; any fragment with a hit is
discarded, the source genome is re-BLASTed, and a replacement is sampled from
the non-phage regions of the same genome.

**Sampling rules:**
- Sample from ALL contigs of length ≥ `max(10 kb, 2 × segment_length)` (so 2k/4k
  → 10 kb floor, 8k → 16 kb floor). Multiple contigs per genome participate.
- Reproducible: per-accession RNG seeded from a stable `hashlib.md5(SEED:acc)`
  hash, so threading order doesn't affect outputs.
- Hits the target count exactly by allocating across genomes that produce
  candidate segments, then redistributing any shortfall.

**Per-split target** is read live from `${PHAGE_DATASETS[i]}/${SPLIT}_segments.fasta`
at job time (`grep -c '^>'`). This guarantees 1:1 balance even if step 03's
counts change in a future rerun.

**Inputs** (from `config.sh`): `GTDB_GENOMES`, `GTDB_SELECTION_DIR`,
`PHAGE_DATASETS`, `BACTERIA_DATASETS`, `SEED`, `THREADS`.

**Outputs:** `${BACTERIA_DATASETS[i]}/{train,dev,test}_segments.fasta` (+ `.tsv`).

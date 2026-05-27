# Step 04 — Segment subsampling

**Status:** ⏳ to port

**Purpose:** Cut fixed-length (2k/4k/8k) non-overlapping segments from phage and
bacterial genomes; balance bacterial segment counts to the phage counts.

**Implements (current, top-level):**
- `subsample_segments.py` — phage (multi-FASTA)
- `subsample_gtdb_segments.py` — bacteria (gzipped nested dirs)

**Inputs (from config.sh):** `INPHARED_FASTA`, `GTDB_GENOMES`,
`SEGMENT_LENGTHS`, `SAMPLE_PER_BP`, `MIN_LENGTH_{2K,4K,8K}`,
`TARGET_SEGMENTS_{TRAIN,DEV,TEST}`, `SEED`

**Decided — multi-contig sampling (alignment-only, no classifier):**
- Sample from **all contigs ≥ a minimum length** (~10 kb; ~16 kb for the 8k set),
  proportional to contig length — NOT just the longest contig (the HQ set is 86%
  fragmented MAGs, so "longest contig" is an arbitrary fragment).
- Preserves all ~15,865 genera (`longest_contig` retention check confirms ~0 genera
  lost at a 10 kb floor). The ~10 kb floor is only tiny-contig hygiene.
- **Plasmids are retained** (valid bacterial negatives) — no contig classification.
- **Prophage filtering happens here, at the SEGMENT level:** after sampling, drop
  any segment overlapping a phage BLAST hit (≥90% id, ≥200 bp; reuse step-03
  results) and backfill to target counts. Keeps the genome/genus (no taxonomic
  skew); fully alignment-based.
- Code action: replace longest-contig-only selection in
  `subsample_gtdb_segments.py` with per-contig sampling over eligible (≥ min-length)
  contigs; fix the docstring (it wrongly says "concatenate"); keep the
  reproducible-RNG and undershoot fixes below.
- **#3 Reproducibility:** bacterial RNG uses `random.Random(seed + hash(acc))`;
  built-in `hash()` is salted per process → not reproducible. Use a stable hash
  (`hashlib`) or rely on `PYTHONHASHSEED=0` (now set in config.sh).
- **Undershoot:** bacterial `total_segments // len(accessions)` divides by all
  accessions incl. not-found/too-short → undershoots target. Count valid genomes
  first (as the phage script does).
- Methods text already corrected: coverage scales with length (2k≈20%, 4k≈40%,
  8k≈80%); 8k uses min-length 10kb.

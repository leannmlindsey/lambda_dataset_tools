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

**Decided — chromosome-only sampling (NOT a flaw):**
- Sample only the chromosome (the longest contig). This is a deliberate
  purity choice: plasmids/secondary contigs concentrate MGEs and phage-like
  elements that would blur the phage/bacteria boundary, and small contigs carry
  more contamination/misassembly. It is justified **only** because step 02
  restricts the bacterial set to closed/isolate genomes (so "longest contig" is
  genuinely the chromosome, not an arbitrary draft fragment).
- Code action: keep the longest-contig selection but (a) fix the docstring,
  which wrongly says "concatenate", to state chromosome-only by design, and
  (b) optionally assert the genome is single-replicon-ish (guard against
  silently sampling one fragment of a draft that slipped through step 02).
- **#3 Reproducibility:** bacterial RNG uses `random.Random(seed + hash(acc))`;
  built-in `hash()` is salted per process → not reproducible. Use a stable hash
  (`hashlib`) or rely on `PYTHONHASHSEED=0` (now set in config.sh).
- **Undershoot:** bacterial `total_segments // len(accessions)` divides by all
  accessions incl. not-found/too-short → undershoots target. Count valid genomes
  first (as the phage script does).
- Methods text already corrected: coverage scales with length (2k≈20%, 4k≈40%,
  8k≈80%); 8k uses min-length 10kb.

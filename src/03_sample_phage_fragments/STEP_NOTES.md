# Step 03 — Phage fragment sampling

**Role in the new workflow:** this step **sets the fragment count** for the
dataset; step 04 (bacterial fragment sampling) will be told to produce the same
per-split totals so the merged dataset is balanced 1:1.

**No filtering / no masking on the phage side.** This is a deliberate choice:
- We will not BLAST phage fragments against anything to "clean" them.
- We will not mask phage regions that look bacteria-like (integrases, AMG
  carriers, etc.) because the boundaries of what is "phage signal" vs "host-
  derived content carried by phage" are unclear and a masking choice could
  remove legitimate phage signal the model should learn.

**Sampling rules (per existing `subsample_segments.py`):**
- 1 fragment per `SAMPLE_PER_BP` bp of genome length (default 10,000), seeded.
- Non-overlapping placements within each genome.
- Min genome length = `MIN_LENGTH_{2K,4K,8K}` from config.

**Inputs** (from `config.sh`): `INPHARED_FASTA`, `PHAGE_DATASETS` (train/dev/test
accessions written by step 01 live in `PHAGE_DATASETS[0]`), `SEED`,
`SAMPLE_PER_BP`, `MIN_LENGTH_*`.

**Outputs:** `${PHAGE_DATASETS[i]}/{train,dev,test}_segments.fasta` (+ `.tsv`
metadata) for the three segment lengths.

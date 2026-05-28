# Step 05 — Prophage check + replacement (per user spec)

**What this step does, end to end:**
1. BLAST every bacterial fragment from step 04 against INPHARED at
   **≥70% identity, ≥200 bp**.
2. For each fragment with a hit → discard it.
3. For each *source genome* of a discarded fragment → BLAST the full genome
   against INPHARED at the same threshold and build per-contig padded masked
   intervals around the hits.
4. Sample exactly **one replacement fragment per discarded fragment** from the
   *same genome*, drawn from un-masked regions.
5. Write `${SPLIT}_segments_checked.fasta` + `.tsv` alongside step 04's
   originals (originals are never modified, so the chain is auditable).

**Phage fragments are not checked** — per user instruction, the boundaries of
"phage signal" are unclear (integrases, AMG-carriers, host-shared genes), and a
threshold-based phage-side filter could remove legitimate phage signal the model
should learn.

**Three slurms (run sequentially):**
- `05a_blast_fragments.slurm` — builds the INPHARED BLAST DB (cached, skipped on
  reruns), BLASTs each of the 9 fragment FASTAs (3 lengths × 3 splits), filters
  to the threshold, and runs `find_affected.py` to emit per-(length×split) hit
  IDs + replacement counts + the union of affected genome accessions.
- `05b_blast_affected_genomes.slurm` — gathers the affected genomes into a
  combined FASTA, BLASTs it against INPHARED at the same threshold, and writes
  `affected_masks.tsv` (padded merged per-contig intervals).
- `05c_resample.slurm` — drops hit fragments and samples replacements from clean
  regions, writing `${SPLIT}_segments_checked.fasta` + `.tsv`. Per-replacement
  RNG seeded by `hashlib(SEED:accession:seg:split:replace)` for reproducibility.

**Thresholds (`config.sh`):**
- `PROPHAGE_CHECK_MIN_IDENTITY=70`
- `PROPHAGE_CHECK_MIN_ALIGN_LEN=200`
- `PROPHAGE_CHECK_MASK_PAD=500`
- `PROPHAGE_CHECK_DIR=${WORK_DIR}/prophage_check`
- `PHAGE_BLASTDB`, `BLAST_MAX_TARGET_SEQS=2000`, `BLAST_EVALUE=1e-5`

**Outputs:**
- `${PROPHAGE_CHECK_DIR}/fragment_raw_${seg}_${split}.tsv` + `fragment_hits_*.tsv`
- `${PROPHAGE_CHECK_DIR}/hit_ids_${seg}_${split}.txt`, `replacement_counts_*.tsv`
- `${PROPHAGE_CHECK_DIR}/affected_accessions.txt`
- `${PROPHAGE_CHECK_DIR}/affected_bacteria.fna`, `affected_blast_*.tsv`,
  `affected_masks.tsv`
- `${BACTERIA_DATASETS[i]}/${SPLIT}_segments_checked.{fasta,tsv}` — what step 06
  reads for the bacterial side of the merged dataset.

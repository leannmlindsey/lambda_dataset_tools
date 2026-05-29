# Step 08 — Bacteria-only FPR dataset (PHROG-sized)

**Role:** a held-out bacterial-only test set for measuring false-positive rate.
Same size as the user's existing PHROG-based phage FNR set (~36,176 rows per
length), built from a **fresh bacterial pool** completely disjoint from the
training selection, and run through the same prophage check / replacement
pipeline as step 05 so the FPR negatives are as clean as the training negatives.

**Three slurms, sequential:**
- `08a_select_and_sample.slurm` — picks `FPR_N_GENOMES` HQ-filtered bacteria that
  are NOT in `${GTDB_SELECTION_DIR}/all_accessions.txt`, then plain-samples
  `FPR_TARGET_FRAGMENTS` fragments per length using step 04's sampler.
- `08b_prophage_check.slurm` — BLAST FPR fragments vs INPHARED at 70%/200, find
  affected genomes, BLAST those full genomes, build per-contig masks, and use
  step 05's `resample_replacements.py` to drop hit fragments + sample
  replacements from clean regions. Reuses the step-05 helpers verbatim.
- `08c_to_phrog.slurm` — convert the checked FPR FASTAs + metadata into
  PHROG-style CSVs:
  `seq_id, start, end, sequence, label`  (`label = 0` for every row).

**Final outputs:**
```
${FPR_DATA_DIR}/
├── 2k/bacteria_segments_2k.csv
├── 4k/bacteria_segments_4k.csv
└── 8k/bacteria_segments_8k.csv
```
`seq_id` = the bacterial CONTIG accession (e.g. `NC_001234.1`) — the direct
analog of the PHROG `MH572405` seq_id. `start`/`end` are 1-based inclusive.

**Config:**
- `FPR_SELECTION_DIR`, `FPR_DATA_DIR`, `FPR_CHECK_DIR`
- `FPR_TARGET_FRAGMENTS=36176`  (PHROG row count)
- `FPR_N_GENOMES=18000`         (fresh-pool size; comfortably yields the target at ~2 per genome)
- `PROPHAGE_CHECK_MIN_IDENTITY=70`, `PROPHAGE_CHECK_MIN_ALIGN_LEN=200`,
  `PROPHAGE_CHECK_MASK_PAD=500`  (same thresholds as step 05).

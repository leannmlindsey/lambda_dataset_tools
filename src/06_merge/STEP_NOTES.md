# Step 06 — Merge into 1:1 CSVs

**Role:** combine the step-03 phage fragments + step-05c checked-bacteria
fragments into one CSV per (length, split). This is the file the model trains /
evaluates on.

**Output structure:**
```
${MERGED_DIR}/
├── 2k/{train,dev,test}.csv
├── 4k/{train,dev,test}.csv
└── 8k/{train,dev,test}.csv
```

**CSV columns:** `segment_id, sequence, label, source`
- `label = 1` for phage, `label = 0` for bacteria
- `source = 'inphared'` (phage) or `'gtdb'` (bacteria)

**Properties:**
- Per file: exactly 50:50 phage/bacteria (counts come straight from step 03 / 05c,
  which are already 1:1 balanced — confirmed in the run-verification table).
- Per file: shuffled with a seeded RNG so reruns are bit-identical.
- Phage fragments are unchanged from step 03 (no per-phage filtering, per spec).
- Bacterial fragments are the *checked* ones from step 05c (any fragment that
  hit INPHARED at ≥70% / ≥200 bp was replaced from the same genome's non-phage
  regions; no genome dropped).

**Inputs (config):** `PHAGE_DATASETS`, `BACTERIA_DATASETS`, `MERGED_DIR`, `SEED`.

**Output validation:** the slurm fails if any of the 9 input FASTAs is missing
or if any CSV ends up with 0 rows. The Python prints
`balance: OK / MISMATCH` per file as a sanity check on the 1:1.

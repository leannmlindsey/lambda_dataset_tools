# Step 07 — Shuffled-nucleotide control

**Role:** generate a negative control test set. For each row of the merged
test CSV, shuffle the `sequence` column's nucleotides in place. Length,
mononucleotide composition (so GC content), and labels are preserved; all
k-mer / positional sequence structure is destroyed.

A model that performs well on this CSV is relying on composition (e.g. GC
content) rather than actual sequence — useful for diagnosing what the model
has actually learned.

**Output structure:**
```
${SHUFFLED_DIR}/
├── 2k/test_shuffled.csv
├── 4k/test_shuffled.csv
└── 8k/test_shuffled.csv
```
Same columns as the merged CSVs (`segment_id, sequence, label, source`); only
the `sequence` column changes.

**Reproducibility:** each row's shuffle uses an RNG seeded from
`hashlib(SEED:segment_id)`, so reruns produce bit-identical output regardless
of row order or threading.

**Per spec:** only the *test* split is shuffled (train/dev untouched).

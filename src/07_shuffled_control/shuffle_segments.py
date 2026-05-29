#!/usr/bin/env python3
"""
Step 07 -- shuffled-nucleotide control.

Reads a merged CSV (from step 06), per-row shuffles the `sequence` column's
nucleotides, and writes a new CSV with the same columns and row order. Length,
mononucleotide composition (so GC content), and labels are preserved; all
k-mer / positional sequence signal is destroyed. Used as a control: a model
that performs well on this CSV is relying on composition rather than sequence.

Each row's shuffle uses an RNG seeded from `hashlib(SEED:segment_id)`, so the
result is reproducible even if the input row order changes.
"""

import argparse
import csv
import hashlib
import random
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Per-row nucleotide shuffle of a merged CSV")
    p.add_argument("--input-csv", required=True)
    p.add_argument("--output-csv", required=True)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def stable_rng(global_seed, segment_id):
    h = hashlib.md5(f"{global_seed}:{segment_id}".encode()).hexdigest()
    return random.Random(int(h[:8], 16))


def main():
    args = parse_args()
    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(args.input_csv) as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            seq = list(row["sequence"])
            stable_rng(args.seed, row["segment_id"]).shuffle(seq)
            row["sequence"] = "".join(seq)
            writer.writerow(row)
            n += 1

    print(f"  shuffled {n:,} rows  ->  {out_path}")


if __name__ == "__main__":
    main()

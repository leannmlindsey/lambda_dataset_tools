#!/usr/bin/env python3
"""
Step 06 -- merge phage + checked-bacteria fragments into one CSV per (length,
split), labeled and shuffled (seeded).

Columns:  segment_id, sequence, label, source
  label  = 1 for phage, 0 for bacteria
  source = 'inphared' (phage) or 'gtdb' (bacteria)

Inputs (per length, split):
  --phage-fasta      step-03 phage fragments
  --bacteria-fasta   step-05 *_segments_checked.fasta  (post prophage check)

The merged FASTA is shuffled with a seeded RNG so reruns are bit-identical.
"""

import argparse
import csv
import random
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Merge phage + checked-bacteria fragments into a labeled CSV")
    p.add_argument("--phage-fasta", required=True)
    p.add_argument("--bacteria-fasta", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def read_fasta(path):
    """Yield (segment_id, sequence_str) tuples in file order."""
    cur_id = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    yield cur_id, "".join(buf)
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id is not None:
            yield cur_id, "".join(buf)


def main():
    args = parse_args()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    n_phage = 0
    for seg_id, seq in read_fasta(args.phage_fasta):
        rows.append([seg_id, seq, 1, "inphared"])
        n_phage += 1
    n_bact = 0
    for seg_id, seq in read_fasta(args.bacteria_fasta):
        rows.append([seg_id, seq, 0, "gtdb"])
        n_bact += 1

    rng = random.Random(args.seed)
    rng.shuffle(rows)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["segment_id", "sequence", "label", "source"])
        w.writerows(rows)

    bal = "OK" if n_phage == n_bact else f"MISMATCH (phage={n_phage} vs bacteria={n_bact})"
    print(f"  phage: {n_phage:,}  bacteria: {n_bact:,}  total: {len(rows):,}  balance: {bal}")
    print(f"  -> {out_path}")


if __name__ == "__main__":
    main()

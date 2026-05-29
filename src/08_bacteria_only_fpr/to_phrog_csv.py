#!/usr/bin/env python3
"""
Step 08c -- convert the checked FPR bacterial fragments into PHROG-style CSVs.

PHROG CSV columns (matches the user's existing phage_segments_2k_1k.csv):
    seq_id, start, end, sequence, label

- seq_id  : the bacterial CONTIG accession (NCBI-style, e.g. NC_001234.1) --
            this is the analog of MH572405 in the PHROG phage CSV.
- start/end: 1-based, inclusive (PHROG convention).
- sequence: the segment nucleotides.
- label   : 0 (bacteria) for every row in this dataset.
"""

import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Convert checked FPR bacterial fragments to PHROG CSV format")
    p.add_argument("--fasta", required=True, help="checked-bacteria FASTA from step 08b's resample")
    p.add_argument("--metadata", required=True,
                   help="matching metadata TSV (segment_id, accession, contig, start, end, length)")
    p.add_argument("--output", required=True)
    p.add_argument("--label", type=int, default=0)
    return p.parse_args()


def read_fasta(path):
    cur_id, buf = None, []
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
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Load metadata: segment_id -> (contig, start, end). The metadata uses 0-based
    # half-open from Python slicing; PHROG uses 1-based inclusive. Convert.
    meta = {}
    with open(args.metadata) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            seg_id, _acc, contig, start, end = parts[0], parts[1], parts[2], int(parts[3]), int(parts[4])
            meta[seg_id] = (contig, start + 1, end)   # 0-based half-open -> 1-based inclusive

    n_written = 0
    n_missing = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["seq_id", "start", "end", "sequence", "label"])
        for seg_id, seq in read_fasta(args.fasta):
            if seg_id not in meta:
                n_missing += 1
                continue
            contig, start_1, end_1 = meta[seg_id]
            w.writerow([contig, start_1, end_1, seq, args.label])
            n_written += 1

    print(f"  wrote {n_written:,} rows to {out}")
    if n_missing:
        print(f"  WARNING: {n_missing} FASTA records had no metadata entry (skipped)")


if __name__ == "__main__":
    main()

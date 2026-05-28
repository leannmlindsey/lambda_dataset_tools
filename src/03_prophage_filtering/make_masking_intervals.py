#!/usr/bin/env python3
"""
Build per-contig masked intervals from the bacteria-vs-phage BLAST hits.

Step 03 output -> step 04 input. The BLAST in this pipeline is bacteria-as-query
(see step 03 design rationale), so bacterial sequence coordinates live in the
QUERY columns, not the subject columns:

  BLAST outfmt 6:
    0 qseqid  1 sseqid  2 pident  3 length  4 qstart  5 qend
    6 sstart  7 send    8 qlen    9 slen    10 evalue 11 bitscore
  Bacterial side (used here): qseqid (col 0), qstart/qend (4/5), qlen (8).
  Phage side (ignored here):  sseqid (col 1), sstart/send (6/7), slen (9).

For every hit meeting the identity/length thresholds, the hit's interval on the
bacterial contig is padded by --pad bp on each side, clamped to [1, qlen], and
merged with overlapping/adjacent intervals per contig. Step 04 then samples
bacterial segments only from the UN-masked remainder, never crossing a masked
region.
"""

import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Build per-contig masked intervals (bacteria-as-query BLAST)")
    p.add_argument("--blast-hits", "-b", required=True, help="BLAST tabular hits (bacteria=query, phage=subject)")
    p.add_argument("--output", "-o", required=True, help="Output masked-intervals TSV")
    p.add_argument("--min-identity", type=float, default=90.0)
    p.add_argument("--min-length", type=int, default=200)
    p.add_argument("--pad", type=int, default=500)
    return p.parse_args()


def merge_intervals(intervals):
    """Merge overlapping/adjacent (gap<=1) 1-based inclusive intervals."""
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def main():
    args = parse_args()
    by_contig = defaultdict(list)
    scanned = 0
    kept = 0
    with open(args.blast_hits) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            scanned += 1
            try:
                contig = parts[0]          # qseqid = bacterial contig
                pident = float(parts[2])
                length = int(parts[3])
                qstart = int(parts[4])     # bacterial coord
                qend = int(parts[5])
                qlen = int(parts[8])       # bacterial contig length
            except (ValueError, IndexError):
                continue
            if pident < args.min_identity or length < args.min_length:
                continue
            s = max(1, min(qstart, qend) - args.pad)
            e = min(qlen, max(qstart, qend) + args.pad)
            by_contig[contig].append((s, e))
            kept += 1

    total_bp = 0
    n_intervals = 0
    with open(args.output, "w") as out:
        out.write("contig_id\tstart\tend\n")
        for contig in sorted(by_contig):
            for s, e in merge_intervals(by_contig[contig]):
                out.write(f"{contig}\t{s}\t{e}\n")
                total_bp += e - s + 1
                n_intervals += 1

    print("=" * 60)
    print("Bacterial masked intervals (bacteria-vs-phage BLAST)")
    print("=" * 60)
    print(f"Hits scanned:                 {scanned:,}")
    print(f"Hits passing >= {args.min_identity}% / >= {args.min_length} bp: {kept:,}")
    print(f"Padding each side:            {args.pad} bp")
    print(f"Contigs with masked regions:  {len(by_contig):,}")
    print(f"Merged masked intervals:      {n_intervals:,}")
    print(f"Total masked bp:              {total_bp:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

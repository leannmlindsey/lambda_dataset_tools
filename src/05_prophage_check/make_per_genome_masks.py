#!/usr/bin/env python3
"""
Build per-contig masked intervals from the affected-genome BLAST hits (step 05b).
Bacteria-as-query: bacterial coords live on the QUERY side (qseqid / qstart /
qend / qlen). Same logic and format as the earlier make_masking_intervals.py.

Output: contig_id <tab> start <tab> end  (1-based, inclusive, merged, padded).
This is what step 05c uses to sample replacement fragments inside un-masked
regions of the affected genomes.
"""

import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Per-contig masked intervals from affected-genome BLAST hits")
    p.add_argument("--blast-hits", "-b", required=True, help="Filtered bacteria-vs-phage hits")
    p.add_argument("--output", "-o", required=True)
    p.add_argument("--min-identity", type=float, default=70.0)
    p.add_argument("--min-length", type=int, default=200)
    p.add_argument("--pad", type=int, default=500)
    return p.parse_args()


def merge(intervals):
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
    kept = 0
    with open(args.blast_hits) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                contig = parts[0]
                pident = float(parts[2])
                length = int(parts[3])
                qstart = int(parts[4])
                qend = int(parts[5])
                qlen = int(parts[8])
            except (ValueError, IndexError):
                continue
            if pident < args.min_identity or length < args.min_length:
                continue
            s = max(1, min(qstart, qend) - args.pad)
            e = min(qlen, max(qstart, qend) + args.pad)
            by_contig[contig].append((s, e))
            kept += 1

    n_intervals = 0
    total_bp = 0
    with open(args.output, "w") as out:
        out.write("contig_id\tstart\tend\n")
        for contig in sorted(by_contig):
            for s, e in merge(by_contig[contig]):
                out.write(f"{contig}\t{s}\t{e}\n")
                n_intervals += 1
                total_bp += e - s + 1

    print(f"Affected-genome masks: hits passing filter = {kept:,}; "
          f"contigs = {len(by_contig):,}; merged intervals = {n_intervals:,}; "
          f"total bp = {total_bp:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

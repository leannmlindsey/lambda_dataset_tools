#!/usr/bin/env python3
"""
Build per-contig masked intervals from phage-vs-bacteria BLAST hits.

Step 03 output -> step 04 input. For every BLAST hit meeting the identity and
alignment-length thresholds, the hit's interval on the bacterial contig (the
subject) is recorded, padded by --pad bp on each side (prophage boundaries
inferred from alignment to known phages are approximate), clamped to the contig,
and merged with overlapping/adjacent intervals. Step 04 then samples bacterial
segments only from the UN-masked remainder, and never rejoins flanks across a
masked region (no chimeric junctions).

Input: BLAST tabular (outfmt 6) with columns
  0 qseqid  1 sseqid  2 pident  3 length  4 qstart  5 qend
  6 sstart  7 send    8 qlen    9 slen    10 evalue 11 bitscore
(sseqid = bacterial contig; sstart/send = coords on it; slen = contig length)

Output: TSV  contig_id <tab> start <tab> end   (1-based, inclusive, merged)

Usage:
    python make_masking_intervals.py \
        --blast-hits phage_vs_bacteria.tsv \
        --output masked_intervals.tsv \
        --min-identity 90 --min-length 200 --pad 500
"""

import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Merge phage BLAST hits into per-contig masked intervals")
    p.add_argument("--blast-hits", "-b", required=True, help="BLAST tabular hits (outfmt 6, see header)")
    p.add_argument("--output", "-o", required=True, help="Output masked-intervals TSV")
    p.add_argument("--min-identity", type=float, default=90.0)
    p.add_argument("--min-length", type=int, default=200)
    p.add_argument("--pad", type=int, default=500, help="bp padding each side of a hit (default: 500)")
    return p.parse_args()


def merge_intervals(intervals):
    """Merge overlapping or adjacent (gap<=1) 1-based inclusive intervals."""
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def main():
    args = parse_args()

    by_contig = defaultdict(list)
    contig_len = {}
    kept = 0
    scanned = 0

    with open(args.blast_hits) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            scanned += 1
            try:
                contig = parts[1]
                pident = float(parts[2])
                length = int(parts[3])
                sstart = int(parts[6])
                send = int(parts[7])
                slen = int(parts[9])
            except (ValueError, IndexError):
                continue

            if pident < args.min_identity or length < args.min_length:
                continue

            start = min(sstart, send)
            end = max(sstart, send)
            # Pad and clamp to the contig.
            start = max(1, start - args.pad)
            end = min(slen, end + args.pad)
            by_contig[contig].append((start, end))
            contig_len[contig] = slen
            kept += 1

    total_masked_bp = 0
    n_intervals = 0
    with open(args.output, "w") as out:
        out.write("contig_id\tstart\tend\n")
        for contig in sorted(by_contig):
            for start, end in merge_intervals(by_contig[contig]):
                out.write(f"{contig}\t{start}\t{end}\n")
                total_masked_bp += end - start + 1
                n_intervals += 1

    print("=" * 60)
    print("Masked intervals from phage BLAST hits")
    print("=" * 60)
    print(f"Hits scanned:                 {scanned:,}")
    print(f"Hits passing >= {args.min_identity}% / >= {args.min_length} bp: {kept:,}")
    print(f"Padding each side:            {args.pad} bp")
    print(f"Contigs with masked regions:  {len(by_contig):,}")
    print(f"Merged masked intervals:      {n_intervals:,}")
    print(f"Total masked bp:              {total_masked_bp:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

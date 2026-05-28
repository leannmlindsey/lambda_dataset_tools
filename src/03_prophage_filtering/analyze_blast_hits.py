#!/usr/bin/env python3
"""
Descriptive analyses of the bacteria-vs-phage BLAST (read-only; the data is not
filtered or modified -- these summaries justify and characterize the masking).

BLAST direction in this pipeline: bacteria-as-query, INPHARED-as-subject. So:
  qseqid = bacterial contig, qstart/qend/qlen = bacterial coords/length
  sseqid = phage,            sstart/send/slen = phage    coords/length

1. THRESHOLD SENSITIVITY -- for identity in {80,90,95} x alignment length in
   {200,500,1000}, report passing hits, bacterial bp that would be masked
   (merged per contig, no padding), and contigs affected.

2. REVERSE CONTAMINATION (phage side) -- at the chosen id/len, what fraction of
   phages have a strong bacterial hit, and how much of each phage is covered
   (merged phage intervals / phage length). Descriptive; no phage-side masking.
"""

import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Threshold-sensitivity + reverse-contamination summaries (bacteria-as-query)")
    p.add_argument("--blast-hits", "-b", required=True)
    p.add_argument("--rev-identity", type=float, default=90.0)
    p.add_argument("--rev-length", type=int, default=200)
    p.add_argument("--output", "-o", default=None)
    return p.parse_args()


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def main():
    args = parse_args()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    # Load into tuples. Indices match the bacteria-as-query convention.
    hits = []  # (contig, phage, pident, length, qstart, qend, sstart, send, slen)
    slen_by_phage = {}
    with open(args.blast_hits) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            try:
                rec = (p[0], p[1], float(p[2]), int(p[3]),
                       int(p[4]), int(p[5]), int(p[6]), int(p[7]), int(p[9]))
            except (ValueError, IndexError):
                continue
            hits.append(rec)
            slen_by_phage[p[1]] = rec[8]

    emit("=" * 70)
    emit("Bacteria-vs-phage BLAST analysis")
    emit("=" * 70)
    emit(f"Total raw hits: {len(hits):,}")
    emit(f"Distinct phages with any hit: {len({h[1] for h in hits}):,}")
    emit(f"Distinct bacterial contigs with any hit: {len({h[0] for h in hits}):,}")
    emit("")

    # ---- 1. Threshold sensitivity (bacterial bp that would be masked) ----
    emit("Threshold sensitivity (bacterial sequence that would be masked):")
    emit(f"{'identity':>9}{'min_len':>9}{'hits':>12}{'contigs':>10}{'masked_bp':>15}")
    emit("-" * 70)
    for ident in (80.0, 90.0, 95.0):
        for mlen in (200, 500, 1000):
            by_contig = defaultdict(list)
            n = 0
            for contig, phage, pid, ln, qs, qe, ss, se, slen in hits:
                if pid >= ident and ln >= mlen:
                    by_contig[contig].append((min(qs, qe), max(qs, qe)))
                    n += 1
            masked_bp = sum(e - s + 1 for c in by_contig
                            for s, e in merge_intervals(by_contig[c]))
            emit(f"{ident:>8.0f}%{mlen:>9,}{n:>12,}{len(by_contig):>10,}{masked_bp:>15,}")
    emit("")

    # ---- 2. Reverse contamination (phage side) ----
    emit(f"Reverse contamination at >= {args.rev_identity}% / >= {args.rev_length} bp "
         f"(host-like content in phages):")
    by_phage = defaultdict(list)
    for contig, phage, pid, ln, qs, qe, ss, se, slen in hits:
        if pid >= args.rev_identity and ln >= args.rev_length:
            by_phage[phage].append((min(ss, se), max(ss, se)))

    n_phage_total = len(slen_by_phage)
    n_phage_hit = len(by_phage)
    coverages = []
    for phage, ivs in by_phage.items():
        covered = sum(e - s + 1 for s, e in merge_intervals(ivs))
        slen = slen_by_phage.get(phage, 0)
        if slen > 0:
            coverages.append(min(100.0, covered / slen * 100))
    coverages.sort()

    def pct(p):
        if not coverages:
            return 0.0
        i = min(len(coverages) - 1, int(p / 100 * len(coverages)))
        return coverages[i]

    frac = (n_phage_hit / n_phage_total * 100) if n_phage_total else 0.0
    emit(f"  Phages with a strong bacterial hit: {n_phage_hit:,} / {n_phage_total:,} ({frac:.1f}%)")
    if coverages:
        emit("  Of those, % of phage length covered by bacterial hits:")
        emit(f"    median={pct(50):.2f}%  mean={sum(coverages)/len(coverages):.2f}%  "
             f"p90={pct(90):.2f}%  max={max(coverages):.2f}%")
    emit("  (Descriptive only -- no phage-side masking is applied.)")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nReport written to {args.output}")


if __name__ == "__main__":
    main()

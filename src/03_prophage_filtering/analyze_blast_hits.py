#!/usr/bin/env python3
"""
Two descriptive analyses of the phage-vs-bacteria BLAST (read-only; neither
filters the data — they justify/characterize the masking choices).

1. THRESHOLD SENSITIVITY — for identity in {80,90,95} x alignment length in
   {200,500,1000}, report passing hits, bacterial bp that would be masked
   (merged per contig, no padding), and contigs affected. Shows whether the
   chosen 90% / 200 bp cut is stable (small moves -> small changes).

2. REVERSE CONTAMINATION (phage side) — at the chosen id/len, what fraction of
   phage genomes have a strong bacterial hit, and how much of each phage is
   covered (merged query intervals / qlen). Quantifies host-like content in the
   positives; we do NOT act on it (documented limitation).

Input: BLAST tabular (outfmt 6) with columns
  0 qseqid 1 sseqid 2 pident 3 length 4 qstart 5 qend
  6 sstart 7 send   8 qlen   9 slen   10 evalue 11 bitscore

Usage:
    python analyze_blast_hits.py --blast-hits raw_hits.tsv \
        [--rev-identity 90 --rev-length 200] [--output report.txt]
"""

import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Threshold-sensitivity + reverse-contamination summaries")
    p.add_argument("--blast-hits", "-b", required=True, help="RAW BLAST tabular hits (outfmt 6)")
    p.add_argument("--rev-identity", type=float, default=90.0)
    p.add_argument("--rev-length", type=int, default=200)
    p.add_argument("--output", "-o", default=None, help="Optional report file")
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

    # Load once into memory as tuples.
    hits = []  # (qseqid, sseqid, pident, length, qstart, qend, sstart, send, qlen)
    qlen = {}
    with open(args.blast_hits) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            try:
                rec = (p[0], p[1], float(p[2]), int(p[3]),
                       int(p[4]), int(p[5]), int(p[6]), int(p[7]), int(p[8]))
            except (ValueError, IndexError):
                continue
            hits.append(rec)
            qlen[p[0]] = rec[8]

    emit("=" * 70)
    emit("Phage-vs-bacteria BLAST analysis")
    emit("=" * 70)
    emit(f"Total raw hits: {len(hits):,}")
    emit(f"Distinct phages with any hit: {len({h[0] for h in hits}):,}")
    emit("")

    # ---- 1. Threshold sensitivity (bacterial masked bp) ----
    emit("Threshold sensitivity (bacterial sequence that would be masked):")
    emit(f"{'identity':>9}{'min_len':>9}{'hits':>12}{'contigs':>10}{'masked_bp':>15}")
    emit("-" * 70)
    for ident in (80.0, 90.0, 95.0):
        for mlen in (200, 500, 1000):
            by_contig = defaultdict(list)
            n = 0
            for q, s, pid, ln, qs, qe, ss, se, ql in hits:
                if pid >= ident and ln >= mlen:
                    by_contig[s].append((min(ss, se), max(ss, se)))
                    n += 1
            masked_bp = sum(e - st + 1 for s in by_contig
                            for st, e in merge_intervals(by_contig[s]))
            emit(f"{ident:>8.0f}%{mlen:>9,}{n:>12,}{len(by_contig):>10,}{masked_bp:>15,}")
    emit("")

    # ---- 2. Reverse contamination (phage side) ----
    emit(f"Reverse contamination at >= {args.rev_identity}% / >= {args.rev_length} bp "
         f"(host-like content in phages):")
    by_phage = defaultdict(list)
    for q, s, pid, ln, qs, qe, ss, se, ql in hits:
        if pid >= args.rev_identity and ln >= args.rev_length:
            by_phage[q].append((min(qs, qe), max(qs, qe)))

    n_phage_total = len(qlen)
    n_phage_hit = len(by_phage)
    coverages = []
    for q, ivs in by_phage.items():
        covered = sum(e - s + 1 for s, e in merge_intervals(ivs))
        ql = qlen.get(q, 0)
        if ql > 0:
            coverages.append(min(100.0, covered / ql * 100))
    coverages.sort()

    def pct(p):
        if not coverages:
            return 0.0
        i = min(len(coverages) - 1, int(p / 100 * len(coverages)))
        return coverages[i]

    frac = (n_phage_hit / n_phage_total * 100) if n_phage_total else 0.0
    emit(f"  Phages with a strong bacterial hit: {n_phage_hit:,} / {n_phage_total:,} ({frac:.1f}%)")
    if coverages:
        mean_cov = sum(coverages) / len(coverages)
        emit(f"  Of those, % of phage length covered by bacterial hits:")
        emit(f"    median={pct(50):.2f}%  mean={mean_cov:.2f}%  p90={pct(90):.2f}%  max={max(coverages):.2f}%")
    emit("  (Descriptive only — no phage-side masking is applied.)")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nReport written to {args.output}")


if __name__ == "__main__":
    main()

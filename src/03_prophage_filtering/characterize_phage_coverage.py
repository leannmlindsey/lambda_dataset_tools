#!/usr/bin/env python3
"""
Characterize the phage bacterial-coverage distribution from the bacteria-vs-phage
BLAST, joined with INPHARED metadata.

For each phage we compute its fraction-of-length covered by strong bacterial
hits and report the coverage distribution + composition of the high-coverage
tail (BACPHLIP lifestyle, family, host). The expected biology: the tail should
be enriched for temperate lifestyle and known temperate-phage taxa that have
integrated copies across many bacterial chromosomes -- i.e., legitimate phages,
not mislabeled non-phage entries (INPHARED is curated).

BLAST direction: bacteria-as-query, phages-as-subject. Phage coords are on the
SUBJECT side:
    1 sseqid = phage accession
    6 sstart / 7 send = phage coords
    9 slen = phage genome length

Outputs an optional report and a tail TSV (accessions sorted by coverage
descending) -- the input you'd feed to a GBK/Pharokka gene-content lookup.
"""

import argparse
from collections import defaultdict

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Per-phage coverage join with INPHARED metadata (bacteria-as-query)")
    p.add_argument("--blast-hits", "-b", required=True)
    p.add_argument("--inphared-metadata", "-m", required=True)
    p.add_argument("--min-identity", type=float, default=90.0)
    p.add_argument("--min-length", type=int, default=200)
    p.add_argument("--coverage-tail", type=float, default=75.0)
    p.add_argument("--output", "-o", default=None)
    p.add_argument("--output-tail", default=None)
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


def find_col(df, keywords):
    for c in df.columns:
        cl = c.lower()
        if any(kw in cl for kw in keywords):
            return c
    return None


def strip_version(acc):
    return str(acc).split(".")[0]


def main():
    args = parse_args()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    # 1. Per-phage bacterial coverage: group by sseqid, merge sstart/send, divide by slen.
    slen_by_phage = {}
    by_phage = defaultdict(list)
    with open(args.blast_hits) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            try:
                phage = p[1]
                pident = float(p[2])
                length = int(p[3])
                sstart = int(p[6])
                send = int(p[7])
                slen = int(p[9])
            except (ValueError, IndexError):
                continue
            slen_by_phage[phage] = slen
            if pident >= args.min_identity and length >= args.min_length:
                by_phage[phage].append((min(sstart, send), max(sstart, send)))

    cov = {}
    for phage, ivs in by_phage.items():
        covered = sum(e - s + 1 for s, e in merge_intervals(ivs))
        slen = slen_by_phage.get(phage, 0)
        if slen > 0:
            cov[phage] = min(100.0, covered / slen * 100)

    # 2. Join with INPHARED metadata.
    meta = pd.read_csv(args.inphared_metadata, sep="\t", low_memory=False)
    acc_col = find_col(meta, ("accession", "genome_id", "genome", "id"))
    life_col = find_col(meta, ("lifestyle", "bacphlip", "temperate", "virulent"))
    fam_col = find_col(meta, ("family",))
    class_col = find_col(meta, ("classification", "taxonomy"))
    host_col = find_col(meta, ("host",))
    if acc_col is None:
        emit("ERROR: cannot identify an accession column in the metadata. Columns:")
        emit(", ".join(meta.columns))
        return

    cov_by_norm = {strip_version(p): c for p, c in cov.items()}
    meta["_acc_norm"] = meta[acc_col].astype(str).apply(strip_version)
    meta["bacterial_coverage_pct"] = meta["_acc_norm"].map(cov_by_norm).fillna(0.0)

    n_total = len(meta)
    n_any_hit = sum(1 for c in cov_by_norm.values() if c > 0)
    n_tail = int((meta["bacterial_coverage_pct"] >= args.coverage_tail).sum())

    emit("=" * 70)
    emit("Phage bacterial-coverage characterization (bacteria-as-query BLAST)")
    emit("=" * 70)
    emit(f"INPHARED entries (metadata):           {n_total:,}")
    emit(f"Strong-hit (>= {args.min_identity}% / >= {args.min_length} bp) phages: "
         f"{n_any_hit:,} ({n_any_hit/n_total*100:.1f}% of INPHARED)")
    emit(f"Tail (>= {args.coverage_tail:.0f}% bacterial coverage):     "
         f"{n_tail:,} ({n_tail/n_total*100:.2f}% of INPHARED)")
    emit("")
    emit(f"Detected metadata columns: accession='{acc_col}', "
         f"lifestyle='{life_col}', family='{fam_col}', class='{class_col}', host='{host_col}'")
    emit("")

    bins = [-0.001, 0.001, 25, 50, 75, 100.01]
    labels = ["0% (no hit)", "0-25%", "25-50%", "50-75%", "75-100%"]
    meta["_cov_bin"] = pd.cut(meta["bacterial_coverage_pct"], bins=bins, labels=labels)
    emit("Phage count per coverage bin:")
    for lab, n in meta["_cov_bin"].value_counts(dropna=False).reindex(labels).fillna(0).astype(int).items():
        emit(f"  {lab:<14}: {n:,}")
    emit("")

    if life_col is not None:
        emit(f"{life_col} distribution per coverage bin:")
        ct = pd.crosstab(meta["_cov_bin"], meta[life_col], dropna=False)
        emit(ct.to_string())
        emit("")
        emit(f"{life_col} in the >= {args.coverage_tail:.0f}% tail "
             f"(expected: temperate-heavy if biology explains it):")
        tail = meta[meta["bacterial_coverage_pct"] >= args.coverage_tail]
        for v, n in tail[life_col].value_counts(dropna=False).items():
            emit(f"  {('(missing)' if pd.isna(v) else str(v))}: {n:,}")
        emit("")

    def top_in_tail(col, k=15):
        if col is None:
            return
        tail = meta[meta["bacterial_coverage_pct"] >= args.coverage_tail]
        if tail.empty:
            return
        emit(f"Top {k} '{col}' values in the tail:")
        for v, n in tail[col].value_counts(dropna=False).head(k).items():
            emit(f"  {('(missing)' if pd.isna(v) else str(v))}: {n:,}")
        emit("")

    top_in_tail(fam_col)
    top_in_tail(host_col)
    top_in_tail(class_col)

    if args.output_tail:
        keep_cols = [c for c in [acc_col, "bacterial_coverage_pct", life_col, fam_col, class_col, host_col] if c]
        tail = meta[meta["bacterial_coverage_pct"] >= args.coverage_tail][keep_cols] \
            .sort_values("bacterial_coverage_pct", ascending=False)
        tail.to_csv(args.output_tail, sep="\t", index=False)
        emit(f"Wrote tail TSV ({len(tail):,} rows) -> {args.output_tail}")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nReport written to {args.output}")


if __name__ == "__main__":
    main()

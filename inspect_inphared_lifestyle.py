#!/usr/bin/env python3
"""
Report the BACPHLIP lifestyle distribution (temperate / virulent / unknown)
across the INPHARED phage source database and, separately, across the 8,684
cluster representatives selected for the LAMBDA phage dataset.

Output is printed to stdout. No files written.

Usage:
    python inspect_inphared_lifestyle.py

If the metadata or vclust paths differ on your system, edit META_TSV and
VCLUST_TSV at the top of the file.
"""

import pandas as pd

META_TSV  = "/projects/bfzj/llindsey1/black_and_white/data/inphared/14Apr2025_data.tsv"
VCLUST_TSV = "/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv"


def find_first_column(df, keywords):
    """Return the first column in df whose name matches any of the keywords (case-insensitive)."""
    for c in df.columns:
        if any(kw in c.lower() for kw in keywords):
            return c
    return None


def main():
    meta = pd.read_csv(META_TSV, sep="\t", low_memory=False)
    print(f"INPHARED metadata: {len(meta)} rows, {len(meta.columns)} columns")
    print(f"Columns: {list(meta.columns)}")

    lifestyle_cols = [
        c for c in meta.columns
        if any(kw in c.lower() for kw in ("lifestyle", "bacphlip", "temperate", "virulent"))
    ]
    print(f"\nLifestyle-like columns: {lifestyle_cols}")

    for col in lifestyle_cols:
        print(f"\n=== {col} distribution (all INPHARED) ===")
        print(meta[col].value_counts(dropna=False).to_string())

    # vclust file: read first, inspect columns
    vclust = pd.read_csv(VCLUST_TSV, sep="\t", low_memory=False)
    print(f"\nvclust file: {len(vclust)} rows")
    print(f"Columns: {list(vclust.columns)}")

    # Heuristic: first column is genome ID, second is cluster ID
    genome_col = vclust.columns[0]
    cluster_col = vclust.columns[1]
    reps = vclust.groupby(cluster_col)[genome_col].first().tolist()
    print(f"\nNumber of clusters (= cluster representatives): {len(reps)}")
    print(f"First 3 cluster reps: {reps[:3]}")

    # Pick the accession column in the metadata
    acc_col = find_first_column(meta, ("accession", "genome_id", "genome", "id"))
    if acc_col is None:
        print("\nERROR: could not auto-detect accession column in metadata.")
        print("Inspect the column list above and edit this script to set acc_col manually.")
        return
    print(f"\nUsing metadata accession column: {acc_col}")

    # Normalize both sides (strip version suffix if present)
    rep_set = {str(r).split(".")[0] for r in reps}
    meta["acc_norm"] = meta[acc_col].astype(str).apply(lambda s: s.split(".")[0])
    rep_meta = meta[meta["acc_norm"].isin(rep_set)]
    print(f"Matched {len(rep_meta)} of {len(reps)} cluster reps in metadata")

    if not lifestyle_cols:
        print("\nNo lifestyle column detected. Inspect the column names listed above.")
        return

    for col in lifestyle_cols:
        print(f"\n=== {col} distribution among LAMBDA cluster representatives ({len(rep_meta)}) ===")
        counts = rep_meta[col].value_counts(dropna=False)
        print(counts.to_string())
        # Percent breakdown
        print("\nPercent:")
        print((counts / len(rep_meta) * 100).round(2).to_string())


if __name__ == "__main__":
    main()

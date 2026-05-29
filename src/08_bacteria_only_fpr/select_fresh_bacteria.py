#!/usr/bin/env python3
"""
Step 08a (part 1) -- select a FRESH bacterial genome pool for the FPR test set.

Same HQ filter as step 02 (CheckM2 completeness >= X%, contamination <= Y%),
but EXCLUDES every accession already in step 02's selection
(${GTDB_SELECTION_DIR}/all_accessions.txt), so the FPR set never overlaps the
training pool. From the remaining HQ representatives, take a random N
(seeded for reproducibility) -- this becomes the FPR source pool.
"""

import argparse
import random
import re

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Select a fresh bacterial pool for FPR (excludes step-02 selection)")
    p.add_argument("--metadata", required=True, help="bac120_metadata.tsv")
    p.add_argument("--exclude-accessions", required=True,
                   help="all_accessions.txt from step 02 -- these are excluded")
    p.add_argument("--output", required=True, help="Output accession list")
    p.add_argument("--n-genomes", type=int, required=True)
    p.add_argument("--min-completeness", type=float, default=95.0)
    p.add_argument("--max-contamination", type=float, default=5.0)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    # Same HQ filter as step 02: representatives + completeness + contamination.
    df = df[df["accession"] == df["gtdb_genome_representative"]].copy()
    df = df[
        (df["checkm2_completeness"] >= args.min_completeness)
        & (df["checkm2_contamination"] <= args.max_contamination)
    ]
    n_hq = len(df)

    with open(args.exclude_accessions) as f:
        excluded = {l.strip() for l in f if l.strip()}
    df = df[~df["accession"].isin(excluded)]
    n_remaining = len(df)

    n_target = min(args.n_genomes, n_remaining)
    rng = random.Random(args.seed)
    selected = rng.sample(df["accession"].tolist(), n_target)
    selected.sort()

    with open(args.output, "w") as f:
        for acc in selected:
            f.write(acc + "\n")

    print("=" * 60)
    print("FPR fresh-bacteria selection")
    print("=" * 60)
    print(f"HQ representatives in metadata:         {n_hq:,}")
    print(f"Excluded (already in step-02 selection): {len(excluded):,}")
    print(f"Remaining pool:                         {n_remaining:,}")
    print(f"Selected for FPR:                       {n_target:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

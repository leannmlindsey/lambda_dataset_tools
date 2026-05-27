#!/usr/bin/env python3
"""
Report how many genera survive different bacterial assembly-contiguity filters.

Decision support for step 02. We sample only the chromosome (longest contig) of
each bacterial genome, which is only well-justified for closed/isolate genomes.
CheckM2 *completeness* is not the same as assembly *contiguity*, so a genome can
pass the >=95% completeness filter and still be a fragmented MAG. This script
quantifies, among the high-quality GTDB representatives, how many genera would be
retained under progressively stricter isolate/contiguity requirements, so the
cutoff is chosen from real numbers rather than guessed.

It is read-only and writes nothing except an optional report file.

Usage:
    python analyze_bacteria_contiguity.py \
        --metadata bac120_metadata.tsv \
        [--min-completeness 95 --max-contamination 5] \
        [--output contiguity_retention_report.txt]
"""

import argparse
import re
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="GTDB bacterial contiguity retention analysis")
    p.add_argument("--metadata", "-m", required=True, help="Path to bac120_metadata.tsv")
    p.add_argument("--min-completeness", type=float, default=95.0)
    p.add_argument("--max-contamination", type=float, default=5.0)
    p.add_argument("--output", "-o", default=None, help="Optional report file")
    return p.parse_args()


def extract_genus(taxonomy):
    m = re.search(r"g__([^;]*)", str(taxonomy))
    return m.group(1) if m else "Unknown"


def is_mag_sag(category):
    """MAG/SAG/environmental (i.e. not an isolate) by ncbi_genome_category text."""
    c = str(category).lower()
    return any(k in c for k in ("metagenome", "single cell", "environmental"))


def main():
    args = parse_args()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)

    # Same representative + quality filter as select_gtdb_representatives.py
    df = df[df["accession"] == df["gtdb_genome_representative"]].copy()
    df = df[
        (df["checkm2_completeness"] >= args.min_completeness)
        & (df["checkm2_contamination"] <= args.max_contamination)
    ].copy()

    needed = ["ncbi_assembly_level", "ncbi_genome_category", "gtdb_taxonomy"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        emit(f"ERROR: metadata is missing expected column(s): {missing}")
        emit("Available columns:")
        emit(", ".join(df.columns))
        sys.exit(1)

    df["genus"] = df["gtdb_taxonomy"].apply(extract_genus)
    level = df["ncbi_assembly_level"].astype(str).str.strip()

    isolate = ~df["ncbi_genome_category"].apply(is_mag_sag)
    complete = level == "Complete Genome"
    complete_or_chrom = level.isin(["Complete Genome", "Chromosome"])

    baseline_genera = set(df["genus"].unique())
    n_base = len(baseline_genera)

    emit("=" * 70)
    emit("Bacterial assembly-contiguity retention analysis")
    emit("=" * 70)
    emit(f"Metadata: {args.metadata}")
    emit(f"Quality filter: completeness >= {args.min_completeness}%, "
         f"contamination <= {args.max_contamination}%")
    emit(f"High-quality representatives: {len(df):,}")
    emit(f"Baseline genera (current pipeline): {n_base:,}")
    emit("")

    emit("ncbi_assembly_level distribution (HQ reps):")
    for val, cnt in level.value_counts(dropna=False).items():
        emit(f"  {val}: {cnt:,}")
    emit("")
    emit("ncbi_genome_category distribution (HQ reps):")
    for val, cnt in df["ncbi_genome_category"].value_counts(dropna=False).items():
        emit(f"  {val if str(val).strip() else '(isolate/none)'}: {cnt:,}")
    emit("")

    policies = [
        ("baseline (all HQ reps)", pd.Series(True, index=df.index)),
        ("isolate only", isolate),
        ("Complete Genome (any source)", complete),
        ("isolate + Complete Genome", isolate & complete),
        ("isolate + Complete/Chromosome", isolate & complete_or_chrom),
    ]

    emit(f"{'policy':<34}{'genomes':>10}{'genera':>9}{'%genera':>9}{'lost':>8}")
    emit("-" * 70)
    lost_examples = {}
    for label, mask in policies:
        sub = df[mask]
        genera = set(sub["genus"].unique())
        pct = (len(genera) / n_base * 100) if n_base else 0.0
        lost = n_base - len(genera)
        emit(f"{label:<34}{len(sub):>10,}{len(genera):>9,}{pct:>8.1f}%{lost:>8,}")
        lost_examples[label] = sorted(baseline_genera - genera)

    # Detail the genera lost by the strictest practical policy.
    strict_label = "isolate + Complete Genome"
    lost = lost_examples[strict_label]
    emit("")
    emit(f"Genera with NO genome under '{strict_label}': {len(lost):,}")
    if lost:
        emit("  (these would be dropped under the strict policy; a tiered policy "
             "would keep them via a draft fallback)")
        emit("  examples: " + ", ".join(lost[:15]) + (" ..." if len(lost) > 15 else ""))

    # Genera retained if we require at least one genome whose longest contig is
    # long enough to sample from (longest_contig >= threshold). This is the
    # diversity cost of a minimum-contig-length floor in sampling (step 04);
    # expected to be ~0 given the median longest contig is ~360 kb. Plasmid
    # removal is handled by geNomad classification (step 03), NOT by length.
    if "longest_contig" in df.columns:
        lc = pd.to_numeric(df["longest_contig"], errors="coerce")
        emit("")
        emit("Genera retained by minimum longest-contig length (sampling-floor cost):")
        emit(f"{'min longest_contig':<22}{'genomes':>10}{'genera':>9}{'%genera':>9}{'lost':>8}")
        emit("-" * 61)
        for thresh in (2000, 4000, 8000, 10000, 16000, 50000, 100000):
            mask = lc >= thresh
            genera = df[mask]["genus"].nunique()
            pct = (genera / n_base * 100) if n_base else 0.0
            emit(f"{('>= ' + format(thresh, ',') + ' bp'):<22}"
                 f"{int(mask.sum()):>10,}{genera:>9,}{pct:>8.1f}%{n_base - genera:>8,}")

    # Contiguity descriptors, to sanity-check fragment sampling feasibility.
    for col in ("contig_count", "n50_contigs", "longest_contig"):
        if col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce")
            emit("")
            emit(f"{col}: median={s.median():.0f}  "
                 f"p90={s.quantile(0.9):.0f}  max={s.max():.0f}")

    if args.output:
        with open(args.output, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"\nReport written to {args.output}")


if __name__ == "__main__":
    main()

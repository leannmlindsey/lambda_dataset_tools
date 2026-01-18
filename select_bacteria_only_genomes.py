#!/usr/bin/env python3
"""
Select bacterial genomes for Bacteria-Only benchmark that do NOT overlap with fine-tuning datasets.

This script:
1. Reads the GTDB metadata to get all high-quality genera
2. Reads the existing fine-tuning accession lists (train/dev/test)
3. Reads the contaminated accessions list (prophage-containing)
4. Identifies genera NOT used in fine-tuning
5. Selects representatives from those genera for Bacteria-Only analysis
"""

import argparse
import random
import re
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select bacterial genomes for Bacteria-Only benchmark (no overlap with fine-tuning)"
    )
    parser.add_argument(
        "--metadata", "-m",
        required=True,
        help="Path to GTDB bac120_metadata.tsv"
    )
    parser.add_argument(
        "--finetuning-dir", "-f",
        required=True,
        help="Directory containing fine-tuning accession lists (train/dev/test_accessions.txt)"
    )
    parser.add_argument(
        "--contaminated",
        default=None,
        help="Path to contaminated_accessions.txt (prophage-containing genomes to exclude)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for Bacteria-Only accession lists"
    )
    parser.add_argument(
        "--n-genomes",
        type=int,
        default=None,
        help="Number of genomes to select (default: all available)"
    )
    parser.add_argument(
        "--min-completeness",
        type=float,
        default=95.0,
        help="Minimum CheckM2 completeness (default: 95.0)"
    )
    parser.add_argument(
        "--max-contamination",
        type=float,
        default=5.0,
        help="Maximum CheckM2 contamination (default: 5.0)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    return parser.parse_args()


def load_accessions(filepath):
    """Load accessions from a file."""
    accessions = set()
    with open(filepath, 'r') as f:
        for line in f:
            acc = line.strip()
            if acc:
                accessions.add(acc)
    return accessions


def extract_taxonomy_level(taxonomy_string, level):
    """Extract a specific taxonomy level from GTDB taxonomy string."""
    pattern = f'{level}__([^;]*)'
    match = re.search(pattern, taxonomy_string)
    return match.group(1) if match else "Unknown"


def main():
    args = parse_args()
    random.seed(args.seed)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    finetuning_dir = Path(args.finetuning_dir)

    print(f"Select Bacteria-Only Genomes (Non-Overlapping)")
    print(f"=" * 60)

    # Load GTDB metadata
    print(f"\nLoading GTDB metadata from {args.metadata}...")
    df = pd.read_csv(args.metadata, sep='\t', low_memory=False)
    print(f"  Total genomes in metadata: {len(df):,}")

    # Filter for representatives
    df_reps = df[df['accession'] == df['gtdb_genome_representative']].copy()
    print(f"  Representatives: {len(df_reps):,}")

    # Filter for quality
    df_hq = df_reps[
        (df_reps['checkm2_completeness'] >= args.min_completeness) &
        (df_reps['checkm2_contamination'] <= args.max_contamination)
    ].copy()
    print(f"  High-quality (>={args.min_completeness}% complete, <={args.max_contamination}% contam): {len(df_hq):,}")

    # Extract genus
    df_hq['genus'] = df_hq['gtdb_taxonomy'].apply(lambda x: extract_taxonomy_level(x, 'g'))
    unique_genera = df_hq['genus'].unique()
    print(f"  Unique genera: {len(unique_genera):,}")

    # Load fine-tuning accessions
    print(f"\nLoading fine-tuning accessions from {finetuning_dir}...")
    finetuning_accessions = set()
    for split in ['train', 'dev', 'test']:
        acc_file = finetuning_dir / f"{split}_accessions.txt"
        if acc_file.exists():
            accs = load_accessions(acc_file)
            finetuning_accessions.update(accs)
            print(f"  {split}: {len(accs):,} accessions")

    print(f"  Total fine-tuning accessions: {len(finetuning_accessions):,}")

    # Load contaminated accessions if provided
    contaminated_accessions = set()
    if args.contaminated and Path(args.contaminated).exists():
        contaminated_accessions = load_accessions(args.contaminated)
        print(f"  Contaminated accessions to exclude: {len(contaminated_accessions):,}")

    # Find genera used in fine-tuning
    finetuning_genera = set()
    for _, row in df_hq.iterrows():
        if row['accession'] in finetuning_accessions:
            finetuning_genera.add(row['genus'])

    print(f"  Genera used in fine-tuning: {len(finetuning_genera):,}")

    # Find available genera (not used in fine-tuning)
    available_genera = set(unique_genera) - finetuning_genera
    print(f"\n  Available genera for Bacteria-Only: {len(available_genera):,}")

    # Filter df_hq to only available genera and exclude contaminated
    df_available = df_hq[
        (df_hq['genus'].isin(available_genera)) &
        (~df_hq['accession'].isin(contaminated_accessions))
    ].copy()
    print(f"  Genomes in available genera (excluding contaminated): {len(df_available):,}")

    # Select one genome per available genus
    selected = []
    for genus in sorted(available_genera):
        genus_genomes = df_available[df_available['genus'] == genus]
        if len(genus_genomes) > 0:
            # Pick one randomly
            chosen = genus_genomes.sample(n=1, random_state=args.seed).iloc[0]
            selected.append({
                'accession': chosen['accession'],
                'genus': genus,
                'gtdb_taxonomy': chosen['gtdb_taxonomy'],
                'genome_size': chosen.get('genome_size', 'NA'),
                'checkm2_completeness': chosen['checkm2_completeness'],
                'checkm2_contamination': chosen['checkm2_contamination']
            })

    print(f"  Representatives selected: {len(selected):,}")

    # Subsample if requested
    if args.n_genomes and args.n_genomes < len(selected):
        random.shuffle(selected)
        selected = selected[:args.n_genomes]
        print(f"  After subsampling: {len(selected):,}")

    # Save outputs
    selected_df = pd.DataFrame(selected)

    # Accessions only
    acc_file = output_dir / "bacteria_only_accessions.txt"
    with open(acc_file, 'w') as f:
        for row in selected:
            f.write(f"{row['accession']}\n")
    print(f"\n  Saved: {acc_file}")

    # With metadata
    meta_file = output_dir / "bacteria_only_metadata.tsv"
    selected_df.to_csv(meta_file, sep='\t', index=False)
    print(f"  Saved: {meta_file}")

    # Summary
    summary_file = output_dir / "bacteria_only_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Bacteria-Only Dataset Selection\n")
        f.write(f"{'=' * 50}\n")
        f.write(f"Source: {args.metadata}\n")
        f.write(f"Fine-tuning dir: {args.finetuning_dir}\n")
        f.write(f"Contaminated file: {args.contaminated}\n")
        f.write(f"Quality filters: completeness >= {args.min_completeness}%, contamination <= {args.max_contamination}%\n")
        f.write(f"Seed: {args.seed}\n")
        f.write(f"\n")
        f.write(f"Total genomes in metadata: {len(df):,}\n")
        f.write(f"High-quality representatives: {len(df_hq):,}\n")
        f.write(f"Total unique genera: {len(unique_genera):,}\n")
        f.write(f"Fine-tuning accessions: {len(finetuning_accessions):,}\n")
        f.write(f"Fine-tuning genera: {len(finetuning_genera):,}\n")
        f.write(f"Available genera: {len(available_genera):,}\n")
        f.write(f"Contaminated excluded: {len(contaminated_accessions):,}\n")
        f.write(f"Selected for Bacteria-Only: {len(selected):,}\n")
        f.write(f"\n")
        f.write(f"NO OVERLAP with fine-tuning (verified by genus)\n")
    print(f"  Saved: {summary_file}")

    print(f"\nDone!")


if __name__ == "__main__":
    main()

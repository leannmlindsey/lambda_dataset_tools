#!/usr/bin/env python3
"""
Select phage genomes for Phage-Only benchmark that do NOT overlap with fine-tuning datasets.

This script:
1. Reads the vclust_info.tsv to get all clusters
2. Reads the existing fine-tuning accession lists (train/dev/test)
3. Identifies clusters NOT used in fine-tuning
4. Selects representatives from those clusters for Phage-Only analysis
"""

import argparse
import random
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select phage genomes for Phage-Only benchmark (no overlap with fine-tuning)"
    )
    parser.add_argument(
        "--vclust", "-v",
        required=True,
        help="Path to vclust_info.tsv"
    )
    parser.add_argument(
        "--finetuning-dir", "-f",
        required=True,
        help="Directory containing fine-tuning accession lists (train/dev/test_accessions.txt)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for Phage-Only accession lists"
    )
    parser.add_argument(
        "--n-genomes",
        type=int,
        default=None,
        help="Number of genomes to select (default: all available)"
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


def main():
    args = parse_args()
    random.seed(args.seed)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    finetuning_dir = Path(args.finetuning_dir)

    print(f"Select Phage-Only Genomes (Non-Overlapping)")
    print(f"=" * 60)

    # Load vclust file
    print(f"\nLoading vclust clusters from {args.vclust}...")
    cluster_to_genomes = defaultdict(list)
    all_genomes = set()

    with open(args.vclust, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                genome = parts[0]
                cluster = int(parts[1])
                cluster_to_genomes[cluster].append(genome)
                all_genomes.add(genome)

    print(f"  Total genomes in vclust: {len(all_genomes):,}")
    print(f"  Total clusters: {len(cluster_to_genomes):,}")

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

    # Find clusters used in fine-tuning
    finetuning_clusters = set()
    for cluster, genomes in cluster_to_genomes.items():
        for genome in genomes:
            if genome in finetuning_accessions:
                finetuning_clusters.add(cluster)
                break

    print(f"  Clusters used in fine-tuning: {len(finetuning_clusters):,}")

    # Find available clusters (not used in fine-tuning)
    available_clusters = set(cluster_to_genomes.keys()) - finetuning_clusters
    print(f"\n  Available clusters for Phage-Only: {len(available_clusters):,}")

    # Select representatives from available clusters
    available_representatives = []
    for cluster in sorted(available_clusters):
        # First genome in cluster is the representative
        rep = cluster_to_genomes[cluster][0]
        available_representatives.append((cluster, rep))

    print(f"  Available representatives: {len(available_representatives):,}")

    # Subsample if requested
    if args.n_genomes and args.n_genomes < len(available_representatives):
        random.shuffle(available_representatives)
        selected = available_representatives[:args.n_genomes]
        print(f"  Selected (subsampled): {len(selected):,}")
    else:
        selected = available_representatives
        print(f"  Selected (all available): {len(selected):,}")

    # Save outputs
    # All selected accessions
    all_file = output_dir / "phage_only_accessions.txt"
    with open(all_file, 'w') as f:
        for cluster, rep in selected:
            f.write(f"{rep}\n")
    print(f"\n  Saved: {all_file}")

    # With cluster info
    cluster_file = output_dir / "phage_only_with_clusters.tsv"
    with open(cluster_file, 'w') as f:
        f.write("accession\tcluster\n")
        for cluster, rep in selected:
            f.write(f"{rep}\t{cluster}\n")
    print(f"  Saved: {cluster_file}")

    # Summary
    summary_file = output_dir / "phage_only_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Phage-Only Dataset Selection\n")
        f.write(f"{'=' * 50}\n")
        f.write(f"Source: {args.vclust}\n")
        f.write(f"Fine-tuning dir: {args.finetuning_dir}\n")
        f.write(f"Seed: {args.seed}\n")
        f.write(f"\n")
        f.write(f"Total genomes in vclust: {len(all_genomes):,}\n")
        f.write(f"Total clusters: {len(cluster_to_genomes):,}\n")
        f.write(f"Fine-tuning accessions: {len(finetuning_accessions):,}\n")
        f.write(f"Fine-tuning clusters: {len(finetuning_clusters):,}\n")
        f.write(f"Available clusters: {len(available_clusters):,}\n")
        f.write(f"Selected for Phage-Only: {len(selected):,}\n")
        f.write(f"\n")
        f.write(f"NO OVERLAP with fine-tuning (verified by cluster)\n")
    print(f"  Saved: {summary_file}")

    print(f"\nDone!")


if __name__ == "__main__":
    main()

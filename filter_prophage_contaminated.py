#!/usr/bin/env python3
"""
Filter out prophage-contaminated bacterial accessions based on BLAST results.

Takes BLAST results (phages vs bacteria) and removes bacterial genomes
that contain prophage sequences from the accession lists.

Input: BLAST results in format 6 (standard tabular)
Output: Filtered accession lists for train/dev/test splits
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter prophage-contaminated bacterial accessions"
    )
    parser.add_argument(
        "--blast-results", "-b",
        required=True,
        help="Combined BLAST results file (format 6 std)"
    )
    parser.add_argument(
        "--input-accessions-dir", "-i",
        required=True,
        help="Directory containing original {train,dev,test}_accessions.txt"
    )
    parser.add_argument(
        "--output-accessions-dir", "-o",
        required=True,
        help="Directory for filtered accession lists"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=90.0,
        help="Minimum percent identity to consider a hit (default: 90)"
    )
    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=200,
        help="Minimum alignment length in bp (default: 200)"
    )
    parser.add_argument(
        "--output-contaminated",
        default=None,
        help="Output file for list of contaminated accessions (optional)"
    )
    parser.add_argument(
        "--output-summary",
        default=None,
        help="Output file for summary statistics (optional)"
    )
    return parser.parse_args()


def extract_accession_from_seqid(seqid):
    """
    Extract genome accession from BLAST subject sequence ID.

    Examples:
        GCA_000005845.2_ASM584v2_genomic -> GCA_000005845.2
        GCF_000005845.2_genomic -> GCF_000005845.2
        GCA_964231105.1_genomic -> GCA_964231105.1
    """
    # Pattern: GCA_XXXXXXXXX.Y or GCF_XXXXXXXXX.Y
    match = re.match(r'(GC[AF]_\d+\.\d+)', seqid)
    if match:
        return match.group(1)

    # Try without version number
    match = re.match(r'(GC[AF]_\d+)', seqid)
    if match:
        return match.group(1)

    # Return None if no match
    return None


def normalize_accession(accession):
    """
    Normalize accession for matching (strip GB_/RS_ prefix).
    """
    if accession.startswith('GB_'):
        return accession[3:]
    elif accession.startswith('RS_'):
        return accession[3:]
    return accession


def load_accessions(filepath):
    """Load accessions from file, preserving original format."""
    accessions = []
    with open(filepath, 'r') as f:
        for line in f:
            acc = line.strip()
            if acc:
                accessions.append(acc)
    return accessions


def main():
    args = parse_args()

    output_dir = Path(args.output_accessions_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Filter Prophage-Contaminated Bacterial Accessions")
    print("=" * 60)
    print(f"BLAST results: {args.blast_results}")
    print(f"Input accessions: {args.input_accessions_dir}")
    print(f"Output accessions: {args.output_accessions_dir}")
    print(f"Min identity: {args.min_identity}%")
    print(f"Min alignment length: {args.min_alignment_length} bp")
    print()

    # Step 1: Parse BLAST results and extract contaminated accessions
    print("Step 1: Parsing BLAST results...")

    contaminated_accessions = set()
    phages_with_hits = set()
    total_hits = 0
    filtered_hits = 0
    hit_stats = defaultdict(int)  # accession -> hit count

    with open(args.blast_results, 'r') as f:
        for line in f:
            if not line.strip():
                continue

            total_hits += 1
            fields = line.strip().split('\t')

            if len(fields) < 12:
                continue

            # Standard BLAST format 6 columns:
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            qseqid = fields[0]      # phage ID
            sseqid = fields[1]      # bacterial contig ID
            pident = float(fields[2])  # percent identity
            length = int(fields[3])    # alignment length

            # Apply filters (alignment length should already be filtered in v2, but double-check)
            if length < args.min_alignment_length:
                continue
            if pident < args.min_identity:
                continue

            filtered_hits += 1

            # Extract bacterial genome accession
            accession = extract_accession_from_seqid(sseqid)
            if accession:
                contaminated_accessions.add(accession)
                hit_stats[accession] += 1
                phages_with_hits.add(qseqid)

    print(f"  Total BLAST hits: {total_hits}")
    print(f"  Hits passing filters: {filtered_hits}")
    print(f"  Unique contaminated accessions: {len(contaminated_accessions)}")
    print(f"  Unique phages with hits: {len(phages_with_hits)}")
    print()

    # Step 2: Filter accession lists for each split
    print("Step 2: Filtering accession lists...")

    splits_stats = {}

    for split in ['train', 'dev', 'test', 'all']:
        input_file = Path(args.input_accessions_dir) / f"{split}_accessions.txt"
        output_file = output_dir / f"{split}_accessions.txt"

        if not input_file.exists():
            print(f"  {split}: skipped (file not found)")
            continue

        # Load original accessions
        original_accessions = load_accessions(input_file)

        # Filter out contaminated accessions
        # Need to handle both prefixed (GB_GCA_...) and unprefixed (GCA_...) formats
        filtered_accessions = []
        removed = 0

        for acc in original_accessions:
            normalized = normalize_accession(acc)
            if normalized in contaminated_accessions:
                removed += 1
            else:
                filtered_accessions.append(acc)

        # Write filtered accessions
        with open(output_file, 'w') as f:
            for acc in filtered_accessions:
                f.write(f"{acc}\n")

        splits_stats[split] = {
            'original': len(original_accessions),
            'filtered': len(filtered_accessions),
            'removed': removed
        }

        print(f"  {split}: {len(original_accessions)} -> {len(filtered_accessions)} (removed {removed})")

    print()

    # Step 3: Write contaminated accessions list
    if args.output_contaminated:
        print(f"Step 3: Writing contaminated accessions to {args.output_contaminated}")
        with open(args.output_contaminated, 'w') as f:
            for acc in sorted(contaminated_accessions):
                f.write(f"{acc}\t{hit_stats[acc]}\n")

    # Step 4: Write summary
    if args.output_summary:
        print(f"Step 4: Writing summary to {args.output_summary}")
        with open(args.output_summary, 'w') as f:
            f.write("Prophage Contamination Filtering Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"BLAST results: {args.blast_results}\n")
            f.write(f"Min identity: {args.min_identity}%\n")
            f.write(f"Min alignment length: {args.min_alignment_length} bp\n\n")
            f.write(f"Total BLAST hits: {total_hits}\n")
            f.write(f"Hits passing filters: {filtered_hits}\n")
            f.write(f"Unique contaminated accessions: {len(contaminated_accessions)}\n")
            f.write(f"Unique phages with hits: {len(phages_with_hits)}\n\n")
            f.write("Accession list filtering:\n")
            for split, stats in splits_stats.items():
                f.write(f"  {split}: {stats['original']} -> {stats['filtered']} (removed {stats['removed']})\n")
            f.write("\n")
            f.write("Top 20 most contaminated bacterial genomes:\n")
            top_contaminated = sorted(hit_stats.items(), key=lambda x: -x[1])[:20]
            for acc, count in top_contaminated:
                f.write(f"  {acc}: {count} hits\n")

    # Print top contaminated genomes
    print()
    print("Top 10 most contaminated bacterial genomes:")
    top_contaminated = sorted(hit_stats.items(), key=lambda x: -x[1])[:10]
    for acc, count in top_contaminated:
        print(f"  {acc}: {count} hits")

    print()
    print("=" * 60)
    print("Filtering complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()

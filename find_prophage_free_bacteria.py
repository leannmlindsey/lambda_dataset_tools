#!/usr/bin/env python3
"""
Find bacterial genomes with ZERO phage BLAST hits.

These are truly prophage-free bacteria suitable for a "bacteria-only" benchmark.

Usage:
    python find_prophage_free_bacteria.py \
        --blast-results /path/to/validation_blast_raw.tsv \
        --contig-map /path/to/contig_to_genome_map.tsv \
        --accessions /path/to/test_accessions.txt \
        --metadata /path/to/bac120_metadata.tsv \
        --output prophage_free_bacteria.txt \
        --output-metadata prophage_free_bacteria_metadata.tsv \
        --n-select 100
"""

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find bacterial genomes with zero phage BLAST hits"
    )
    parser.add_argument(
        "--blast-results", "-b",
        required=True,
        help="Raw BLAST results (validation_blast_raw.tsv)"
    )
    parser.add_argument(
        "--contig-map", "-c",
        required=True,
        help="Contig to genome mapping file"
    )
    parser.add_argument(
        "--accessions", "-a",
        required=True,
        help="Bacterial accessions to consider (e.g., test_accessions.txt)"
    )
    parser.add_argument(
        "--metadata", "-m",
        required=True,
        help="GTDB metadata file (bac120_metadata.tsv)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output file with prophage-free accessions"
    )
    parser.add_argument(
        "--output-metadata",
        help="Output file with metadata for selected genomes"
    )
    parser.add_argument(
        "--n-select", "-n",
        type=int,
        default=100,
        help="Number of genomes to select (default: 100)"
    )
    parser.add_argument(
        "--balance-taxonomy",
        action="store_true",
        help="Try to balance selection across phyla"
    )
    return parser.parse_args()


def load_contig_map(contig_map_file):
    """Load contig to genome mapping."""
    contig_to_genome = {}
    with open(contig_map_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_to_genome[parts[0]] = parts[1]
    return contig_to_genome


def load_accessions(accessions_file):
    """Load accession list."""
    accessions = set()
    with open(accessions_file, 'r') as f:
        for line in f:
            acc = line.strip()
            if acc:
                accessions.add(acc)
    return accessions


def find_genomes_with_hits(blast_file, contig_to_genome):
    """Find all genomes that have ANY BLAST hits."""
    genomes_with_hits = set()

    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Column 2 is the subject (bacterial contig)
                contig = parts[1]
                if contig in contig_to_genome:
                    genomes_with_hits.add(contig_to_genome[contig])

    return genomes_with_hits


def load_metadata(metadata_file, accessions_of_interest):
    """Load GTDB metadata for accessions of interest."""
    metadata = {}

    with open(metadata_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find important columns
        acc_idx = header.index('accession') if 'accession' in header else 0
        tax_idx = header.index('gtdb_taxonomy') if 'gtdb_taxonomy' in header else None
        comp_idx = header.index('checkm2_completeness') if 'checkm2_completeness' in header else None
        cont_idx = header.index('checkm2_contamination') if 'checkm2_contamination' in header else None
        size_idx = header.index('genome_size') if 'genome_size' in header else None

        for line in f:
            parts = line.strip().split('\t')
            acc = parts[acc_idx]

            # GTDB accessions have RS_ or GB_ prefix, our accessions might not
            acc_clean = acc.replace('RS_', '').replace('GB_', '')

            if acc in accessions_of_interest or acc_clean in accessions_of_interest:
                taxonomy = parts[tax_idx] if tax_idx else ""

                # Parse taxonomy to get phylum
                phylum = ""
                if taxonomy:
                    tax_parts = taxonomy.split(';')
                    for part in tax_parts:
                        if part.startswith('p__'):
                            phylum = part[3:]
                            break

                metadata[acc_clean] = {
                    'accession': acc_clean,
                    'full_accession': acc,
                    'taxonomy': taxonomy,
                    'phylum': phylum,
                    'completeness': float(parts[comp_idx]) if comp_idx and parts[comp_idx] else 0,
                    'contamination': float(parts[cont_idx]) if cont_idx and parts[cont_idx] else 0,
                    'genome_size': int(parts[size_idx]) if size_idx and parts[size_idx] else 0
                }

    return metadata


def select_balanced(prophage_free, metadata, n_select):
    """Select genomes balanced across phyla."""
    # Group by phylum
    by_phylum = defaultdict(list)
    for acc in prophage_free:
        if acc in metadata:
            phylum = metadata[acc]['phylum'] or 'Unknown'
            by_phylum[phylum].append(acc)

    print(f"\nPhyla distribution of prophage-free genomes:")
    for phylum, accs in sorted(by_phylum.items(), key=lambda x: -len(x[1])):
        print(f"  {phylum}: {len(accs)}")

    # Select proportionally from each phylum
    selected = []
    total_available = len(prophage_free)

    # First pass: allocate proportionally (minimum 1 per phylum if possible)
    allocations = {}
    for phylum, accs in by_phylum.items():
        proportion = len(accs) / total_available
        allocated = max(1, int(n_select * proportion))
        allocations[phylum] = min(allocated, len(accs))

    # Adjust to hit target
    total_allocated = sum(allocations.values())
    while total_allocated < n_select:
        # Add one more from the largest available phylum
        for phylum in sorted(by_phylum.keys(), key=lambda x: -len(by_phylum[x])):
            if allocations[phylum] < len(by_phylum[phylum]):
                allocations[phylum] += 1
                total_allocated += 1
                break
        else:
            break  # No more available

    while total_allocated > n_select:
        # Remove one from the smallest allocation (if > 1)
        for phylum in sorted(allocations.keys(), key=lambda x: allocations[x]):
            if allocations[phylum] > 1:
                allocations[phylum] -= 1
                total_allocated -= 1
                break
        else:
            break

    # Select from each phylum
    import random
    random.seed(42)

    for phylum, n in allocations.items():
        accs = by_phylum[phylum]
        # Sort by quality (completeness - contamination) and pick top n
        accs_sorted = sorted(accs,
                            key=lambda x: metadata.get(x, {}).get('completeness', 0) -
                                         metadata.get(x, {}).get('contamination', 0),
                            reverse=True)
        selected.extend(accs_sorted[:n])

    return selected[:n_select]


def main():
    args = parse_args()

    print("=" * 60)
    print("Find Prophage-Free Bacteria")
    print("=" * 60)
    print()

    # Load contig mapping
    print(f"Loading contig mapping from {args.contig_map}...")
    contig_to_genome = load_contig_map(args.contig_map)
    print(f"  Loaded {len(contig_to_genome):,} contig mappings")

    # Load accessions to consider
    print(f"\nLoading accessions from {args.accessions}...")
    all_accessions = load_accessions(args.accessions)
    print(f"  Loaded {len(all_accessions):,} accessions")

    # Find genomes with ANY BLAST hits
    print(f"\nScanning BLAST results from {args.blast_results}...")
    genomes_with_hits = find_genomes_with_hits(args.blast_results, contig_to_genome)
    print(f"  Found {len(genomes_with_hits):,} genomes with at least one phage hit")

    # Find prophage-free genomes (no hits at all)
    prophage_free = all_accessions - genomes_with_hits
    print(f"\nProphage-free genomes (ZERO hits): {len(prophage_free):,}")
    print(f"  ({len(prophage_free)/len(all_accessions)*100:.1f}% of test set)")

    if len(prophage_free) == 0:
        print("\nERROR: No prophage-free genomes found!")
        return

    if len(prophage_free) < args.n_select:
        print(f"\nWARNING: Only {len(prophage_free)} prophage-free genomes available, "
              f"requested {args.n_select}")
        args.n_select = len(prophage_free)

    # Load metadata for selection
    print(f"\nLoading metadata from {args.metadata}...")
    metadata = load_metadata(args.metadata, prophage_free)
    print(f"  Loaded metadata for {len(metadata):,} genomes")

    # Select genomes
    if args.balance_taxonomy:
        print(f"\nSelecting {args.n_select} balanced across phyla...")
        selected = select_balanced(prophage_free, metadata, args.n_select)
    else:
        # Just take first n (or random)
        import random
        random.seed(42)
        selected = list(prophage_free)
        random.shuffle(selected)
        selected = selected[:args.n_select]

    print(f"\nSelected {len(selected)} prophage-free genomes")

    # Write output accessions
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        for acc in sorted(selected):
            f.write(f"{acc}\n")
    print(f"\nWrote accessions to {args.output}")

    # Write metadata if requested
    if args.output_metadata:
        with open(args.output_metadata, 'w') as f:
            f.write("accession\tphylum\tcompleteness\tcontamination\tgenome_size\ttaxonomy\n")
            for acc in sorted(selected):
                if acc in metadata:
                    m = metadata[acc]
                    f.write(f"{acc}\t{m['phylum']}\t{m['completeness']:.2f}\t"
                           f"{m['contamination']:.2f}\t{m['genome_size']}\t{m['taxonomy']}\n")
                else:
                    f.write(f"{acc}\t\t\t\t\t\n")
        print(f"Wrote metadata to {args.output_metadata}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total test accessions:     {len(all_accessions):,}")
    print(f"Genomes with phage hits:   {len(genomes_with_hits):,} ({len(genomes_with_hits)/len(all_accessions)*100:.1f}%)")
    print(f"Prophage-free genomes:     {len(prophage_free):,} ({len(prophage_free)/len(all_accessions)*100:.1f}%)")
    print(f"Selected for benchmark:    {len(selected):,}")

    if args.balance_taxonomy and metadata:
        # Show phylum distribution of selected
        phylum_counts = defaultdict(int)
        for acc in selected:
            if acc in metadata:
                phylum_counts[metadata[acc]['phylum'] or 'Unknown'] += 1

        print(f"\nPhylum distribution of selected {len(selected)} genomes:")
        for phylum, count in sorted(phylum_counts.items(), key=lambda x: -x[1]):
            print(f"  {phylum}: {count}")


if __name__ == "__main__":
    main()

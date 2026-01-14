#!/usr/bin/env python3
"""
Shuffle nucleotides within each CDS of phage genomes.

This preserves:
- Gene order and boundaries
- GC content (same nucleotides, just reordered)
- Genome length

This destroys:
- Actual sequence (no similarity to original)
- Codon structure
- Any sequence-based features

Input: GenBank files with CDS annotations
Output: FASTA files with shuffled sequences
"""

import argparse
import random
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Shuffle nucleotides within each CDS"
    )
    parser.add_argument(
        "--input-dir", "-i",
        required=True,
        help="Directory containing GenBank files"
    )
    parser.add_argument(
        "--accession-list", "-a",
        required=True,
        help="File with list of accessions to process (one per line)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for shuffled FASTA files"
    )
    parser.add_argument(
        "--shuffle-intergenic",
        action="store_true",
        help="Also shuffle intergenic regions (default: only shuffle CDS)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--suffix",
        default="_genome.gbk",
        help="Suffix for GenBank files (default: _genome.gbk)"
    )
    return parser.parse_args()


def shuffle_sequence(seq):
    """Shuffle nucleotides in a sequence."""
    seq_list = list(str(seq))
    random.shuffle(seq_list)
    return ''.join(seq_list)


def get_cds_regions(record):
    """Extract CDS regions from a GenBank record."""
    cds_regions = []

    for feature in record.features:
        if feature.type == "CDS":
            # Get start and end (0-based)
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand  # 1 or -1
            cds_regions.append((start, end, strand))

    # Sort by start position
    cds_regions.sort(key=lambda x: x[0])

    return cds_regions


def shuffle_genome(record, shuffle_intergenic=False, seed=None):
    """
    Shuffle nucleotides within each CDS of a genome.

    Args:
        record: BioPython SeqRecord
        shuffle_intergenic: If True, also shuffle intergenic regions
        seed: Random seed for this genome

    Returns:
        New sequence with shuffled CDS regions
    """
    if seed is not None:
        random.seed(seed)

    seq = str(record.seq)
    seq_length = len(seq)

    # Get CDS regions
    cds_regions = get_cds_regions(record)

    if not cds_regions:
        # No CDS annotations - shuffle entire sequence or return as-is
        if shuffle_intergenic:
            return shuffle_sequence(seq)
        else:
            return seq

    # Build new sequence
    new_seq = list(seq)  # Mutable copy

    # Track which positions are in CDS
    in_cds = [False] * seq_length

    # Shuffle each CDS region
    for start, end, strand in cds_regions:
        # Mark positions as in CDS
        for i in range(start, end):
            in_cds[i] = True

        # Extract and shuffle CDS sequence
        cds_seq = seq[start:end]
        shuffled_cds = shuffle_sequence(cds_seq)

        # Replace in new sequence
        for i, nt in enumerate(shuffled_cds):
            new_seq[start + i] = nt

    # Optionally shuffle intergenic regions
    if shuffle_intergenic:
        # Find intergenic stretches
        i = 0
        while i < seq_length:
            if not in_cds[i]:
                # Find end of intergenic region
                j = i
                while j < seq_length and not in_cds[j]:
                    j += 1
                # Shuffle this intergenic region
                intergenic_seq = seq[i:j]
                shuffled_intergenic = shuffle_sequence(intergenic_seq)
                for k, nt in enumerate(shuffled_intergenic):
                    new_seq[i + k] = nt
                i = j
            else:
                i += 1

    return ''.join(new_seq)


def main():
    args = parse_args()
    random.seed(args.seed)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load accession list
    with open(args.accession_list, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]

    print(f"Processing {len(accessions)} genomes...")
    print(f"  Input dir: {input_dir}")
    print(f"  Output dir: {output_dir}")
    print(f"  Shuffle intergenic: {args.shuffle_intergenic}")
    print(f"  Seed: {args.seed}")
    print()

    processed = 0
    skipped = 0
    errors = []

    for i, accession in enumerate(accessions):
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i + 1}/{len(accessions)}")

        # Find GenBank file
        gbk_file = input_dir / f"{accession}{args.suffix}"

        if not gbk_file.exists():
            # Try without suffix
            gbk_file = input_dir / f"{accession}.gbk"

        if not gbk_file.exists():
            errors.append(f"{accession}: GenBank file not found")
            skipped += 1
            continue

        try:
            # Parse GenBank file
            record = SeqIO.read(gbk_file, "genbank")

            # Generate unique seed for this genome (reproducible)
            genome_seed = args.seed + hash(accession) % 10000

            # Shuffle the genome
            shuffled_seq = shuffle_genome(
                record,
                shuffle_intergenic=args.shuffle_intergenic,
                seed=genome_seed
            )

            # Create output record
            shuffled_record = SeqRecord(
                Seq(shuffled_seq),
                id=f"{accession}_shuffled",
                description=f"Shuffled CDS nucleotides from {accession}"
            )

            # Write FASTA
            output_file = output_dir / f"{accession}_shuffled.fasta"
            SeqIO.write(shuffled_record, output_file, "fasta")

            processed += 1

        except Exception as e:
            errors.append(f"{accession}: {str(e)}")
            skipped += 1

    print(f"\nDone!")
    print(f"  Processed: {processed}")
    print(f"  Skipped: {skipped}")

    if errors:
        error_file = output_dir / "errors.log"
        with open(error_file, 'w') as f:
            for error in errors:
                f.write(f"{error}\n")
        print(f"  Errors saved to: {error_file}")

    # Save summary
    summary_file = output_dir / "shuffle_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Shuffled CDS Nucleotides\n")
        f.write(f"{'='*40}\n")
        f.write(f"Input dir: {input_dir}\n")
        f.write(f"Accession list: {args.accession_list}\n")
        f.write(f"Total accessions: {len(accessions)}\n")
        f.write(f"Processed: {processed}\n")
        f.write(f"Skipped: {skipped}\n")
        f.write(f"Shuffle intergenic: {args.shuffle_intergenic}\n")
        f.write(f"Seed: {args.seed}\n")

    print(f"  Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()

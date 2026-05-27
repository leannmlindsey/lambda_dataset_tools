#!/usr/bin/env python3
"""
Remove fragment-level leakage between splits.

Given a class's per-split segment FASTAs (train/dev/test), this flags any
dev or test fragment that is near-identical to ANY train fragment and writes
"dereplicated" dev/test FASTAs with those fragments removed. The original
(full) FASTAs are left untouched so both versions are available.

Leakage is judged from a BLASTn search of the query (dev or test) fragments
against a database built from the train fragments. A query fragment is
considered leaked if it has any HSP with:
    percent identity        >= --min-identity   (default 95)
    query coverage per HSP   >= --min-coverage   (default 50)

Rationale: phages are mosaic and bacteria share conserved genes (rRNA,
housekeeping operons), so two fragments from genomes in *different* train/test
clusters can still be near-identical over a gene-sized window. Genome/cluster
splitting does not catch this; the model only ever sees fragments, so leakage
must be removed at the fragment level.

This script only does the filtering. The BLAST search itself is run by
dereplicate_fragments.slurm, which calls this script with the resulting
tabular output. Expected BLAST format (see the slurm wrapper):

    -outfmt "6 qseqid sseqid pident length qlen qcovhsp evalue bitscore"

Usage:
    python filter_leaked_fragments.py \
        --query-fasta dev_segments.fasta \
        --blast-output dev_vs_train.tsv \
        --output-fasta dev_segments.fasta            # in the *_derep dir \
        --output-leaked dev_leaked.tsv \
        --output-report dev_dereplication_report.txt \
        --min-identity 95 --min-coverage 50
"""

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove dev/test fragments that are near-identical to train fragments"
    )
    parser.add_argument(
        "--query-fasta", "-q",
        required=True,
        help="FASTA of the split being dereplicated (dev or test segments)"
    )
    parser.add_argument(
        "--blast-output", "-b",
        required=True,
        help='BLAST tabular output of query vs train fragments, '
             'outfmt "6 qseqid sseqid pident length qlen qcovhsp evalue bitscore"'
    )
    parser.add_argument(
        "--output-fasta", "-o",
        required=True,
        help="Output FASTA with leaked fragments removed (dereplicated split)"
    )
    parser.add_argument(
        "--output-leaked",
        required=True,
        help="Output TSV listing removed fragments and their best train hit"
    )
    parser.add_argument(
        "--output-report",
        required=True,
        help="Output text summary (counts + identity/coverage distribution)"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=95.0,
        help="Minimum percent identity to count as leakage (default: 95)"
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=50.0,
        help="Minimum query coverage per HSP (%%) to count as leakage (default: 50)"
    )
    return parser.parse_args()


def read_fasta_ids_in_order(fasta_file):
    """Return (ordered_ids, id_to_record_lines) so output preserves input order."""
    ids = []
    records = {}  # id -> list of raw lines (header + sequence lines) for fast re-emit
    current_id = None
    current_lines = []

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = current_lines
                current_id = line[1:].split()[0]
                ids.append(current_id)
                current_lines = [line]
            else:
                current_lines.append(line)
        if current_id is not None:
            records[current_id] = current_lines

    return ids, records


def find_leaked(blast_output, min_identity, min_coverage):
    """Scan BLAST tabular output; return {qseqid: (sseqid, pident, qcovhsp)} for best leaking hit.

    Expected columns:
      0 qseqid  1 sseqid  2 pident  3 length  4 qlen  5 qcovhsp  6 evalue  7 bitscore
    """
    leaked = {}  # qseqid -> (best_sseqid, best_pident, best_qcovhsp)

    with open(blast_output, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            try:
                pident = float(parts[2])
                qcovhsp = float(parts[5])
            except ValueError:
                continue

            if pident >= min_identity and qcovhsp >= min_coverage:
                prev = leaked.get(qseqid)
                # Keep the strongest leaking hit (by identity, then coverage) for the report.
                if prev is None or (pident, qcovhsp) > (prev[1], prev[2]):
                    leaked[qseqid] = (sseqid, pident, qcovhsp)

    return leaked


def main():
    args = parse_args()

    print("=" * 60)
    print("Fragment-level dereplication")
    print("=" * 60)
    print(f"Query FASTA:    {args.query_fasta}")
    print(f"BLAST output:   {args.blast_output}")
    print(f"Leakage cutoff: >= {args.min_identity}% identity AND "
          f">= {args.min_coverage}% query coverage")
    print()

    ids, records = read_fasta_ids_in_order(args.query_fasta)
    n_total = len(ids)
    print(f"Query fragments: {n_total:,}")

    leaked = find_leaked(args.blast_output, args.min_identity, args.min_coverage)
    n_leaked = sum(1 for qid in ids if qid in leaked)
    n_kept = n_total - n_leaked
    print(f"Leaked (removed): {n_leaked:,}")
    print(f"Kept:             {n_kept:,}")

    # Write dereplicated FASTA (preserve original order, drop leaked)
    out_fasta = Path(args.output_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as f:
        for qid in ids:
            if qid in leaked:
                continue
            f.writelines(records[qid])

    # Write leaked list
    out_leaked = Path(args.output_leaked)
    out_leaked.parent.mkdir(parents=True, exist_ok=True)
    with open(out_leaked, "w") as f:
        f.write("query_fragment\tbest_train_hit\tpercent_identity\tquery_coverage\n")
        for qid in ids:
            if qid in leaked:
                sseqid, pident, qcov = leaked[qid]
                f.write(f"{qid}\t{sseqid}\t{pident:.2f}\t{qcov:.2f}\n")

    # Identity distribution of the removed fragments (for the report)
    id_bins = {"95-96": 0, "96-97": 0, "97-98": 0, "98-99": 0, "99-100": 0}
    for qid in ids:
        if qid in leaked:
            p = leaked[qid][1]
            if p < 96:
                id_bins["95-96"] += 1
            elif p < 97:
                id_bins["96-97"] += 1
            elif p < 98:
                id_bins["97-98"] += 1
            elif p < 99:
                id_bins["98-99"] += 1
            else:
                id_bins["99-100"] += 1

    out_report = Path(args.output_report)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    pct = (n_leaked / n_total * 100) if n_total else 0.0
    with open(out_report, "w") as f:
        f.write("Fragment-level Dereplication Report\n")
        f.write("=" * 40 + "\n")
        f.write(f"Query FASTA:        {args.query_fasta}\n")
        f.write(f"Min identity:       {args.min_identity}%\n")
        f.write(f"Min query coverage: {args.min_coverage}%\n")
        f.write("\n")
        f.write(f"Total fragments:    {n_total}\n")
        f.write(f"Leaked (removed):   {n_leaked} ({pct:.2f}%)\n")
        f.write(f"Kept:               {n_kept}\n")
        f.write("\n")
        f.write("Identity distribution of removed fragments:\n")
        for label, count in id_bins.items():
            f.write(f"  {label}%: {count}\n")

    print(f"\nWrote dereplicated FASTA -> {out_fasta}")
    print(f"Wrote leaked list        -> {out_leaked}")
    print(f"Wrote report             -> {out_report}")
    print("\nDone!")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Read the (length x split) filtered BLAST hit tables from step 05a + the matching
fragment metadata TSVs from step 04, and emit:

  - hit_ids_${seg}_${split}.txt        : segment IDs that hit (one per line)
  - replacement_counts_${seg}_${split}.tsv : accession <tab> count_to_replace
  - affected_accessions.txt            : union of accessions needing full-genome BLAST

Step 05b reads `affected_accessions.txt` to gather + re-BLAST those genomes.
Step 05c reads the per-(length x split) hit IDs and replacement counts.

Filtered hit table is BLAST outfmt 6: qseqid (= fragment segment_id) is col 0.
Fragment metadata TSV columns: segment_id, accession, contig, start, end, length.
"""

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Identify hit fragments + their source genomes")
    p.add_argument("--hits-dir", required=True, help="Dir containing fragment_hits_${seg}_${split}.tsv")
    p.add_argument("--fragment-dirs", required=True,
                   help="Comma-separated list of BACTERIA_DATASETS dirs (one per length, in order 2k,4k,8k)")
    p.add_argument("--seg-lengths", default="2000,4000,8000")
    p.add_argument("--splits", default="train,dev,test")
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_hit_segment_ids(hits_tsv):
    """qseqid (col 0) of each row -> set of hit fragment IDs."""
    ids = set()
    with open(hits_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                ids.add(parts[0])
    return ids


def load_metadata(meta_tsv):
    """segment_id -> accession (skip header)."""
    seg2acc = {}
    with open(meta_tsv) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                seg2acc[parts[0]] = parts[1]
    return seg2acc


def main():
    args = parse_args()
    seg_lengths = [int(x) for x in args.seg_lengths.split(",")]
    splits = args.splits.split(",")
    bact_dirs = args.fragment_dirs.split(",")
    assert len(bact_dirs) == len(seg_lengths), "fragment-dirs must match seg-lengths order"
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    affected_union = set()
    summary = []

    for seg, bact_dir in zip(seg_lengths, bact_dirs):
        for split in splits:
            hits_tsv = Path(args.hits_dir) / f"fragment_hits_{seg}_{split}.tsv"
            meta_tsv = Path(bact_dir) / f"{split}_segments.tsv"
            hit_ids = load_hit_segment_ids(hits_tsv) if hits_tsv.exists() else set()
            seg2acc = load_metadata(meta_tsv)

            hit_ids_in_meta = hit_ids & set(seg2acc)
            if len(hit_ids_in_meta) != len(hit_ids):
                # Some hits had no metadata row (shouldn't happen) -- log it
                print(f"  WARNING: {seg}/{split}: {len(hit_ids) - len(hit_ids_in_meta)} hit IDs "
                      "not found in metadata")

            # per-accession replacement counts (one fragment to replace per hit)
            per_acc = defaultdict(int)
            for sid in hit_ids_in_meta:
                per_acc[seg2acc[sid]] += 1

            # Write outputs for this (seg, split)
            with open(out_dir / f"hit_ids_{seg}_{split}.txt", "w") as fh:
                for sid in sorted(hit_ids_in_meta):
                    fh.write(sid + "\n")
            with open(out_dir / f"replacement_counts_{seg}_{split}.tsv", "w") as fh:
                fh.write("accession\tcount_to_replace\n")
                for acc, n in sorted(per_acc.items()):
                    fh.write(f"{acc}\t{n}\n")

            affected_union.update(per_acc)
            summary.append((seg, split, len(hit_ids_in_meta), len(per_acc)))

    # Union of affected accessions across all (seg x split)
    with open(out_dir / "affected_accessions.txt", "w") as fh:
        for acc in sorted(affected_union):
            fh.write(acc + "\n")

    print("=" * 60)
    print("Find affected fragments + genomes")
    print("=" * 60)
    print(f"{'seg':>6}{'split':>8}{'hit_fragments':>16}{'affected_genomes':>20}")
    for seg, split, n_hits, n_aff in summary:
        print(f"{seg:>6}{split:>8}{n_hits:>16,}{n_aff:>20,}")
    print(f"\nUnion of affected genomes: {len(affected_union):,} -> "
          f"{out_dir / 'affected_accessions.txt'}")


if __name__ == "__main__":
    main()

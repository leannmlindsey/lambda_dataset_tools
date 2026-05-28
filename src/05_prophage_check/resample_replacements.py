#!/usr/bin/env python3
"""
Step 05c -- final assembly.

For each (segment length, split):
  1. Read the original step-04 bacterial fragments + metadata.
  2. Drop fragments whose ID is in the hit list (they hit INPHARED at the
     prophage-check threshold).
  3. For each affected genome (one that contributed >=1 hit fragment in this
     length x split), sample EXACTLY `count_to_replace` replacement fragments
     from that genome, avoiding all masked regions in step 05b's mask file.
  4. Write the merged set (kept + replacements) to a NEW
     `${SPLIT}_segments_checked.fasta` + `.tsv` next to the originals. Step 04's
     outputs are not touched, so the chain is auditable.

Reproducible: per-replacement RNG seeded from a stable hash of
`SEED:accession:seg:split:replace`.
"""

import argparse
import gzip
import hashlib
import random
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description="Drop hit fragments + sample replacements from clean regions")
    p.add_argument("--gtdb-dir", required=True)
    p.add_argument("--fragment-dirs", required=True, help="Comma-separated BACTERIA_DATASETS")
    p.add_argument("--seg-lengths", default="2000,4000,8000")
    p.add_argument("--splits", default="train,dev,test")
    p.add_argument("--hits-input-dir", required=True,
                   help="Dir containing hit_ids_${seg}_${split}.txt + replacement_counts_${seg}_${split}.tsv")
    p.add_argument("--masks", required=True,
                   help="Per-contig masks (contig_id, start, end, 1-based inclusive) from step 05b")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---- GTDB accession -> file path (matches the rest of the pipeline) ----
def normalize_acc(acc):
    if acc.startswith("GB_"): return acc[3:]
    if acc.startswith("RS_"): return acc[3:]
    return acc


def accession_to_path(gtdb_dir, accession):
    accession = normalize_acc(accession)
    parts = accession.split("_")
    if len(parts) != 2:
        return None, None
    prefix, num_v = parts
    num = num_v.split(".")[0].zfill(9)
    fna = Path(gtdb_dir) / "database" / prefix / num[0:3] / num[3:6] / num[6:9] / f"{accession}_genomic.fna.gz"
    return fna, accession


def stable_seed(global_seed, *parts):
    key = ":".join([str(global_seed)] + [str(p) for p in parts])
    return int(hashlib.md5(key.encode()).hexdigest()[:8], 16)


# ---- allowed intervals = contig minus masked, each >= seg_len ----
def allowed_intervals(contig_len, masked_1based, seg_len):
    ms = sorted((s - 1, e) for s, e in masked_1based)
    merged = []
    for s, e in ms:
        s = max(0, s); e = min(contig_len, e)
        if e <= s: continue
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    allowed = []
    prev = 0
    for s, e in merged:
        if s > prev:
            allowed.append((prev, s))
        prev = max(prev, e)
    if prev < contig_len:
        allowed.append((prev, contig_len))
    return [(a, b) for a, b in allowed if b - a >= seg_len]


def place_n(slots, n, seg_len, rng, max_attempts=2000):
    weights = [(b - a - seg_len + 1) for _, a, b in slots]
    total_w = sum(w for w in weights if w > 0)
    if total_w <= 0:
        return []
    placed = []
    placed_by_contig = defaultdict(list)
    for _ in range(n):
        for _attempt in range(max_attempts):
            r = rng.randrange(total_w)
            cum = 0; idx = 0
            for i, w in enumerate(weights):
                if w <= 0: continue
                cum += w
                if r < cum:
                    idx = i; break
            cid, a, b = slots[idx]
            s = rng.randint(a, b - seg_len)
            e = s + seg_len
            if all(e <= ps or s >= pe for ps, pe in placed_by_contig[cid]):
                placed.append((cid, s, e))
                placed_by_contig[cid].append((s, e))
                break
        else:
            break
    return placed


def load_masks(path):
    """contig_id -> list of (start, end) 1-based inclusive."""
    m = defaultdict(list)
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            try:
                m[parts[0]].append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    return m


def read_genome(gtdb_dir, accession):
    """Return {contig_id: seq_str}, handling GCA<->GCF fallback."""
    fna, norm = accession_to_path(gtdb_dir, accession)
    if fna is None or not fna.exists():
        alt = norm.replace("GCA_", "GCF_") if norm.startswith("GCA_") else norm.replace("GCF_", "GCA_")
        fna, norm = accession_to_path(gtdb_dir, alt)
        if fna is None or not fna.exists():
            return None
    with gzip.open(fna, "rt") as fh:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}


def load_metadata(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                rows.append({"segment_id": parts[0], "accession": parts[1],
                             "contig": parts[2], "start": int(parts[3]),
                             "end": int(parts[4]), "length": int(parts[5])})
    return rows


def load_fasta(path):
    """{segment_id: seq_str} (first whitespace-token id)."""
    seqs = {}
    cur_id = None; buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if cur_id is not None:
            seqs[cur_id] = "".join(buf)
    return seqs


def write_fasta(path, segments):
    with open(path, "w") as fh:
        for seg_id, seq in segments:
            fh.write(f">{seg_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_metadata(path, records):
    with open(path, "w") as fh:
        fh.write("segment_id\taccession\tcontig\tstart\tend\tlength\n")
        for r in records:
            fh.write(f"{r['segment_id']}\t{r['accession']}\t{r['contig']}\t"
                     f"{r['start']}\t{r['end']}\t{r['length']}\n")


def main():
    args = parse_args()
    seg_lengths = [int(x) for x in args.seg_lengths.split(",")]
    splits = args.splits.split(",")
    bact_dirs = args.fragment_dirs.split(",")
    assert len(bact_dirs) == len(seg_lengths)

    masks = load_masks(args.masks)
    print(f"Loaded masks for {len(masks):,} contigs.")

    hits_dir = Path(args.hits_input_dir)

    for seg_len, bact_dir in zip(seg_lengths, bact_dirs):
        min_contig = max(10000, 2 * seg_len)
        for split in splits:
            print(f"\n--- {seg_len} bp / {split} ---")
            fasta = Path(bact_dir) / f"{split}_segments.fasta"
            meta = Path(bact_dir) / f"{split}_segments.tsv"
            hit_ids_file = hits_dir / f"hit_ids_{seg_len}_{split}.txt"
            counts_file = hits_dir / f"replacement_counts_{seg_len}_{split}.tsv"

            hit_ids = set()
            if hit_ids_file.exists():
                with open(hit_ids_file) as f:
                    hit_ids = {l.strip() for l in f if l.strip()}
            print(f"  hit fragments to drop: {len(hit_ids):,}")

            seqs = load_fasta(fasta)
            meta_rows = load_metadata(meta)

            # Kept = those not in hit_ids
            kept_segments = [(r["segment_id"], seqs[r["segment_id"]]) for r in meta_rows
                             if r["segment_id"] not in hit_ids and r["segment_id"] in seqs]
            kept_records = [r for r in meta_rows if r["segment_id"] not in hit_ids]
            print(f"  kept (no hit): {len(kept_records):,}")

            # Replacements: per accession, sample count_to_replace from clean regions
            counts = {}
            if counts_file.exists():
                with open(counts_file) as f:
                    f.readline()
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) == 2:
                            counts[parts[0]] = int(parts[1])
            print(f"  affected genomes needing replacements: {len(counts):,} "
                  f"(total replacements: {sum(counts.values()):,})")

            replaced = 0
            for acc, n_needed in sorted(counts.items()):
                contigs = read_genome(args.gtdb_dir, acc)
                if not contigs:
                    print(f"    WARN: genome not found for {acc} (needed {n_needed})")
                    continue
                slots = []
                for cid, seq in contigs.items():
                    if len(seq) < min_contig:
                        continue
                    for a, b in allowed_intervals(len(seq), masks.get(cid, []), seg_len):
                        slots.append((cid, a, b))
                if not slots:
                    print(f"    WARN: no clean slots in {acc} (needed {n_needed})")
                    continue
                rng = random.Random(stable_seed(args.seed, acc, seg_len, split, "replace"))
                placed = place_n(slots, n_needed, seg_len, rng)
                for cid, s, e in placed:
                    seg_id = f"{acc}_{cid}_{s}_{e}_repl"
                    seq = contigs[cid][s:e]
                    kept_segments.append((seg_id, seq))
                    kept_records.append({"segment_id": seg_id, "accession": acc,
                                         "contig": cid, "start": s, "end": e,
                                         "length": e - s})
                    replaced += 1
            print(f"  replacements actually placed: {replaced:,}")

            # Write the "_checked" outputs alongside step 04's originals
            out_fasta = Path(bact_dir) / f"{split}_segments_checked.fasta"
            out_meta  = Path(bact_dir) / f"{split}_segments_checked.tsv"
            write_fasta(out_fasta, kept_segments)
            write_metadata(out_meta, kept_records)
            print(f"  -> {out_fasta}  ({len(kept_segments):,} fragments total)")


if __name__ == "__main__":
    main()

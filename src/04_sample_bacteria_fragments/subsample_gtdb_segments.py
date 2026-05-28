#!/usr/bin/env python3
"""
Subsample bacterial segments from GTDB genomes — mask-aware and multi-contig.

Rewrite addressing the pipeline review + the masking design:
- Samples from ALL contigs >= --min-contig-length (not just the longest), so
  plasmids and secondary contigs are included. Placement is proportional to the
  available (un-masked) length.
- Excludes phage-similar regions: reads the per-contig masked intervals produced
  by step 03 and samples only from the UN-masked remainder. It never rejoins
  flanks across a masked region (no chimeric junctions) — a segment must fit
  entirely within one un-masked stretch.
- Reproducible: each genome's RNG is seeded from a STABLE hash (hashlib) of the
  accession, not Python's salted built-in hash(). The global allocation uses its
  own seeded RNG over a sorted accession list, so results are independent of
  thread completion order.
- Hits --total-segments by allocating across genomes that actually yield
  candidate segments, then redistributing any shortfall to genomes with spare
  capacity (fixes the old undershoot from dividing by ALL accessions).

Masked intervals TSV (from step 03): contig_id <tab> start <tab> end
  (1-based, inclusive). A contig with no row is fully un-masked.
"""

import argparse
import gzip
import hashlib
import math
import random
from collections import defaultdict
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description="Mask-aware multi-contig bacterial segment subsampling")
    p.add_argument("--gtdb-dir", "-g", required=True)
    p.add_argument("--accession-list", "-a", required=True)
    p.add_argument("--masked-intervals", "-m", default=None,
                   help="Per-contig masked intervals TSV from step 03 (optional; if omitted, no masking)")
    p.add_argument("--output", "-o", required=True)
    p.add_argument("--output-metadata", default=None)
    p.add_argument("--segment-length", type=int, default=2000)
    p.add_argument("--min-contig-length", type=int, default=10000,
                   help="Only sample from contigs at least this long (tiny-contig hygiene)")
    p.add_argument("--total-segments", type=int, required=True)
    p.add_argument("--threads", "-t", type=int, default=8)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def normalize_accession(acc):
    if acc.startswith("GB_"):
        return acc[3:]
    if acc.startswith("RS_"):
        return acc[3:]
    return acc


def accession_to_path(gtdb_dir, accession):
    accession = normalize_accession(accession)
    parts = accession.split("_")
    if len(parts) != 2:
        return None, None
    prefix, number_version = parts
    number = number_version.split(".")[0].zfill(9)
    fna = Path(gtdb_dir) / "database" / prefix / number[0:3] / number[3:6] / number[6:9] / f"{accession}_genomic.fna.gz"
    return fna, accession


def stable_seed(global_seed, accession):
    """Deterministic per-accession seed (independent of PYTHONHASHSEED)."""
    h = hashlib.md5(f"{global_seed}:{accession}".encode()).hexdigest()
    return int(h[:8], 16)


def load_masked_intervals(path):
    """contig_id -> list of (start, end) 1-based inclusive."""
    masked = defaultdict(list)
    if not path:
        return masked
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                masked[parts[0]].append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    return masked


def allowed_intervals(contig_len, masked, min_len):
    """Un-masked sub-intervals of a contig, 0-based half-open, each >= min_len.

    masked: list of (start, end) 1-based inclusive.
    """
    # to 0-based half-open, merge
    ms = sorted((s - 1, e) for s, e in masked)
    merged = []
    for s, e in ms:
        s = max(0, s)
        e = min(contig_len, e)
        if e <= s:
            continue
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
    return [(a, b) for a, b in allowed if b - a >= min_len]


def place_segments(slots, n, seg_len, rng, max_attempts=1000):
    """Place up to n non-overlapping segments across allowed slots.

    slots: list of (contig_id, astart, aend) 0-based half-open, each >= seg_len.
    Returns list of (contig_id, start, end) in deterministic insertion order.
    """
    weights = [aend - astart - seg_len + 1 for _, astart, aend in slots]
    total_w = sum(w for w in weights if w > 0)
    if total_w <= 0:
        return []
    placed = []
    placed_by_contig = defaultdict(list)
    for _ in range(n):
        for _attempt in range(max_attempts):
            r = rng.randrange(total_w)
            cum = 0
            idx = 0
            for i, w in enumerate(weights):
                if w <= 0:
                    continue
                cum += w
                if r < cum:
                    idx = i
                    break
            cid, astart, aend = slots[idx]
            start = rng.randint(astart, aend - seg_len)
            end = start + seg_len
            if all(end <= ps or start >= pe for ps, pe in placed_by_contig[cid]):
                placed.append((cid, start, end))
                placed_by_contig[cid].append((start, end))
                break
        else:
            break
    return placed


def process_genome(task):
    accession, gtdb_dir, masked, seg_len, min_contig_len, cap, global_seed = task
    fna, norm = accession_to_path(gtdb_dir, accession)
    if fna is None:
        return (accession, [], "invalid")
    if not fna.exists():
        # try the other prefix (GCA <-> GCF)
        alt = norm.replace("GCA_", "GCF_") if norm.startswith("GCA_") else norm.replace("GCF_", "GCA_")
        fna, norm = accession_to_path(gtdb_dir, alt)
        if fna is None or not fna.exists():
            return (accession, [], "not_found")
    try:
        with gzip.open(fna, "rt") as fh:
            contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
    except Exception as e:
        return (accession, [], f"error:{e}")

    slots = []
    contig_seq = {}
    for cid, seq in contigs.items():
        if len(seq) < min_contig_len:
            continue
        for astart, aend in allowed_intervals(len(seq), masked.get(cid, []), seg_len):
            slots.append((cid, astart, aend))
        contig_seq[cid] = seq
    if not slots:
        return (accession, [], "no_eligible")

    rng = random.Random(stable_seed(global_seed, accession))
    placed = place_segments(slots, cap, seg_len, rng)
    candidates = []
    for i, (cid, start, end) in enumerate(placed):
        seg_id = f"{accession}_{cid}_{start}_{end}"
        candidates.append((seg_id, contig_seq[cid][start:end], cid, start, end))
    return (accession, candidates, "success")


def main():
    args = parse_args()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(args.accession_list) as f:
        accessions = [l.strip() for l in f if l.strip()]

    masked = load_masked_intervals(args.masked_intervals)
    # Candidate cap per genome: enough that the pool comfortably exceeds the target.
    cap = max(2, math.ceil(args.total_segments / max(1, len(accessions))) + 3)

    print("GTDB mask-aware subsampling")
    print("=" * 50)
    print(f"  accessions:        {len(accessions)}")
    print(f"  segment length:    {args.segment_length}")
    print(f"  min contig length: {args.min_contig_length}")
    print(f"  masked contigs:    {len(masked)}")
    print(f"  target segments:   {args.total_segments}")
    print(f"  per-genome cap:    {cap}")
    print(f"  seed:              {args.seed}")

    tasks = [(acc, args.gtdb_dir, masked, args.segment_length,
              args.min_contig_length, cap, args.seed) for acc in accessions]

    candidates_by_acc = {}
    stats = defaultdict(int)
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_genome, t): t[0] for t in tasks}
        done = 0
        for fut in as_completed(futures):
            acc, cands, status = fut.result()
            stats[status] += 1
            if cands:
                candidates_by_acc[acc] = cands
            done += 1
            if done % 1000 == 0:
                print(f"  processed {done}/{len(accessions)}")

    # Deterministic allocation to hit the target exactly.
    valid = sorted(candidates_by_acc)
    alloc_rng = random.Random(args.seed)
    alloc_rng.shuffle(valid)
    n_valid = len(valid)
    if n_valid == 0:
        print("ERROR: no genomes produced candidate segments.")
        raise SystemExit(1)
    base, extra = divmod(args.total_segments, n_valid)
    chosen = {}
    shortfall = 0
    for i, acc in enumerate(valid):
        want = base + (1 if i < extra else 0)
        avail = len(candidates_by_acc[acc])
        chosen[acc] = min(want, avail)
        shortfall += want - chosen[acc]
    # redistribute shortfall to genomes with spare candidates
    while shortfall > 0:
        progressed = False
        for acc in valid:
            if shortfall == 0:
                break
            if chosen[acc] < len(candidates_by_acc[acc]):
                chosen[acc] += 1
                shortfall -= 1
                progressed = True
        if not progressed:
            break

    # Emit (sorted accession order for stable output).
    segments = []
    for acc in sorted(chosen):
        segments.extend(candidates_by_acc[acc][:chosen[acc]])

    with open(out, "w") as fh:
        for seg_id, seq, *_ in segments:
            fh.write(f">{seg_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    if args.output_metadata:
        with open(args.output_metadata, "w") as fh:
            fh.write("segment_id\taccession\tcontig\tstart\tend\tlength\n")
            for seg_id, seq, cid, start, end in segments:
                acc = seg_id.rsplit("_", 3)[0] if "_" in seg_id else seg_id
                fh.write(f"{seg_id}\t{acc}\t{cid}\t{start}\t{end}\t{end - start}\n")

    print("\nSummary:")
    print(f"  success genomes:   {stats['success']}")
    print(f"  not found:         {stats['not_found']}")
    print(f"  no eligible seq:   {stats['no_eligible']}")
    print(f"  genomes used:      {sum(1 for a in chosen if chosen[a] > 0)}")
    print(f"  total segments:    {len(segments)} (target {args.total_segments})")
    if shortfall > 0:
        print(f"  WARNING: short by {shortfall} (insufficient eligible sequence)")


if __name__ == "__main__":
    main()

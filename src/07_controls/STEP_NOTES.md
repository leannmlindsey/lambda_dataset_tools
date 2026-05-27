# Step 07 — Composition controls

**Status:** ⏳ to port

**Purpose:** Controls that hold sequence composition fixed while destroying
sequence-specific signal, to test whether a model relies on composition.

**Implements (current, top-level):** `shuffle_segments.py`,
`shuffle_test_segments.slurm`, `create_gc_control_csvs.sh`

**Pending fixes (review #6):**
- The current full within-segment shuffle preserves only **mononucleotide**
  (GC) composition and destroys all k-mer/codon structure, so it does not
  isolate GC. Add a **dinucleotide-preserving** shuffle (Altschul–Erikson /
  `uShuffle`) to separate "uses GC" from "uses k-mer/codon motifs".
- Rename outputs so "GC control" is not overclaimed (e.g. `mono_shuffle` vs
  `dinuc_shuffle`).
- Run controls on both the full and dereplicated test sets.

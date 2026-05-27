# Step 09 — Package for release

**Status:** ⏳ to port

**Purpose:** Build the tar.gz archives + metadata bundle for Zenodo.

**Implements (current, top-level):** `create_lambda_final.sh`

**Pending fixes / additions:**
- Include **both** the full and dereplicated ("hard") merged datasets, and the
  step-05 `*_dereplication_report.txt` + `*_leaked.tsv` audit files.
- Include the pinned `environment.yml` and a versions manifest (vclust,
  pharokka, bakta, checkm2, blast, genomad).
- Write `manifests/` md5 checksums of all released files so a rerun can be
  verified bit-identical.
- Bundle the corrected `METHODS_DATA_SELECTION.md` and the resolved source
  counts (GTDB release, INPHARED rep count) once verified.

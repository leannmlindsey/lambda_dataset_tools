#!/bin/bash
#SBATCH --account=bfzj-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --gpus-per-node=1
#SBATCH --job-name=09_package
#SBATCH --output=09_package_%j.out
#SBATCH --error=09_package_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g

# Step 09: assemble the Zenodo / HuggingFace-ready release.
#  - copies the train/val/test, shuffled-control, and FPR CSVs into RELEASE_DIR
#  - copies metadata (accession lists, pipeline git commit)
#  - writes a top-level README.md (Zenodo+HF dataset card)
#  - writes checksums.md5 for every data file
# Notes:
#  * Internally we used "dev"; the release renames to "val" to match HF convention.
#  * Genome-wide evaluation datasets stay on biowulf -- not built by this pipeline.

set -euo pipefail
source "${LAMBDA_CONFIG:-/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/src/config.sh}"
source ~/.bashrc
conda activate "${CONDA_ENV}"

echo "Step 09: package release -> ${RELEASE_DIR}  (start $(date))"
rm -rf "${RELEASE_DIR}"
mkdir -p "${RELEASE_DIR}"/{train_val_test,shuffled_controls,fpr_test,metadata/phage_accessions,metadata/bacteria_accessions}

# ---- train/val/test (dev -> val on copy) ----
for LEN in 2k 4k 8k; do
    OUT="${RELEASE_DIR}/train_val_test/${LEN}"
    mkdir -p "${OUT}"
    cp "${MERGED_DIR}/${LEN}/train.csv" "${OUT}/train.csv"
    cp "${MERGED_DIR}/${LEN}/dev.csv"   "${OUT}/val.csv"
    cp "${MERGED_DIR}/${LEN}/test.csv"  "${OUT}/test.csv"
done

# ---- shuffled controls ----
for LEN in 2k 4k 8k; do
    OUT="${RELEASE_DIR}/shuffled_controls/${LEN}"
    mkdir -p "${OUT}"
    cp "${SHUFFLED_DIR}/${LEN}/test_shuffled.csv" "${OUT}/test_shuffled.csv"
done

# ---- FPR ----
for LEN in 2k 4k 8k; do
    OUT="${RELEASE_DIR}/fpr_test/${LEN}"
    mkdir -p "${OUT}"
    cp "${FPR_DATA_DIR}/${LEN}/bacteria_segments_${LEN}.csv" "${OUT}/bacteria_segments_${LEN}.csv"
done

# ---- metadata: accession lists ----
PHAGE_ACC="${PHAGE_DATASETS[0]}"             # step 01 wrote train/dev/test_accessions.txt here
cp "${PHAGE_ACC}/train_accessions.txt" "${RELEASE_DIR}/metadata/phage_accessions/train.txt"
cp "${PHAGE_ACC}/dev_accessions.txt"   "${RELEASE_DIR}/metadata/phage_accessions/val.txt"
cp "${PHAGE_ACC}/test_accessions.txt"  "${RELEASE_DIR}/metadata/phage_accessions/test.txt"
cp "${GTDB_SELECTION_DIR}/train_accessions.txt" "${RELEASE_DIR}/metadata/bacteria_accessions/train.txt"
cp "${GTDB_SELECTION_DIR}/dev_accessions.txt"   "${RELEASE_DIR}/metadata/bacteria_accessions/val.txt"
cp "${GTDB_SELECTION_DIR}/test_accessions.txt"  "${RELEASE_DIR}/metadata/bacteria_accessions/test.txt"
cp "${FPR_SELECTION_DIR}/fpr_accessions.txt"    "${RELEASE_DIR}/metadata/fpr_bacteria_accessions.txt"

# ---- pipeline_version.txt = git commit ----
( cd "${PIPELINE_ROOT}" && git rev-parse HEAD 2>/dev/null || echo "(no git)" ) > \
    "${RELEASE_DIR}/metadata/pipeline_version.txt"

# ---- top-level README (Zenodo + HF dataset card) ----
cat > "${RELEASE_DIR}/README.md" <<EOF
---
license: cc-by-4.0
language:
  - en
tags:
  - genomics
  - bacteriophage
  - prophage-detection
  - dna-language-model
  - benchmark
size_categories:
  - 100K<n<1M
---

# LAMBDA — LAnguage Model Bacteriophage Detection Assessment (${RELEASE_VERSION})

Phage vs bacteria sequence-classification benchmark for genomic language models,
balanced 1:1 at three fragment lengths (2 kb, 4 kb, 8 kb).

## Directory layout

\`\`\`
LAMBDA_${RELEASE_VERSION}/
├── train_val_test/{2k,4k,8k}/{train,val,test}.csv   # main fine-tuning datasets
├── shuffled_controls/{2k,4k,8k}/test_shuffled.csv   # composition-preserved control
├── fpr_test/{2k,4k,8k}/bacteria_segments_*.csv      # bacteria-only PHROG-sized FPR
├── metadata/
│   ├── phage_accessions/{train,val,test}.txt
│   ├── bacteria_accessions/{train,val,test}.txt
│   ├── fpr_bacteria_accessions.txt
│   └── pipeline_version.txt
├── checksums.md5
└── README.md
\`\`\`

## CSV schemas

**train_val_test** and **shuffled_controls** — \`segment_id, sequence, label, source\`
- \`label\` : 1 = phage, 0 = bacteria
- \`source\`: \`inphared\` (phage) or \`gtdb\` (bacteria)

**fpr_test** — \`seq_id, start, end, sequence, label\` (matches the existing
PHROG phage-FNR CSV format)
- \`seq_id\` : NCBI-style contig accession (e.g. \`NC_001234.1\`)
- \`start\`/\`end\` : 1-based inclusive
- \`label\` : 0 for every row (bacteria-only)

## How the data was built (one line per step)

1. **Phage selection** — INPHARED vclust 95% ANI reps, one random rep / cluster, 80/10/10 split by cluster.
2. **Bacteria selection** — GTDB HQ reps (CheckM2 ≥95% complete, ≤5% contam), one random genome / genus, 80/10/10 split by genus (15,865 genera).
3. **Phage fragments** — sample 2k/4k/8k at ~1 per 10 kb of genome length.
4. **Bacterial fragments** — sample 2k/4k/8k from all contigs ≥10 kb, matching the per-split phage count for a 1:1 balance.
5. **Prophage check** — BLAST every bacterial fragment vs INPHARED at ≥70% identity, ≥200 bp; for each hit drop the fragment, BLAST the source genome, and resample a replacement from the same genome's non-phage regions (no genome dropped).
6. **Merge** — combine + shuffle into the 1:1 CSVs (label 1 phage / 0 bacteria).
7. **Shuffled control** — per-row nucleotide shuffle of each test CSV (composition preserved, sequence destroyed).
8. **FPR set** — fresh HQ bacteria pool disjoint from selection, same prophage check, output in PHROG CSV format.

Reproducibility: every step is seeded (single global \`SEED=42\` driving stable
per-accession hashes). Pipeline commit recorded in \`metadata/pipeline_version.txt\`.

## Citation

(Add citation when published.)
EOF

# ---- checksums ----
echo "Computing MD5 checksums..."
( cd "${RELEASE_DIR}" && find . -type f ! -name 'checksums.md5' -print0 \
    | sort -z | xargs -0 md5sum > checksums.md5 )

# ---- summary ----
echo ""
echo "=========================================="
echo "Release ready at: ${RELEASE_DIR}"
echo "=========================================="
echo ""
echo "Top-level structure:"
( cd "${RELEASE_DIR}" && find . -maxdepth 2 -mindepth 1 -type d | sort )
echo ""
echo "Total size: $(du -sh "${RELEASE_DIR}" | cut -f1)"
echo "Total files: $(find "${RELEASE_DIR}" -type f | wc -l)"
echo "Done: $(date)"

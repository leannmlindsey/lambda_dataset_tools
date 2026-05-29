# Step 09 — Package release

**Role:** assemble a Zenodo / HuggingFace-ready release directory at
`${RELEASE_DIR}` (default `${DATA_DIR}/LAMBDA_v1`).

**Renamed for the release:** internal `dev` → published `val` (HF convention).

**Layout:**
```
LAMBDA_v1/
├── README.md                              # dataset card (Zenodo + HF compatible)
├── train_val_test/{2k,4k,8k}/{train,val,test}.csv
├── shuffled_controls/{2k,4k,8k}/test_shuffled.csv
├── fpr_test/{2k,4k,8k}/bacteria_segments_*.csv
├── metadata/
│   ├── phage_accessions/{train,val,test}.txt
│   ├── bacteria_accessions/{train,val,test}.txt
│   ├── fpr_bacteria_accessions.txt
│   └── pipeline_version.txt              # git commit of this pipeline
└── checksums.md5                         # md5 of every data file
```

Genome-wide eval datasets stay on biowulf — they're not built by this pipeline
and aren't part of this release.

**Configurable:** `RELEASE_VERSION` (default `v1`) and `RELEASE_DIR`.

# LAMBDA Dataset Construction Tools

Tools for constructing the LAMBDA (LAnguage Model Bacteriophage Detection Assessment) benchmark datasets for evaluating genomic language models on prophage detection tasks.

## Overview

This repository contains scripts for:
- Downloading and processing INPHARED phage genome database
- Selecting and filtering GTDB bacterial representative genomes
- **Removing prophage-contaminated bacterial genomes** (BLAST-based filtering)
- Creating balanced phage/bacterial datasets with cluster-aware train/dev/test splits
- Subsampling genome segments for language model training (2k, 4k, 8k lengths)
- Generating GC-content control sequences (shuffled nucleotides)
- Creating taxonomically balanced annotation datasets (Phage-Only, Bacteria-Only)

> **See [`PATHS.md`](PATHS.md) for a comprehensive reference of all data paths and file locations.**

## Quick Start

```bash
# 1. Create phage segments (already done)
sbatch subsample_inphared.slurm

# 2. Create prophage-filtered bacterial segments (recommended)
sbatch blast_phage_vs_gtdb_selected_v3.slurm  # BLAST phages vs bacteria
sbatch filter_and_subsample_gtdb_v3.slurm     # Filter and subsample

# 3. Validate prophage removal
sbatch validate_prophage_removal_v3.slurm

# 4. Create GC-content controls
sbatch shuffle_test_segments.slurm

# 5. Create annotation datasets for benchmarking
sbatch select_bacteria_only_balanced.slurm    # Select 100 balanced bacteria
sbatch run_bakta_bacteria.slurm               # Annotate with Bakta
```

## Installation

### Prerequisites
- Python 3.8+
- Conda or Mamba package manager
- BLAST+ (for prophage filtering)
- Bakta (for bacterial annotation)

### Create Conda Environment

```bash
conda create -n lambda_tools python=3.10 -y
conda activate lambda_tools

conda install -c conda-forge -c bioconda \
    pandas \
    biopython \
    matplotlib \
    numpy \
    blast \
    -y
```

---

## Scripts Reference

### Data Download & Selection

| Script | Description |
|--------|-------------|
| `download_inphared.py` | Download INPHARED phage genome database |
| `extract_representatives.py` | Extract INPHARED cluster representatives, split into train/dev/test |
| `select_gtdb_representatives.py` | Select GTDB bacterial genomes with quality/taxonomic filtering |

### Segment Subsampling

| Script | Description |
|--------|-------------|
| `subsample_segments.py` | Subsample segments from INPHARED genomes |
| `subsample_gtdb_segments.py` | Subsample segments from GTDB genomes (handles nested dirs) |
| `subsample_inphared.slurm` | SLURM job for INPHARED 2k segments |
| `subsample_inphared_4k.slurm` | SLURM job for INPHARED 4k segments |
| `subsample_inphared_8k.slurm` | SLURM job for INPHARED 8k segments |
| `subsample_gtdb.slurm` | SLURM job for GTDB 2k segments (original, with prophage) |
| `subsample_gtdb_4k.slurm` | SLURM job for GTDB 4k segments (original) |
| `subsample_gtdb_8k.slurm` | SLURM job for GTDB 8k segments (original) |

### Prophage Filtering (v3 - Recommended)

| Script | Description |
|--------|-------------|
| `blast_phage_vs_gtdb_selected_v3.slurm` | BLAST INPHARED vs GTDB (max_target_seqs=2000) |
| `filter_and_subsample_gtdb_v3.slurm` | Filter contaminated genomes, subsample 2k/4k/8k |
| `validate_prophage_removal_v3.slurm` | Validate no prophage sequences remain |
| `create_contig_to_genome_map.py` | Map BLAST contig IDs to genome accessions |
| `create_contig_map.slurm` | SLURM job for contig mapping |

### BLAST Database Creation

| Script | Description |
|--------|-------------|
| `create_phage_blastdb.slurm` | Create BLAST database from INPHARED phages |
| `create_gtdb_blastdb.slurm` | Create BLAST database from all GTDB bacteria |
| `create_gtdb_selected_blastdb.slurm` | Create BLAST database from selected GTDB (~16k) |

### GC-Content Control

| Script | Description |
|--------|-------------|
| `shuffle_segments.py` | Shuffle nucleotides within segments (preserves GC) |
| `shuffle_test_segments.slurm` | Create shuffled test segments for GC control |

### Benchmark-Specific Datasets

| Script | Description |
|--------|-------------|
| `select_bacteria_only_balanced.py` | Select 100 taxonomically balanced bacteria |
| `select_bacteria_only_balanced.slurm` | SLURM job for balanced selection |
| `test_bakta_single.slurm` | Test Bakta on single genome before full run |
| `run_bakta_bacteria.slurm` | Run Bakta annotation on 100 bacteria |

### Dataset Merging

| Script | Description |
|--------|-------------|
| `merge_and_shuffle.py` | Merge phage/bacteria segments with labels |
| `merge_datasets.slurm` | SLURM job for merging all segment lengths |

### Analysis & Evaluation

| Script | Description |
|--------|-------------|
| `analyze_blast_filtering.sh` | Analyze BLAST filtering steps in detail |
| `investigate_missed_hits.sh` | Investigate remaining prophage hits |
| `calculate_prophage_metrics.py` | Calculate TP/TN/FP/FN, Precision, Recall, F1 |
| `compare_accessions_v2.sh` | Compare original vs contaminated accessions |

### Clustering (Optional)

| Script | Description |
|--------|-------------|
| `run_mash_distances.slurm` | Compute MASH all-vs-all distances |
| `cluster_mash_distances.py` | Cluster genomes based on MASH distances |

---

## Workflow

### Step 1: Create Phage Dataset

```bash
# Extract cluster representatives from vclust output
python extract_representatives.py \
    --clusters /path/to/vclust_info.tsv \
    --output ./inphared_dataset \
    --split-ratio 80:10:10 \
    --seed 42

# Subsample segments at 2k, 4k, 8k
sbatch subsample_inphared.slurm
sbatch subsample_inphared_4k.slurm
sbatch subsample_inphared_8k.slurm
```

### Step 2: Create Prophage-Filtered Bacterial Dataset (Recommended)

```bash
# 1. BLAST all phages against selected bacteria (max_target_seqs=2000)
sbatch blast_phage_vs_gtdb_selected_v3.slurm

# 2. Filter contaminated genomes and subsample segments
sbatch filter_and_subsample_gtdb_v3.slurm

# 3. Validate prophage removal
sbatch validate_prophage_removal_v3.slurm
```

**Filtering Parameters:**
- Min alignment length: ≥200 bp
- Min identity: ≥90%
- max_target_seqs: 2000

**Results (v3):**
- Original genomes: 15,865
- Contaminated (removed): 971 (6.1%)
- Clean (retained): 14,894 (93.9%)

### Step 3: Create GC-Content Control Dataset

```bash
sbatch shuffle_test_segments.slurm
```

This shuffles nucleotides within each test segment, preserving GC content but destroying sequence patterns.

### Step 4: Create Annotation Datasets for Benchmarking

**Phage-Only:** Uses test set genomes (869) with Pharokka annotations (already available).

**Bacteria-Only:**
```bash
# Select 100 taxonomically balanced bacteria
sbatch select_bacteria_only_balanced.slurm

# Test Bakta on one genome first
conda activate bakta_env
sbatch test_bakta_single.slurm

# Run Bakta on all 100
sbatch run_bakta_bacteria.slurm
```

---

## Dataset Structure

```
# Phage segments (INPHARED)
inphared_dataset/                    # 2k segments
├── train_accessions.txt
├── dev_accessions.txt
├── test_accessions.txt
├── train_segments.fasta
├── dev_segments.fasta
├── test_segments.fasta
└── train_segments.tsv

inphared_dataset_4k/                 # 4k segments
inphared_dataset_8k/                 # 8k segments

# Bacteria segments - FILTERED (prophage-free) ✓ RECOMMENDED
gtdb_dataset_filtered_v3/            # 2k segments
├── train_accessions.txt
├── dev_accessions.txt
├── test_accessions.txt
├── contaminated_accessions.txt      # Removed genomes
├── train_segments.fasta
├── dev_segments.fasta
└── test_segments.fasta

gtdb_dataset_filtered_v3_4k/         # 4k segments
gtdb_dataset_filtered_v3_8k/         # 8k segments

# GC-content control (shuffled)
shuffled_gc_control/
├── inphared_2k_test_shuffled.fasta
├── gtdb_2k_test_shuffled.fasta
└── ...

# Bacteria-only annotation dataset
bacteria_only_dataset/
├── bacteria_only_100_accessions.txt
└── bacteria_only_100_metadata.tsv

bacteria_only_annotations/           # Bakta output
└── {accession}/
    ├── {accession}.gff3
    ├── {accession}.faa
    └── ...
```

---

## Key Parameters

### Segment Subsampling
- `--segment-length`: Length of each segment (2000, 4000, or 8000 bp)
- `--min-length`: Minimum genome length to include (default: 5000 bp)
- `--sample-per-bp`: Take 1 sample per N bp (default: 10000)
- `--total-segments`: Target total segments (for balancing datasets)

### Quality Filtering (GTDB)
- `--min-completeness`: Minimum CheckM2 completeness (default: 95%)
- `--max-contamination`: Maximum CheckM2 contamination (default: 5%)

### Prophage Filtering
- Min alignment length: ≥200 bp
- Min identity: ≥90%
- BLAST max_target_seqs: 2000

### Data Splitting
- `--split-ratio`: Train:dev:test ratio (default: 80:10:10)
- `--seed`: Random seed for reproducibility (default: 42)

---

## Data Leakage Prevention

- **Phage**: Split by vclust cluster (95% ANI) - no cluster spans multiple splits
- **Bacteria**: Split by GTDB genus - no genus spans multiple splits
- **Prophage removal**: BLAST-based filtering ensures bacterial segments don't contain phage sequences

---

## Dataset Statistics

| Dataset | Train | Dev | Test | Genomes Used |
|---------|------:|----:|-----:|-------------:|
| **Phage (INPHARED)** | 27,661 | 3,703 | 3,427 | 7,044 |
| **Bacteria (filtered v3)** | 27,000 | 3,400 | 3,400 | 14,894 |

---

## File Paths

See [`PATHS.md`](PATHS.md) for a comprehensive reference of all data paths, including:
- Source data locations (INPHARED, GTDB)
- BLAST database paths
- Generated dataset directories
- Intermediate/working files

---

## Citation

If you use these tools or the LAMBDA benchmark, please cite:

```
[Citation to be added]
```

## License

[License to be added]

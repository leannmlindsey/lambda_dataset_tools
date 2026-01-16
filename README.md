# LAMBDA Dataset Construction Tools

Tools for constructing the LAMBDA ( LAnguage Model Bacteriophage Detection Assessment) benchmark datasets for prophage detection in bacterial hosts using genomic language models.

## Overview

This repository contains scripts for:
- Downloading and processing INPHARED phage genome database
- Selecting and filtering GTDB bacterial representative genomes
- Creating balanced phage/bacterial datasets with cluster-aware train/dev/test splits
- Subsampling genome segments for language model training (2k, 4k, 8k lengths)
- Merging and shuffling phage/bacteria datasets for classification
- Generating shuffled control sequences
- Creating BLAST databases for prophage detection/removal
- Calculating prophage prediction metrics

## Installation

### Prerequisites
- Python 3.8+
- Conda or Mamba package manager

### Create Conda Environment

```bash
# Create new environment
conda create -n lambda_tools python=3.10 -y

# Activate environment
conda activate lambda_tools

# Install required packages
conda install -c conda-forge -c bioconda \
    pandas \
    biopython \
    matplotlib \
    numpy \
    blast \
    -y

# Or install via pip (except BLAST which requires conda/bioconda)
pip install pandas biopython matplotlib numpy
```

### Verify Installation

```bash
python -c "import pandas; import Bio; import matplotlib; import numpy; print('All packages installed successfully!')"
```

## Scripts

### Data Download

| Script | Description |
|--------|-------------|
| `download_inphared.py` | Download INPHARED phage genome database |
| `download_inphared.slurm` | SLURM job script for downloading INPHARED |

### Dataset Selection

| Script | Description |
|--------|-------------|
| `extract_representatives.py` | Extract INPHARED cluster representatives and split into train/dev/test |
| `select_gtdb_representatives.py` | Select GTDB bacterial genomes with quality filtering and taxonomic subsampling |

### Segment Processing

| Script | Description |
|--------|-------------|
| `subsample_segments.py` | Subsample segments from INPHARED genomes (multi-FASTA input) |
| `subsample_gtdb_segments.py` | Subsample segments from GTDB genomes (handles nested dirs and gzipped files) |
| `subsample_inphared.slurm` | SLURM job for INPHARED 2k segment subsampling (train/dev/test) |
| `subsample_gtdb.slurm` | SLURM job for GTDB 2k segment subsampling (train/dev/test) |
| `subsample_inphared_4k.slurm` | SLURM job for INPHARED 4k segment subsampling |
| `subsample_inphared_8k.slurm` | SLURM job for INPHARED 8k segment subsampling |
| `subsample_gtdb_4k.slurm` | SLURM job for GTDB 4k segment subsampling |
| `subsample_gtdb_8k.slurm` | SLURM job for GTDB 8k segment subsampling |

### Dataset Merging

| Script | Description |
|--------|-------------|
| `merge_and_shuffle.py` | Merge phage/bacteria segments and shuffle with labels |
| `merge_datasets.slurm` | SLURM job for merging all segment lengths (2k, 4k, 8k) |

### Control Sequences

| Script | Description |
|--------|-------------|
| `shuffle_cds_nucleotides.py` | Shuffle nucleotides within CDS regions to create control sequences |
| `shuffle_cds.slurm` | SLURM job for creating shuffled control sequences |

### BLAST Database Creation

| Script | Description |
|--------|-------------|
| `create_phage_blastdb.slurm` | Create BLAST database from all INPHARED phage genomes |
| `create_gtdb_blastdb.slurm` | Create BLAST database from all GTDB bacterial genomes (~136k) |
| `create_gtdb_selected_blastdb.slurm` | Create BLAST database from selected GTDB genomes (~16k) |

### Evaluation

| Script | Description |
|--------|-------------|
| `calculate_prophage_metrics.py` | Calculate nucleotide-level TP/TN/FP/FN, Precision, Recall, F1, MCC |
| `plot_processing_time.py` | Generate bar plots of processing times |

### Clustering (Optional)

| Script | Description |
|--------|-------------|
| `run_mash_distances.slurm` | SLURM job for computing MASH all-vs-all distances |
| `cluster_mash_distances.py` | Cluster genomes based on MASH distances |

## Usage

### 1. Download INPHARED Database

```bash
python download_inphared.py \
    --output ./inphared_data \
    --category all \
    --decompress
```

### 2. Create Phage Dataset (INPHARED)

```bash
# Extract cluster representatives and split
python extract_representatives.py \
    --clusters vclust_info.tsv \
    --output ./inphared_dataset \
    --split-ratio 80:10:10 \
    --seed 42
```

### 3. Create Bacterial Dataset (GTDB)

```bash
# Select high-quality representatives, one per genus
python select_gtdb_representatives.py \
    --metadata bac120_metadata.tsv \
    --output ./gtdb_dataset \
    --min-completeness 95 \
    --max-contamination 5 \
    --seed 42
```

### 4. Subsample Segments

Use the SLURM scripts to subsample segments for all splits (train/dev/test):

```bash
# INPHARED phage segments (2k, 4k, or 8k)
sbatch subsample_inphared.slurm      # 2k segments
sbatch subsample_inphared_4k.slurm   # 4k segments
sbatch subsample_inphared_8k.slurm   # 8k segments

# GTDB bacterial segments (matched to phage counts)
# First, update TOTAL_SEGMENTS_* in the scripts to match INPHARED output
sbatch subsample_gtdb.slurm          # 2k segments
sbatch subsample_gtdb_4k.slurm       # 4k segments
sbatch subsample_gtdb_8k.slurm       # 8k segments
```

Or run manually:

```bash
# Phage segments (1 per 10kb)
python subsample_segments.py \
    --input-fasta 14Apr2025_genomes.fa \
    --accession-list ./inphared_dataset/train_accessions.txt \
    --output ./inphared_dataset/train_segments.fasta \
    --output-metadata ./inphared_dataset/train_segments.tsv \
    --segment-length 2000 \
    --min-length 5000 \
    --sample-per-bp 10000 \
    --seed 42

# Bacterial segments (matched to phage count using --total-segments)
python subsample_gtdb_segments.py \
    --gtdb-dir /path/to/gtdb/fna/gtdb_genomes_reps_r226 \
    --accession-list ./gtdb_dataset/train_accessions.txt \
    --output ./gtdb_dataset/train_segments.fasta \
    --output-metadata ./gtdb_dataset/train_segments.tsv \
    --total-segments 27000 \
    --threads 16 \
    --seed 42
```

### 5. Merge and Shuffle Datasets

Combine phage and bacteria segments into labeled datasets for classification:

```bash
# Merge all segment lengths (2k, 4k, 8k) for train/dev/test
sbatch merge_datasets.slurm

# Or run manually for a single dataset
python merge_and_shuffle.py \
    --phage-fasta ./inphared_dataset/train_segments.fasta \
    --bacteria-fasta ./gtdb_dataset/train_segments.fasta \
    --output-fasta ./merged_datasets/2k/train_merged.fasta \
    --output-labels ./merged_datasets/2k/train_labels.tsv \
    --seed 42
```

### 6. Create Shuffled Controls

Create negative control sequences by shuffling nucleotides within CDS regions:

```bash
sbatch shuffle_cds.slurm

# Or run manually
python shuffle_cds_nucleotides.py \
    --input-dir ./pharokka_annotations \
    --accession-list ./inphared_dataset/all_representatives.txt \
    --output-dir ./shuffled_genomes \
    --shuffle-intergenic \
    --seed 42
```

### 7. Create BLAST Databases

Create BLAST databases for sequence similarity searches:

```bash
# INPHARED phage database
sbatch create_phage_blastdb.slurm

# GTDB bacterial database (selected genomes only, ~60-80 GB)
sbatch create_gtdb_selected_blastdb.slurm

# GTDB bacterial database (all representatives, ~1+ TB)
sbatch create_gtdb_blastdb.slurm
```

### 8. Remove Prophages from Bacterial Dataset (Optional)

To create a "prophage-free" bacterial dataset, BLAST INPHARED sequences against the selected GTDB database and remove any bacterial genomes with significant hits:

```bash
# BLAST all INPHARED phages against selected GTDB bacteria
blastn \
    -db /path/to/gtdb/blastdb_selected/gtdb_bacteria_selected \
    -query /path/to/inphared/14Apr2025_genomes.fa \
    -out prophage_hits.txt \
    -outfmt 6 \
    -evalue 1e-10 \
    -num_threads 16

# Parse results to identify bacterial genomes with prophage hits
# Remove those accessions from the GTDB accession lists
# Re-run subsampling on the filtered dataset
```

### 9. Evaluate Predictions

```bash
python calculate_prophage_metrics.py \
    --ground-truth ground_truth.csv \
    --predictions predictions.csv \
    --fasta-dir ./genome_fastas \
    --output results.csv
```

## Dataset Structure

After running the pipeline, you will have:

```
inphared_dataset/                 # 2k segments (default)
├── train_accessions.txt          # Training genome accessions
├── dev_accessions.txt            # Development genome accessions
├── test_accessions.txt           # Test genome accessions
├── all_representatives.txt       # All selected accessions
├── train_segments.fasta          # 2000nt training segments
├── dev_segments.fasta            # 2000nt development segments
├── test_segments.fasta           # 2000nt test segments
├── train_segments.tsv            # Segment metadata
└── summary.txt                   # Dataset statistics

inphared_dataset_4k/              # 4k segments
├── train_segments.fasta
├── dev_segments.fasta
└── test_segments.fasta

inphared_dataset_8k/              # 8k segments
├── train_segments.fasta
├── dev_segments.fasta
└── test_segments.fasta

gtdb_dataset/                     # 2k segments (default)
├── train_accessions.txt
├── dev_accessions.txt
├── test_accessions.txt
├── all_accessions.txt
├── train_segments.fasta
├── dev_segments.fasta
├── test_segments.fasta
├── all_metadata.tsv              # Full metadata for selected genomes
└── summary.txt

gtdb_dataset_4k/                  # 4k segments
gtdb_dataset_8k/                  # 8k segments

merged_datasets/                  # Combined phage + bacteria
├── 2k/
│   ├── train_merged.fasta        # Shuffled phage + bacteria segments
│   ├── train_labels.tsv          # Labels (segment_id, label, source)
│   ├── dev_merged.fasta
│   ├── dev_labels.tsv
│   ├── test_merged.fasta
│   └── test_labels.tsv
├── 4k/
└── 8k/

shuffled_genomes/                 # Control sequences
├── <accession>_shuffled.fasta    # Shuffled CDS nucleotides
└── shuffle_summary.txt

blastdb/                          # BLAST databases
├── inphared_phages.*             # INPHARED phage database
├── gtdb_bacteria.*               # All GTDB bacteria (if created)
└── gtdb_bacteria_selected.*      # Selected GTDB bacteria
```

## Key Parameters

### Segment Subsampling
- `--segment-length`: Length of each segment (default: 2000 bp)
- `--min-length`: Minimum genome length to include (default: 5000 bp)
- `--sample-per-bp`: Take 1 sample per N bp (default: 10000)
- `--total-segments`: Target total segments (for balancing datasets)

### Quality Filtering (GTDB only)
- `--min-completeness`: Minimum CheckM2 completeness (default: 95%)
- `--max-contamination`: Maximum CheckM2 contamination (default: 5%)

### Data Splitting
- `--split-ratio`: Train:dev:test ratio (default: 80:10:10)
- `--seed`: Random seed for reproducibility (default: 42)

## Data Leakage Prevention

- **Phage**: Split by vclust cluster (95% ANI) - no cluster spans multiple splits
- **Bacteria**: Split by GTDB genus - no genus spans multiple splits

## GTDB File Structure

GTDB genomes are stored in a nested directory structure with gzipped files:

```
gtdb_genomes_reps_r226/
└── database/
    ├── GCA/
    │   └── XXX/XXX/XXX/
    │       └── GCA_XXXXXXXXX.1_genomic.fna.gz
    └── GCF/
        └── XXX/XXX/XXX/
            └── GCF_XXXXXXXXX.1_genomic.fna.gz
```

The `subsample_gtdb_segments.py` script handles:
- Converting accessions to file paths (e.g., `GB_GCA_964231105.1` → `database/GCA/964/231/105/GCA_964231105.1_genomic.fna.gz`)
- Stripping `GB_` and `RS_` prefixes from GTDB accessions
- Reading gzipped files directly without decompression to disk
- Parallel processing with configurable thread count

## Citation

If you use these tools or the LAMBDA benchmark, please cite:

```
[Citation to be added]
```

## License

[License to be added]

## Contact

[Contact information to be added]

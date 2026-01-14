# LAMBDA Dataset Construction Tools

Tools for constructing the LAMBDA (Large-scale Assessment of Model-Based Detection Algorithms) benchmark datasets for prophage detection in genomic language models.

## Overview

This repository contains scripts for:
- Downloading and processing INPHARED phage genome database
- Selecting and filtering GTDB bacterial representative genomes
- Creating balanced phage/bacterial datasets with cluster-aware train/dev/test splits
- Subsampling genome segments for language model training
- Generating shuffled control sequences
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
    -y

# Or install via pip
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
| `subsample_segments.py` | Subsample 2000nt segments from genomes (supports balanced mode) |
| `shuffle_cds_nucleotides.py` | Shuffle nucleotides within CDS regions to create control sequences |

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

# Bacterial segments (matched to phage count)
python subsample_segments.py \
    --input-fasta gtdb_genomes.fa \
    --accession-list ./gtdb_dataset/train_accessions.txt \
    --output ./gtdb_dataset/train_segments.fasta \
    --total-segments 20000 \
    --seed 42
```

### 5. Create Shuffled Controls

```bash
python shuffle_cds_nucleotides.py \
    --input-dir ./pharokka_annotations \
    --accession-list ./inphared_dataset/all_representatives.txt \
    --output-dir ./shuffled_genomes \
    --shuffle-intergenic \
    --seed 42
```

### 6. Evaluate Predictions

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
inphared_dataset/
├── train_accessions.txt      # Training genome accessions
├── dev_accessions.txt        # Development genome accessions
├── test_accessions.txt       # Test genome accessions
├── train_segments.fasta      # 2000nt training segments
├── dev_segments.fasta        # 2000nt development segments
├── test_segments.fasta       # 2000nt test segments
└── summary.txt               # Dataset statistics

gtdb_dataset/
├── train_accessions.txt
├── dev_accessions.txt
├── test_accessions.txt
├── train_segments.fasta
├── dev_segments.fasta
├── test_segments.fasta
├── all_metadata.tsv          # Full metadata for selected genomes
└── summary.txt
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

## Citation

If you use these tools or the LAMBDA benchmark, please cite:

```
[Citation to be added]
```

## License

[License to be added]

## Contact

[Contact information to be added]

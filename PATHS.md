# LAMBDA Project - Important Paths Reference

> **Note on path prefixes:**
> - `/projects/bfzj/` - Shared project storage (persistent)
> - `/work/hdd/bfzj/` - High-capacity work storage (for large intermediate files)
> - `/u/llindsey1/` - User home/scratch directories

---

## 1. Source Data

### INPHARED (Phage Data)
| Description | Path |
|-------------|------|
| Base directory | `/u/llindsey1/llindsey/black_and_white/data/inphared/` |
| Genomes FASTA | `/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa` |
| vclust clustering | `/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv` |
| Pharokka annotations | `/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka/` |

### GTDB (Bacteria Data)
| Description | Path |
|-------------|------|
| Base directory | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/` |
| Metadata file | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/metadata/bac120_metadata.tsv` |
| Representative genomes | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226/` |

### External Tools
| Description | Path |
|-------------|------|
| Bakta database | `/u/llindsey1/scratch/bakta_db/db` |

---

## 2. BLAST Databases

| Description | Path | Created by |
|-------------|------|------------|
| INPHARED phages | `/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages` | `create_phage_blastdb.slurm` |
| GTDB selected bacteria | `/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected` | `create_gtdb_selected_blastdb.slurm` |
| GTDB full bacteria | `/work/hdd/bfzj/llindsey1/gtdb_blastdb/gtdb_bacteria` | `create_gtdb_blastdb.slurm` |

---

## 3. Scripts Directory

| Description | Path |
|-------------|------|
| Lambda dataset tools | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/` |

---

## 4. Generated Datasets

### 4.1 Phage Segments (INPHARED)
| Length | Directory | Created by |
|--------|-----------|------------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset/` | `subsample_inphared.slurm` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset_4k/` | `subsample_inphared_4k.slurm` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset_8k/` | `subsample_inphared_8k.slurm` |

**Contents per directory:**
- `train_accessions.txt`, `dev_accessions.txt`, `test_accessions.txt`
- `train_segments.fasta`, `dev_segments.fasta`, `test_segments.fasta`
- `train_segments.tsv`, `dev_segments.tsv`, `test_segments.tsv`

### 4.2 Bacteria Segments - Original (with prophage contamination)
| Length | Directory | Created by |
|--------|-----------|------------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset/` | `subsample_gtdb.slurm` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_4k/` | `subsample_gtdb_4k.slurm` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_8k/` | `subsample_gtdb_8k.slurm` |

> ⚠️ **NOT RECOMMENDED** - Contains prophage sequences. Use filtered v3 instead.

### 4.3 Bacteria Segments - Filtered v3 (prophage-free) ✓ RECOMMENDED
| Length | Directory | Created by |
|--------|-----------|------------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3/` | `filter_and_subsample_gtdb_v3.slurm` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3_4k/` | `filter_and_subsample_gtdb_v3.slurm` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3_8k/` | `filter_and_subsample_gtdb_v3.slurm` |

**Additional files in 2k directory:**
- `contaminated_accessions.txt` - Genomes removed due to prophage hits
- `all_accessions.txt` - All clean genome accessions

### 4.4 GC-Content Control (Shuffled Test Segments)
| Description | Path | Created by |
|-------------|------|------------|
| Output directory | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/shuffled_gc_control/` | `shuffle_test_segments.slurm` |

**Output files:**
- `inphared_2k_test_shuffled.fasta`, `inphared_4k_test_shuffled.fasta`, `inphared_8k_test_shuffled.fasta`
- `gtdb_2k_test_shuffled.fasta`, `gtdb_4k_test_shuffled.fasta`, `gtdb_8k_test_shuffled.fasta`

### 4.5 Merged Datasets
| Description | Path | Created by |
|-------------|------|------------|
| Output directory | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/merged_datasets/` | `merge_datasets.slurm` |

---

## 5. Benchmark-Specific Datasets

### 5.1 Phage-Only (Full Genomes with Annotations)
| Description | Path |
|-------------|------|
| Accession list | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset/test_accessions.txt` |
| Full genomes | `/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa` |
| Pharokka annotations | `/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka/` |

**Uses:** 869 test set genomes with functional annotations

### 5.2 Bacteria-Only (100 Taxonomically Balanced with Annotations)
| Description | Path | Created by |
|-------------|------|------------|
| Accession list | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/bacteria_only_dataset/bacteria_only_100_accessions.txt` | `select_bacteria_only_balanced.slurm` |
| Metadata | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/bacteria_only_dataset/bacteria_only_100_metadata.tsv` | `select_bacteria_only_balanced.slurm` |
| Full genomes | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226/` | - |
| Bakta annotations | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/bacteria_only_annotations/` | `run_bakta_bacteria.slurm` |

**Uses:** 100 taxonomically balanced genomes from test set

---

## 6. Intermediate/Working Files

### 6.1 BLAST Results
| Description | Path | Created by |
|-------------|------|------------|
| Prophage BLAST v3 results | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v3/` | `blast_phage_vs_gtdb_selected_v3.slurm` |
| Combined filtered results | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v3/blast_filtered_combined.tsv` | `blast_phage_vs_gtdb_selected_v3.slurm` |
| Contig-to-genome map | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contig_to_genome_map.tsv` | `create_contig_map.slurm` |

### 6.2 Validation Results
| Description | Path | Created by |
|-------------|------|------------|
| Validation directory | `/work/hdd/bfzj/llindsey1/prophage_validation_v3/` | `validate_prophage_removal_v3.slurm` |
| Validation BLAST results | `/work/hdd/bfzj/llindsey1/prophage_validation_v3/blast_results_validation.tsv` | `validate_prophage_removal_v3.slurm` |

---

## 7. Key Files Reference

### INPHARED vclust_info.tsv
- **Path:** `/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv`
- Contains cluster assignments for all phage genomes
- First genome in each cluster is the representative
- Used by `extract_representatives.py` to create train/dev/test splits
- **Total genomes:** 28,665 → **8,683 clusters**

### GTDB bac120_metadata.tsv
- **Path:** `/projects/bfzj/llindsey1/black_and_white/data/gtdb/metadata/bac120_metadata.tsv`
- Contains metadata for all bacterial genomes
- Includes taxonomy, CheckM2 quality scores, genome size
- Used by `select_gtdb_representatives.py` to select high-quality representatives
- **Selected representatives:** 15,865 → **14,894 after prophage filtering**

---

## 8. Quick Reference - Common Paths

```bash
# Scripts
SCRIPT_DIR="/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools"

# Source data
INPHARED_FASTA="/u/llindsey1/llindsey/black_and_white/data/inphared/14Apr2025_genomes.fa"
GTDB_GENOMES="/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226"
GTDB_METADATA="/projects/bfzj/llindsey1/black_and_white/data/gtdb/metadata/bac120_metadata.tsv"

# BLAST databases
PHAGE_BLASTDB="/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages"
BACTERIA_BLASTDB="/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected"

# Recommended datasets (prophage-free)
PHAGE_DATASET="${SCRIPT_DIR}/inphared_dataset"
BACTERIA_DATASET="${SCRIPT_DIR}/gtdb_dataset_filtered_v3"

# Bakta
BAKTA_DB="/u/llindsey1/scratch/bakta_db/db"
```

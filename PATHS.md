# LAMBDA Project - Important Paths Reference

## Source Data

### INPHARED (Phage Data)
| Description | Path |
|-------------|------|
| Base directory | `/projects/bfzj/llindsey1/black_and_white/data/inphared/` |
| Genomes FASTA | `/projects/bfzj/llindsey1/black_and_white/data/inphared/14Apr2025_genomes.fa` |
| vclust clustering | `/projects/bfzj/llindsey1/black_and_white/data/inphared/vclust_info.tsv` |
| Pharokka annotations | `/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka/` |
| BLAST database | `/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages` |

### GTDB (Bacteria Data)
| Description | Path |
|-------------|------|
| Base directory | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/` |
| Taxonomy files | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/taxonomy/` |
| Genome FNA files | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/` |
| Representative genomes | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226/` |
| Metadata file | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/taxonomy/bac120_metadata.tsv` |

---

## BLAST Databases

| Description | Path |
|-------------|------|
| INPHARED phages | `/u/llindsey1/llindsey/black_and_white/data/inphared/blastdb/inphared_phages` |
| GTDB selected bacteria | `/work/hdd/bfzj/llindsey1/gtdb_blastdb_selected/gtdb_bacteria_selected` |
| GTDB full bacteria | `/work/hdd/bfzj/llindsey1/gtdb_blastdb/gtdb_bacteria` |

---

## Scripts Directory

| Description | Path |
|-------------|------|
| Lambda dataset tools | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/` |

---

## Generated Datasets

### Phage Segments (INPHARED)
| Length | Path |
|--------|------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset/` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset_4k/` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset_8k/` |

### Bacteria Segments - Original (with prophage contamination)
| Length | Path |
|--------|------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset/` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_4k/` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_8k/` |

### Bacteria Segments - Filtered v3 (prophage-free)
| Length | Path |
|--------|------|
| 2k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3/` |
| 4k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3_4k/` |
| 8k | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3_8k/` |

### GC-Content Control (Shuffled)
| Description | Path |
|-------------|------|
| Output directory | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/shuffled_gc_control/` |

### Phage-Only Dataset (Full Genomes with Annotations)
| Description | Path |
|-------------|------|
| Accession list | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/inphared_dataset/test_accessions.txt` |
| Full genomes | `/projects/bfzj/llindsey1/black_and_white/data/inphared/14Apr2025_genomes.fa` |
| Pharokka annotations | `/projects/bfzj/llindsey1/black_and_white/data/inphared/pharokka/` |

### Bacteria-Only Dataset (Full Genomes with Annotations)
| Description | Path |
|-------------|------|
| Accession list | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3/test_accessions.txt` |
| Full genomes | `/projects/bfzj/llindsey1/black_and_white/data/gtdb/fna/gtdb_genomes_reps_r226/` |
| Bakta annotations | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/bacteria_only_annotations/` |
| Bakta database | `/projects/bfzj/llindsey1/bakta_db/db` (UPDATE if different) |

---

## Intermediate/Working Files

### BLAST Results
| Description | Path |
|-------------|------|
| Prophage BLAST v2 | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/` |
| Prophage BLAST v3 | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v3/` |
| Contig-to-genome map | `/work/hdd/bfzj/llindsey1/prophage_blast_results_v2/contig_to_genome_map.tsv` |

### Accession Lists
| Description | Path |
|-------------|------|
| Clean bacterial accessions | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3/train_accessions.txt` |
| Contaminated accessions | `/projects/bfzj/llindsey1/black_and_white/scripts/lambda_dataset_tools/gtdb_dataset_filtered_v3/contaminated_accessions.txt` |

---

## Key Files in Source Data

### INPHARED vclust_info.tsv
- Contains cluster assignments for all phage genomes
- First genome in each cluster is the representative
- Used by `extract_representatives.py` to create train/dev/test splits

### GTDB bac120_metadata.tsv
- Contains metadata for all bacterial genomes
- Includes taxonomy, CheckM2 quality scores, genome size
- Used by `select_gtdb_representatives.py` to select high-quality representatives

---

## Notes

- Paths starting with `/projects/` are on the shared project storage
- Paths starting with `/work/hdd/` are on the high-capacity work storage
- Paths starting with `/u/` are in user home directories

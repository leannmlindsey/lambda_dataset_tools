# Unmapped Reads BLAST Analysis Report

**Date:** January 18, 2026
**Analyst:** [Your name]
**Data from:** IgA/IgG viral capture experiment

---

## Overview

Analyzed unmapped reads from the viral capture experiment that did not map to the assembled viral genomes. Goal: identify what these sequences are.

## Methods

- **Input:** 500,000 subsampled reads (10,000 from each of 50 FASTQ files)
- **Databases searched:**
  - INPHARED phage database (~18,000 phage genomes)
  - GTDB bacterial database (~16,000 selected bacterial genomes)
- **BLAST parameters:** E-value threshold 1e-10, megablast

## Key Results

| Category | Reads | Percent |
|----------|------:|--------:|
| **Phage only** | 15,576 | 3.1% |
| **Bacteria only** | 13,467 | 2.6% |
| **Both phage & bacteria** | 109,331 | 21.8% |
| **No hits** | 361,626 | 72.3% |
| **Total** | 500,000 | 100% |

## Interpretation

### 1. Reads hitting both phage AND bacteria (21.8%)
These likely represent **prophage sequences** - viral DNA integrated into bacterial genomes. The IgA/IgG method may have captured:
- Induced prophages still carrying bacterial flanking sequences
- Reads spanning prophage-bacteria boundaries
- Temperate phages with homology to both databases

### 2. Phage-only reads (3.1%)
These are likely **true phage sequences** not present in bacterial genomes - possibly:
- Lytic phages
- Novel phages worth investigating further

### 3. Bacteria-only reads (2.6%)
Possible explanations:
- Bacterial contamination in the sample
- Bacteria-derived sequences in viral particles (generalized transduction)
- Bacterial sequences not homologous to known phages

### 4. No hits (72.3%)
Could be:
- **Eukaryotic viruses** (INPHARED only contains prokaryotic phages)
- **Novel/divergent sequences** not in databases
- **Host contamination** (human sequences not fully removed)
- **Low-quality sequences**

## Top Hits

### Top Phages
| Hits | Phage Accession |
|-----:|-----------------|
| 34,655 | PP129558 |
| 26,775 | MZ092003 |
| 24,607 | LT985568 |
| 22,796 | MK170160 |
| 12,806 | LT993247 |

### Top Bacteria
| Hits | Bacterial Contig |
|-----:|------------------|
| 20,412 | CABHON010000015.1 |
| 15,481 | CABHON010000009.1 |
| 15,080 | CABHON010000005.1 |
| 8,765 | CALMMC010000168.1 |
| 8,375 | NZ_CP046176.1 |

*Note: The CABHON contigs are likely from the same bacterial genome - worth identifying this organism.*

## Recommendations

1. **For the "no hits" reads (72.3%):**
   - BLAST against NCBI nt database to identify eukaryotic viruses, human contamination, etc.
   - Analysis submitted to Biowulf

2. **For the phage-only reads (3.1%):**
   - These may represent novel phages worth assembling/characterizing

3. **Identify the dominant bacteria:**
   - The CABHON genome appears repeatedly - identify this organism
   - May indicate a specific bacterial contaminant or host

4. **Consider running on full dataset:**
   - Current analysis used 10k reads per sample
   - Full analysis would give exact counts

## Output Files

All results saved to: `/work/hdd/bfzj/llindsey1/unmapped_blast_results/`

| File | Description |
|------|-------------|
| `summary.txt` | Overall summary statistics |
| `phage_blast_results.tsv` | All phage BLAST hits |
| `bacteria_blast_results.tsv` | All bacteria BLAST hits |
| `phage_only_read_ids.txt` | Read IDs hitting only phage |
| `bacteria_only_read_ids.txt` | Read IDs hitting only bacteria |
| `both_read_ids.txt` | Read IDs hitting both |
| `top_bacteria_analysis.txt` | Detailed bacteria hit analysis |

## Next Steps

- [ ] Identify CABHON bacterial genome
- [ ] Run nt BLAST on no-hit reads (Biowulf)
- [ ] Run full dataset analysis (optional)
- [ ] Investigate top phage-only sequences

---

*Analysis performed using INPHARED and GTDB BLAST databases*

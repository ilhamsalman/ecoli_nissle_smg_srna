# Workflow

## Overview

This workflow describes the stepwise analysis of RNA-seq data from
*Escherichia coli* Nissle 1917 under simulated microgravity (SMG) versus normal gravity (NG).

The pipeline is organized into modular stages, each of which can be executed independently.

---

## 00_pipeline — Read processing and alignment

### Purpose

* Trim reads
* Align reads to reference genome
* Infer library strandedness

### Dependencies

* fastp
* bowtie2
* samtools
* RSeQC

### Inputs

* Raw paired-end FASTQ files
* Reference genome (FASTA)

### Scripts

```bash
bash scripts/00_pipeline/02_trim_fastq.sh
bash scripts/00_pipeline/03_align_bowtie2.sh
bash scripts/00_pipeline/04_infer_strandedness.sh
```

### Outputs

* Trimmed FASTQ files
* Sorted BAM files
* Strandedness report

---

## 01_annogesic — sRNA prediction

### Purpose

* Generate coverage tracks (WIG)
* Perform genome-wide sRNA detection

### Dependencies

* ANNOgesic
* bedtools
* samtools
* Python

### Inputs

* Sorted BAM files
* Reference genome

### Scripts

```bash
bash scripts/01_annogesic/make_wigs.sh
bash scripts/01_annogesic/annogesic.sh
```

### Outputs

* Predicted sRNA annotations (GFF, tables)
* Results in `results/sRNA_detection/01_annogesic/`

---

## 02_filter — Filtering predicted sRNAs

### Purpose

* Apply filtering criteria to ANNOgesic predictions

### Dependencies

* Python 3
* pandas (if used)

### Scripts

```bash
python scripts/02_filter/filter_annogesic_srna.py
```

### Outputs

* Filtered sRNA annotations

---

## 03_shared_srna — Shared sRNA identification

### Purpose

* Identify shared sRNAs across conditions
* Convert annotations to SAF format

### Dependencies

* Python 3

### Scripts

```bash
python scripts/03_shared_srna/make_saf_from_annogesic.py
python scripts/03_shared_srna/shared_srna_from_saf.py
```

### Outputs

* SAF files
* Shared sRNA tables

---

## 04_srna_de — Differential expression analysis

### Purpose

* Quantify expression of sRNAs
* Identify differentially expressed sRNAs

### Dependencies

* featureCounts
* R
* limma
* edgeR

### Scripts

```bash
bash scripts/04_srna_de/01_run_featurecounts_filtered.sh
Rscript scripts/04_srna_de/02_run_voom_phase_tables.R
Rscript scripts/04_srna_de/03_compare_sig_exp_sta_by_SAF.R
```

### Outputs

* Count matrices
* Differential expression tables

---

## 05_annotate_known_srna — Annotation against known sRNAs

### Purpose

* Compare predicted sRNAs to known databases
* Annotate known vs novel sRNAs

### Dependencies

* Python 3
* bedtools
* BLAST+ (blastn)
* R

### Scripts

```bash
python scripts/05_annotate_known_srna/01_make_bed_from_annogesic.py
bash scripts/05_annotate_known_srna/02_extract_fasta.sh
bash scripts/05_annotate_known_srna/03_blast_vs_bsrd.sh
python scripts/05_annotate_known_srna/04_filter_annotate_escherichia.py
Rscript scripts/05_annotate_known_srna/05_intersect_known_srna.R
```

### Outputs

* Annotated sRNA tables
* Known vs novel classification

---

## 06_srna_deg — sRNA–gene interaction analysis

### Purpose

* Identify overlaps between sRNAs and differentially expressed genes
* Detect antisense relationships

### Dependencies

* Python 3
* bedtools

### Scripts

```bash
python scripts/06_srna_deg/01_sig_antisense_srna_to_bed.py
bash scripts/06_srna_deg/02_deg_tables_to_bed.sh
bash scripts/06_srna_deg/03_bedtools_antisense_overlap.sh
bash scripts/06_srna_deg/04_overlap_expression_patterns.sh
```

### Outputs

* Overlap tables
* Antisense interaction results

---

## 07_intarna — RNA–RNA interaction prediction

### Purpose

* Predict interactions between sRNAs and mRNA targets

### Dependencies

* Python 3
* IntaRNA

### Scripts

```bash
python scripts/07_intarna/01_extract_sequences.py
bash scripts/07_intarna/02_run_intarna.sh
bash scripts/07_intarna/03_combine_results.sh
python scripts/07_intarna/04_filter_results.py
```

### Outputs

* Interaction prediction tables
* Filtered candidate interactions

---

## 08_conservation — Conservation analysis

### Purpose

* Assess conservation of sRNAs across *E. coli* strains

### Dependencies

* BLAST+ (blastn)

### Scripts

```bash
bash scripts/08_conservation/01_download_genomes.sh
bash scripts/08_conservation/02_blast_and_summarize.sh
```

### Outputs

* Conservation summaries
* BLAST results

---

## 09_figures — Visualization

### Purpose

* Generate publication-ready figures

### Dependencies

* R
* Python (matplotlib/seaborn if used)

### Scripts

```bash
Rscript scripts/09_figures/make_heatmap_sig_from_counts.R
Rscript scripts/09_figures/make_pca_from_counts.R
Rscript scripts/09_figures/make_volcano_from_voom_csv.R
Rscript scripts/09_figures/make_srna_type_distribution.R
python scripts/09_figures/make_srna_venn_pubstyle.py
```

### Outputs

* Heatmaps
* PCA plots
* Volcano plots
* Venn diagrams

---


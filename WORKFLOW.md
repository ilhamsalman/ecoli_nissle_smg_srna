# Workflow

## Overview

This workflow describes the stepwise analysis of RNA-seq data from
*Escherichia coli* Nissle 1917 under simulated microgravity (SMG) versus normal gravity (NG).

The pipeline is organized into numbered analysis stages. Scripts may require inputs produced by earlier stages; review script configuration and input paths before execution. The commands below identify the scripts present in this repository.

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
bash scripts/00_pipeline/01_download_fastq.sh
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

* Filtered sRNA annotations in `results/sRNA_detection/02_filter/`

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

* SAF files and shared-locus tables in `results/sRNA_detection/03_shared_filtered/`

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
Rscript scripts/04_srna_de/03_compare_candidates_exp_sta_by_SAF.R
```

### Outputs

* Count matrices
* Differential expression tables, including complete results in `results/sRNA_detection/04_srna_de_filtered/`

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

## 06_intarna — RNA–RNA interaction prediction

### Purpose

* Predict interactions between sRNAs and mRNA targets

### Dependencies

* Python 3
* IntaRNA

### Scripts

```bash
python scripts/06_intarna/01_extract_sequences.py
bash scripts/06_intarna/02_run_intarna.sh
bash scripts/06_intarna/03_combine_results.sh
python scripts/06_intarna/04_filter_candidates.py
```

### Outputs

* Interaction prediction tables and filtered candidate interactions in `results/sRNA_detection/06_intarna/`

---

## 07_conservation — Conservation analysis

### Purpose

* Assess conservation of sRNAs across *E. coli* strains

### Dependencies

* BLAST+ (blastn)

### Scripts

```bash
bash scripts/07_conservation/conservation_from_selected_srnas.sh
```

### Outputs

* Conservation summaries
* BLAST results

---

## 08_figures — Visualization

### Purpose

* Generate publication-ready figures

### Dependencies

* R
* Python (matplotlib/seaborn if used)

### Scripts

```bash
Rscript scripts/08_figures/make_volcano_from_voom_csv.R
Rscript scripts/08_figures/make_srna_type_distribution.R
Rscript scripts/08_figures/make_cis_element_class_distribution.R
python scripts/08_figures/make_srna_venn_pubstyle.py
```

### Outputs

* Cis-regulatory element class-distribution plots
* Volcano plots
* sRNA-type distribution plot
* Venn diagram

---


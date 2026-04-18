# *Escherichia coli* Nissle 1917 — sRNA responses under simulated microgravity

This repository contains scripts, metadata, and curated analysis outputs for the identification and characterization of small RNAs (sRNAs) in *Escherichia coli* Nissle 1917 under simulated microgravity (SMG) compared with normal gravity (NG).

The analysis focuses on:

* Genome-wide sRNA detection using **ANNOgesic**
* Condition-specific comparison across growth phases (exponential vs stationary)
* Prediction of regulatory RNA elements (riboswitches and RNA thermometers)

---

## Dataset

**GEO accession:** GSE147272
*Transcriptional profiling of probiotic* *E. coli* Nissle 1917 *under simulated microgravity*

### RNA-seq runs (SRA)

**Simulated microgravity (SMG)**

* Exponential: SRR11357510, SRR11357511
* Stationary: SRR11357512, SRR11357513

**Normal gravity (NG)**

* Exponential: SRR11357514, SRR11357515
* Stationary: SRR11357516, SRR11357517

---

## Sample mapping

| sample_id | SRR         | growth_phase | gravity |
| --------- | ----------- | ------------ | ------- |
| MG_exp_1  | SRR11357510 | exponential  | MG      |
| MG_exp_2  | SRR11357511 | exponential  | MG      |
| MG_sta_1  | SRR11357512 | stationary   | MG      |
| MG_sta_2  | SRR11357513 | stationary   | MG      |
| NG_exp_1  | SRR11357514 | exponential  | NG      |
| NG_exp_2  | SRR11357515 | exponential  | NG      |
| NG_sta_1  | SRR11357516 | stationary   | NG      |
| NG_sta_2  | SRR11357517 | stationary   | NG      |

---

## Repository structure

```
scripts/     → pipeline scripts organized by analysis stage
results/     → curated outputs and final results
meta/        → sample metadata and input tables
```

---

## Workflow

The analysis is organized into numbered stages under `scripts/`.

### 00 — Data acquisition and preprocessing

* Download FASTQ from SRA
  `scripts/00_pipeline/01_download_fastq.sh`
* Read trimming (fastp)
  `scripts/00_pipeline/02_trim_fastq.sh`
* Alignment (Bowtie2)
  `scripts/00_pipeline/03_align_bowtie2.sh`
* Strandedness inference (RSeQC)
  `scripts/00_pipeline/04_infer_strandedness.sh`

---

### 01 — ANNOgesic preparation and annotation

* WIG generation
* ANNOgesic execution for sRNA prediction
  `scripts/01_annogesic/`

---

### 02 — Filtering

* Filtering ANNOgesic sRNA predictions
  `scripts/02_filter/filter_annogesic_srna.py`

---

### 03 — Shared sRNA detection

* Generate SAF files
* Identify shared sRNAs across conditions
  `scripts/03_shared_srna/`

---

### 04 — Differential expression analysis

* FeatureCounts quantification
* Voom/limma differential expression
  `scripts/04_srna_de/`

---

### 05 — Known sRNA annotation

* BLAST against BSRD
* Annotation and filtering
  `scripts/05_annotate_known_srna/`


---

### 06 — sRNA–mRNA interaction prediction

* IntaRNA analysis
  `scripts/06_intarna/`

---

### 07 — Conservation analysis

* Genome collection and BLAST analysis
  `scripts/07_conservation/`

---

###  Visualization

* volcano plots, Venn diagrams
  `scripts/figures/`

---

## Reference genome

* NCBI accession: **CP007799.1**
* Used for alignment and annotation steps

---


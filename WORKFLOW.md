# Workflow

## Scope and execution notes

This document describes the repository's **current numbered scripts and committed output locations** for computational analysis of *Escherichia coli* Nissle 1917 RNA-seq data under simulated microgravity (SMG) and normal gravity (NG). Exponential and stationary phases are processed separately where applicable.

Run commands from the **repository root**. They are stage-level entry points, not a guarantee that every step will run on a fresh clone: the pipeline also requires appropriate local FASTQ files, reference genome/annotation files, software environments, and inputs from earlier stages. Check the input paths configured inside each script before execution.

The first eight analysis stages already existed; the revised workflow separates cis-element–DEG comparison into **stage 08** and moves figure generation to **stage 09**.

## Stage 00 — `scripts/00_pipeline/`: read processing and alignment

**Purpose:** Retrieve paired-end reads, trim reads with fastp, align to the EcN reference using Bowtie2, and infer library strandedness using RSeQC.

**Main tools:** SRA tools, fastp, Bowtie2, samtools, RSeQC.

```bash
bash scripts/00_pipeline/01_download_fastq.sh
bash scripts/00_pipeline/02_trim_fastq.sh
bash scripts/00_pipeline/03_align_bowtie2.sh
bash scripts/00_pipeline/04_infer_strandedness.sh
```

**Inputs:** GSE147272-associated SRA runs and EcN reference FASTA. **Outputs:** trimmed reads, aligned/sorted BAM files, strandedness assessment. Downstream feature counting uses reverse-stranded libraries.

## Stage 01 — `scripts/01_annogesic/`: ANNOgesic prediction

**Purpose:** Generate WIG coverage tracks and identify candidate sRNAs and additional regulatory RNA-element classes through ANNOgesic/Rfam.

**Main tools:** ANNOgesic, bedtools, samtools, Python.

```bash
bash scripts/01_annogesic/make_wigs.sh
bash scripts/01_annogesic/annogesic.sh
```

**Outputs:** sRNA annotations under `results/sRNA_detection/01_annogesic/`; curated Rfam-based cis-element exports are available in `results/riboswitch_thermometer/exp/` and `sta/`.

## Stage 02 — `scripts/02_filter/`: sRNA filtering

**Purpose:** Retain candidate sRNA loci passing transcript-length and coverage/expression criteria; filtering uses ANNOgesic's reported `End - Start` length convention (50–500 nt) and a phase-specific coverage threshold.

```bash
python scripts/02_filter/filter_annogesic_srna.py
```

**Outputs:** `results/sRNA_detection/02_filter/`.

## Stage 03 — `scripts/03_shared_srna/`: phase comparison

**Purpose:** Generate SAF-formatted loci and identify shared sRNA genomic regions across growth phases.

```bash
python scripts/03_shared_srna/make_saf_from_annogesic.py
python scripts/03_shared_srna/shared_srna_from_saf.py
```

**Outputs:** `results/sRNA_detection/03_shared_filtered/`.

**Interpretation:** The filtered ANNOgesic sets contain 3,752 exponential and 519 stationary sRNA loci, including 24 loci meeting the specified reciprocal-overlap criterion. This comparison is distinct from the 79 nominal expression-change records (78 distinct loci).

## Stage 04 — `scripts/04_srna_de/`: sRNA differential expression

**Purpose:** Count sRNA loci, filter low-count features, and compare SMG with NG separately within each phase using edgeR filtering, TMM normalisation, and limma–voom, with Benjamini–Hochberg FDR adjustment.

**Main tools:** featureCounts, R, limma, edgeR.

```bash
bash scripts/04_srna_de/01_run_featurecounts_filtered.sh
Rscript scripts/04_srna_de/02_run_voom_phase_tables.R
Rscript scripts/04_srna_de/03_compare_candidates_exp_sta_by_SAF.R
```

**Outputs:** `results/sRNA_detection/04_srna_de_filtered/`, including `voom_exp_s2_all.csv`, `voom_sta_s2_all.csv`, and phase-specific candidate tables.

**Interpretation:** 697 exponential and 317 stationary loci were tested; 44 and 35 candidate records met **nominal P < 0.05**, respectively. No candidate sRNA remained significant after multiple-testing correction. Do not describe these as FDR-significant DE sRNAs.

## Stage 05 — `scripts/05_annotate_known_srna/`: known sRNA comparison

**Purpose:** Compare predicted sRNA sequences against BSRD and annotate potential known-sRNA matches.

**Main tools:** Python, BLAST+, bedtools, R.

```bash
python scripts/05_annotate_known_srna/01_make_bed_from_annogesic.py
bash scripts/05_annotate_known_srna/02_extract_fasta.sh
bash scripts/05_annotate_known_srna/03_blast_vs_bsrd.sh
python scripts/05_annotate_known_srna/04_filter_annotate_escherichia.py
Rscript scripts/05_annotate_known_srna/05_intersect_known_srna.R
```

**Outputs:** `results/sRNA_detection/05_known_srna/`. Sequence matches to GcvB/RydB are homology-based candidate annotations.

## Stage 06 — `scripts/06_intarna/`: selected sRNA–mRNA interaction analysis

**Purpose:** Extract selected intergenic sRNA sequences and phase-specific curated DEG target sequences; run IntaRNA; combine and filter predicted interactions.

**Main tools:** Python 3, Biopython, pandas, IntaRNA; a reference FASTA and matching GenBank annotation are needed for sequence extraction.

**Selected inputs:**

- `results/sRNA_detection/06_intarna/selected_exp_ids.txt`
- `results/sRNA_detection/06_intarna/selected_sta_ids.txt`
- `results/sRNA_detection/06_intarna/curated_mrna_tables/TableS1_EXP_curated.tsv`
- `results/sRNA_detection/06_intarna/curated_mrna_tables/TableS2_STA_curated.tsv`

The selected spreadsheets `selected_exp_srnas.xlsx` and `selected_sta_srnas.xlsx` are also retained under `results/sRNA_detection/06_intarna/`.

### 06a. Extract selected sRNAs

```bash
python scripts/06_intarna/01a_extract_srna_sequences.py
```

This script expects the reference at `ref/EcN_CP007799.1.fa` and phase-specific BED inputs under `results/sRNA_detection/05_known_srna/`. Confirm or adjust those paths locally.

**Outputs:** `results/sRNA_detection/06_intarna/sequences/srnas/exp/` and `sta/`.

### 06b. Extract curated mRNA targets

The `01b` script requires explicit paths to the reference FASTA and GenBank annotation. Supply the actual local file paths; the placeholders below are not repository files.

```bash
python scripts/06_intarna/01b_extract_mrna_sequences.py \
  --table results/sRNA_detection/06_intarna/curated_mrna_tables/TableS1_EXP_curated.tsv \
  --genome /path/to/EcN_reference.fa \
  --annotation /path/to/EcN_annotation.gb \
  --phase EXP \
  --outdir results/sRNA_detection/06_intarna/sequences/mrnas/deg_exp

python scripts/06_intarna/01b_extract_mrna_sequences.py \
  --table results/sRNA_detection/06_intarna/curated_mrna_tables/TableS2_STA_curated.tsv \
  --genome /path/to/EcN_reference.fa \
  --annotation /path/to/EcN_annotation.gb \
  --phase STA \
  --outdir results/sRNA_detection/06_intarna/sequences/mrnas/deg_sta
```

**Target windows:** Default 200 nt upstream of the CDS start plus the first 500 nt of the CDS (check script options if changing these values). Duplicate gene symbols are disambiguated with locus tags, e.g. `kfiA__ECOLIN_25915.fasta`.

### 06c. Run IntaRNA

```bash
bash scripts/06_intarna/02_run_intarna.sh --threads 4
```

The script pairs selected sRNAs only with **same-phase DEG targets**. Output CSVs are written to `results/sRNA_detection/06_intarna/intarna_exp/` and `intarna_sta/`.

### 06d. Combine predictions and organise the final tables

```bash
bash scripts/06_intarna/03_combine_results.sh
```

**Important current-path detail:** This script writes combined CSV/TSV files into `results/sRNA_detection/06_intarna/` **itself**. The repository's curated final copies reside in `results/sRNA_detection/06_intarna/combined_results/`. After regeneration, move the combined files into that directory, or update the script output directory before running it. Do not assume the existing script writes there automatically. The filtering script checks the root location first, so remove or relocate stale root-level copies before filtering the curated `combined_results/` tables.

### 06e. Filter and prioritise interactions

```bash
python scripts/06_intarna/04_filter_interactions.py
```

The filtering script accepts combined results either at the `06_intarna/` root or inside `combined_results/`. It applies an energy criterion of **ΔG ≤ −10 kcal/mol** and records filtered/prioritised predicted interactions.

**Outputs:** `results/sRNA_detection/06_intarna/filtered_interactions/` (including `exp_filtered_interactions.tsv`, `sta_filtered_interactions.tsv`, locus-level TSV files, and summary files). These are **in silico** interaction predictions, not experimentally validated targets.

## Stage 07 — `scripts/07_conservation/`: cross-strain conservation

**Purpose:** BLAST selected sRNA sequences against *E. coli* reference strains and summarise sequence conservation at the configured identity/coverage thresholds.

```bash
bash scripts/07_conservation/conservation_from_selected_srnas.sh \
  results/sRNA_detection/06_intarna/selected_exp_srnas.xlsx \
  results/sRNA_detection/06_intarna/selected_sta_srnas.xlsx
```

The script reads selected-sRNA spreadsheets (or paths passed as its first two arguments), expects locally available genome/BLAST inputs, and generates tables under `results/sRNA_detection/07_conservation/`. Check its local configuration before running. The manuscript discusses K-12 MG1655, O157:H7 Sakai and CFT073.

## Stage 08 — `scripts/08_cis_deg/`: cis-element–DEG association

**Purpose:** Compare ANNOgesic `Associated_CDS` locus tags with phase-specific curated DEG entries. This is a **locus-tag association**, not a test of coordinate overlap or causal RNA regulation.

**Inputs:**

- `results/riboswitch_thermometer/exp/annotation_RNA_thermometer.csv`
- `results/riboswitch_thermometer/exp/annotation_riboswitch.csv`
- `results/riboswitch_thermometer/sta/annotation_RNA_thermometer.csv`
- `results/riboswitch_thermometer/sta/annotation_riboswitch.csv`
- `results/sRNA_detection/06_intarna/curated_mrna_tables/TableS1_EXP_curated.tsv`
- `results/sRNA_detection/06_intarna/curated_mrna_tables/TableS2_STA_curated.tsv`
- Phase-specific curated DEG mRNA FASTA headers when symbol-to-locus resolution is required.

```bash
python scripts/08_cis_deg/compare_cis_elements_to_degs.py
```

**Outputs:**

- `results/riboswitch_thermometer/deg_comparison/all_predictions_with_DEG_flags.tsv`
- `results/riboswitch_thermometer/deg_comparison/same_phase_DEG_associations.tsv`
- `results/riboswitch_thermometer/deg_comparison/deg_mapping_audit.tsv`
- `results/riboswitch_thermometer/deg_comparison/summary.tsv`

**Expected phase-labelled module records:** 43 exponential/14 stationary RNA-thermometer-module predictions and 27 exponential/11 stationary riboswitch-module predictions. There are **two same-phase DEG associations**, involving *rcnA* (exponential; TPP class) and *manA* (stationary; Phe_leader). The results do not establish genomic overlap or functional regulation.

## Stage 09 — `scripts/09_figures/`: reproducible figures

**Purpose:** Regenerate figures from the relevant analysis tables.

```bash
# Figure 5, both ANNOgesic cis-element panels (directly runnable from repo root)
Rscript scripts/09_figures/make_cis_element_class_distribution.R

# Shared-sRNA Venn diagram (counts are configured in this script)
python scripts/09_figures/make_srna_venn_pubstyle.py

# sRNA-type distributions: this script requires exponential/stationary input tables
Rscript scripts/09_figures/make_srna_type_distribution.R \
  results/sRNA_detection/02_filter/srna_exponential.filtered.tsv \
  results/sRNA_detection/02_filter/srna_stationary.filtered.tsv \
  results/figures/srna_type/srna_type_distribution

# Volcano plots: each run requires an all-results CSV and output prefix
Rscript scripts/09_figures/make_volcano_from_voom_csv.R \
  results/sRNA_detection/04_srna_de_filtered/voom_exp_s2_all.csv \
  results/figures/volcano/volcano_exp_s2
Rscript scripts/09_figures/make_volcano_from_voom_csv.R \
  results/sRNA_detection/04_srna_de_filtered/voom_sta_s2_all.csv \
  results/figures/volcano/volcano_sta_s2
```

**Figure output directories:** `results/figures/cis_elements/`, `shared_srna/`, `srna_type/`, and `volcano/`.

**Figure 5 sources and checks:** Figure 5 is generated from the **GitHub-tracked** phase-specific ANNOgesic CSV exports. The script checks the expected module-record counts and writes `results/figures/cis_elements/figure_5_plot_counts.csv` as a count audit. The ANNOgesic modules include other RNA classes, notably Phe_leader (leader RNA), suhB (sRNA), and drz-agam-1 (ribozyme); these should not be reported as established thermometers/riboswitches.

## Interpretation and reproducibility safeguards

- SMG-versus-NG log2 fold changes and nominal P-value candidate calls must not be presented as FDR-significant results.
- The 79 phase-specific nominal sRNA candidate records represent 78 genomic loci because one locus appears in both phases.
- IntaRNA free-energy ranking, predicted RNA structures, conservation matches, and `Associated_CDS` matches are predictions/annotations requiring experimental verification.
- Preserve original ANNOgesic phase provenance (`exp` and `sta`) and the original prediction tables when regenerating downstream plots and comparisons.
- Review the local reference-input paths and package versions before rerunning a stage; scripts with positional/required command-line inputs cannot be run with no arguments.

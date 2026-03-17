# 04_srna_de — sRNA differential expression analysis

This folder quantifies **filtered ANNOgesic-predicted sRNAs** and performs differential expression (DE)
between MG and NG conditions separately for exponential and stationary phases.

Workflow:
1) featureCounts quantification using SAF annotations
2) limma–voom DE analysis
3) Optional overlap of significant loci between phases using SAF coordinates

---

## Inputs

### SAF annotation files (filtered)
- results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf
- results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf

NOTE: SAF contig names must match BAM headers (e.g., CP007799.1).

### BAM files
Exponential phase:
- annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/MG_exp_1.bam
- annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/MG_exp_2.bam
- annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/NG_exp_1.bam
- annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/NG_exp_2.bam

Stationary phase:
- annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/MG_sta_1.bam
- annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/MG_sta_2.bam
- annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/NG_sta_1.bam
- annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/NG_sta_2.bam

### Filtered ANNOgesic tables (for sRNA type annotation)
- results/sRNA_detection/02_filter/srna_exponential.filtered.tsv
- results/sRNA_detection/02_filter/srna_stationary.filtered.tsv

---

## Step 1 — Quantify sRNAs with featureCounts

Run from project root:

bash scripts/04_srna_de/01_run_featurecounts_filtered.sh

Outputs:
- results/sRNA_detection/04_srna_de_filtered/counts_exp_s1.txt
- results/sRNA_detection/04_srna_de_filtered/counts_exp_s2.txt
- results/sRNA_detection/04_srna_de_filtered/counts_sta_s1.txt
- results/sRNA_detection/04_srna_de_filtered/counts_sta_s2.txt
(+ corresponding *.summary files)

---

## Step 2 — Differential expression with limma–voom

Run DE separately for each phase (recommended stranded mode: s1):

Rscript scripts/04_srna_de/02_run_voom_phase_tables.R \
  results/sRNA_detection/04_srna_de_filtered/counts_exp_s1.txt \
  results/sRNA_detection/02_filter/srna_exponential.filtered.tsv \
  exp__ exponential \
  results/sRNA_detection/04_srna_de_filtered/voom_exp_s1

Rscript scripts/04_srna_de/02_run_voom_phase_tables.R \
  results/sRNA_detection/04_srna_de_filtered/counts_sta_s1.txt \
  results/sRNA_detection/02_filter/srna_stationary.filtered.tsv \
  sta__ stationary \
  results/sRNA_detection/04_srna_de_filtered/voom_sta_s1

Outputs per run:
- <prefix>_all.csv          full results table
- <prefix>_significant.csv  significant loci (FDR ≤ 0.05)

LogFC is reported as MG vs NG (positive = up in MG).

---

## Step 3 — Overlap of significant loci between phases (optional)

# Usage:
# Rscript 03_compare_sig_exp_sta_by_SAF.R <exp_sig> <sta_sig> <exp_saf> <sta_saf> <min_recip> <out_csv>
#
# Example:
# Rscript scripts/04_srna_de/03_compare_sig_exp_sta_by_SAF.R \
#   results/sRNA_detection/04_srna_de_filtered/voom_exp_s1_significant.csv \
#   results/sRNA_detection/04_srna_de_filtered/voom_sta_s1_significant.csv \
#   results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf \
#   results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf \
#   0.5 \
#   results/sRNA_detection/04_srna_de_filtered/shared_sig_exp_sta_minrecip_0.5_shared_overlap.csv

Outputs:
- results/sRNA_detection/04_srna_de_filtered/voom_exp_s1_withSAF.csv
- results/sRNA_detection/04_srna_de_filtered/voom_sta_s1_withSAF.csv
- results/sRNA_detection/04_srna_de_filtered/shared_sig_exp_sta_by_SAF_overlap.csv

Set overlap stringency

# 04_srna_de — sRNA expression analysis

This module quantifies filtered ANNOgesic-predicted sRNAs and performs
expression analysis between MG and NG conditions for exponential and
stationary phases.

## Workflow

1) featureCounts quantification  
2) limma–voom analysis  
3) overlap of candidate loci between phases  

---

## Full pipeline (run from project root)

```bash
# Step 1 — Quantification
bash scripts/04_srna_de/01_run_featurecounts_filtered.sh

# Step 2 — limma-voom analysis (exponential)
Rscript scripts/04_srna_de/02_run_voom_phase_tables.R \
  results/sRNA_detection/04_srna_de_filtered/counts_exp_s2.txt \
  results/sRNA_detection/02_filter/srna_exponential.filtered.tsv \
  exp__ exponential \
  results/sRNA_detection/04_srna_de_filtered/voom_exp_s2

# Step 2 — limma-voom analysis (stationary)
Rscript scripts/04_srna_de/02_run_voom_phase_tables.R \
  results/sRNA_detection/04_srna_de_filtered/counts_sta_s2.txt \
  results/sRNA_detection/02_filter/srna_stationary.filtered.tsv \
  sta__ stationary \
  results/sRNA_detection/04_srna_de_filtered/voom_sta_s2

# Step 3 — overlap between phases
Rscript scripts/04_srna_de/03_compare_candidates_exp_sta_by_SAF.R \
  results/sRNA_detection/04_srna_de_filtered/voom_exp_s2_candidates.csv \
  results/sRNA_detection/04_srna_de_filtered/voom_sta_s2_candidates.csv \
  results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf \
  results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf \
  0.5 \
  results/sRNA_detection/04_srna_de_filtered/shared_candidates_exp_sta_minrecip_0.5.csv
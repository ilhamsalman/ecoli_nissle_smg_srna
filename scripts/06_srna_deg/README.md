# sRNA–DEG Antisense Overlap Analysis

This module performs coordinate-based antisense overlap analysis between
significant antisense sRNAs and differentially expressed genes (DEGs)
in exponential and stationary phases.

All commands must be executed from the project root:

cd ~/projects/ecoli_nissle_smg_srna


---------------------------------------------------------------------

## STEP 1 — Extract Significant Antisense sRNAs to BED

Script:
scripts/06_srna_deg/01_sig_antisense_srna_to_bed.py

Run:

python3 scripts/06_srna_deg/01_sig_antisense_srna_to_bed.py \
  --exp_sig results/sRNA_detection/04_srna_de_filtered/voom_exp_s1_significant.csv \
  --sta_sig results/sRNA_detection/04_srna_de_filtered/voom_sta_s1_significant.csv \
  --exp_saf results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf \
  --sta_saf results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf \
  --outdir results/sRNA_detection/06_srna_deg/outputs

Outputs:

results/sRNA_detection/06_srna_deg/outputs/
├── sig_antisense_srna_exp.bed
└── sig_antisense_srna_sta.bed

---------------------------------------------------------------------

## STEP 2 — Convert DEG Tables (with Coordinates) to BED6

Input tables (must contain chromosome, start, end, strand):

results/sRNA_detection/06_srna_deg/inputs/
├── DEG_exp_with_coordinates.tsv
└── DEG_sta_with_coordinates.tsv

Run:

bash scripts/06_srna_deg/02_deg_tables_to_bed.sh

This generates:

results/sRNA_detection/06_srna_deg/inputs/
├── DEG_exp_with_fc.bed
└── DEG_sta_with_fc.bed



---------------------------------------------------------------------

## STEP 3 — Antisense Coordinate Overlap

Script:
scripts/06_srna_deg/03_bedtools_antisense_overlap.sh

Run:

bash scripts/06_srna_deg/03_bedtools_antisense_overlap.sh

Internally executes:

bedtools intersect -S -wa -wb

Outputs:

results/sRNA_detection/06_srna_deg/outputs/
├── antisense_deg_exp_overlap.bed
└── antisense_deg_sta_overlap.bed



Strand specificity is enforced using -S (opposite strand only).

---------------------------------------------------------------------

## STEP 4 — Expression Pattern Classification

Script:
scripts/06_srna_deg/04_overlap_expression_patterns.sh

Run:

bash scripts/06_srna_deg/04_overlap_expression_patterns.sh

Outputs:

results/sRNA_detection/06_srna_deg/outputs/
├── antisense_deg_exp_analysis.tsv
└── antisense_deg_sta_analysis.tsv


---------------------------------------------------------------------



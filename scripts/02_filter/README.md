# sRNA ANNOgesic Filtering

## Purpose
This step defines a **high-confidence sRNA candidate set** prior to downstream analyses
(shared-locus comparison, featureCounts quantification, and differential expression).

Filtering reduces low-confidence transcript fragments and improves statistical robustness.
Filtering is applied **per growth phase** to avoid bias from different coverage distributions.

---

## Script
filter_annogesic_srna.py

---

## Inputs
ANNOgesic sRNA tables (CSV or TSV), typically located in:

results/sRNA_detection/01_annogesic/

Example files:
- srna_exponential.csv
- srna_stationary.csv

Required columns:
- Name
- Start
- End
- Strand
- sRNA_length
- Best_avg_coverage
- Avg_coverage:* sample columns

---

## Filtering Criteria
1. Length filter:
   retain loci between 50–500 nt

2. Phase-specific coverage threshold:
   minimum coverage = 25th percentile of Best_avg_coverage
   locus must meet this threshold in ≥2 samples

These criteria define a high-confidence descriptive set.
Differential expression filtering is performed later using filterByExpr.

---

## Usage
Run from project root:

python3 scripts/02_filter/filter_annogesic_srna.py \
  --in_tsv results/sRNA_detection/01_annogesic/srna_exponential.csv \
  --out_tsv results/sRNA_detection/02_filter/srna_exponential.filtered.tsv

python3 scripts/02_filter/filter_annogesic_srna.py \
  --in_tsv results/sRNA_detection/01_annogesic/srna_stationary.csv \
  --out_tsv results/sRNA_detection/02_filter/srna_stationary.filtered.tsv

---

## Outputs
Written to:
results/sRNA_detection/02_filter/

Files:
- srna_exponential.filtered.tsv
- srna_stationary.filtered.tsv

Each run prints a summary such as:
[filter] input=X kept=Y min_len=50 max_len=500 min_cov(q=0.25)=...

---


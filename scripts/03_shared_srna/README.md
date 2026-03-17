# Shared sRNA Loci Detection

This folder contains the workflow used to identify **shared sRNA loci**
between exponential and stationary growth phases using filtered
ANNOgesic predictions.

The analysis converts filtered ANNOgesic tables into SAF format and
detects shared loci based on reciprocal genomic overlap.

---

## Inputs

Filtered ANNOgesic sRNA tables:

results/sRNA_detection/02_filter/srna_exponential.filtered.tsv  
results/sRNA_detection/02_filter/srna_stationary.filtered.tsv

Reference genome contig used for counting:

CP007799.1

NOTE:
The contig name must match BAM headers when SAF files are used for
featureCounts.

---

## Step 1 — Convert filtered ANNOgesic tables to SAF format

Run from project root:

python scripts/03_shared_srna/make_saf_from_annogesic.py \
  --infile results/sRNA_detection/02_filter/srna_exponential.filtered.tsv \
  --outfile results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf \
  --phase_label exp \
  --chr_value CP007799.1

python scripts/03_shared_srna/make_saf_from_annogesic.py \
  --infile results/sRNA_detection/02_filter/srna_stationary.filtered.tsv \
  --outfile results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf \
  --phase_label sta \
  --chr_value CP007799.1

Outputs:

results/sRNA_detection/03_shared_filtered/
  srna_exp.filtered.saf
  srna_sta.filtered.saf

---

## Step 2 — Identify shared sRNA loci

Default overlap threshold (moderate stringency):

python scripts/03_shared_srna/shared_srna_from_saf.py \
  --exp_saf results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf \
  --sta_saf results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf \
  --min_recip 0.5 \
  --out_prefix results/sRNA_detection/03_shared_filtered/shared

Optional stricter overlap:

--min_recip 0.8

Outputs:

shared_loci_minrecip_0.5.csv  
shared_loci_minrecip_0.8.csv

---

## Interpretation

Shared loci represent genomic regions where filtered sRNA predictions
from both growth phases overlap above a defined reciprocal threshold.

0.5 overlap captures broader shared regions.  
0.8 overlap highlights loci with stable boundaries across phases.

---


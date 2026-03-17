# 01_pipeline — RNA-seq preprocessing and alignment

This folder contains the core preprocessing steps used to generate alignments (BAM) and strand-specific coverage tracks (WIG) for downstream ANNOgesic analysis.

## Inputs
- Sample sheet:
  - `meta/samples_nissle_GSE147272.tsv` (columns: `sample_id`, `SRR`, `growth_phase`, `gravity`)
- Reference files (expected in `ref/`):
  - `EcN_CP007799.1.fa` (reference FASTA used for Bowtie2 alignment)
  - Bowtie2 index files: `EcN_CP007799.1*.bt2`
  - `EcN_ASM71459v1.gtf` (used for strandedness inference)

## Scripts (run order)

bash scripts/01_pipeline/01_download_fastq.sh

bash scripts/01_pipeline/02_trim_fastq.sh

bash scripts/01_pipeline/03_align_bowtie2.sh

bash scripts/01_pipeline/04_infer_strandedness.sh


NOTE: Raw FASTQ, trimmed FASTQ, BAM, WIG and other large intermediate files are not included in this repository.

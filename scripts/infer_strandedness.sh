#!/usr/bin/env bash
# Diagnostic: Infer RNA-seq library strandedness
#
# Tool  : RSeQC infer_experiment.py
# Input : aligned BAM file(s) + genome annotation (GTF/BED)
# Output: strandedness inference summary (stdout)
#
# Notes :
#   - Used to determine whether the library is unstranded, FR, or RF
#   - In this study, libraries were inferred as strand-specific, reverse-stranded (RF)
#     with ~97.7% consistency
#   - This step justifies strand handling in bam_to_wig_RF.sh

# Example command used:
infer_experiment.py \
  -i alignment/bam/SRR11357510.sorted.bam \
  -r meta/annotation.gtf

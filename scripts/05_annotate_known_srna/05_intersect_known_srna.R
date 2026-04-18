#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Intersect candidate sRNAs with BLAST-annotated known sRNAs
## - Candidate tables: produced by 04_srna_de (voom/limma)
## - Known tables: produced by BLAST filtering (04_filter_annotate_escherichia.py)
## Output: phase-specific tables of known candidate sRNAs
## ------------------------------------------------------------

options(stringsAsFactors = FALSE)

## ---------- INPUTS ----------
exp_candidates_file <- "results/sRNA_detection/04_srna_de_filtered/voom_exp_s2_candidates.csv"
sta_candidates_file <- "results/sRNA_detection/04_srna_de_filtered/voom_sta_s2_candidates.csv"

exp_known_file <- "results/sRNA_detection/05_known_srna/known_sRNAs_exponential_Escherichia.csv"
sta_known_file <- "results/sRNA_detection/05_known_srna/known_sRNAs_stationary_Escherichia.csv"

## ---------- OUTPUTS ----------
out_exp <- "results/sRNA_detection/05_known_srna/exp_known_candidate_sRNAs.csv"
out_sta <- "results/sRNA_detection/05_known_srna/sta_known_candidate_sRNAs.csv"

## ---------- HELPERS ----------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

normalize_candidate_id <- function(x) {
  # exp__sRNA_00044  -> sRNA_00044
  # sta__sRNA_00044  -> sRNA_00044
  # exp_sRNA_00044   -> sRNA_00044
  x <- trimws(as.character(x))
  x <- sub("^(exp|sta)(__|_)+", "", x)
  x
}

normalize_qseqid <- function(x) {
  # sRNA_00044::CP007799.1:14067-14327(-) -> sRNA_00044
  x <- trimws(as.character(x))
  x <- sub("::.*$", "", x)
  x
}

## ---------- CHECK FILES ----------
for (f in c(exp_candidates_file, sta_candidates_file, exp_known_file, sta_known_file)) {
  stop_if_missing(f)
}

## ---------- READ DATA ----------
exp_candidates <- read.csv(exp_candidates_file, check.names = FALSE)
sta_candidates <- read.csv(sta_candidates_file, check.names = FALSE)

exp_known <- read.csv(exp_known_file, check.names = FALSE)
sta_known <- read.csv(sta_known_file, check.names = FALSE)

## ---------- SANITY CHECKS ----------
if (!("Geneid" %in% colnames(exp_candidates))) stop("EXP candidate table missing column 'Geneid'")
if (!("Geneid" %in% colnames(sta_candidates))) stop("STA candidate table missing column 'Geneid'")

if (!("qseqid" %in% colnames(exp_known))) stop("EXP known table missing column 'qseqid'")
if (!("qseqid" %in% colnames(sta_known))) stop("STA known table missing column 'qseqid'")

## ---------- NORMALIZE IDS ----------
exp_candidates$srna_id <- normalize_candidate_id(exp_candidates$Geneid)
sta_candidates$srna_id <- normalize_candidate_id(sta_candidates$Geneid)

exp_known$srna_id <- normalize_qseqid(exp_known$qseqid)
sta_known$srna_id <- normalize_qseqid(sta_known$qseqid)

## ---------- DIAGNOSTICS ----------
cat("EXP candidate sRNAs:", nrow(exp_candidates), "\n")
cat("STA candidate sRNAs:", nrow(sta_candidates), "\n")
cat("EXP known hits (rows):", nrow(exp_known), "\n")
cat("STA known hits (rows):", nrow(sta_known), "\n")

cat("Unique EXP candidate srna_id:", length(unique(exp_candidates$srna_id)), "\n")
cat("Unique STA candidate srna_id:", length(unique(sta_candidates$srna_id)), "\n")
cat("Unique EXP known srna_id:", length(unique(exp_known$srna_id)), "\n")
cat("Unique STA known srna_id:", length(unique(sta_known$srna_id)), "\n")

## ---------- INTERSECT + MERGE ----------
exp_hits <- intersect(unique(exp_candidates$srna_id), unique(exp_known$srna_id))
sta_hits <- intersect(unique(sta_candidates$srna_id), unique(sta_known$srna_id))

cat("Known ∩ candidate (exponential):", length(exp_hits), "\n")
cat("Known ∩ candidate (stationary): ", length(sta_hits), "\n")

exp_known_candidates <- merge(exp_candidates, exp_known, by = "srna_id")
sta_known_candidates <- merge(sta_candidates, sta_known, by = "srna_id")

## Prefer best hit per srna_id (if duplicates exist)
## Sort by higher bitscore, then lower evalue
if (nrow(exp_known_candidates) > 0 && all(c("bitscore","evalue") %in% colnames(exp_known_candidates))) {
  exp_known_candidates <- exp_known_candidates[order(-exp_known_candidates$bitscore, exp_known_candidates$evalue), ]
  exp_known_candidates <- exp_known_candidates[!duplicated(exp_known_candidates$srna_id), ]
}
if (nrow(sta_known_candidates) > 0 && all(c("bitscore","evalue") %in% colnames(sta_known_candidates))) {
  sta_known_candidates <- sta_known_candidates[order(-sta_known_candidates$bitscore, sta_known_candidates$evalue), ]
  sta_known_candidates <- sta_known_candidates[!duplicated(sta_known_candidates$srna_id), ]
}

## Order by raw P-value if present
if ("P.Value" %in% colnames(exp_known_candidates)) {
  exp_known_candidates <- exp_known_candidates[order(exp_known_candidates$P.Value), ]
}
if ("P.Value" %in% colnames(sta_known_candidates)) {
  sta_known_candidates <- sta_known_candidates[order(sta_known_candidates$P.Value), ]
}

## ---------- WRITE OUTPUT ----------
write.csv(exp_known_candidates, out_exp, row.names = FALSE)
write.csv(sta_known_candidates, out_sta, row.names = FALSE)

cat("\nWrote:\n  ", out_exp, "\n  ", out_sta, "\n", sep="")

## ---------- REPORT ----------
if (nrow(exp_known_candidates) > 0) {
  cat("\nExponential: known candidate sRNAs (top rows):\n")
  keep_cols <- intersect(
    c("srna_id","known_sRNA_name","organism","pident","length","evalue","bitscore",
      "logFC","P.Value","direction"),
    colnames(exp_known_candidates)
  )
  print(head(exp_known_candidates[, keep_cols, drop = FALSE], 20))
} else {
  cat("\nExponential: no overlap between known and candidate sRNAs.\n")
}

if (nrow(sta_known_candidates) > 0) {
  cat("\nStationary: known candidate sRNAs (top rows):\n")
  keep_cols <- intersect(
    c("srna_id","known_sRNA_name","organism","pident","length","evalue","bitscore",
      "logFC","P.Value","direction"),
    colnames(sta_known_candidates)
  )
  print(head(sta_known_candidates[, keep_cols, drop = FALSE], 20))
} else {
  cat("\nStationary: no overlap between known and candidate sRNAs.\n")
}

cat("\nDone.\n")
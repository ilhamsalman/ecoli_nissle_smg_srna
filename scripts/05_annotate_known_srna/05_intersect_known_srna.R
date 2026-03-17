#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Intersect DE-significant sRNAs with BLAST-annotated known sRNAs
## - DE tables: produced by 04_srna_de (voom/limma)
## - Known tables: produced by BLAST filtering (04_filter_annotate_escherichia.py)
## Output: phase-specific tables of known + significant sRNAs
## ------------------------------------------------------------

options(stringsAsFactors = FALSE)

## ---------- INPUTS ----------
exp_sig_file <- "results/sRNA_detection/04_srna_de_filtered/voom_exp_s1_significant.csv"
sta_sig_file <- "results/sRNA_detection/04_srna_de_filtered/voom_sta_s1_significant.csv"

exp_known_file <- "results/sRNA_detection/05_known_srna/known_sRNAs_exponential_Escherichia.csv"
sta_known_file <- "results/sRNA_detection/05_known_srna/known_sRNAs_stationary_Escherichia.csv"

## ---------- OUTPUTS ----------
out_exp <- "results/sRNA_detection/05_known_srna/exp_known_SMG_responsive_sRNAs.csv"
out_sta <- "results/sRNA_detection/05_known_srna/sta_known_SMG_responsive_sRNAs.csv"

## ---------- HELPERS ----------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

normalize_de_id <- function(x) {
  # exp__sRNA_00044  -> sRNA_00044
  # sta__sRNA_00044  -> sRNA_00044
  # exp_sRNA_00044   -> sRNA_00044  (handles older variants)
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
for (f in c(exp_sig_file, sta_sig_file, exp_known_file, sta_known_file)) stop_if_missing(f)

## ---------- READ DATA ----------
exp_sig <- read.csv(exp_sig_file, check.names = FALSE)
sta_sig <- read.csv(sta_sig_file, check.names = FALSE)

exp_known <- read.csv(exp_known_file, check.names = FALSE)
sta_known <- read.csv(sta_known_file, check.names = FALSE)

## ---------- SANITY CHECKS ----------
if (!("Geneid" %in% colnames(exp_sig))) stop("EXP sig table missing column 'Geneid'")
if (!("Geneid" %in% colnames(sta_sig))) stop("STA sig table missing column 'Geneid'")

if (!("qseqid" %in% colnames(exp_known))) stop("EXP known table missing column 'qseqid'")
if (!("qseqid" %in% colnames(sta_known))) stop("STA known table missing column 'qseqid'")

## ---------- NORMALIZE IDS ----------
exp_sig$srna_id <- normalize_de_id(exp_sig$Geneid)
sta_sig$srna_id <- normalize_de_id(sta_sig$Geneid)

exp_known$srna_id <- normalize_qseqid(exp_known$qseqid)
sta_known$srna_id <- normalize_qseqid(sta_known$qseqid)

## Quick diagnostics (prints to console)
cat("EXP significant sRNAs:", nrow(exp_sig), "\n")
cat("STA significant sRNAs:", nrow(sta_sig), "\n")
cat("EXP known hits (rows):", nrow(exp_known), "\n")
cat("STA known hits (rows):", nrow(sta_known), "\n")

cat("Unique EXP sig srna_id:", length(unique(exp_sig$srna_id)), "\n")
cat("Unique STA sig srna_id:", length(unique(sta_sig$srna_id)), "\n")
cat("Unique EXP known srna_id:", length(unique(exp_known$srna_id)), "\n")
cat("Unique STA known srna_id:", length(unique(sta_known$srna_id)), "\n")

## ---------- INTERSECT + MERGE ----------
exp_hits <- intersect(unique(exp_sig$srna_id), unique(exp_known$srna_id))
sta_hits <- intersect(unique(sta_sig$srna_id), unique(sta_known$srna_id))

cat("Known ∩ significant (exponential):", length(exp_hits), "\n")
cat("Known ∩ significant (stationary): ", length(sta_hits), "\n")

# Merge (keep only those that intersect)
exp_known_sig <- merge(exp_sig, exp_known, by = "srna_id")
sta_known_sig <- merge(sta_sig, sta_known, by = "srna_id")

# Prefer best hit per srna_id (if duplicates exist)
# Sort: higher bitscore, lower evalue
if (nrow(exp_known_sig) > 0 && all(c("bitscore","evalue") %in% colnames(exp_known_sig))) {
  exp_known_sig <- exp_known_sig[order(-exp_known_sig$bitscore, exp_known_sig$evalue), ]
  exp_known_sig <- exp_known_sig[!duplicated(exp_known_sig$srna_id), ]
}
if (nrow(sta_known_sig) > 0 && all(c("bitscore","evalue") %in% colnames(sta_known_sig))) {
  sta_known_sig <- sta_known_sig[order(-sta_known_sig$bitscore, sta_known_sig$evalue), ]
  sta_known_sig <- sta_known_sig[!duplicated(sta_known_sig$srna_id), ]
}

# Order by adjusted p-value if present
if ("adj.P.Val" %in% colnames(exp_known_sig)) exp_known_sig <- exp_known_sig[order(exp_known_sig$adj.P.Val), ]
if ("adj.P.Val" %in% colnames(sta_known_sig)) sta_known_sig <- sta_known_sig[order(sta_known_sig$adj.P.Val), ]

## ---------- WRITE OUTPUT ----------
write.csv(exp_known_sig, out_exp, row.names = FALSE)
write.csv(sta_known_sig, out_sta, row.names = FALSE)

cat("\nWrote:\n  ", out_exp, "\n  ", out_sta, "\n", sep="")

## ---------- REPORT (compact) ----------
if (nrow(exp_known_sig) > 0) {
  cat("\nExponential: known & significant sRNAs (top rows):\n")
  keep_cols <- intersect(c("srna_id","known_sRNA_name","organism","pident","length","evalue","bitscore",
                           "logFC","adj.P.Val","direction"),
                         colnames(exp_known_sig))
  print(head(exp_known_sig[, keep_cols, drop=FALSE], 20))
} else {
  cat("\nExponential: no overlap between known and significant sRNAs.\n")
}

if (nrow(sta_known_sig) > 0) {
  cat("\nStationary: known & significant sRNAs (top rows):\n")
  keep_cols <- intersect(c("srna_id","known_sRNA_name","organism","pident","length","evalue","bitscore",
                           "logFC","adj.P.Val","direction"),
                         colnames(sta_known_sig))
  print(head(sta_known_sig[, keep_cols, drop=FALSE], 20))
} else {
  cat("\nStationary: no overlap between known and significant sRNAs.\n")
}

cat("\nDone.\n")

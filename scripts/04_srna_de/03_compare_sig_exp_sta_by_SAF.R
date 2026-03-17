#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(
    "Usage:\n",
    "  Rscript 03_compare_sig_exp_sta_by_SAF.R <exp_sig> <sta_sig> <exp_saf> <sta_saf> <min_recip> <out_csv>\n"
  )
}

exp_sig_file <- args[1]
sta_sig_file <- args[2]
exp_saf_file <- args[3]
sta_saf_file <- args[4]
min_recip    <- as.numeric(args[5])
out_csv      <- args[6]

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}
for (f in c(exp_sig_file, sta_sig_file, exp_saf_file, sta_saf_file)) stop_if_missing(f)

# Load inputs
exp_sig <- read_csv(exp_sig_file, show_col_types = FALSE)
sta_sig <- read_csv(sta_sig_file, show_col_types = FALSE)

exp_saf <- read_tsv(exp_saf_file, show_col_types = FALSE)
sta_saf <- read_tsv(sta_saf_file, show_col_types = FALSE)

# Validate columns
if (!("Geneid" %in% names(exp_sig))) stop("EXP sig table missing column: Geneid")
if (!("Geneid" %in% names(sta_sig))) stop("STA sig table missing column: Geneid")

req_saf_cols <- c("GeneID","Chr","Start","End","Strand")
if (!all(req_saf_cols %in% names(exp_saf))) stop("EXP SAF missing cols: ", paste(setdiff(req_saf_cols, names(exp_saf)), collapse=", "))
if (!all(req_saf_cols %in% names(sta_saf))) stop("STA SAF missing cols: ", paste(setdiff(req_saf_cols, names(sta_saf)), collapse=", "))

# Attach SAF coords to significant tables
exp_saf2 <- exp_saf %>%
  rename(Geneid = GeneID) %>%
  mutate(Start = as.integer(Start), End = as.integer(End), Strand = as.character(Strand), Chr = as.character(Chr))

sta_saf2 <- sta_saf %>%
  rename(Geneid = GeneID) %>%
  mutate(Start = as.integer(Start), End = as.integer(End), Strand = as.character(Strand), Chr = as.character(Chr))

exp_sig2 <- exp_sig %>%
  left_join(exp_saf2 %>% select(Geneid, Chr, Start, End, Strand), by = "Geneid")

sta_sig2 <- sta_sig %>%
  left_join(sta_saf2 %>% select(Geneid, Chr, Start, End, Strand), by = "Geneid")

exp_ok <- exp_sig2 %>% filter(!is.na(Chr), !is.na(Start), !is.na(End), !is.na(Strand))
sta_ok <- sta_sig2 %>% filter(!is.na(Chr), !is.na(Start), !is.na(End), !is.na(Strand))

cat("EXP significant rows:", nrow(exp_sig2), "\n")
cat("STA significant rows:", nrow(sta_sig2), "\n")
cat("EXP missing coords:", sum(is.na(exp_sig2$Chr)), "\n")
cat("STA missing coords:", sum(is.na(sta_sig2$Chr)), "\n")

# Reciprocal-overlap pairing (Chr+Strand, overlap, recip>=min_recip)
pairs <- exp_ok %>%
  inner_join(sta_ok, by = c("Chr","Strand"), suffix = c("_exp","_sta"), relationship = "many-to-many") %>%
  filter(.data$Start_exp <= .data$End_sta, .data$Start_sta <= .data$End_exp) %>%
  mutate(
    overlap_start = pmax(.data$Start_exp, .data$Start_sta),
    overlap_end   = pmin(.data$End_exp,   .data$End_sta),
    overlap_len   = .data$overlap_end - .data$overlap_start + 1,
    len_exp       = .data$End_exp - .data$Start_exp + 1,
    len_sta       = .data$End_sta - .data$Start_sta + 1,
    recip_exp     = .data$overlap_len / .data$len_exp,
    recip_sta     = .data$overlap_len / .data$len_sta
  ) %>%
  filter(.data$overlap_len > 0) %>%
  filter(.data$recip_exp >= min_recip, .data$recip_sta >= min_recip) %>%
  arrange(desc(.data$overlap_len))

cat("Reciprocal-overlap pairs (min_recip=", min_recip, "): ", nrow(pairs), "\n", sep="")

# Write ONLY the overlap table (single output)
write_csv(pairs, out_csv)
cat("Wrote:\n  ", out_csv, "\n", sep="")

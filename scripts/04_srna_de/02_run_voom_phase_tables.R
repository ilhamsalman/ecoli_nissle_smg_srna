#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage:\n",
    "  Rscript run_voom_phase_tables.R <counts_file> <annogesic_table> <phase_prefix> <phase_label> <out_prefix>\n\n",
    "Example:\n",
    "  Rscript run_voom_phase_tables.R \\\n",
    "    results/sRNA_detection/03_srna_de/counts_exp_s1.txt \\\n",
    "    results/sRNA_detection/01_annogesic/srna_exponential.csv \\\n",
    "    exp__ exponential \\\n",
    "    results/sRNA_detection/03_srna_de/voom_exp_s1\n"
  )
}

counts_file  <- args[1]
anno_file    <- args[2]
phase_prefix <- args[3]   # exp__ or sta__
phase_label  <- args[4]   # exponential / stationary
out_prefix   <- args[5]   # e.g. results/.../voom_exp_s1  (no extension)

out_all_csv <- paste0(out_prefix, "_all.csv")
out_sig_csv <- paste0(out_prefix, "_significant.csv")

# ----------------------------
# helpers
# ----------------------------
normalize_type <- function(x) {
  x <- tolower(as.character(x))
  x[is.na(x) | x == ""] <- "other"
  ifelse(grepl("intergenic", x), "Intergenic",
         ifelse(grepl("antisense", x), "Antisense", "Other"))
}

# ----------------------------
# 1) Read featureCounts output
# ----------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

# featureCounts standard columns: Geneid, Chr, Start, End, Strand, Length, <samples...>
if (!("Geneid" %in% colnames(fc))) stop("counts file missing Geneid column: ", counts_file)

sample_cols <- colnames(fc)[7:ncol(fc)]
cat("Sample columns:\n")
print(data.frame(column = sample_cols))

# Parse MG/NG from filenames (robust to full paths)
sample_names <- basename(sample_cols)
group <- ifelse(grepl("^MG", sample_names), "MG",
                ifelse(grepl("^NG", sample_names), "NG", NA))
if (any(is.na(group))) {
  stop("Could not infer MG/NG for these sample columns:\n",
       paste(sample_cols[is.na(group)], collapse = "\n"))
}
group <- factor(group, levels = c("MG", "NG"))

cat("\nInferred groups:\n")
print(data.frame(sample = sample_names, group = group))
cat("\nGroup counts:\n")
print(table(group))

cts <- fc[, 7:ncol(fc)]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid
colnames(cts) <- sample_names

design <- model.matrix(~ group)

cat("\nDone.\n")
cat("Total candidates (from counts): ", nrow(cts), "\n", sep="")

# ----------------------------
# 2) edgeR + filterByExpr
# ----------------------------
y <- DGEList(counts = cts, group = group)

keep <- filterByExpr(y, design)
cat("Passed filterByExpr (tested): ", sum(keep), "\n", sep="")

y_f <- y[keep, , keep.lib.sizes = FALSE]
y_f <- calcNormFactors(y_f)

# ----------------------------
# 3) voom-limma DE (only tested features)
# ----------------------------
v <- voom(y_f, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# coef "groupNG" is NG - MG ; flip to MG - NG (SMG effect)
tab <- topTable(fit, coef = "groupNG", number = Inf, sort.by = "P")
tab$logFC <- -tab$logFC
tab$Geneid <- rownames(tab)

# ----------------------------
# 4) Build full table (including filtered-out rows)
# ----------------------------
all_df <- data.frame(
  Geneid = rownames(cts),
  phase = phase_label,
  passed_filter = keep,
  stringsAsFactors = FALSE
)

tab$core_p <- tab$P.Value
tab$core_fdr <- tab$adj.P.Val

tab2 <- tab[, c("Geneid","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
all_df <- merge(all_df, tab2, by = "Geneid", all.x = TRUE)

# Direction (SMG effect) — you said “all FDR-significant loci”, so |logFC| cutoff = 0.
all_df$direction <- "Not significant"
all_df$direction[!is.na(all_df$adj.P.Val) & all_df$adj.P.Val <= 0.05 & all_df$logFC > 0] <- "Up (SMG)"
all_df$direction[!is.na(all_df$adj.P.Val) & all_df$adj.P.Val <= 0.05 & all_df$logFC < 0] <- "Down (SMG)"

# For filtered-out rows, make it explicit
all_df$direction[all_df$passed_filter == FALSE] <- "Not tested (filtered out)"

# ----------------------------
# 5) Read ANNOgesic table + add sRNA_types
# ----------------------------
anno <- tryCatch(
  read.delim(anno_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) read.csv(anno_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
)

if (!("Name" %in% colnames(anno))) {
  stop("ANNOgesic table must contain column 'Name'. Found: ",
       paste(colnames(anno), collapse = ", "))
}

type_col <- NULL
if ("sRNA_types" %in% colnames(anno)) type_col <- "sRNA_types"
if (is.null(type_col) && "sRNA_type" %in% colnames(anno)) type_col <- "sRNA_type"
if (is.null(type_col)) {
  anno$sRNA_types <- NA
  type_col <- "sRNA_types"
}

anno$Geneid <- paste0(phase_prefix, trimws(anno$Name))
anno_keep <- unique(anno[, c("Geneid", type_col)])
colnames(anno_keep)[colnames(anno_keep) == type_col] <- "sRNA_types"

all_df <- merge(all_df, anno_keep, by = "Geneid", all.x = TRUE)
all_df$sRNA_types <- normalize_type(all_df$sRNA_types)

# Order columns nicely
all_df <- all_df[, c("Geneid","phase","sRNA_types","passed_filter",
                     "logFC","AveExpr","t","P.Value","adj.P.Val","B","direction")]

# ----------------------------
# 6) Write outputs
# ----------------------------
write.csv(all_df, out_all_csv, row.names = FALSE)

sig <- all_df[!is.na(all_df$adj.P.Val) & all_df$adj.P.Val <= 0.05 & all_df$passed_filter == TRUE, ]
write.csv(sig, out_sig_csv, row.names = FALSE)

cat("Significant (FDR<=0.05, |logFC|>=0): ", nrow(sig), "\n", sep="")
cat("Wrote:\n  ", out_all_csv, "\n  ", out_sig_csv, "\n", sep="")

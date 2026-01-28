suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript run_voom_phase.R <counts_file> <annogesic_table> <phase_prefix> <phase_label> <out_all_csv> <out_sig_csv> <out_volcano_png>\n",
       "Example:\n",
       "  Rscript run_voom_phase.R counts_exp_frac.txt srna_exponential.csv exp__ exponential out_all.csv out_sig.csv volcano.png")
}

counts_file   <- args[1]
anno_file     <- args[2]
phase_prefix  <- args[3]   # e.g. "exp__" or "sta__" (must match Geneid prefix in counts)
phase_label   <- args[4]   # e.g. "exponential" or "stationary" (for a column in output)
out_all_csv   <- args[5]
out_sig_csv   <- args[6]
out_volcano   <- args[7]

# ----------------------------
# 1) Read featureCounts output
# ----------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

cts <- fc[, 7:ncol(fc)]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid
colnames(cts) <- colnames(fc)[7:ncol(fc)]

# Expect 4 samples in order: MG, MG, NG, NG
group <- factor(c("MG", "MG", "NG", "NG"))
design <- model.matrix(~ group)

# ----------------------------
# 2) edgeR object + filtering
# ----------------------------
y <- DGEList(counts = cts, group = group)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# ----------------------------
# 3) voom-limma DE
# ----------------------------
v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# limma coef "groupNG" is NG - MG
tab <- topTable(fit, coef = "groupNG", number = Inf, sort.by = "P")

# Convert to SMG effect = MG - NG (flip sign)
tab$logFC <- -tab$logFC

tab$Geneid <- rownames(tab)
tab$phase <- phase_label

tab$direction <- ifelse(tab$adj.P.Val < 0.05 & abs(tab$logFC) >= 1,
                        ifelse(tab$logFC > 0, "Up (SMG)", "Down (SMG)"),
                        "NotSig")

# ----------------------------
# 4) Read ANNOgesic table and merge sRNA type
# ----------------------------
# Your ANNOgesic files are tab-delimited but named .csv
anno <- tryCatch(
  read.delim(anno_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) read.csv(anno_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
)

# We need at least Name and sRNA_types. If sRNA_types missing, we'll still run.
if (!("Name" %in% colnames(anno))) {
  stop("ANNOgesic table must contain a column named 'Name'. Columns found: ",
       paste(colnames(anno), collapse = ", "))
}

anno$Geneid <- paste0(phase_prefix, anno$Name)

# Pick type column if present
type_col <- NULL
if ("sRNA_types" %in% colnames(anno)) type_col <- "sRNA_types"
if (is.null(type_col) && "sRNA_type" %in% colnames(anno)) type_col <- "sRNA_type"
if (is.null(type_col)) {
  # create empty column so downstream is consistent
  anno$sRNA_types <- NA
  type_col <- "sRNA_types"
}

anno_keep <- unique(anno[, c("Geneid", type_col)])

tab2 <- merge(tab, anno_keep, by = "Geneid", all.x = TRUE)
colnames(tab2)[colnames(tab2) == type_col] <- "sRNA_types"

# Order columns nicely
tab2 <- tab2[, c("Geneid", "phase", "sRNA_types",
                 "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "direction")]

# ----------------------------
# 5) Write outputs
# ----------------------------
write.csv(tab2, out_all_csv, row.names = FALSE)

sig <- tab2[tab2$adj.P.Val < 0.05 & abs(tab2$logFC) >= 1, ]
write.csv(sig, out_sig_csv, row.names = FALSE)

cat("Done.\n", sep = "")
cat("Significant (FDR<0.05 & |logFC|>=1): ", nrow(sig), "\n", sep = "")

# ----------------------------
# 6) Volcano plot (SMG: MG - NG)
# ----------------------------
p <- ggplot(tab2, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.7, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = paste0("Volcano: SMG effect (MG - NG) - ", phase_label),
    x = "log2FC (MG - NG)",
    y = "-log10(P-value)"
  ) +
  theme_bw()

ggsave(out_volcano, p, width = 7, height = 6, dpi = 300)

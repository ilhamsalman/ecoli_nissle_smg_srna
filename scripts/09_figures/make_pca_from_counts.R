#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage:\n",
    "  Rscript scripts/5.figures/make_pca_from_counts.R <counts_file> <out_prefix> <phase_label> [top_var]\n\n",
    "Arguments:\n",
    "  counts_file : featureCounts output (counts_*.txt)\n",
    "  out_prefix  : prefix for outputs (no extension)\n",
    "  phase_label : e.g. exponential / stationary (used only for filenames/log)\n",
    "  top_var     : optional, number of most-variable features to use (default 1000; use 0 for all)\n\n",
    "Example:\n",
    "  Rscript scripts/5.figures/make_pca_from_counts.R \\\n",
    "    results/sRNA_detection/03_srna_de/counts_exp_s1.txt \\\n",
    "    results/figures/pca_exp_s1 \\\n",
    "    exponential \\\n",
    "    1000\n"
  )
}

counts_file <- args[1]
out_prefix  <- args[2]
phase_label <- args[3]
top_var     <- if (length(args) >= 4) as.integer(args[4]) else 1000

if (!file.exists(counts_file)) stop("Missing counts file: ", counts_file)
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Read featureCounts matrix
# ----------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
if (!("Geneid" %in% colnames(fc))) stop("counts file missing Geneid column: ", counts_file)

sample_cols  <- colnames(fc)[7:ncol(fc)]
sample_names <- basename(sample_cols)

cts <- fc[, 7:ncol(fc), drop = FALSE]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid
colnames(cts) <- sample_names

# infer MG/NG from filenames
group <- ifelse(grepl("^MG", sample_names), "MG",
         ifelse(grepl("^NG", sample_names), "NG", NA))
if (any(is.na(group))) {
  stop("Could not infer MG/NG for these sample columns:\n",
       paste(sample_names[is.na(group)], collapse = "\n"))
}
group <- factor(group, levels = c("MG", "NG"))

# ----------------------------
# 2) Filter low-expression (for stable PCA)
# ----------------------------
design <- model.matrix(~ group)
y <- DGEList(counts = cts, group = group)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# logCPM
logcpm <- cpm(y, log = TRUE, prior.count = 2)

# ----------------------------
# 3) Choose features for PCA
# ----------------------------
if (!is.na(top_var) && top_var > 0 && nrow(logcpm) > top_var) {
  rv <- apply(logcpm, 1, var)
  idx <- order(rv, decreasing = TRUE)[seq_len(top_var)]
  mat <- logcpm[idx, , drop = FALSE]
} else {
  mat <- logcpm
}

# PCA on samples: transpose to samples x features
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", pct[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct[2])

df <- data.frame(
  sample = rownames(pca$x),
  group  = group[match(rownames(pca$x), sample_names)],
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  stringsAsFactors = FALSE
)

# ----------------------------
# 4) Plot (NO text labels; legend only)
# ----------------------------
p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 6, alpha = 0.95) +
  labs(x = pc1_lab, y = pc2_lab, color = "Group") +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.15),
    legend.position = "right"
  )

# Make it wide so PC1 visually stretches vs PC2
png(paste0(out_prefix, ".png"), width = 2400, height = 1400, res = 200)
print(p)
dev.off()

pdf(paste0(out_prefix, ".pdf"), width = 10, height = 5.5)
print(p)
dev.off()

cat("OK\n")
cat("Counts:  ", counts_file, "\n", sep = "")
cat("Features used for PCA: ", nrow(mat), "\n", sep = "")
cat("Wrote:\n  ", out_prefix, ".png\n  ", out_prefix, ".pdf\n", sep = "")


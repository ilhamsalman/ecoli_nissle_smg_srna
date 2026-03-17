#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage:\n",
    "  Rscript scripts/5.figures/make_heatmap_sig_from_counts.R <counts_file> <voom_all_csv> <out_prefix> <topN>\n\n",
    "Example:\n",
    "  Rscript scripts/5.figures/make_heatmap_sig_from_counts.R \\\n",
    "    results/sRNA_detection/03_srna_de/counts_exp_s1.txt \\\n",
    "    results/sRNA_detection/03_srna_de/voom_exp_s1_all.csv \\\n",
    "    results/figures/heatmap_exp_s1_sig_top40 \\\n",
    "    40\n"
  )
}

counts_file <- args[1]
voom_all_csv <- args[2]
out_prefix <- args[3]
topN <- as.integer(args[4])

if (!file.exists(counts_file)) stop("Missing counts file: ", counts_file)
if (!file.exists(voom_all_csv)) stop("Missing voom table: ", voom_all_csv)

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Read counts (featureCounts)
# ----------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
if (!("Geneid" %in% colnames(fc))) stop("counts file missing Geneid column: ", counts_file)

cts <- fc[, 7:ncol(fc), drop = FALSE]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid

sample_cols <- colnames(fc)[7:ncol(fc)]
sample_names <- basename(sample_cols)
colnames(cts) <- sample_names

# infer MG/NG from filenames
group <- ifelse(grepl("^MG", sample_names), "MG",
                ifelse(grepl("^NG", sample_names), "NG", NA))
if (any(is.na(group))) {
  stop("Could not infer MG/NG for:\n", paste(sample_names[is.na(group)], collapse="\n"))
}
group <- factor(group, levels = c("MG","NG"))

# ----------------------------
# 2) Read voom_all table, pick significant only
# ----------------------------
voom <- read.csv(voom_all_csv, stringsAsFactors = FALSE, check.names = FALSE)

need <- c("Geneid","adj.P.Val","passed_filter")
missing <- setdiff(need, colnames(voom))
if (length(missing) > 0) stop("voom_all_csv missing cols: ", paste(missing, collapse=", "))

sig <- voom[!is.na(voom$adj.P.Val) & voom$adj.P.Val <= 0.05 & voom$passed_filter == TRUE, ]
if (nrow(sig) == 0) {
  stop("No significant loci found (FDR<=0.05) in: ", voom_all_csv,
       "\nHeatmap cannot be generated because the matrix would be empty.")
}

sig <- sig[order(sig$adj.P.Val, decreasing = FALSE), ]
if (!is.na(topN) && topN > 0 && nrow(sig) > topN) sig <- sig[1:topN, ]

sig_ids <- unique(sig$Geneid)
sig_ids <- sig_ids[sig_ids %in% rownames(cts)]
if (length(sig_ids) == 0) stop("Significant Geneids not found in counts matrix (Geneid mismatch).")

cts_sig <- cts[sig_ids, , drop = FALSE]

# ----------------------------
# 3) logCPM transform + row Z-score for plotting
# ----------------------------
y <- DGEList(counts = cts_sig, group = group)
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 2)

# Row Z-score (center/scale each sRNA across samples)
zmat <- t(scale(t(logcpm)))
zmat[is.na(zmat)] <- 0

ann_col <- data.frame(Group = group)
rownames(ann_col) <- colnames(zmat)

# ----------------------------
# 4) Colors: "previous" style (classic blue->white->red)
# ----------------------------
hm_cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

# Helper to save reliably
save_heatmap <- function(file, device = c("png","pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    # cairo tends to be more reliable on some Linux setups
    grDevices::png(file, width = 1800, height = 1600, res = 300, type = "cairo-png")
  } else {
    grDevices::pdf(file, width = 7.5, height = 7.5)
  }

  pheatmap(
    zmat,
    color = hm_cols,
    annotation_col = ann_col,
    fontsize = 9,
    fontsize_row = 6,
    main = "",
    border_color = NA
  )
  grDevices::dev.off()
}

save_heatmap(paste0(out_prefix, ".png"), "png")
save_heatmap(paste0(out_prefix, ".pdf"), "pdf")

cat("OK\n")
cat("Counts file: ", counts_file, "\n", sep="")
cat("Voom table:  ", voom_all_csv, "\n", sep="")
cat("Sig used:    ", length(sig_ids), " (after topN)\n", sep="")
cat("Wrote:\n  ", out_prefix, ".png\n  ", out_prefix, ".pdf\n", sep="")

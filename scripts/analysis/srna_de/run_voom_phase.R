suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript run_voom_phase.R <counts_file> <phase_label> <annogesic_csv> <out_all_csv> <out_sig_csv> <out_volcano_png>")
}

counts_file   <- args[1]
phase_label   <- args[2]   # "exp" or "sta" (or "exponential"/"stationary" if your IDs match)
annogesic_csv <- args[3]
out_all_csv   <- args[4]
out_sig_csv   <- args[5]
out_volcano   <- args[6]

# ---- read featureCounts output ----
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

cts <- fc[, 7:ncol(fc)]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid
colnames(cts) <- colnames(fc)[7:ncol(fc)]

# Expect 4 samples: MG, MG, NG, NG (in that exact order)
group <- factor(c("MG", "MG", "NG", "NG"))
design <- model.matrix(~ group)

# ---- voom-limma pipeline ----
y <- DGEList(counts = cts)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# coef groupNG is (NG - MG)
tab <- topTable(fit, coef = "groupNG", number = Inf, sort.by = "P")
tab$Geneid <- rownames(tab)

# ---- SMG effect (MG - NG) ----
# Convert (NG - MG) into (MG - NG) by flipping sign
tab$logFC_SMG <- -tab$logFC

tab$SMG_direction <- "NotSig"
tab$SMG_direction[tab$adj.P.Val < 0.05 & tab$logFC_SMG >=  1] <- "Up (SMG)"
tab$SMG_direction[tab$adj.P.Val < 0.05 & tab$logFC_SMG <= -1] <- "Down (SMG)"

# ---- attach ANNOgesic annotation: sRNA_types (antisense/intergenic/...) ----
anno <- read.delim(annogesic_csv, stringsAsFactors = FALSE, check.names = FALSE)

# FeatureCounts Geneid looks like exp__sRNA_XXXXX or sta__sRNA_XXXXX
# Make join key srna_id = sRNA_XXXXX
tab$srna_id <- sub("^(exp|sta)__", "", tab$Geneid)

# Ensure required cols exist
if (!("Name" %in% colnames(anno))) stop("ANNOgesic file must contain column: Name")
if (!("sRNA_types" %in% colnames(anno))) stop("ANNOgesic file must contain column: sRNA_types")

anno_sub <- unique(anno[, c("Name", "sRNA_types")])
colnames(anno_sub) <- c("srna_id", "sRNA_types")

tab2 <- merge(tab, anno_sub, by = "srna_id", all.x = TRUE, sort = FALSE)

# add phase label column (nice for combined downstream)
tab2$phase <- phase_label

# ---- write outputs ----
# full table
write.csv(tab2, out_all_csv, row.names = FALSE)

# significant table (FDR<0.05 and abs(logFC_SMG)>=1)
sig <- tab2[tab2$adj.P.Val < 0.05 & abs(tab2$logFC_SMG) >= 1, ]
write.csv(sig, out_sig_csv, row.names = FALSE)

cat("Done.\n", sep = "")
cat("Total tested:", nrow(tab2), "\n", sep = " ")
cat("Significant (FDR<0.05 & |logFC_SMG|>=1):", nrow(sig), "\n", sep = " ")

# ---- volcano plot (SMG effect) ----
p <- ggplot(tab2, aes(x = logFC_SMG, y = -log10(P.Value), color = SMG_direction)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = paste0("Volcano: SMG (MG - NG) effect — ", phase_label),
    x = "log2FC (MG - NG)",
    y = "-log10(p-value)"
  ) +
  theme_bw()

ggsave(out_volcano, p, width = 7, height = 6, dpi = 200)

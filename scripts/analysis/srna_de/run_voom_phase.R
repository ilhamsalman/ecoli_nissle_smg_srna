suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_voom_phase.R <counts_file> <phase_label> <out_all_csv> <out_sig_csv> <out_volcano_png>")
}

counts_file <- args[1]
phase_label <- args[2]
out_all_csv <- args[3]
out_sig_csv <- args[4]
out_volcano <- args[5]

# Read featureCounts output
fc <- read.delim(counts_file, comment.char="#", check.names=FALSE)

cts <- fc[, 7:ncol(fc)]
cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
cts <- as.matrix(cts)
rownames(cts) <- fc$Geneid
colnames(cts) <- colnames(fc)[7:ncol(fc)]

# Expect 4 samples in order: MG, MG, NG, NG
if (ncol(cts) != 4) stop("This script expects exactly 4 samples (MG, MG, NG, NG). Found: ", ncol(cts))

group <- factor(c("MG","MG","NG","NG"))
design <- model.matrix(~ group)

y <- DGEList(counts=cts)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

tab <- topTable(fit, coef="groupNG", number=Inf, sort.by="P")
tab$Geneid <- rownames(tab)
tab$phase <- phase_label

# direction labels (based on thresholds)
tab$direction <- "NotSig"
tab$direction[tab$adj.P.Val < 0.05 & tab$logFC >= 1]  <- "Up (NG)"
tab$direction[tab$adj.P.Val < 0.05 & tab$logFC <= -1] <- "Down (NG)"

# Write full table
tab_out <- tab[, c("Geneid","phase","logFC","AveExpr","t","P.Value","adj.P.Val","B","direction")]
write.csv(tab_out, out_all_csv, row.names=FALSE)

# Significant subset (FDR<0.05 & |logFC|>=1)
sig <- tab_out[tab_out$adj.P.Val < 0.05 & abs(tab_out$logFC) >= 1, ]
write.csv(sig, out_sig_csv, row.names=FALSE)

# Volcano plot
p <- ggplot(tab_out, aes(x=logFC, y=-log10(adj.P.Val), color=direction)) +
  geom_point(alpha=0.6, size=1.2) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(
    title=paste0("sRNA differential expression (", phase_label, "): NG vs MG"),
    x="log2 fold change (NG vs MG)",
    y="-log10(FDR)"
  ) +
  theme_bw()

ggsave(out_volcano, p, width=6, height=5)

cat("Done.\n")
cat("Significant (FDR<0.05 & |logFC|>=1): ", nrow(sig), "\n", sep="")

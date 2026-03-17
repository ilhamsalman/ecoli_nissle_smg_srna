suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage:\n",
    "  Rscript make_volcano_from_voom_csv.R <voom_all_csv> <out_prefix>\n\n",
    "Example:\n",
    "  Rscript make_volcano_from_voom_csv.R results/sRNA_detection/03_srna_de/voom_exp_s1_all.csv ",
    "results/figures/volcano_exp_s1\n"
  )
}

in_csv <- args[1]
out_prefix <- args[2]

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

df <- read_csv(in_csv, show_col_types = FALSE)

# keep only tested rows
df <- df %>% filter(!is.na(adj.P.Val), !is.na(logFC))

df <- df %>%
  mutate(volcano_group = case_when(
    adj.P.Val <= 0.05 & logFC > 0 ~ "Up (SMG)",
    adj.P.Val <= 0.05 & logFC < 0 ~ "Down (SMG)",
    TRUE ~ "Not significant"
  ))

p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = volcano_group)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(
      "Up (SMG)" = "#D73027",
      "Down (SMG)" = "#4575B4",
      "Not significant" = "grey70"
    ),
    name = "SMG response"
  ) +
  labs(x = "log2 FC (MG - NG)", y = "-log10(FDR)") +
  theme_bw() +
  theme(plot.title = element_blank(), legend.position = "right")

ggsave(paste0(out_prefix, ".png"), p, width = 7, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = 6)

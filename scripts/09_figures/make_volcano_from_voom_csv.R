suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage:\n",
    "  Rscript make_volcano_from_voom_csv.R <voom_all_csv> <out_prefix>\n"
  )
}

in_csv <- args[1]
out_prefix <- args[2]

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

df <- read_csv(in_csv, show_col_types = FALSE)

df <- df %>%
  filter(passed_filter == TRUE, !is.na(P.Value), !is.na(logFC)) %>%
  mutate(
    plot_p = pmax(P.Value, 1e-300),
    volcano_group = case_when(
      P.Value <= 0.05 & logFC > 0 ~ "Up",
      P.Value <= 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

p <- ggplot(df, aes(x = logFC, y = -log10(plot_p), color = volcano_group)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(
      "Down" = "#4575B4",
      "Not significant" = "grey70",
      "Up" = "#D73027"
    ),
    breaks = c("Down", "Not significant", "Up"),
    name = "SMG response"
  ) +
  labs(
    x = expression(log[2] * "FC (SMG - NG)"),
    y = expression(-log[10] * "(P-value)")
  ) +
  theme_bw() +
  theme(plot.title = element_blank(), legend.position = "right")

ggsave(paste0(out_prefix, ".png"), p, width = 7, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = 6)
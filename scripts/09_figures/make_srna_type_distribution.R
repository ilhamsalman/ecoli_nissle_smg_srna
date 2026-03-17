#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage:\n",
    "  Rscript make_srna_type_distribution.R <exp_table> <sta_table> <out_prefix>\n\n",
    "Example:\n",
    "  Rscript scripts/09_figures/make_srna_type_distribution.R \\\n",
    "    results/sRNA_detection/02_filter/srna_exponential.filtered.tsv \\\n",
    "    results/sRNA_detection/02_filter/srna_stationary.filtered.tsv \\\n",
    "    results/figures/predicted_srna_candidates_type_distribution_filtered\n"
  )
}

exp_file   <- args[1]
sta_file   <- args[2]
out_prefix <- args[3]

if (!file.exists(exp_file)) stop("Missing: ", exp_file)
if (!file.exists(sta_file)) stop("Missing: ", sta_file)
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# Robust reader: try TSV first, then CSV
read_table_auto <- function(path) {
  tryCatch(
    read_tsv(path, show_col_types = FALSE),
    error = function(e) read_csv(path, show_col_types = FALSE)
  )
}

exp <- read_table_auto(exp_file)
sta <- read_table_auto(sta_file)

if (!("sRNA_types" %in% names(exp))) stop("Missing column 'sRNA_types' in: ", exp_file)
if (!("sRNA_types" %in% names(sta))) stop("Missing column 'sRNA_types' in: ", sta_file)

normalize_type <- function(x) {
  x <- tolower(as.character(x))
  x[is.na(x) | x == ""] <- "Other"
  ifelse(grepl("intergenic", x), "Intergenic",
         ifelse(grepl("antisense", x), "Antisense", "Other"))
}

exp2 <- exp %>%
  mutate(Phase = "Exponential",
         Type  = normalize_type(.data$sRNA_types))

sta2 <- sta %>%
  mutate(Phase = "Stationary",
         Type  = normalize_type(.data$sRNA_types))

df <- bind_rows(exp2, sta2) %>%
  count(Phase, Type, name = "Count") %>%
  filter(Type %in% c("Intergenic", "Antisense")) %>%
  mutate(
    Type  = factor(Type, levels = c("Intergenic", "Antisense")),
    Phase = factor(Phase, levels = c("Stationary", "Exponential"))
  )

p <- ggplot(df, aes(x = Type, y = Count, fill = Phase)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65,
           color = "black",
           linewidth = 0.3) +
  geom_text(aes(label = format(Count, big.mark=",")),
            position = position_dodge(width = 0.75),
            vjust = -0.4,
            size = 4) +
  labs(x = NULL, y = "Number of predicted sRNA candidates", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top") +
  expand_limits(y = max(df$Count) * 1.10)

ggsave(paste0(out_prefix, ".pdf"), p, width = 5.2, height = 4.8)
ggsave(paste0(out_prefix, ".png"), p, width = 5.2, height = 4.8, dpi = 600)

# Print counts table to stdout for logging
print(df)

cat("Wrote:\n  ", out_prefix, ".pdf\n  ", out_prefix, ".png\n", sep = "")

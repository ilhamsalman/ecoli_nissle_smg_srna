suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# -----------------
# Input files
# -----------------
exp_file <- "srna_exponential.csv"
sta_file <- "srna_stationary.csv"

# -----------------
# Read ANNOgesic output (TAB-delimited)
# -----------------
exp <- read_tsv(exp_file, show_col_types = FALSE)
sta <- read_tsv(sta_file, show_col_types = FALSE)

# -----------------
# Normalize sRNA types
# -----------------
normalize_type <- function(x) {
  x <- tolower(x)
  ifelse(grepl("intergenic", x), "Intergenic",
         ifelse(grepl("antisense", x), "Antisense", "Other"))
}

exp2 <- exp %>%
  mutate(
    Phase = "Exponential",
    Type = ifelse(is.na(sRNA_types), "Other", normalize_type(sRNA_types))
  )

sta2 <- sta %>%
  mutate(
    Phase = "Stationary",
    Type = ifelse(is.na(sRNA_types), "Other", normalize_type(sRNA_types))
  )

# -----------------
# Count categories
# -----------------
df <- bind_rows(exp2, sta2) %>%
  count(Phase, Type, name = "Count") %>%
  filter(Type %in% c("Intergenic", "Antisense")) %>%
  mutate(
    Type = factor(Type, levels = c("Intergenic", "Antisense")),
    Phase = factor(Phase, levels = c("Stationary", "Exponential"))
  )

# -----------------
# Plot (NO TITLE)
# -----------------
p <- ggplot(df, aes(x = Type, y = Count, fill = Phase)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "Stationary" = "#E69F00",
      "Exponential" = "#0072B2"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of sRNAs",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

# -----------------
# Save
# -----------------
ggsave("srna_intergenic_antisense_distribution.pdf",
       p, width = 5.2, height = 4.8)
ggsave("srna_intergenic_antisense_distribution.png",
       p, width = 5.2, height = 4.8, dpi = 600)

# -----------------
# Print counts (for Results)
# -----------------
print(df)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
})

cat("Working directory:", getwd(), "\n")

if (!dir.exists("results")) {
  stop("Run this script from the project root directory: ecoli_nissle_smg_srna")
}

out_dir <- "results/figures/cis_elements"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Panel A: RNA thermometer classes
# -----------------------------
thermo <- tribble(
  ~Class,            ~Exponential, ~Stationary,
  "Phe_leader",      8,            26,
  "hsp17",           2,            4,
  "ROSE_2",          1,            4,
  "FourU",           1,            3,
  "cspA",            1,            2,
  "Lambda_thermo",   0,            2,
  "ROSE_3",          1,            1,
  "PrfA",            0,            1
)

thermo_order <- thermo %>%
  arrange(Stationary, Exponential) %>%
  pull(Class)

thermo_long <- thermo %>%
  pivot_longer(
    cols = c("Exponential", "Stationary"),
    names_to = "Phase",
    values_to = "Count"
  ) %>%
  mutate(
    Class = factor(Class, levels = thermo_order),
    Phase = factor(Phase, levels = c("Exponential", "Stationary"))
  )

p_thermo <- ggplot(thermo_long, aes(x = Count, y = Class, fill = Phase)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62
  ) +
  geom_text(
    data = thermo_long %>% filter(Count > 0),
    aes(label = Count),
    position = position_dodge(width = 0.72),
    hjust = -0.15,
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(0, 30),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_fill_manual(
    values = c(
      "Exponential" = "#d95f5f",
      "Stationary"  = "#8c2d2d"
    )
  ) +
  labs(
    x = "Number of predicted RNA thermometers",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

# -----------------------------
# Panel B: Riboswitch classes
# IMPORTANT:
# Replace these class counts with your final verified values from Tables S7-S8.
# The class names below reflect biological classes and remove the placeholder bucket.
# -----------------------------
ribo <- tribble(
  ~Class,         ~Exponential, ~Stationary,
  "suhB",         3,            11,
  "TPP",          2,            3,
  "SAM",          0,            2,
  "FMN",          0,            1,
  "Lysine",       0,            1,
  "Mg_sensor",    1,            1,
  "glmS",         0,            1,
  "PreQ1",        0,            1,
  "yybP-ykoY",    1,            1,
  "speF",         1,            1,
  "serC",         1,            1
)

ribo_order <- ribo %>%
  arrange(Stationary, Exponential) %>%
  pull(Class)

ribo_long <- ribo %>%
  pivot_longer(
    cols = c("Exponential", "Stationary"),
    names_to = "Phase",
    values_to = "Count"
  ) %>%
  mutate(
    Class = factor(Class, levels = ribo_order),
    Phase = factor(Phase, levels = c("Exponential", "Stationary"))
  )

p_ribo <- ggplot(ribo_long, aes(x = Count, y = Class, fill = Phase)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62
  ) +
  geom_text(
    data = ribo_long %>% filter(Count > 0),
    aes(label = Count),
    position = position_dodge(width = 0.72),
    hjust = -0.15,
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(0, 14),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_fill_manual(
    values = c(
      "Exponential" = "#4c9ed9",
      "Stationary"  = "#1f5a85"
    )
  ) +
  labs(
    x = "Number of predicted riboswitches",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

# -----------------------------
# Save outputs: PNG + PDF only
# -----------------------------
ggsave(
  filename = file.path(out_dir, "panel_A_rna_thermometer_classes.png"),
  plot = p_thermo,
  width = 7.2,
  height = 4.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "panel_A_rna_thermometer_classes.pdf"),
  plot = p_thermo,
  width = 7.2,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "panel_B_riboswitch_classes.png"),
  plot = p_ribo,
  width = 7.2,
  height = 4.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "panel_B_riboswitch_classes.pdf"),
  plot = p_ribo,
  width = 7.2,
  height = 4.8,
  device = cairo_pdf
)

message("Done. Figures saved to: ", out_dir)

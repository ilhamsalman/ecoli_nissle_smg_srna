#!/usr/bin/env Rscript
# Figure 5: Rfam family distributions from GitHub-tracked ANNOgesic CSV exports.
# Run from the repository root:
#   Rscript scripts/08_figures/make_cis_element_class_distribution.R
# There are NO hard-coded family counts and NO dependencies on local annogesic_* directories.
# Phase provenance was verified independently: annogesic_exp = exp WIGs;
# annogesic_sta = sta WIGs. Therefore exp has 43 thermometer/27 riboswitch
# module records; sta has 14 thermometer/11 riboswitch module records.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source_dir <- "results/riboswitch_thermometer"
out_dir <- "results/figures/cis_elements"
if (!dir.exists(source_dir)) {
  stop("Run from the repository root: missing ", source_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_annogesic <- function(phase, filename, expected) {
  path <- file.path(source_dir, phase, filename)
  if (!file.exists(path)) stop("Missing source file: ", path)
  # Despite the .csv extension, ANNOgesic writes tab-delimited tables.
  x <- read.delim(path, sep = "\t", header = TRUE, comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
  if (!"Rfam_name" %in% names(x)) stop("Missing Rfam_name column: ", path)
  if (nrow(x) != expected) {
    stop("Unexpected record count in ", path, ": ", nrow(x),
         "; expected ", expected, ". Check the GitHub-tracked ANNOgesic CSV exports.")
  }
  if (anyNA(x$Rfam_name) || any(!nzchar(trimws(x$Rfam_name)))) {
    stop("Missing or blank Rfam_name found in ", path)
  }
  data.frame(Class = trimws(x$Rfam_name),
             Phase = if (phase == "exp") "Exponential" else "Stationary",
             stringsAsFactors = FALSE)
}

make_counts <- function(filename, exp_expected, sta_expected) {
  bind_rows(
    read_annogesic("exp", filename, exp_expected),
    read_annogesic("sta", filename, sta_expected)
  ) %>%
    count(Class, Phase, name = "Count") %>%
    complete(Class, Phase = c("Exponential", "Stationary"),
             fill = list(Count = 0L))
}

thermo_counts <- make_counts("annotation_RNA_thermometer.csv", 43L, 14L)
ribo_counts   <- make_counts("annotation_riboswitch.csv", 27L, 11L)

# Retain exact Rfam names, including SAM_alpha and SAM-SAH as separate families.
# suhB and drz-agam-1 are kept as module predictions, not reclassified
# as biologically confirmed riboswitches.
plot_counts <- function(counts, axis_label, palette) {
  order <- counts %>%
    pivot_wider(names_from = Phase, values_from = Count) %>%
    arrange(Stationary, Exponential) %>%
    pull(Class)
  dat <- counts %>%
    mutate(Class = factor(Class, levels = order),
           Phase = factor(Phase, levels = c("Exponential", "Stationary")))

  ggplot(dat, aes(x = Count, y = Class, fill = Phase)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(data = dat %>% filter(Count > 0), aes(label = Count),
              position = position_dodge(width = 0.72), hjust = -0.15,
              size = 3.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.17))) +
    scale_fill_manual(values = palette) +
    labs(x = axis_label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.title = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10))
}

p_thermo <- plot_counts(
  thermo_counts,
  "Number of ANNOgesic thermometer-module predictions",
  c("Exponential" = "#d95f5f", "Stationary" = "#8c2d2d")
)
p_ribo <- plot_counts(
  ribo_counts,
  "Number of ANNOgesic riboswitch-module predictions",
  c("Exponential" = "#4c9ed9", "Stationary" = "#1f5a85")
)

ggsave(file.path(out_dir, "panel_A_rna_thermometer_classes.png"),
       p_thermo, width = 7.2, height = 4.8, dpi = 600)
ggsave(file.path(out_dir, "panel_A_rna_thermometer_classes.pdf"),
       p_thermo, width = 7.2, height = 4.8)
ggsave(file.path(out_dir, "panel_B_riboswitch_classes.png"),
       p_ribo, width = 7.2, height = 5.8, dpi = 600)
ggsave(file.path(out_dir, "panel_B_riboswitch_classes.pdf"),
       p_ribo, width = 7.2, height = 5.8)

# Save the values actually plotted as an audit trail.
write.csv(bind_rows(mutate(thermo_counts, Panel = "A"),
                    mutate(ribo_counts, Panel = "B")) %>%
            select(Panel, Class, Phase, Count),
          file.path(out_dir, "figure_5_plot_counts.csv"), row.names = FALSE)

cat("Figure 5 regenerated using GitHub-tracked ANNOgesic CSV exports.\n")
cat("Thermometer-module records: exp = 43, sta = 14\n")
cat("Riboswitch-module records:  exp = 27, sta = 11\n")
cat("Saved figure files and figure_5_plot_counts.csv in: ", out_dir, "\n", sep = "")

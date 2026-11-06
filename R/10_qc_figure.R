# =============================================================================
# 10_qc_figure.R  --  QC figures (before/after) from the depositor metadata
# -----------------------------------------------------------------------------
# Input : data/raw/<GSE>/<GSM>_metaData.rds   (all called barcodes, pre-QC)
#         results/tables/01_cell_numbers.csv
# Output: results/figures/10_qc_*.png
#
# The pre-QC object is not kept in memory (it would double the memory), so the
# "before" distributions are re-read from the small per-sample metadata files.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
library(Matrix)

log_msg("STEP 10 -- QC figures")
samples <- PARAMS$samples
qc <- PARAMS$qc

md <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples$sample[i]
  d <- readRDS(file.path(RAW_DIR, samples$dataset[i], paste0(s, "_metaData.rds")))
  data.frame(sample = s, dataset = samples$dataset[i], group = samples$group[i],
             nGene = d$nGene, nUMI = d$nUMI, percent.mt = d$percent.mito)
}) %>% bind_rows()

long <- md %>% select(sample, dataset, nGene, nUMI) %>%
  pivot_longer(c(nGene, nUMI), names_to = "metric", values_to = "value") %>%
  filter(is.finite(value), value > 0)

p1 <- ggplot(long, aes(sample, value, fill = dataset)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_y_log10() + labs(x = NULL, y = NULL, title = "Depth per cell before filtering") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p1, "10_qc_violin", 12, 4.5)

p1b <- ggplot(md, aes(sample, percent.mt, fill = dataset)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_hline(yintercept = qc$max_mito_pct, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 40)) +
  labs(x = NULL, y = "mitochondrial %", title = "Mitochondrial content (dashed = cutoff)") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p1b, "10_qc_mito", 9, 4)

cn <- read.csv(file.path(TAB_DIR, "01_cell_numbers.csv"))
cn_long <- cn %>% select(sample, cells_called, pass_qc, cells_final) %>%
  pivot_longer(-sample, names_to = "stage", values_to = "n") %>%
  mutate(stage = factor(stage, levels = c("cells_called", "pass_qc", "cells_final"),
                        labels = c("called", "after QC", "after doublets")))
p2 <- ggplot(cn_long, aes(sample, n, fill = stage)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(called = "grey75", `after QC` = "#4E79A7",
                               `after doublets` = "#E15759")) +
  labs(x = NULL, y = "cells", fill = NULL, title = "Cells retained per library") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p2, "10_qc_cell_numbers", 9, 4.5)

p3 <- ggplot(md, aes(nGene, percent.mt, colour = dataset)) +
  geom_point(size = 0.15, alpha = 0.25) +
  geom_hline(yintercept = qc$max_mito_pct, linetype = "dashed") +
  geom_vline(xintercept = c(qc$min_genes_per_cell, qc$max_genes_per_cell), linetype = "dashed") +
  scale_x_log10() + labs(x = "genes per cell", y = "mitochondrial %",
                         title = "QC thresholds") + theme_pub()
save_fig(p3, "10_qc_scatter", 6, 5)

log_msg("STEP 10 done")

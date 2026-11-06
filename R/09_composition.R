# =============================================================================
# 09_composition.R  --  DC infiltration and subtype composition
# -----------------------------------------------------------------------------
# Input : results/objects/03_annotated.rds  (all cells, for DC fraction)
#         results/objects/05_dc_annotated.rds (DC subtypes)
# Output: results/tables/09_*.csv, results/figures/09_*.png
#
# Only two patients (P7, P8) have paired tumour/normal tissue, so the comparison
# is descriptive; no hypothesis test is reported for n=2 pairs.  This is stated
# explicitly in the figure caption and in docs/limitations.md.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(ggplot2); library(tidyr)})

log_msg("STEP 09 -- composition")
all <- readRDS(file.path(OBJ_DIR, "03_annotated.rds"))
dc  <- readRDS(file.path(OBJ_DIR, "05_dc_annotated.rds"))

# ---- 1. DC fraction of all cells, per sample ---------------------------
comp <- all@meta.data %>% count(sample, group, patient, cell_type) %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
dc_frac <- comp %>% filter(cell_type == "DC") %>%
  select(sample, group, patient, n_dc = n, frac_dc = frac)
write.csv(dc_frac, file.path(TAB_DIR, "09_dc_fraction_by_sample.csv"), row.names = FALSE)

pd <- position_dodge(0.2)
p1 <- ggplot(dc_frac, aes(group, 100 * frac_dc, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.5) +
  geom_point(size = 2, position = pd) +
  geom_line(aes(group = patient), colour = "grey40", linewidth = 0.4, position = pd) +
  scale_fill_manual(values = c(tumour = "#D62728", normal = "#1F77B4")) +
  labs(x = NULL, y = "DC (% of all cells)", title = "DC infiltration") +
  theme_pub() + theme(legend.position = "none")
save_fig(p1, "09_dc_infiltration", 5, 5)

# ---- 2. subtype composition within DC ----------------------------------
sub <- dc@meta.data %>% count(sample, group, patient, dc_subtype) %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(sub, file.path(TAB_DIR, "09_dc_subtype_fraction.csv"), row.names = FALSE)

p2 <- ggplot(sub, aes(sample, frac, fill = dc_subtype)) +
  geom_col(width = 0.8) + scale_fill_manual(values = celltype_pal(levels(sub$dc_subtype))) +
  geom_text(data = sub %>% group_by(sample) %>% summarise(n = sum(n), .groups = "drop"),
            aes(sample, 1.02, label = paste0("n=", n)), inherit.aes = FALSE, size = 2.6) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08)) +
  labs(x = NULL, y = "Fraction of DC", title = "DC subtype composition per library") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p2, "09_dc_subtype_composition", 8, 5)

# pooled composition across all DCs (library-level counts are too small to
# support per-sample proportions; see docs/limitations.md)
pooled <- dc@meta.data %>% count(dc_subtype) %>% mutate(frac = n / sum(n))
write.csv(pooled, file.path(TAB_DIR, "09_dc_subtype_pooled.csv"), row.names = FALSE)
p2b <- ggplot(pooled, aes(reorder(dc_subtype, -frac), frac, fill = dc_subtype)) +
  geom_col(width = 0.7) + scale_fill_manual(values = celltype_pal(pooled$dc_subtype)) +
  geom_text(aes(label = paste0(n, "\n(", round(100 * frac, 1), "%)")), vjust = -0.2, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(pooled$frac) * 1.2)) +
  labs(x = NULL, y = "Fraction of all DC", title = "Pooled DC subtype composition") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p2b, "09_dc_subtype_pooled", 6, 4.5)

p3 <- ggplot(sub, aes(dc_subtype, frac, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6) +
  geom_point(position = position_dodge(0.6), size = 1.8) +
  scale_fill_manual(values = c(tumour = "#D62728", normal = "#1F77B4")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Fraction of DC", title = "Subtype fraction: tumour vs normal (n=2 pairs)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p3, "09_dc_subtype_tumor_vs_normal", 8, 5)

# ---- 3. cDC1 / cDC2 balance --------------------------------------------
bal <- sub %>% select(sample, group, patient, dc_subtype, frac) %>%
  pivot_wider(names_from = dc_subtype, values_from = frac, values_fill = 0)
if (all(c("cDC1", "cDC2") %in% colnames(bal))) {
  bal$cDC1_to_cDC2 <- bal$cDC1 / pmax(bal$cDC2, 1e-9)
  write.csv(bal, file.path(TAB_DIR, "09_cDC1_cDC2_balance.csv"), row.names = FALSE)
  p4 <- ggplot(bal, aes(group, cDC1_to_cDC2, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.5) +
    geom_point(size = 2) + scale_fill_manual(values = c(tumour = "#D62728", normal = "#1F77B4")) +
    labs(x = NULL, y = "cDC1 / cDC2 ratio", title = "cDC1 versus cDC2 balance") +
    theme_pub() + theme(legend.position = "none")
  save_fig(p4, "09_cdc1_cdc2_ratio", 5, 5)
}
log_msg("STEP 09 done")

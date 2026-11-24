# Stage 09: library-level descriptive composition from stages 03 and 05.
source("R/00_setup.R")
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr)})

log_msg("STEP 09 -- composition")
major <- read.csv(file.path(TAB_DIR, "03_celltype_by_sample.csv")) %>%
  rename(cell_type = Var1, sample = Var2, n = Freq) %>%
  sample_metadata()
dc <- read.csv(file.path(TAB_DIR, "05_dc_subtype_by_sample.csv")) %>%
  sample_metadata()
stopifnot(setequal(unique(major$sample), PARAMS$samples$sample),
          !anyDuplicated(major[c("sample", "cell_type")]),
          !anyDuplicated(dc[c("sample", "dc_subtype")]))

dc_frac <- major %>%
  group_by(sample, dataset, group, patient) %>%
  summarise(n_atlas = sum(n), n_dc = sum(n[cell_type %in% DC_LABELS]), .groups = "drop") %>%
  left_join(dc %>% group_by(sample) %>% summarise(n_retained_dc = sum(n), .groups = "drop"),
            by = "sample") %>%
  mutate(n_retained_dc = replace_na(n_retained_dc, 0),
         frac_dc = n_dc / n_atlas, frac_retained_dc = n_retained_dc / n_atlas)
stopifnot(all(dc_frac$n_atlas > 0), all(dc_frac$n_retained_dc <= dc_frac$n_dc))
write.csv(dc_frac, file.path(TAB_DIR, "09_dc_fraction_by_sample.csv"), row.names = FALSE)

pooled <- dc %>% count(dc_subtype, wt = n) %>% mutate(frac = n / sum(n))
pal <- celltype_pal(pooled$dc_subtype)
sub <- expand_grid(sample = dc_frac$sample, dc_subtype = pooled$dc_subtype) %>%
  left_join(select(dc, sample, dc_subtype, n), by = c("sample", "dc_subtype")) %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(select(dc_frac, sample, n_retained_dc), by = "sample") %>%
  sample_metadata() %>%
  mutate(frac = if_else(n_retained_dc > 0, n / n_retained_dc, NA_real_))
write.csv(sub, file.path(TAB_DIR, "09_dc_subtype_fraction.csv"), row.names = FALSE)
write.csv(pooled, file.path(TAB_DIR, "09_dc_subtype_pooled.csv"), row.names = FALSE)

paired <- dc_frac %>% group_by(patient) %>% filter(n_distinct(group) == 2L) %>% ungroup()
paired_summary <- paired %>%
  select(patient, group, frac_dc, frac_retained_dc) %>%
  pivot_wider(names_from = group, values_from = c(frac_dc, frac_retained_dc)) %>%
  mutate(candidate_tumour_minus_normal = frac_dc_tumour - frac_dc_normal,
         retained_tumour_minus_normal = frac_retained_dc_tumour - frac_retained_dc_normal)
write.csv(paired_summary, file.path(TAB_DIR, "09_paired_tissue_descriptive.csv"), row.names = FALSE)

paired_plot <- function(d, y, ylabel) {
  ggplot(d, aes(group, .data[[y]], group = patient, colour = patient)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    labs(x = NULL, y = ylabel, colour = "Donor") + theme_pub()
}

paired_long <- paired %>%
  pivot_longer(c(frac_dc, frac_retained_dc), names_to = "definition", values_to = "fraction") %>%
  mutate(definition = factor(definition, c("frac_dc", "frac_retained_dc"),
                             c("DC candidates", "Retained DCs")))
p1 <- paired_plot(paired_long, "fraction", "Fraction of atlas cells") +
  facet_wrap(~ definition) +
  labs(title = "Paired tumour and adjacent tissue", subtitle = "Two donors, descriptive")
save_fig(p1, "09_dc_infiltration", 7, 4.2)

p2 <- ggplot(sub, aes(sample, frac, fill = dc_subtype)) +
  geom_col(width = 0.8, na.rm = TRUE) +
  scale_fill_manual(values = pal) +
  geom_text(data = dc_frac, aes(sample, 1.025, label = n_retained_dc),
            inherit.aes = FALSE, size = 2.7) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.08)) +
  labs(x = NULL, y = "Fraction of retained DCs", fill = "Subtype",
       title = "DC composition by library") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p2, "09_dc_subtype_composition", 8, 5)

p2b <- ggplot(pooled, aes(reorder(dc_subtype, -frac), frac, fill = dc_subtype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = pal) +
  geom_text(aes(label = paste0(n, "\n", scales::percent(frac, accuracy = 0.1))),
            vjust = -0.2, size = 3) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, max(pooled$frac) * 1.2)) +
  labs(x = NULL, y = "Fraction of retained DCs", title = "Pooled DC composition") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p2b, "09_dc_subtype_pooled", 6, 4.5)

p3 <- paired_plot(filter(sub, patient %in% paired$patient), "frac", "Fraction of retained DCs") +
  facet_wrap(~ dc_subtype, ncol = 3, scales = "free_y") +
  labs(title = "Subtype composition within matched tissues", subtitle = "Two donors, descriptive")
save_fig(p3, "09_dc_subtype_tumor_vs_normal", 8, 5)

bal <- sub %>%
  select(sample, dataset, group, patient, dc_subtype, frac) %>%
  pivot_wider(names_from = dc_subtype, values_from = frac) %>%
  mutate(cDC1_to_cDC2 = if_else(cDC2 > 0, cDC1 / cDC2, NA_real_))
write.csv(bal, file.path(TAB_DIR, "09_cDC1_cDC2_balance.csv"), row.names = FALSE)
p4 <- ggplot(bal, aes(group, cDC1_to_cDC2, colour = dataset)) +
  geom_point(position = position_jitter(width = 0.1, height = 0, seed = SEED),
             size = 2.3, na.rm = TRUE) +
  labs(x = NULL, y = "cDC1 / cDC2", colour = "GEO series", title = "DC balance by library") +
  theme_pub()
save_fig(p4, "09_cdc1_cdc2_ratio", 5.5, 4.5)
log_msg("STEP 09 done")

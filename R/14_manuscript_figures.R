# Stage 14: publication figures using the reviewed atlas estimates.
started <- proc.time()[["elapsed"]]
source("R/00_setup.R")
source("R/pseudobulk.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

log_msg("STEP 14 -- manuscript figures")
cfg <- PARAMS$pseudobulk
read_table <- function(name) read.csv(file.path(TAB_DIR, paste0(name, ".csv")))
count_label <- function(n) format(n, big.mark = ",", trim = TRUE, scientific = FALSE)
p_label <- function(p) sprintf("%.3g", p)
label_lookup <- function(labels) {
  force(labels)
  function(breaks) unname(labels[names(breaks)])
}
neg_log <- function(p) -log10(pmax(p, .Machine$double.xmin))
series_colours <- c(GSE172577 = "#117B75", GSE193766 = "#8D4F91")
migdc_colour <- "#117B75"
cdc2_colour <- "#B05C40"
neutral_colour <- "#A5B0B7"
program_order <- c("Maturation", "NF-kB feedback", "Regulatory expression", "MHC class II")

manuscript_theme <- function(size = 10) {
  theme_pub(size) +
    theme(legend.position = "bottom",
          plot.title.position = "plot",
          plot.title = element_text(face = "bold", colour = "#203541"),
          plot.subtitle = element_text(colour = "#43545C", size = rel(0.9),
                                       margin = margin(b = 5)),
          plot.caption = element_text(hjust = 0, colour = "#43545C"),
          strip.text = element_text(face = "bold"),
          plot.tag = element_text(face = "bold", colour = "#203541"))
}
save_composite <- function(p, filename, title, width, height) {
  p <- p + plot_annotation(title = title, tag_levels = "A",
                             theme = theme(plot.title = element_text(face = "bold", size = 16,
                                                                    colour = "#203541")))
  save_fig(p, filename, width = width, height = height, dpi = 300)
  log_msg("  Saved ", filename, " (PNG and PDF)")
}

# Fig 1: atlas and retained DC composition.
library_labels <- setNames(PARAMS$samples$patient, PARAMS$samples$sample)
library_groups <- c("GSE172577\nPrimary tumour", "GSE193766\nPrimary tumour",
                    "GSE193766\nAdjacent tissue")
label_libraries <- function(rows) {
  sample_metadata(rows) %>%
    mutate(sample = factor(sample, levels = PARAMS$samples$sample),
           tissue = if_else(group == "tumour", "Primary tumour", "Adjacent tissue"),
           library_group = factor(paste(dataset, tissue, sep = "\n"), levels = library_groups))
}
major <- read_table("03_celltype_by_sample") %>%
  rename(cell_type = Var1, sample = Var2, n = Freq) %>%
  label_libraries()
major_order <- major %>% count(cell_type, wt = n, sort = TRUE)
major <- major %>%
  group_by(sample) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = rev(major_order$cell_type)))
atlas_composition <- ggplot(major, aes(sample, cell_type, fill = fraction)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  facet_grid(~ library_group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = library_labels) +
  scale_fill_gradient(low = "#F2F5F4", high = "#155F58", transform = "sqrt",
                      breaks = c(0, 0.05, 0.20, 0.50), labels = scales::label_percent(),
                      name = "Library share", guide = guide_colourbar(barwidth = grid::unit(38, "mm"))) +
  labs(title = "Major cell types by library",
       subtitle = paste(count_label(sum(major$n)), "cells across",
                        n_distinct(major$cell_type), "annotated types"),
       x = "Patient library", y = NULL) +
  manuscript_theme(9) + theme(panel.grid = element_blank())

dc_types <- c("cDC1", "cDC2", "LC-like", "migDC", "pDC", "Proliferating DC")
dc_pooled <- read_table("09_dc_subtype_pooled")
dc_legend <- setNames(paste0(dc_pooled$dc_subtype, " (", dc_pooled$n, ")"), dc_pooled$dc_subtype)
dc_by_library <- read_table("05_dc_subtype_by_sample") %>%
  label_libraries() %>%
  mutate(dc_subtype = factor(dc_subtype, levels = dc_types))
dc_composition <- ggplot(dc_by_library, aes(sample, n, fill = dc_subtype)) +
  geom_col(width = 0.72) +
  facet_grid(~ library_group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = library_labels) +
  scale_fill_manual(values = celltype_pal(dc_types), labels = dc_legend,
                    guide = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(title = "Retained DC subtypes by library",
       subtitle = paste(count_label(sum(dc_pooled$n)), "retained DCs in 10 libraries"),
       x = "Patient library", y = "Retained DCs", fill = "Subtype (pooled n)") +
  manuscript_theme(9) + theme(panel.grid.major.x = element_blank())

dc <- read_atlas_object("05_dc_annotated.rds")
marker_sets <- setNames(list(c("CLEC9A", "XCR1", "CADM1"), c("CD1C", "CLEC10A", "FCER1A"),
                            c("CD207", "CD1A", "CD1E"), c("LAMP3", "CCR7", "FSCN1"),
                            c("LILRA4", "CLEC4C", "TCF4"), c("MKI67", "TOP2A", "SPP1")), dc_types)
marker_genes <- unlist(marker_sets, use.names = FALSE)
marker_groups <- setNames(rep(names(marker_sets), lengths(marker_sets)), marker_genes)
expression <- as.matrix(SeuratObject::LayerData(dc, assay = "RNA", layer = "data", features = marker_genes))
stopifnot(setequal(rownames(expression), marker_genes))
subtypes <- dc@meta.data[colnames(expression), "dc_subtype"]
marker_values <- bind_rows(lapply(dc_types, function(subtype) {
  values <- expression[, subtypes == subtype, drop = FALSE]
  data.frame(gene = rownames(values), dc_subtype = subtype,
             mean_expression = log1p(rowMeans(expm1(values))), detected = rowMeans(values > 0))
})) %>%
  group_by(gene) %>%
  mutate(scaled_expression = as.numeric(scale(mean_expression))) %>%
  ungroup() %>%
  mutate(marker_group = factor(marker_groups[gene], levels = dc_types),
         dc_subtype = factor(dc_subtype, levels = rev(dc_types)),
         gene = factor(gene, levels = marker_genes))
marker_plot <- ggplot(marker_values, aes(gene, dc_subtype)) +
  geom_point(aes(colour = scaled_expression, size = detected)) +
  facet_grid(~ marker_group, scales = "free_x", space = "free_x",
             labeller = as_labeller(c(cDC1 = "cDC1", cDC2 = "cDC2", "LC-like" = "LC-like",
                                      migDC = "migDC", pDC = "pDC", "Proliferating DC" = "Cycling"))) +
  scale_colour_gradient2(low = "#D6DDE0", mid = "#F4F5F2", high = "#126961",
                         name = "Gene-scaled mean") +
  scale_size_continuous(range = c(0, 4.3), limits = c(0, 1), labels = scales::label_percent(),
                        name = "Detected") +
  labs(title = "Core annotation markers", x = NULL, y = NULL) +
  manuscript_theme(9) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid.major.x = element_blank())

tissue_levels <- c("normal", "tumour")
tissue_labels <- c("Adjacent", "Primary tumour")
fraction_labels <- c(frac_dc = "DC candidates", frac_retained_dc = "Retained DCs")
fractions <- read_table("09_dc_fraction_by_sample") %>%
  sample_metadata() %>%
  pivot_longer(all_of(names(fraction_labels)), names_to = "measure", values_to = "fraction") %>%
  mutate(tissue = factor(group, levels = tissue_levels, labels = tissue_labels))
matched <- read_table("09_paired_tissue_descriptive") %>%
  pivot_longer(starts_with("frac_"), names_pattern = "^(frac_.*)_(normal|tumour)$",
               names_to = c("measure", "group"), values_to = "fraction") %>%
  mutate(tissue = factor(group, levels = tissue_levels, labels = tissue_labels))
infiltration_plot <- ggplot(fractions, aes(tissue, fraction)) +
  geom_line(data = matched, aes(group = patient), colour = "#7C7688", linewidth = 0.65) +
  geom_point(data = filter(fractions, !patient %in% matched$patient),
             aes(colour = dataset), size = 2.1,
             position = position_jitter(width = 0.12, height = 0, seed = SEED)) +
  geom_point(data = filter(fractions, patient %in% matched$patient), aes(colour = dataset), size = 2.4) +
  facet_wrap(~ measure, labeller = as_labeller(fraction_labels)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  scale_colour_manual(values = series_colours) +
  labs(title = "DC infiltration by tissue",
       subtitle = "Lines connect matched P7 and P8 tissue libraries",
       x = NULL, y = "Share of atlas cells", colour = "GEO series") +
  manuscript_theme(9)
save_composite((atlas_composition | dc_composition) / (marker_plot | infiltration_plot) +
                 plot_layout(widths = c(1.32, 1), heights = c(1.3, 1)),
               "Fig1_discovery_atlas", "Discovery atlas", 13, 10.3)

# Fig 2: paired differential expression and robustness.
de <- read_table("11_de_migDC_vs_cDC2")
volcano_genes <- c("LAMP3", "CCR7", "CCL19", "CCL22", "CD274", "IDO1", "CD200", "PDCD1LG2",
                   "CD83", "BIRC3", "SOCS2", "CD1C", "CLEC10A", "CD300E")
volcano <- de %>%
  mutate(log_fdr = neg_log(fdr),
         direction = case_when(fdr >= cfg$fdr ~ "FDR at or above threshold",
                               log2_fc > 0 ~ "Higher in migDC", TRUE ~ "Higher in cDC2"))
missing_volcano <- setdiff(volcano_genes, volcano$gene)
volcano_plot <- ggplot(volcano, aes(log2_fc, log_fdr)) +
  geom_hline(yintercept = neg_log(cfg$fdr), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.35) +
  geom_point(aes(colour = direction), size = 1, alpha = 0.53) +
  ggrepel::geom_text_repel(data = filter(volcano, gene %in% volcano_genes), aes(label = gene),
                           seed = SEED, size = 2.65, max.overlaps = Inf, min.segment.length = 0,
                           colour = "#203541", box.padding = 0.28) +
  scale_colour_manual(values = c("Higher in migDC" = migdc_colour, "Higher in cDC2" = cdc2_colour,
                                 "FDR at or above threshold" = neutral_colour),
                      breaks = c("Higher in migDC", "Higher in cDC2", "FDR at or above threshold"),
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Paired migDC vs cDC2",
       subtitle = paste0(count_label(nrow(de)), " genes. ", count_label(sum(de$fdr < cfg$fdr)),
                         " with FDR < ", cfg$fdr),
       x = "Paired log2 fold change", y = "-log10(edgeR FDR)", colour = NULL,
       caption = if (length(missing_volcano)) paste(paste(missing_volcano, collapse = ", "),
                                                     "fell below the expression filter.")) +
  manuscript_theme()

matched_tests <- de %>%
  filter(is.finite(cell_level_log2_fc), is.finite(cell_level_p_adj)) %>%
  mutate(wilcox_hit = cell_level_p_adj < cfg$fdr, edger_hit = fdr < cfg$fdr,
         call = case_when(wilcox_hit & edger_hit ~ "Both",
                          wilcox_hit ~ "Wilcoxon only", edger_hit ~ "edgeR only", TRUE ~ "Neither"))
method_axes <- c("Adjusted significance (-log10)", "Effect (log2 FC)")
method_points <- bind_rows(
  transmute(matched_tests, metric = method_axes[1], x = neg_log(cell_level_p_adj), y = neg_log(fdr), call),
  transmute(matched_tests, metric = method_axes[2], x = cell_level_log2_fc, y = log2_fc, call)
) %>%
  mutate(metric = factor(metric, levels = method_axes),
         call = factor(call, levels = c("Neither", "Both", "Wilcoxon only", "edgeR only"))) %>%
  arrange(call)
significance_lines <- data.frame(metric = factor(method_axes[1], levels = method_axes),
                                  line = neg_log(cfg$fdr))
method_plot <- ggplot(method_points, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_vline(data = significance_lines, aes(xintercept = line), colour = "grey50", linewidth = 0.35) +
  geom_hline(data = significance_lines, aes(yintercept = line), colour = "grey50", linewidth = 0.35) +
  geom_point(aes(colour = call), size = 0.9, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free", nrow = 1) +
  scale_colour_manual(values = c("Neither" = "#BBC4C9", "Both" = "#4B7387",
                                 "Wilcoxon only" = "#D99151", "edgeR only" = "#0E766C")) +
  labs(title = "Cell-level and donor-paired tests",
       subtitle = paste0(count_label(nrow(matched_tests)), " shared genes. ",
                         sum(matched_tests$call == "Wilcoxon only"), " Wilcoxon-only and ",
                         sum(matched_tests$call == "edgeR only"), " edgeR-only hits."),
       x = "Cell-level Wilcoxon", y = "Donor-paired edgeR", colour = NULL) +
  manuscript_theme()

sensitivity <- read_table("11_sensitivity")
robustness <- bind_rows(
  transmute(read_table("11_leave_one_donor_out"), label = paste("Omit", omitted_donor),
            category = "Leave one donor out", n_donors = n_remaining, rho = rho_with_primary),
  transmute(filter(sensitivity, analysis %in% paste0("min_", cfg$sensitivity_min_cells)),
            label = paste(min_cells, "cells"), category = "Cell minimum", n_donors, rho = rho_with_primary),
  transmute(filter(sensitivity, analysis == "normalization"), label = norm,
            category = "Normalization", n_donors, rho = rho_with_primary),
  transmute(filter(sensitivity, analysis == "strict_filter_min10"),
            label = paste0("Strict filter (", min_count, "/", min_total_count, ")"),
            category = "Expression filter", n_donors, rho = rho_with_primary),
  transmute(read_table("11_marker_gate_sensitivity"), label = paste0("Marker gate (", min_cells, " cells)"),
            category = "Marker gate", n_donors, rho = rho_with_primary)
) %>%
  mutate(label = paste0(label, "  (n = ", n_donors, ")"),
         label = factor(label, levels = rev(label)))
min_rho <- min(robustness$rho, na.rm = TRUE)
rho_labels <- with(robustness, setNames(ifelse(is.finite(rho), sprintf("%.5f", rho), ""),
                                       as.character(label)))
robustness_plot <- ggplot(robustness, aes(rho, label)) +
  geom_vline(xintercept = 1, colour = "grey65", linetype = "dashed") +
  geom_point(data = filter(robustness, is.finite(rho)), aes(colour = category), size = 2.5) +
  geom_text(data = filter(robustness, !is.finite(rho)),
            aes(x = (min_rho + 1) / 2, label = "Insufficient donors"), colour = "grey43", size = 2.8) +
  scale_x_continuous(limits = c(min_rho, 1), expand = expansion(mult = 0.06),
                     labels = scales::label_number(accuracy = 0.01)) +
  scale_y_discrete(limits = levels(robustness$label),
                   sec.axis = dup_axis(name = NULL, labels = label_lookup(rho_labels))) +
  scale_colour_manual(values = c("Leave one donor out" = "#31667D", "Cell minimum" = "#0D786C",
                                 "Normalization" = "#A25F41", "Expression filter" = "#82638D",
                                 "Marker gate" = "#74683C")) +
  labs(title = "Robustness of paired gene effects", x = "Spearman rho with the primary model",
       y = NULL, colour = NULL) +
  manuscript_theme(9) +
  theme(panel.grid.major.y = element_line(colour = "grey95"),
        axis.text.y.right = element_text(hjust = 0), axis.ticks.y.right = element_blank())
save_composite(((volcano_plot | method_plot) + plot_layout(widths = c(1, 1.34))) /
                 robustness_plot + plot_layout(heights = c(1.12, 1)),
               "Fig2_pseudobulk_de_and_robustness", "Paired differential expression and robustness", 13, 8.9)


# Fig 3: held-out programs and inflammatory pathway direction.
program_points <- read_table("11_program_donor_effects") %>%
  filter(analysis == "primary") %>%
  mutate(program = factor(program, levels = rev(program_order)))
program_summary <- read_table("11_program_summary") %>%
  filter(analysis == "primary") %>%
  mutate(program = factor(program, levels = rev(program_order)))
program_labels <- with(program_summary,
  setNames(sprintf("Mean %+.3f  |  Paired t P %s\nSign-flip P %s  |  fry FDR %s",
                   mean_effect, p_label(p_paired_t), p_label(p_sign_flip), p_label(fdr_fry)),
           as.character(program)))
program_plot <- ggplot(program_points, aes(effect, program)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_errorbar(data = program_summary, aes(y = program, xmin = lower_95_t, xmax = upper_95_t),
                inherit.aes = FALSE, orientation = "y", width = 0.13, colour = "#657786") +
  geom_point(aes(colour = dataset), size = 2.6,
             position = position_jitter(width = 0, height = 0.07, seed = SEED)) +
  geom_point(data = program_summary, aes(x = mean_effect), shape = 23, fill = "white", size = 3.4) +
  scale_y_discrete(sec.axis = dup_axis(name = NULL, labels = label_lookup(program_labels))) +
  scale_colour_manual(values = series_colours) +
  labs(title = "Held-out expression programs",
       subtitle = "Five paired tumour patients. Diamonds and bars show means and 95% t intervals.",
       x = "Paired log2 CPM difference (migDC - cDC2)", y = NULL, colour = "GEO series") +
  manuscript_theme() +
  theme(panel.grid.major.y = element_blank(), axis.text.y.right = element_text(hjust = 0),
        axis.ticks.y.right = element_blank())

camera <- read_table("11_paired_hallmark_camera") %>% arrange(FDR, PValue)
camera_core <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ALLOGRAFT_REJECTION",
                 "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
camera_display <- camera %>%
  filter(pathway %in% union(head(camera$pathway, 8), camera_core)) %>%
  mutate(label = tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", pathway)))),
         label = if_else(pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF-alpha via NF-kB", label),
         label = if_else(pathway %in% camera_core, paste0(label, " (q=", p_label(FDR), ")"), label),
         score = if_else(Direction == "Up", neg_log(FDR), -neg_log(FDR)),
         evidence = case_when(FDR >= cfg$fdr ~ "FDR at or above threshold",
                              Direction == "Up" ~ "Favours migDC", TRUE ~ "Favours cDC2"),
         label = factor(label, levels = rev(label)))
camera_plot <- ggplot(camera_display, aes(score, label, fill = evidence)) +
  geom_vline(xintercept = c(-1, 1) * neg_log(cfg$fdr), colour = "grey60", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_col(width = 0.67) +
  scale_fill_manual(values = c("Favours cDC2" = cdc2_colour, "Favours migDC" = migdc_colour,
                               "FDR at or above threshold" = neutral_colour),
                    guide = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Hallmark competitive enrichment",
       subtitle = "CAMERA on the gene universe outside the assignment set",
       x = "Signed -log10(CAMERA FDR)", y = NULL, fill = NULL) +
  manuscript_theme(9) + theme(panel.grid.major.y = element_blank())

primary_md <- read.csv(file.path(TAB_DIR, "11_primary_pseudobulk_libraries.csv"), row.names = "pseudobulk")
raw_counts <- as.matrix(read.delim(gzfile(file.path(TAB_DIR, "11_primary_pseudobulk_counts.tsv.gz")),
                                    row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(raw_counts), rownames(primary_md)))
primary_md$patient <- factor(primary_md$patient)
primary_md$dc_subtype <- factor(primary_md$dc_subtype, levels = c("cDC2", "migDC"))
primary_y <- fit_paired_bulk(list(counts = raw_counts, metadata = primary_md))$y
stopifnot(setequal(rownames(primary_y), de$gene))
mechanism_sets <- list("NF-kB feedback and survival" = c("BIRC3", "RELB", "NFKBIA", "TNFAIP3", "TRAF1", "TNIP1"),
                       "cDC2-associated myeloid genes" = c("TREM1", "OLR1", "CD300E", "THBS1", "F13A1", "CLEC10A"))
mechanism_genes <- unlist(mechanism_sets, use.names = FALSE)
mechanism_group <- setNames(rep(names(mechanism_sets), lengths(mechanism_sets)), mechanism_genes)
# Use the primary effective library sizes for all descriptive genes, including
# those whose expression was too sparse to enter the edgeR model.
effective_lib <- with(primary_y$samples, lib.size * norm.factors)
mechanism_cpm <- edgeR::cpm(raw_counts[mechanism_genes, , drop = FALSE],
                            lib.size = effective_lib, log = TRUE, prior.count = 1)
hi <- primary_md$dc_subtype == "migDC"
lo <- primary_md$dc_subtype == "cDC2"
stopifnot(identical(primary_md$patient[hi], primary_md$patient[lo]))
mechanism_diff <- mechanism_cpm[, hi, drop = FALSE] - mechanism_cpm[, lo, drop = FALSE]
colnames(mechanism_diff) <- as.character(primary_md$patient[hi])
mechanism_labels <- setNames(ifelse(mechanism_genes %in% de$gene, mechanism_genes,
                                    paste0(mechanism_genes, "*")), mechanism_genes)
mechanism_data <- gene_effects_long(mechanism_diff, mechanism_genes) %>%
  mutate(block = factor(mechanism_group[gene], levels = names(mechanism_sets)),
         gene = factor(gene, levels = rev(mechanism_genes)))
effect_limit <- max(abs(mechanism_data$log2_cpm_difference))
mechanism_plot <- ggplot(mechanism_data, aes(patient, gene, fill = log2_cpm_difference)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_wrap(~ block, scales = "free_y", space = "free_y", ncol = 1) +
  scale_y_discrete(labels = mechanism_labels) +
  scale_fill_gradient2(low = cdc2_colour, mid = "#F5F4F0", high = migdc_colour,
                       limits = c(-effect_limit, effect_limit), name = "migDC - cDC2") +
  labs(title = "Paired inflammatory gene expression",
       x = "Tumour patient", y = NULL,
       caption = "Colour shows TMM log2 CPM differences. * Below the DE expression filter.") +
  manuscript_theme(9) + theme(panel.grid = element_blank())
save_composite(program_plot / ((camera_plot | mechanism_plot) + plot_layout(widths = c(1.06, 1))) +
                 plot_layout(heights = c(0.82, 1.28)),
               "Fig3_programs_and_inflammation_paradox", "Expression programs and inflammatory context", 13, 9.3)

# Fig 4: lineage comparisons and patient-level TME covariation.
lineages <- c(primary = "migDC - cDC2", migDC_vs_cDC1 = "migDC - cDC1",
              migDC_vs_LClike = "migDC - LC-like", cDC1_vs_cDC2 = "cDC1 - cDC2",
              LClike_vs_cDC2 = "LC-like - cDC2")
lineage_genes <- c("IDO1", "CD274", "PDCD1LG2", "CD200", "CD83", "TNFAIP3")
lineage_data <- read_table("11_panel_donor_effects") %>%
  filter(analysis %in% names(lineages), gene %in% lineage_genes)
lineage_missing <- expand_grid(analysis = names(lineages), gene = lineage_genes) %>%
  anti_join(distinct(lineage_data, analysis, gene), by = c("analysis", "gene"))
lineage_levels <- function(rows) {
  mutate(rows, analysis = factor(analysis, levels = names(lineages), labels = unname(lineages)),
         gene = factor(gene, levels = rev(lineage_genes)))
}
lineage_data <- lineage_levels(lineage_data)
lineage_missing <- lineage_levels(lineage_missing)
lineage_mean <- lineage_data %>%
  group_by(analysis, gene) %>%
  summarise(log2_cpm_difference = mean(log2_cpm_difference), .groups = "drop")
lineage_plot <- ggplot(lineage_data, aes(log2_cpm_difference, gene)) +
  geom_vline(xintercept = 0, colour = "grey63", linetype = "dashed") +
  geom_point(aes(colour = dataset), size = 1.9,
             position = position_jitter(width = 0, height = 0.10, seed = SEED)) +
  geom_point(data = lineage_mean, shape = 23, fill = "white", size = 2.8) +
  geom_text(data = lineage_missing, aes(x = 0, y = gene), label = "Filtered", inherit.aes = FALSE,
            colour = "grey48", hjust = -0.12, size = 2.5) +
  facet_wrap(~ analysis, nrow = 1) +
  scale_colour_manual(values = series_colours) +
  labs(title = "Maturation and regulatory genes across DC lineages",
       subtitle = "Dots are tumour patients. Diamonds show paired means.",
       x = "Paired log2 CPM difference", y = NULL, colour = "GEO series") +
  manuscript_theme(9) + theme(panel.grid.major.y = element_blank())

tme <- read_table("13_tme_ligand_receptor_correlations")
profiles <- read_table("13_tme_donor_profiles")
tme_choices <- tribble(
  ~predictor,                       ~response,                   ~axis_label,
  "dc_CCL17",                       "treg_CCR4",                 "Retained DC CCL17 / Treg CCR4",
  "dc_CCL19",                       "tcell_CCR7",                "Retained DC CCL19 / T cells CCR7",
  "dc_CCL22",                       "treg_CCR4",                 "Retained DC CCL22 / Treg CCR4",
  "migdc_of_dc",                    "treg_of_tnk",               "migDC share / Treg share of T and NK",
  "migdc_CD274",                    "tcell_PDCD1",               "migDC CD274 / T cells PDCD1",
  "dc_IDO1",                        "treg_CTLA4",                "Retained DC IDO1 / Treg CTLA4",
  "primary_maturation",             "tcell_inhibitory_module",   "Maturation / T-cell inhibitory module",
  "primary_nfkb_feedback",          "tcell_inhibitory_module",   "NF-kB feedback / T-cell inhibitory module",
  "primary_regulatory_expression",  "tcell_inhibitory_module",   "Regulatory expression / T-cell inhibitory module",
  "primary_regulatory_expression",  "treg_CCR4",                 "Regulatory expression / Treg CCR4"
)
tme_display <- left_join(tme_choices, tme, by = c("predictor", "response")) %>%
  mutate(axis_label = factor(axis_label, levels = rev(axis_label)),
         highlighted = predictor %in% c("dc_CCL17", "dc_CCL19"))
tme_fdr <- setNames(paste("q", p_label(tme_display$fdr)), as.character(tme_display$axis_label))
covariation_plot <- ggplot(tme_display, aes(rho, axis_label)) +
  geom_vline(xintercept = 0, colour = "grey63", linewidth = 0.4) +
  geom_col(aes(fill = highlighted), width = 0.64) +
  scale_fill_manual(values = c("FALSE" = "#8196A4", "TRUE" = migdc_colour), guide = "none") +
  scale_y_discrete(sec.axis = dup_axis(name = NULL, labels = label_lookup(tme_fdr))) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  labs(title = "Patient-level TME covariation",
       subtitle = "Selected ligand, composition and paired-program axes",
       x = "Spearman rho", y = NULL) +
  manuscript_theme(9) +
  theme(panel.grid.major.y = element_blank(), axis.text.y.right = element_text(hjust = 0),
        axis.ticks.y.right = element_blank())

scatter_axes <- tme_display[tme_display$highlighted, ]
scatter_names <- c("CCL17 and Treg CCR4", "CCL19 and T cells CCR7")
scatter_labels <- sprintf("%s\nrho = %.2f  |  BH FDR = %s  |  n = %d",
                          scatter_names, scatter_axes$rho, p_label(scatter_axes$fdr), scatter_axes$n_donors)
tme_points <- bind_rows(lapply(seq_len(nrow(scatter_axes)), function(i) {
  transmute(profiles, dataset, panel = scatter_labels[i],
            ligand = .data[[scatter_axes$predictor[i]]], receptor = .data[[scatter_axes$response[i]]])
})) %>%
  mutate(panel = factor(panel, levels = scatter_labels))
scatter_plot <- ggplot(tme_points, aes(ligand, receptor)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey54", linewidth = 0.6) +
  geom_point(aes(colour = dataset), size = 2.5) +
  facet_wrap(~ panel, ncol = 1, scales = "free") +
  scale_colour_manual(values = series_colours) +
  labs(title = "Chemokine and receptor expression",
       subtitle = "Tumour patients and pooled linear trends",
       x = "Retained DC ligand (log2 CPM)", y = "Receptor (log2 CPM)", colour = "GEO series") +
  manuscript_theme(9)
save_composite(lineage_plot / ((covariation_plot | scatter_plot) + plot_layout(widths = c(1.35, 1))) +
                 plot_layout(heights = c(0.92, 1.32)),
               "Fig4_lineage_and_tme_context", "DC lineage and tumour microenvironment context", 13.2, 9.1)


# Fig 5: external gating, replication, and shared-gene program effects.
external_classification <- read_table("12_classification")
external_counts <- external_classification %>%
  filter(gate_min_detected == 2L, subtype %in% c("migDC", "cDC1", "cDC2")) %>%
  group_by(author_cell_type, subtype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(subtype = factor(subtype, levels = c("migDC", "cDC1", "cDC2")))
mac_cdc2_share <- with(external_counts, sum(n[author_cell_type == "Macrophages" & subtype == "cDC2"]) /
                        sum(n[subtype == "cDC2"]))
apc_inventory <- read_table("12_cohort_inventory")
author_colours <- c("Dendritic.cells" = "#376681", "Macrophages" = "#B05C40")
apc_plot <- ggplot(external_counts, aes(subtype, n, fill = author_cell_type)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.69) +
  geom_text(aes(label = n), position = position_dodge(width = 0.76), vjust = -0.45, size = 3.1) +
  scale_fill_manual(values = author_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Author labels and gated DC subtypes",
       subtitle = sprintf("%s APCs in %d oral tumours\n%.1f%% of gated cDC2 from Macrophages",
                          count_label(sum(apc_inventory$n_myeloid_apc)), nrow(apc_inventory),
                          100 * mac_cdc2_share),
       x = "2-of-3 gate", y = "Cells", fill = "Author label") +
  manuscript_theme(9) + theme(panel.grid.major.x = element_blank())

external_focus <- c("primary", "hybrid_author_migDC_vs_cd163neg_cDC2")
external_headings <- setNames(c("Myeloid APC pool", "Author migDC / CD163-negative cDC2"), external_focus)
replication_summary <- read_table("12_cross_cohort_summary") %>%
  filter(analysis %in% external_focus) %>%
  mutate(panel = sprintf("%s (n = %d)\nrho = %.3f  |  %.1f%% of discovery hits concordant",
                         external_headings[analysis], n_donors, rho, 100 * concordant_hit_rate))
replication_panels <- replication_summary$panel[match(external_focus, replication_summary$analysis)]
replication_data <- read_table("12_cross_cohort_genes") %>%
  filter(analysis %in% external_focus, !annotation_gene) %>%
  left_join(select(replication_summary, analysis, panel), by = "analysis") %>%
  mutate(panel = factor(panel, levels = replication_panels),
         replicated = concordant & fdr_discovery < cfg$fdr & fdr_external < cfg$fdr) %>%
  arrange(replicated)
replication_plot <- ggplot(replication_data, aes(log2_fc_discovery, log2_fc_external)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_point(aes(colour = replicated), alpha = 0.48, size = 0.9) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = c("FALSE" = "#ADB9C0", "TRUE" = migdc_colour),
                      labels = c("Other shared genes", "Concordant with FDR < 0.05 in both")) +
  labs(title = "Cross-cohort migDC vs cDC2 replication",
       subtitle = "Shared genes outside the assignment gene set",
       x = "Discovery log2 fold change", y = "GSE181919 log2 fold change", colour = NULL) +
  manuscript_theme(9)

cohort_models <- c(primary = "Primary, 10 cells", min_5 = "Five-cell minimum",
                   cd163_neg_min5 = "CD163-negative, five cells")
shared <- read_table("12_shared_program_donor_effects") %>%
  filter(analysis %in% names(cohort_models))
joint <- read_table("12_combined_program_summary") %>%
  filter(analysis %in% names(cohort_models), status == "fitted")
cohort_summary <- bind_rows(lapply(split(shared, list(shared$analysis, shared$cohort), drop = TRUE), function(rows) {
  cbind(analysis = rows$analysis[1], cohort = rows$cohort[1], summarize_programs(rows))
}))
model_labels <- joint %>%
  distinct(analysis, n_discovery, n_external, n_donors) %>%
  mutate(panel = sprintf("%s\nDiscovery %d  |  GSE181919 %d  |  Combined %d",
                         cohort_models[analysis], n_discovery, n_external, n_donors))
model_order <- model_labels$panel[match(names(cohort_models), model_labels$analysis)]
cohort_offsets <- c(discovery = 0.25, GSE181919 = 0.04, combined = -0.17)
cohort_labels <- c(discovery = "Discovery", GSE181919 = "GSE181919", combined = "Combined")
cohort_colours <- c(Discovery = "#356782", GSE181919 = "#B05C40", Combined = "#117B75")
comparison_points <- bind_rows(
  cohort_summary,
  mutate(joint, cohort = "combined")
) %>%
  left_join(select(model_labels, analysis, panel), by = "analysis") %>%
  mutate(panel = factor(panel, levels = model_order),
         row = match(program, rev(program_order)) + cohort_offsets[cohort],
         cohort = factor(cohort, levels = names(cohort_labels), labels = unname(cohort_labels)))
joint_labels <- joint %>%
  left_join(select(model_labels, analysis, panel), by = "analysis") %>%
  mutate(panel = factor(panel, levels = model_order),
         row = match(program, rev(program_order)) - 0.40,
         label = sprintf("Flip P %s  |  t P %s  |  Fisher FDR %s",
                         p_label(p_sign_flip), p_label(p_paired_t), p_label(fdr_fry_meta)))
label_x <- min(comparison_points$lower_95_t, na.rm = TRUE)
cohort_program_plot <- ggplot(comparison_points, aes(mean_effect, row, colour = cohort)) +
  geom_hline(yintercept = seq_len(length(program_order) - 1) + 0.5, colour = "grey91") +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_errorbar(aes(xmin = lower_95_t, xmax = upper_95_t), orientation = "y",
                width = 0.095, linewidth = 0.53) +
  geom_point(aes(shape = cohort), size = 2.55) +
  geom_text(data = joint_labels, aes(x = label_x, y = row, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.25, colour = "#405461") +
  facet_wrap(~ panel, nrow = 1) +
  scale_y_continuous(breaks = seq_along(program_order), labels = rev(program_order),
                     limits = c(0.48, length(program_order) + 0.39)) +
  scale_colour_manual(values = cohort_colours) +
  scale_shape_manual(values = c(Discovery = 1, GSE181919 = 15, Combined = 18)) +
  labs(title = "Shared-gene program effects across cohorts",
       subtitle = "Means and 95% t intervals. P and FDR labels are for the combined cohort.",
       x = "Paired log2 CPM difference (migDC - cDC2)", y = NULL, colour = NULL, shape = NULL) +
  manuscript_theme(9) + theme(panel.grid.major.y = element_blank())
save_composite(((apc_plot | replication_plot) + plot_layout(widths = c(1, 2.13))) /
                 cohort_program_plot + plot_layout(heights = c(1, 1.15)),
               "Fig5_external_validation", "Independent external validation", 13.7, 10.1)

log_msg("STEP 14 done. Elapsed seconds ", sprintf("%.2f", proc.time()[["elapsed"]] - started))

# Patient-level tumour microenvironment covariation across immune populations.
source("R/00_setup.R")
source("R/pseudobulk.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

log_msg("STEP 13 -- tumour microenvironment covariation")
cfg <- PARAMS$pseudobulk
samples <- select(PARAMS$samples, sample, patient, dataset, group)
tnk_types <- c("T cells", "Treg", "NK cells")
immune_types <- c(tnk_types, "B cells", "Plasma cells", "Macrophages", "DC", "pDC",
                  "Mast cells", "Neutrophils")
profile_types <- c(dc = "Retained DC", migdc = "migDC", treg = "Treg",
                   tcell = "T cells", mac = "Macrophages")
migdc_min_cells <- min(cfg$sensitivity_min_cells)
sensitivity_id <- paste0("min_", migdc_min_cells)

atlas <- read_atlas_object("03_annotated.rds")
atlas_md <- atlas@meta.data
stopifnot(!anyNA(atlas_md$cell_type), setequal(atlas_md$sample, samples$sample))
atlas_sizes <- atlas_md %>%
  group_by(sample) %>%
  summarise(n_total = n(),
            n_immune = sum(cell_type %in% immune_types),
            n_tnk = sum(cell_type %in% tnk_types),
            n_treg = sum(cell_type == "Treg"),
            n_tcell = sum(cell_type == "T cells"),
            n_macrophages = sum(cell_type == "Macrophages"), .groups = "drop")

# Match the counts layer by cell name before aggregating each tumour population.
atlas_counts <- SeuratObject::LayerData(atlas, assay = "RNA", layer = "counts")
immune_md <- atlas_md[colnames(atlas_counts), , drop = FALSE]
immune_cells <- immune_md$group == "tumour" &
  immune_md$cell_type %in% profile_types[c("treg", "tcell", "mac")]
immune_md <- mutate(immune_md[immune_cells, , drop = FALSE], dc_subtype = cell_type)
immune_bulk <- aggregate_umi(atlas_counts[, immune_cells, drop = FALSE], immune_md)
rm(atlas, atlas_counts, immune_md)
invisible(gc())

dc <- read_atlas_object("05_dc_annotated.rds")
dc_counts <- SeuratObject::LayerData(dc, assay = "RNA", layer = "counts")
dc_md <- dc@meta.data[colnames(dc_counts), , drop = FALSE]
stopifnot(!anyNA(dc_md$dc_subtype),
          identical(as.character(dc_md$sample), as.character(atlas_md[rownames(dc_md), "sample"])))
dc_sizes <- dc_md %>%
  group_by(sample) %>%
  summarise(n_dc_retained = n(),
            n_migdc = sum(dc_subtype == "migDC"),
            n_cdc2 = sum(dc_subtype == "cDC2"),
            n_cdc1 = sum(dc_subtype == "cDC1"),
            n_lclike = sum(dc_subtype == "LC-like"), .groups = "drop")

proportion <- function(numerator, denominator) {
  if_else(denominator > 0L, numerator / denominator, NA_real_)
}
composition <- samples %>%
  left_join(atlas_sizes, by = "sample") %>%
  left_join(dc_sizes, by = "sample") %>%
  mutate(across(starts_with("n_"), ~ coalesce(.x, 0L)),
         migdc_of_dc = proportion(n_migdc, n_dc_retained),
         migdc_of_immune = proportion(n_migdc, n_immune),
         cdc2_of_dc = proportion(n_cdc2, n_dc_retained),
         cdc1_of_dc = proportion(n_cdc1, n_dc_retained),
         lclike_of_dc = proportion(n_lclike, n_dc_retained),
         treg_of_tnk = proportion(n_treg, n_tnk),
         treg_of_immune = proportion(n_treg, n_immune),
         tcell_of_immune = proportion(n_tcell, n_immune),
         mac_of_immune = proportion(n_macrophages, n_immune))
log_msg("  ", sum(composition$n_total), " atlas cells and ", sum(composition$n_dc_retained),
        " retained DCs across ", nrow(composition), " libraries")

dc_bulk <- aggregate_umi(dc_counts, dc_md)
dc_total_bulk <- aggregate_umi(dc_counts, mutate(dc_md, dc_subtype = profile_types[["dc"]]))
ligand_genes <- c("CCL22", "CCL17", "CCL19", "CD274", "PDCD1LG2", "IDO1", "CD80", "CD86", "CD200")
inhibitory_genes <- c("PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2")
profile_genes <- list(dc = ligand_genes, migdc = ligand_genes, treg = c("CCR4", "FOXP3", "CTLA4"),
                      tcell = c("CCR7", inhibitory_genes), mac = ligand_genes)
profile_bulks <- list(dc = dc_total_bulk, migdc = dc_bulk, treg = immune_bulk,
                      tcell = immune_bulk, mac = immune_bulk)

log_cpm_profile <- function(bulk, cell_type, genes, prefix, min_cells) {
  md <- subset(bulk$metadata, group == "tumour" & dc_subtype == cell_type & n_cells >= min_cells)
  if (!nrow(md)) stop("No available tumour pseudobulks for ", cell_type, ".")
  missing <- setdiff(genes, rownames(bulk$counts))
  if (length(missing)) stop(cell_type, " counts lack genes: ", paste(missing, collapse = ", "))
  # Estimate TMM factors on the complete transcriptome within this population.
  y <- edgeR::DGEList(as.matrix(bulk$counts[, rownames(md), drop = FALSE]), samples = md)
  y <- edgeR::normLibSizes(y, method = "TMM")
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  values <- t(log_cpm[genes, , drop = FALSE])
  colnames(values) <- paste(prefix, genes, sep = "_")
  log_msg("  ", cell_type, " expression: ", nrow(md), " tumour libraries")
  data.frame(sample = md$sample, values, row.names = NULL, check.names = FALSE)
}

profiles <- filter(composition, group == "tumour")
for (prefix in names(profile_types)) {
  features <- log_cpm_profile(profile_bulks[[prefix]], profile_types[[prefix]],
                              profile_genes[[prefix]], prefix,
                              min_cells = if (prefix == "migdc") migdc_min_cells else 1L)
  profiles <- left_join(profiles, features, by = "sample")
}
profiles$tcell_inhibitory_module <- rowMeans(profiles[paste0("tcell_", inhibitory_genes)])

# Reuse exported donor program scores from 11_donor_pseudobulk.R when present;
# compute sensitivity-threshold scores directly if only primary scores were exported.
program_keys <- c("Maturation" = "maturation", "Regulatory expression" = "regulatory_expression",
                  "NF-kB feedback" = "nfkb_feedback")
programs <- read_dc_programs()
stopifnot(all(names(program_keys) %in% names(programs)))
programs <- programs[names(program_keys)]
program_thresholds <- setNames(c(cfg$primary_min_cells, migdc_min_cells),
                               c("primary", sensitivity_id))
published_effects <- read.csv(file.path(TAB_DIR, "11_program_donor_effects.csv"))

load_program_effects <- function(analysis) {
  eligibility <- paired_cell_counts(dc_bulk, samples, "migDC", "cDC2", program_thresholds[[analysis]])
  used <- subset(eligibility, group == "tumour" & eligible)
  effects <- filter(published_effects, .data$analysis == .env$analysis, program %in% names(programs))
  if (!nrow(effects) && analysis == sensitivity_id) {
    stopifnot(nrow(used) >= cfg$min_donors)
    log_msg("  Computing ", analysis, " donor program scores for ",
            nrow(used), " matched tumour donors")
    fit <- fit_paired_bulk(paired_bulk(dc_bulk, used, "migDC", "cDC2"))
    effects <- program_differences(paired_effects(fit$y, "migDC", "cDC2"), programs,
                                   cfg$minimum_program_coverage) %>%
      left_join(select(used, patient, dataset), by = "patient")
  }
  expected <- expand_grid(patient = used$patient, program = names(programs)) %>%
    left_join(select(used, patient, dataset), by = "patient")
  keys <- c("patient", "dataset", "program")
  stopifnot(nrow(effects) == nrow(expected), !anyDuplicated(effects[keys]),
            nrow(anti_join(expected, effects, by = keys)) == 0L)
  transmute(effects, patient, dataset, n_genes, n_target_genes,
            feature = paste(.env$analysis, program_keys[program], sep = "_"), effect)
}

program_effects <- bind_rows(lapply(names(program_thresholds), load_program_effects))
coverage_gaps <- program_effects %>%
  distinct(feature, n_genes, n_target_genes) %>%
  filter(n_genes < ceiling(cfg$minimum_program_coverage * n_target_genes))
coverage_labels <- with(coverage_gaps, setNames(sprintf("Coverage %d/%d", n_genes, n_target_genes), feature))
for (feature in names(coverage_labels)) {
  log_msg("  ", feature, " has ", tolower(coverage_labels[[feature]]),
          " and requires ", scales::percent(cfg$minimum_program_coverage), " coverage.")
}
program_profiles <- pivot_wider(program_effects, id_cols = c(patient, dataset),
                                names_from = feature, values_from = effect)
profiles <- left_join(profiles, program_profiles, by = c("patient", "dataset"))
stopifnot(!anyDuplicated(profiles$patient))

# Include the prespecified composition, ligand/receptor and paired-program axes.
direct_pairs <- bind_rows(
  "Chemotaxis and infiltration" = tribble(
    ~predictor,           ~response,
    "migdc_of_dc",        "treg_of_tnk",
    "migdc_of_immune",    "treg_of_immune",
    "migdc_of_immune",    "treg_of_tnk",
    "dc_CCL22",           "treg_CCR4",
    "dc_CCL17",           "treg_CCR4",
    "migdc_CCL22",        "treg_CCR4",
    "migdc_CCL17",        "treg_CCR4",
    "dc_CCL19",           "tcell_CCR7",
    "migdc_CCL19",        "tcell_CCR7"
  ),
  "Checkpoint and T-cell state" = tribble(
    ~predictor,           ~response,
    "dc_CD274",           "tcell_PDCD1",
    "migdc_CD274",        "tcell_PDCD1",
    "dc_PDCD1LG2",        "tcell_PDCD1",
    "migdc_PDCD1LG2",     "tcell_PDCD1",
    "dc_IDO1",            "treg_FOXP3",
    "dc_IDO1",            "treg_CTLA4",
    "migdc_IDO1",         "treg_FOXP3",
    "migdc_IDO1",         "treg_CTLA4",
    "dc_CD80",           "treg_CTLA4",
    "dc_CD86",           "treg_CTLA4",
    "migdc_CD80",        "treg_CTLA4",
    "migdc_CD86",        "treg_CTLA4",
    "dc_CD200",          "tcell_inhibitory_module"
  ),
  .id = "axis"
) %>%
  mutate(analysis = "tumour")
program_pairs <- expand_grid(
  analysis = names(program_thresholds), program = names(program_keys),
  response = c("treg_of_tnk", "treg_CCR4", "treg_FOXP3", "tcell_inhibitory_module")
) %>%
  transmute(axis = "Paired DC programs", analysis,
            predictor = paste(analysis, program_keys[program], sep = "_"), response)
pairs <- bind_rows(direct_pairs, program_pairs)
stopifnot(!anyDuplicated(pairs[c("predictor", "response")]),
          all(c(pairs$predictor, pairs$response) %in% names(profiles)))

spearman <- function(predictor, response) {
  x <- profiles[[predictor]]
  y <- profiles[[response]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  # The asymptotic test needs at least three complete pairs and two varying profiles.
  n_donors <- length(x)
  if (n_donors < 3L) {
    return(data.frame(n_donors, rho = NA_real_, p_value = NA_real_, status = "insufficient_donors"))
  }
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(n_donors, rho = NA_real_, p_value = NA_real_, status = "constant_feature"))
  }
  fit <- cor.test(x, y, method = "spearman", alternative = "two.sided", exact = FALSE)
  data.frame(n_donors, rho = unname(fit$estimate), p_value = fit$p.value, status = "tested")
}

feature_labels <- c(migdc_of_dc = "migDC / retained DC", migdc_of_immune = "migDC / immune",
                    treg_of_tnk = "Treg / T and NK cells", treg_of_immune = "Treg / immune",
                    tcell_inhibitory_module = "T-cell inhibitory module")
for (prefix in names(profile_genes)) {
  genes <- profile_genes[[prefix]]
  feature_labels[paste(prefix, genes, sep = "_")] <- paste(profile_types[[prefix]], genes)
}
for (analysis in names(program_thresholds)) {
  feature_labels[paste(analysis, program_keys, sep = "_")] <- names(program_keys)
}

correlations <- bind_cols(pairs, bind_rows(Map(spearman, pairs$predictor, pairs$response))) %>%
  mutate(predictor_label = unname(feature_labels[predictor]),
         response_label = unname(feature_labels[response]),
         status = if_else(predictor %in% coverage_gaps$feature, "insufficient_program_coverage", status),
         fdr = p.adjust(p_value, "BH")) %>%
  select(axis, analysis, predictor, response, predictor_label, response_label,
         n_donors, rho, p_value, fdr, status)

# Plot all tests with finite estimates and four prespecified patient-level examples.
rho_data <- correlations %>%
  mutate(panel = if_else(axis == "Paired DC programs", paste0(axis, " (", analysis, ")"), axis),
         panel = factor(panel, levels = unique(panel)),
         pair_label = paste(predictor_label, "vs", response_label),
         pair_id = paste(predictor, response, sep = "/"),
         pair_id = factor(pair_id, levels = rev(pair_id)))
rho_labels <- setNames(rho_data$pair_label, rho_data$pair_id)
missing_axes <- filter(rho_data, status != "tested") %>%
  mutate(note = case_when(
    status == "insufficient_program_coverage" ~ unname(coverage_labels[predictor]),
    status == "constant_feature" ~ "Constant profile",
    TRUE ~ "Insufficient pairs"
  ))
rho_plot <- ggplot(rho_data, aes(rho, pair_id)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_col(fill = "#476B81", width = 0.68, na.rm = TRUE) +
  geom_point(data = filter(rho_data, rho == 0), colour = "#476B81", size = 1.4) +
  geom_text(data = missing_axes, aes(x = 0.05, label = note), hjust = 0, colour = "grey40", size = 2.6) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  scale_y_discrete(labels = function(pair_id) rho_labels[pair_id]) +
  labs(x = "Spearman rho", y = NULL, title = "Covariation across tumour patients") +
  theme_pub(10) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold", hjust = 0))

scatter_plot <- function(x, y, title, x_label, y_label, percent = FALSE) {
  result <- filter(correlations, predictor == .env$x, response == .env$y)
  points <- profiles %>%
    transmute(dataset, x_value = .data[[x]], y_value = .data[[y]]) %>%
    filter(is.finite(x_value), is.finite(y_value))
  p <- ggplot(points, aes(x_value, y_value)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey55", linewidth = 0.55) +
    geom_point(aes(colour = dataset), size = 2.6) +
    scale_colour_manual(values = c(GSE172577 = "#117B75", GSE193766 = "#8D4F91")) +
    labs(title = title, subtitle = sprintf("Spearman rho = %.2f | n = %d | BH FDR = %.3g",
                                         result$rho, result$n_donors, result$fdr),
         x = x_label, y = y_label, colour = "GEO series") +
    theme_pub(10)
  if (percent) {
    p <- p + scale_x_continuous(labels = scales::label_percent()) +
      scale_y_continuous(labels = scales::label_percent())
  }
  p
}

scatter_plots <- list(
  scatter_plot("migdc_of_dc", "treg_of_tnk", "Cell composition",
               "migDC / retained DC", "Treg / T and NK cells", percent = TRUE),
  scatter_plot("dc_CCL17", "treg_CCR4", "Chemotaxis",
               "Retained DC CCL17 (log2 CPM)", "Treg CCR4 (log2 CPM)"),
  scatter_plot("migdc_CD274", "tcell_PDCD1", "Checkpoint and T-cell state",
               "migDC CD274 (log2 CPM)", "T cells PDCD1 (log2 CPM)"),
  scatter_plot("primary_regulatory_expression", "tcell_inhibitory_module",
               "Paired regulatory program (primary)",
               "Regulatory expression (paired log2 CPM difference)",
               "T-cell inhibitory module\n(mean log2 CPM)")
)
figure <- rho_plot + wrap_plots(scatter_plots, ncol = 1) +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(title = "Tumour microenvironment covariation in OSCC",
                  caption = "Each point is a tumour patient. Grey lines show pooled linear trends.",
                  tag_levels = "A")
figure <- figure & theme(legend.position = "bottom")

write.csv(composition, file.path(TAB_DIR, "13_tme_cell_proportions.csv"), row.names = FALSE)
write.csv(profiles, file.path(TAB_DIR, "13_tme_donor_profiles.csv"), row.names = FALSE)
write.csv(correlations, file.path(TAB_DIR, "13_tme_ligand_receptor_correlations.csv"), row.names = FALSE)
save_fig(figure, "13_tme_covariation", width = 14, height = 14)
log_msg("STEP 13 done. ", nrow(profiles), " tumour profiles, ",
        sum(correlations$status == "tested"), "/", nrow(correlations), " axes tested, ",
        sum(correlations$fdr < cfg$fdr, na.rm = TRUE), " with BH FDR < ", cfg$fdr)

# Donor-paired pseudobulk differential expression and sensitivity analyses.
source("R/00_setup.R")
source("R/pseudobulk.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

log_msg("STEP 11 -- donor-paired expression")
cfg <- PARAMS$pseudobulk
programs <- read_dc_programs()
gates <- read.csv(file.path("config", "external_dc_gates.csv"))
annotation <- read.csv(file.path("config", "dc_programs.csv"))
assignment_genes <- union(gates$gene, annotation$gene[annotation$role == "annotation"])
stopifnot(!any(gates$gene %in% unlist(programs)))
panel <- c("LAMP3", "CCR7", "CCL19", "CCL22", "CD274", "IDO1", "CD1C", "CLEC10A",
           "CLEC9A", "XCR1", "CD83", "CD86", "PDCD1LG2", "CD200", "NFKBIA", "TNFAIP3")
object_dir <- Sys.getenv("OSCC_OBJECT_DIR", unset = OBJ_DIR)
input_name <- "05_dc_annotated.rds"
input_file <- atlas_object_path(input_name, object_dir)
dc <- read_atlas_object(input_name, dirname(input_file))
counts <- SeuratObject::LayerData(dc, assay = "RNA", layer = "counts")
md <- dc@meta.data[colnames(counts), , drop = FALSE]
bulk <- aggregate_umi(counts, md)
samples <- PARAMS$samples

contrasts <- bind_rows(
  data.frame(analysis = "primary", positive = "migDC", negative = "cDC2",
             min_cells = cfg$primary_min_cells, norm = "TMM"),
  data.frame(analysis = paste0("min_", cfg$sensitivity_min_cells), positive = "migDC",
             negative = "cDC2", min_cells = cfg$sensitivity_min_cells, norm = "TMM"),
  data.frame(analysis = "normalization", positive = "migDC", negative = "cDC2",
             min_cells = cfg$primary_min_cells, norm = "TMMwsp"),
  data.frame(analysis = "strict_filter_min10", positive = "migDC", negative = "cDC2",
             min_cells = cfg$primary_min_cells, norm = "TMM",
             min_count = cfg$strict_min_count, min_total_count = cfg$strict_min_total_count),
  data.frame(analysis = c("migDC_vs_cDC1", "cDC1_vs_cDC2"),
             positive = c("migDC", "cDC1"), negative = c("cDC1", "cDC2"),
             min_cells = cfg$secondary_min_cells, norm = "TMM"),
  data.frame(analysis = c("migDC_vs_LClike", "LClike_vs_cDC2"),
             positive = c("migDC", "LC-like"), negative = c("LC-like", "cDC2"),
             min_cells = cfg$primary_min_cells, norm = "TMM")
) %>%
  mutate(min_count = coalesce(min_count, cfg$min_count),
         min_total_count = coalesce(min_total_count, cfg$min_total_count))

run_contrast <- function(spec) {
  eligibility <- paired_cell_counts(bulk, samples, spec$positive, spec$negative, spec$min_cells)
  used <- subset(eligibility, group == "tumour" & eligible)
  if (nrow(used) < cfg$min_donors) {
    return(list(spec = spec, eligibility = eligibility, used = used, fit = NULL))
  }
  pair <- paired_bulk(bulk, used, spec$positive, spec$negative)
  fit <- fit_paired_bulk(pair, norm = spec$norm,
                        min_count = spec$min_count, min_total_count = spec$min_total_count)
  list(spec = spec, eligibility = eligibility, used = used, fit = fit)
}
runs <- setNames(lapply(seq_len(nrow(contrasts)), function(i) run_contrast(contrasts[i, ])),
                 contrasts$analysis)
primary <- runs$primary
if (is.null(primary$fit)) stop("Primary contrast has fewer than ", cfg$min_donors, " eligible donors.")

leave_out <- lapply(primary$used$patient, function(omitted) {
  used <- subset(primary$used, patient != omitted)
  pair <- paired_bulk(bulk, used, "migDC", "cDC2")
  fit_paired_bulk(pair)$de
})
names(leave_out) <- primary$used$patient
primary_de <- primary$fit$de
loo <- vapply(leave_out, function(de) de$log2_fc[match(primary_de$gene, de$gene)],
              numeric(nrow(primary_de)))
colnames(loo) <- paste0("without_", names(leave_out))
primary_de$n_loo_tested <- rowSums(!is.na(loo))
primary_de$n_loo_same_sign <- rowSums(sign(loo) == sign(primary_de$log2_fc), na.rm = TRUE)
primary_de <- cbind(primary_de, loo)

wilcox <- read.csv(file.path(TAB_DIR, "05_de_migDC_vs_cDC2.csv"))
at <- match(primary_de$gene, wilcox$gene)
primary_de$cell_level_log2_fc <- wilcox$avg_log2FC[at]
primary_de$cell_level_p_adj <- wilcox$p_val_adj[at]

summary <- bind_rows(lapply(runs, function(run) {
  de <- run$fit$de
  same_contrast <- run$spec$positive == "migDC" && run$spec$negative == "cDC2"
  common <- if (same_contrast) match(primary_de$gene, de$gene) else rep(NA_integer_, nrow(primary_de))
  joined <- !is.na(common)
  sig <- de$fdr < cfg$fdr
  cbind(run$spec, data.frame(
    n_donors = nrow(run$used), donors = paste(run$used$patient, collapse = " "),
    status = if (is.null(de)) "insufficient_donors" else "fitted",
    n_genes = if (is.null(de)) NA_integer_ else nrow(de),
    n_fdr = if (is.null(de)) NA_integer_ else sum(sig),
    n_fdr_positive = if (is.null(de)) NA_integer_ else sum(sig & de$log2_fc > 0),
    n_fdr_negative = if (is.null(de)) NA_integer_ else sum(sig & de$log2_fc < 0),
    n_shared_with_primary = if (same_contrast && !is.null(de)) sum(joined) else NA_integer_,
    rho_with_primary = if (sum(joined) > 1L) {
      cor(primary_de$log2_fc[joined], de$log2_fc[common[joined]], method = "spearman")
    } else NA_real_
  ))
}))
loo_summary <- bind_rows(lapply(names(leave_out), function(donor) {
  de <- leave_out[[donor]]
  at <- match(primary_de$gene, de$gene)
  common <- !is.na(at)
  data.frame(omitted_donor = donor, n_remaining = nrow(primary$used) - 1L,
             n_genes = nrow(de), n_fdr = sum(de$fdr < cfg$fdr), n_shared = sum(common),
             rho_with_primary = cor(primary_de$log2_fc[common], de$log2_fc[at[common]],
                                    method = "spearman"))
}))

focal <- runs[c("primary", "migDC_vs_cDC1", "cDC1_vs_cDC2", "migDC_vs_LClike", "LClike_vs_cDC2")]
focal <- Filter(function(run) !is.null(run$fit), focal)
donor_effects <- lapply(focal, function(run) {
  paired_effects(run$fit$y, run$spec$positive, run$spec$negative)
})
program_effects <- bind_rows(lapply(names(focal), function(id) {
  out <- program_differences(donor_effects[[id]], programs, cfg$minimum_program_coverage)
  out$analysis <- id
  left_join(out, select(focal[[id]]$used, patient, dataset), by = "patient")
}))
program_fry <- bind_rows(lapply(focal, function(run) fry_programs(run$fit$y, programs)),
                         .id = "analysis")
program_summary <- bind_rows(lapply(split(program_effects, program_effects$analysis), function(effects) {
  cbind(analysis = effects$analysis[1], summarize_programs(effects))
})) %>%
  left_join(select(program_fry, -n_genes_fry), by = c("analysis", "program"))
panel_effects <- bind_rows(lapply(names(focal), function(id) {
  gene_effects_long(donor_effects[[id]], panel) %>%
    mutate(analysis = id) %>%
    left_join(select(focal[[id]]$used, patient, dataset), by = "patient")
}))
program_gene_effects <- gene_effects_long(donor_effects$primary, unlist(programs))
program_by_cohort <- program_effects %>%
  filter(analysis == "primary") %>%
  group_by(dataset, program) %>%
  summarise(n_donors = sum(is.finite(effect)), mean_effect = mean(effect),
            positive_donors = sum(effect > 0, na.rm = TRUE), .groups = "drop")
panel_fits <- bind_rows(lapply(runs, function(run) {
  de <- run$fit$de
  at <- match(panel, de$gene)
  data.frame(analysis = run$spec$analysis, gene = panel,
             log2_fc = if (is.null(de)) NA_real_ else de$log2_fc[at],
             fdr = if (is.null(de)) NA_real_ else de$fdr[at])
}))
focus_detection <- bind_rows(lapply(intersect(c("IDO1", "LAMP3", "CD274", "PDCD1LG2"), rownames(counts)),
                                   function(gene) {
  mutate(md, detected = as.numeric(counts[gene, ]) > 0) %>%
    filter(group == "tumour", patient %in% primary$used$patient,
           dc_subtype %in% c("cDC1", "cDC2", "migDC")) %>%
    group_by(patient, dataset, dc_subtype) %>%
    summarise(n = n(), n_detected = sum(detected), fraction_detected = mean(detected),
              .groups = "drop") %>%
    mutate(gene = gene)
}))

# Match subtype cell numbers within each eligible donor. These repeated subsamples
# quantify sensitivity to cell-number imbalance, not technical or biological sampling CIs.
cells_by_pair <- split(rownames(md), paste(md$sample, md$dc_subtype, sep = "/"))
pair_keys <- unlist(lapply(primary$used$sample, function(s) paste(s, c("migDC", "cDC2"), sep = "/")))
balanced_n <- rep(pmin(primary$used$n_positive, primary$used$n_negative), each = 2L)
primary_genes <- rownames(primary$fit$y)
set.seed(SEED)
balanced <- bind_rows(lapply(seq_len(cfg$balanced_repeats), function(iteration) {
  cells <- unlist(Map(function(key, n) sample(cells_by_pair[[key]], n), pair_keys, balanced_n))
  selected <- aggregate_umi(counts[primary_genes, cells, drop = FALSE], md[cells, , drop = FALSE])
  paired <- paired_bulk(selected, primary$used, "migDC", "cDC2")
  y <- edgeR::normLibSizes(edgeR::DGEList(as.matrix(paired$counts), samples = paired$metadata))
  effects <- paired_effects(y, "migDC", "cDC2")
  scores <- program_differences(effects, programs, cfg$minimum_program_coverage) %>%
    group_by(program) %>% summarise(effect = mean(effect), .groups = "drop")
  genes <- intersect(panel, rownames(effects))
  bind_rows(scores, data.frame(program = genes, effect = rowMeans(effects[genes, , drop = FALSE]))) %>%
    mutate(iteration = iteration)
}))
balanced_summary <- balanced %>%
  group_by(program) %>%
  summarise(n_runs = sum(!is.na(effect)),
            median = median(effect, na.rm = TRUE),
            lower_025 = quantile(effect, 0.025, na.rm = TRUE),
            upper_975 = quantile(effect, 0.975, na.rm = TRUE),
            positive_runs = sum(effect > 0, na.rm = TRUE), .groups = "drop")

# Competitive Hallmark tests use the same paired design and exclude assignment
# genes from the tested universe. Keep TMM library sizes from the main fit.
gmt <- readLines(file.path("data", "genesets", "h.all.symbols.gmt"))
entries <- strsplit(gmt, "\t", fixed = TRUE)
hallmark <- setNames(lapply(entries, function(term) term[-c(1L, 2L)]),
                     vapply(entries, "[", character(1), 1L))
path_y <- primary$fit$y[!primary_genes %in% assignment_genes, , keep.lib.sizes = TRUE]
path_design <- model.matrix(~ patient + dc_subtype, path_y$samples)
sets <- lapply(hallmark, function(genes) which(rownames(path_y) %in% genes))
sets <- sets[lengths(sets) >= PARAMS$enrich$min_gs]
camera <- limma::camera(path_y, sets, design = path_design, contrast = ncol(path_design))
pathways <- data.frame(pathway = rownames(camera), camera,
                       original_set_size = lengths(hallmark)[rownames(camera)],
                       context = "migDC vs cDC2", row.names = NULL)

# Sensitivity to annotation uses the same gates as the external analysis.
# Restrict the alternative comparison to the conventional and migratory pool.
stopifnot(!any(c("CD274", "IDO1") %in% gates$gene))
gate_programs <- programs
gate_programs[["Regulatory expression (incl. CD274/IDO1)"]] <-
  c(programs[["Regulatory expression"]], "CD274", "IDO1")
gated <- mutate(md, marker_gate = classify_dc(counts, gates, n_detected = 2L))
label_agreement <- count(gated, dataset, group, dc_subtype, marker_gate, name = "n")
gate_md <- subset(gated, dc_subtype %in% unique(gates$subtype) & marker_gate %in% unique(gates$subtype))
gate_md$dc_subtype <- gate_md$marker_gate
gate_bulk <- aggregate_umi(counts[, rownames(gate_md), drop = FALSE], gate_md)
gate_eligibility <- bind_rows(lapply(c(cfg$primary_min_cells, cfg$secondary_min_cells), function(n) {
  paired_cell_counts(gate_bulk, samples, "migDC", "cDC2", n)
}))
gate_runs <- lapply(split(gate_eligibility, gate_eligibility$min_cells), function(eligibility) {
  used <- subset(eligibility, group == "tumour" & eligible)
  result <- data.frame(min_cells = eligibility$min_cells[1], n_donors = nrow(used),
                       n_genes = NA_integer_, n_shared = NA_integer_, rho_with_primary = NA_real_)
  if (nrow(used) < cfg$min_donors) return(list(summary = result))
  fit <- fit_paired_bulk(paired_bulk(gate_bulk, used, "migDC", "cDC2"))
  nongate <- subset(primary_de, !gene %in% assignment_genes)
  idx <- match(nongate$gene, fit$de$gene)
  ok <- !is.na(idx)
  result$n_genes <- nrow(fit$de)
  result$n_shared <- sum(ok)
  result$rho_with_primary <- cor(nongate$log2_fc[ok], fit$de$log2_fc[idx[ok]], method = "spearman")
  scores <- program_differences(paired_effects(fit$y, "migDC", "cDC2"),
                                gate_programs, cfg$minimum_program_coverage)
  program_results <- summarize_programs(scores) %>%
    left_join(fry_programs(fit$y, gate_programs), by = "program") %>%
    mutate(min_cells = result$min_cells)
  list(summary = result, programs = program_results)
})

write_result <- function(x, suffix) {
  write.csv(x, file.path(TAB_DIR, paste0("11_", suffix, ".csv")), row.names = FALSE)
}
write_result(bind_rows(lapply(runs, function(run) mutate(run$eligibility, analysis = run$spec$analysis))),
             "donor_eligibility")
write_result(primary_de, "de_migDC_vs_cDC2")
write_result(data.frame(file = input_name, md5 = unname(tools::md5sum(input_file)),
                        n_cells = ncol(counts), n_genes = nrow(counts)), "input_snapshot")
primary_pair <- paired_bulk(bulk, primary$used, "migDC", "cDC2")
primary_library <- data.frame(pseudobulk = rownames(primary_pair$metadata), primary_pair$metadata,
                              umi = Matrix::colSums(primary_pair$counts), row.names = NULL)
write_result(primary_library, "primary_pseudobulk_libraries")
con <- gzfile(file.path(TAB_DIR, "11_primary_pseudobulk_counts.tsv.gz"), "wt")
tryCatch(write.table(as.matrix(primary_pair$counts), con, sep = "\t",
                     quote = FALSE, col.names = NA), finally = close(con))
for (id in setdiff(names(focal), "primary")) write_result(focal[[id]]$fit$de, paste0("de_", id))
write_result(summary, "sensitivity")
write_result(loo_summary, "leave_one_donor_out")
write_result(panel_fits, "panel_model_results")
write_result(panel_effects, "panel_donor_effects")
write_result(program_gene_effects, "program_gene_donor_effects")
write_result(program_effects, "program_donor_effects")
write_result(program_fry, "program_fry")
write_result(program_summary, "program_summary")
write_result(program_by_cohort, "program_by_cohort")
write_result(focus_detection, "gene_detection")
write_result(balanced_summary, "balanced_cell_sensitivity")
write_result(pathways, "paired_hallmark_camera")
write_result(label_agreement, "label_agreement")
write_result(gate_eligibility, "marker_gate_eligibility")
write_result(bind_rows(lapply(gate_runs, "[[", "summary")), "marker_gate_sensitivity")
write_result(bind_rows(lapply(gate_runs, "[[", "programs")), "marker_gate_programs")
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(TAB_DIR, "11_sessionInfo.txt"))

plot_scores <- filter(program_effects, analysis == "primary", is.finite(effect))
plot_means <- filter(program_summary, analysis == "primary", is.finite(mean_effect)) %>%
  mutate(stat_label = sprintf("Mean %+.3f\nPaired t p = %.3g   fry FDR = %.3g",
                              mean_effect, p_paired_t, fdr_fry))
stat_labels <- setNames(plot_means$stat_label, plot_means$program)
p1 <- ggplot(plot_scores, aes(x = effect, y = program)) +
  geom_vline(xintercept = 0, colour = "grey75", linetype = "dashed") +
  geom_errorbar(data = plot_means, aes(y = program, xmin = lower_95_t, xmax = upper_95_t),
                inherit.aes = FALSE, orientation = "y", width = 0.12, colour = "#78838B") +
  geom_point(aes(colour = dataset), size = 2.9, alpha = 0.85,
             position = position_jitter(width = 0, height = 0.08, seed = SEED)) +
  geom_point(data = plot_means, aes(x = mean_effect),
             shape = 23, size = 3.6, fill = "white", colour = "black") +
  scale_colour_manual(values = c(GSE172577 = "#117B75", GSE193766 = "#8D4F91")) +
  scale_y_discrete(sec.axis = dup_axis(name = NULL, labels = function(program) stat_labels[program])) +
  labs(x = "Paired mean log2 CPM difference (migDC - cDC2)", y = NULL,
       colour = "GEO series", title = "Donor-paired DC programs",
       subtitle = paste(nrow(primary$used), "matched tumour donors. Mean and 95% t intervals shown.")) +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
                      axis.text.y.right = element_text(hjust = 0, size = 9),
                      axis.ticks.y.right = element_blank())
save_fig(p1, "11_paired_dc_programs", 10, 4.1)

plot_panel <- filter(panel_effects, gene %in% c("IDO1", "CD274", "PDCD1LG2", "CD200", "CD83", "TNFAIP3")) %>%
  mutate(analysis = factor(analysis,
                           levels = c("primary", "migDC_vs_cDC1", "migDC_vs_LClike",
                                      "cDC1_vs_cDC2", "LClike_vs_cDC2"),
                           labels = c("migDC - cDC2", "migDC - cDC1", "migDC - LC-like",
                                      "cDC1 - cDC2", "LC-like - cDC2")))
p2 <- ggplot(plot_panel, aes(log2_cpm_difference, gene)) +
  geom_vline(xintercept = 0, colour = "grey75", linetype = "dashed") +
  geom_point(aes(colour = dataset), size = 2.4, alpha = 0.8,
             position = position_jitter(width = 0, height = 0.08, seed = SEED)) +
  facet_wrap(~ analysis) +
  scale_colour_manual(values = c(GSE172577 = "#117B75", GSE193766 = "#8D4F91")) +
  labs(x = "Paired log2 CPM difference per donor", y = NULL, colour = "GEO series",
       title = "Gene expression across DC states") +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_fig(p2, "11_dc_lineage_comparators", 10, 4.5)
log_msg("STEP 11 done")
print(summary, row.names = FALSE)
print(program_summary, row.names = FALSE)

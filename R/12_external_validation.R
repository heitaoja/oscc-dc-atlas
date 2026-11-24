# External cohort validation in primary oral-cavity tumours from GSE181919.
source("R/00_setup.R")
source("R/pseudobulk.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

log_msg("STEP 12 -- GSE181919 oral-cavity validation")
cfg <- PARAMS$pseudobulk
external_dir <- file.path(RAW_DIR, "GSE181919")
inputs <- file.path(external_dir, c("GSE181919_Barcode_metadata.txt.gz", "GSE181919_UMI_counts.txt.gz"))
if (!file.exists(inputs[1])) stop("Missing GEO barcode metadata. See data/README.md.")
programs <- read_dc_programs()
gates <- read.csv(file.path("config", "external_dc_gates.csv"))
stopifnot(!anyDuplicated(gates), !any(gates$gene %in% unlist(programs)))
stopifnot(!any(c("CD274", "IDO1") %in% gates$gene))
gate_programs <- programs
gate_programs[["Regulatory expression (incl. CD274/IDO1)"]] <-
  c(programs[["Regulatory expression"]], "CD274", "IDO1")

write_result <- function(x, suffix) {
  write.csv(x, file.path(TAB_DIR, paste0("12_", suffix, ".csv")), row.names = FALSE)
}
authors <- read.delim(gzfile(inputs[1]), row.names = NULL)
stopifnot(!anyDuplicated(authors$row.names),
          all(c("patient.id", "sample.id", "tissue.type", "subsite", "cell.type") %in% names(authors)))
oral <- subset(authors, tissue.type == "CA" & subsite == "OC")
myeloid <- subset(oral, cell.type %in% c("Dendritic.cells", "Macrophages"))
sample_table <- unique(data.frame(sample = oral$sample.id, patient = paste0("Choi_", oral$patient.id),
                                  dataset = "GSE181919", group = "tumour"))
sample_table <- sample_table[order(sample_table$sample), ]
stopifnot(!anyDuplicated(sample_table$patient), !anyDuplicated(sample_table$sample))
myeloid_counts <- myeloid %>%
  group_by(sample.id) %>%
  summarise(n_author_dc = sum(cell.type == "Dendritic.cells"),
            n_author_macrophages = sum(cell.type == "Macrophages"),
            n_myeloid_apc = n(), .groups = "drop")
inventory <- sample_table %>%
  left_join(count(oral, sample.id, name = "n_cells"), by = c("sample" = "sample.id")) %>%
  left_join(myeloid_counts, by = c("sample" = "sample.id")) %>%
  mutate(across(c(n_author_dc, n_author_macrophages, n_myeloid_apc), ~ replace_na(.x, 0L)))
write_result(inventory, "cohort_inventory")
input_info <- data.frame(file = basename(inputs), bytes = file.size(inputs),
                         md5 = unname(tools::md5sum(inputs)), available = file.exists(inputs))
write_result(input_info, "inputs")
log_msg("  ", sum(inventory$n_cells), " cells across ", nrow(inventory), " oral tumours")
log_msg("  Myeloid APCs: ", sum(inventory$n_author_dc), " author DCs and ",
        sum(inventory$n_author_macrophages), " author macrophages")
if ("--inventory-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no")
if (!file.exists(inputs[2])) stop("Missing GEO UMI counts. Only the metadata inventory was generated.")
md <- sample_table[match(myeloid$sample.id, sample_table$sample), , drop = FALSE]
md$author_cell_type <- myeloid$cell.type
rownames(md) <- myeloid$row.names

# Read only relevant columns in bounded chunks. The source is a wide TSV whose
# decompressed size makes a dense read of all 54,239 cells impractical.
read_dc_columns <- function(path, barcode_ids, source_md5 = unname(tools::md5sum(path))) {
  canonical <- function(x) sub("\\.([0-9]+)$", "-\\1", x)
  ids <- canonical(barcode_ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Expected unique, nonempty myeloid APC barcodes.")
  }
  source_bytes <- file.size(path)
  if (is.na(source_bytes) || is.na(source_md5)) stop("Could not fingerprint GEO counts: ", path)
  cache_path <- file.path(OBJ_DIR, "12_gse181919_oral_myeloid_counts.rds")
  valid_counts <- function(counts) {
    inherits(counts, "dgCMatrix") && isTRUE(methods::validObject(counts, test = TRUE)) &&
      identical(colnames(counts), barcode_ids) && nrow(counts) > 0L &&
      !is.null(rownames(counts)) && !anyNA(rownames(counts)) &&
      all(nzchar(rownames(counts))) && !anyDuplicated(rownames(counts)) &&
      all(is.finite(counts@x) & counts@x >= 0 & counts@x == floor(counts@x))
  }
  if (file.exists(cache_path)) {
    cache <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.list(cache) && identical(cache$bytes, source_bytes) && identical(cache$md5, source_md5) &&
        identical(cache$barcodes, barcode_ids) && is.character(cache$source_barcodes) &&
        identical(canonical(cache$source_barcodes), ids) && valid_counts(cache$counts)) {
      log_msg("  Loaded ", ncol(cache$counts), " myeloid APC columns from ", basename(cache_path))
      return(cache$counts)
    }
    log_msg("  Rebuilding myeloid APC cache after a fingerprint or content mismatch.")
  }

  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L)
  first_gene <- readLines(con, n = 1L)
  if (!length(header_line) || !length(first_gene)) stop("GEO counts matrix has no gene rows.")
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  fields <- length(strsplit(first_gene, "\t", fixed = TRUE)[[1]])
  if (length(header) == fields) header <- header[-1L]
  header <- sub('^"(.*)"$', "\\1", header)
  header_ids <- canonical(header)
  at <- match(ids, header_ids)
  if (length(header) != fields - 1L || anyNA(at) || anyDuplicated(at) || anyDuplicated(header_ids)) {
    stop("Cell barcodes do not match uniquely between the UMI matrix and author metadata.")
  }
  log_msg("  Matched ", length(at), " myeloid APC barcodes in ", length(header), " GEO columns. ",
          sum(header[at] != barcode_ids), " matched after suffix normalization.")

  select_columns <- c(1L, at + 1L)
  expected_columns <- paste0("V", select_columns)
  chunk_rows <- 256L
  parts <- list()
  lines <- c(first_gene, readLines(con, n = chunk_rows - 1L))
  while (length(lines)) {
    chunk <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                               sep = "\t", select = select_columns, colClasses = list(character = 1L),
                               showProgress = FALSE, nThread = 2L)
    if (nrow(chunk) != length(lines) || !identical(names(chunk), expected_columns)) {
      stop("Unexpected number of rows or columns while reading GEO counts.")
    }
    parts[[length(parts) + 1L]] <- chunk
    lines <- readLines(con, n = chunk_rows)
  }
  tab <- data.table::rbindlist(parts)
  genes <- tab[[1L]]
  storage <- as.matrix(tab[, -1L])
  dimnames(storage) <- list(genes, barcode_ids)
  if (!is.numeric(storage)) stop("GEO UMI columns are not numeric.")
  counts <- methods::as(Matrix::Matrix(storage, sparse = TRUE), "generalMatrix")
  if (!valid_counts(counts)) {
    stop("Expected unique gene IDs and non-negative integer UMI counts.")
  }
  if (!identical(file.size(path), source_bytes) || !identical(unname(tools::md5sum(path)), source_md5)) {
    stop("The GEO counts file changed during extraction.")
  }
  temporary <- tempfile("12_gse181919_counts_", tmpdir = dirname(cache_path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(list(counts = counts, bytes = source_bytes, md5 = source_md5, barcodes = barcode_ids,
               source_barcodes = header[at]), temporary)
  if (!file.rename(temporary, cache_path)) stop("Could not publish validated myeloid APC cache: ", cache_path)
  log_msg("  Cached ", nrow(counts), " genes by ", ncol(counts), " myeloid APCs.")
  counts
}

log_msg("  Selecting ", nrow(myeloid), " myeloid APCs in ", nrow(sample_table), " oral tumours")
counts <- read_dc_columns(inputs[2], rownames(md), source_md5 = input_info$md5[2L])

fit_validation <- function(bulk, min_cells, positive, negative) {
  eligibility <- paired_cell_counts(bulk, sample_table, positive, negative, min_cells)
  chosen <- subset(eligibility, eligible)
  fit <- if (nrow(chosen) >= cfg$min_donors) {
    fit_paired_bulk(paired_bulk(bulk, chosen, positive, negative))
  } else NULL
  list(eligibility = eligibility, chosen = chosen, fit = fit)
}

main_labels <- classify_dc(counts, gates, n_detected = 2L)
loose_labels <- classify_dc(counts, gates, n_detected = 1L)
if (!"CD163" %in% rownames(counts)) stop("GEO UMI matrix is missing the CD163 purification marker.")
author_dc <- md$author_cell_type == "Dendritic.cells"
cd163_negative <- counts["CD163", ] == 0
cell_pools <- list(
  myeloid_apc = rep(TRUE, ncol(counts)),
  cd163_negative = cd163_negative,
  author_dc = author_dc,
  hybrid = (author_dc & main_labels == "migDC") | (cd163_negative & main_labels == "cDC2")
)
pool_bulks <- lapply(cell_pools, function(included) {
  cells <- included & main_labels %in% gates$subtype
  annotated <- md[cells, , drop = FALSE]
  annotated$dc_subtype <- main_labels[cells]
  if (any(cells)) aggregate_umi(counts[, cells, drop = FALSE], annotated) else {
    list(metadata = transform(annotated, n_cells = integer()))
  }
})

hybrid_id <- "hybrid_author_migDC_vs_cd163neg_cDC2"
specs <- tribble(
  ~analysis,                       ~pool,            ~min_cells,                     ~positive, ~negative,
  "primary",                       "myeloid_apc",    cfg$primary_min_cells,          "migDC",   "cDC2",
  "min_5",                         "myeloid_apc",    cfg$secondary_min_cells,        "migDC",   "cDC2",
  "cd163_neg_min10",               "cd163_negative", cfg$primary_min_cells,          "migDC",   "cDC2",
  "cd163_neg_min5",                "cd163_negative", cfg$secondary_min_cells,        "migDC",   "cDC2",
  hybrid_id,                       "hybrid",         cfg$primary_min_cells,          "migDC",   "cDC2",
  "author_dc_only_min10",          "author_dc",      cfg$primary_min_cells,          "migDC",   "cDC2",
  "author_dc_only_min5",           "author_dc",      cfg$secondary_min_cells,        "migDC",   "cDC2",
  "author_dc_migDC_vs_cDC1_min3",   "author_dc",      min(cfg$sensitivity_min_cells), "migDC",   "cDC1",
  "migDC_vs_cDC1",                 "myeloid_apc",    cfg$secondary_min_cells,        "migDC",   "cDC1",
  "cDC1_vs_cDC2",                  "myeloid_apc",    cfg$secondary_min_cells,        "cDC1",    "cDC2"
) %>%
  mutate(gate = 2L, .after = pool) %>%
  mutate(discovery_contrast = paste(positive, negative, sep = "_vs_"))
runs <- setNames(lapply(seq_len(nrow(specs)), function(i) {
  spec <- specs[i, ]
  fit_validation(pool_bulks[[spec$pool]], spec$min_cells, spec$positive, spec$negative)
}), specs$analysis)

classification <- bind_rows(
  mutate(md, subtype = main_labels, gate_min_detected = 2L),
  mutate(md, subtype = loose_labels, gate_min_detected = 1L)
) %>%
  count(sample, patient, author_cell_type, gate_min_detected, subtype, name = "n")
write_result(classification, "classification")
write_result(bind_rows(lapply(names(runs), function(id) {
  mutate(runs[[id]]$eligibility, analysis = id, pool = specs$pool[match(id, specs$analysis)])
})), "donor_eligibility")

validation_ids <- c("primary", "min_5")
combined_ids <- specs$analysis[specs$discovery_contrast == "migDC_vs_cDC2" & specs$pool != "author_dc"]
fitted_ids <- names(Filter(function(run) !is.null(run$fit), runs))
selected_id <- head(intersect(validation_ids, fitted_ids), 1L)
sensitivity <- bind_rows(lapply(names(runs), function(id) {
  run <- runs[[id]]
  de <- run$fit$de
  cbind(specs[specs$analysis == id, ], data.frame(
    n_donors = nrow(run$chosen), donors = paste(run$chosen$patient, collapse = " "),
    status = if (is.null(de)) "insufficient_donors" else "fitted",
    n_genes = if (is.null(de)) NA_integer_ else nrow(de),
    n_fdr = if (is.null(de)) NA_integer_ else sum(de$fdr < cfg$fdr),
    selected_for_validation = id %in% selected_id
  ))
}))
write_result(sensitivity, "sensitivity")
print(sensitivity, row.names = FALSE)
if (length(selected_id)) {
  log_msg("  Selected ", selected_id, " for independent validation.")
} else {
  log_msg("  Both preset validation thresholds fall below the minimum of ",
          cfg$min_donors, " matched donors.")
}

discovery_contrasts <- unique(specs$discovery_contrast)
discovery_de <- setNames(lapply(discovery_contrasts, function(contrast) {
  read.csv(file.path(TAB_DIR, paste0("11_de_", contrast, ".csv")))
}), discovery_contrasts)
discovery_programs <- read.csv(file.path(TAB_DIR, "11_program_gene_donor_effects.csv"))
stopifnot(!anyDuplicated(discovery_programs[c("gene", "patient")]))
discovery_programs <- as.matrix(xtabs(log2_cpm_difference ~ gene + patient, data = discovery_programs))
discovery_effects <- read.csv(file.path(TAB_DIR, "11_panel_donor_effects.csv"))
discovery_fry <- read.csv(file.path(TAB_DIR, "11_program_fry.csv")) %>%
  filter(analysis == "primary") %>%
  select(program, p_fry_discovery = p_fry, direction_fry_discovery = direction_fry)
annotation_genes <- union(gates$gene, read.csv(file.path("config", "dc_programs.csv")) %>%
                            filter(role == "annotation") %>% pull(gene))

cross_cohort <- list()
replication <- list()
program_results <- list()
program_fry <- list()
program_shared <- list()
panel_results <- list()
empty_de <- discovery_de[["migDC_vs_cDC2"]][0, c("gene", "log2_fc", "log_cpm", "ql_f", "p_value", "fdr")]
for (id in specs$analysis) {
  de <- runs[[id]]$fit$de
  suffix <- paste0("de_", id)
  if (!is.null(de) || file.exists(file.path(TAB_DIR, paste0("12_", suffix, ".csv")))) {
    write_result(if (is.null(de)) empty_de else de, suffix)
  }
}

for (id in fitted_ids) {
  run <- runs[[id]]
  spec <- specs[specs$analysis == id, ]
  effects <- paired_effects(run$fit$y, spec$positive, spec$negative)
  panel_genes <- intersect(unique(discovery_effects$gene), rownames(effects))
  panel_results[[id]] <- mutate(gene_effects_long(effects, panel_genes), analysis = id)
  program_results[[id]] <- mutate(
    program_differences(effects, gate_programs, cfg$minimum_program_coverage), analysis = id)
  program_fry[[id]] <- mutate(fry_programs(run$fit$y, gate_programs), analysis = id, .before = 1)

  discovery <- discovery_de[[spec$discovery_contrast]]
  joined <- inner_join(select(discovery, gene, log2_fc, fdr), select(run$fit$de, gene, log2_fc, fdr),
                       by = "gene", suffix = c("_discovery", "_external")) %>%
    mutate(annotation_gene = gene %in% annotation_genes,
           concordant = sign(log2_fc_discovery) == sign(log2_fc_external))
  compare <- subset(joined, !annotation_gene)
  hits <- subset(compare, fdr_discovery < cfg$fdr)
  replication[[id]] <- data.frame(
    analysis = id, n_shared = nrow(compare),
    rho = if (nrow(compare) > 1L) {
      cor(compare$log2_fc_discovery, compare$log2_fc_external, method = "spearman")
    } else NA_real_,
    concordant_shared = sum(compare$concordant),
    concordance_rate = if (nrow(compare)) mean(compare$concordant) else NA_real_,
    n_discovery_hits = nrow(hits), concordant_hits = sum(hits$concordant),
    concordant_hit_rate = if (nrow(hits)) mean(hits$concordant) else NA_real_,
    external_fdr_hits = sum(hits$concordant & hits$fdr_external < cfg$fdr)
  )

  cross_cohort[[id]] <- mutate(joined, analysis = id)
  if (id %in% combined_ids) {
    # The discovery annotation used CD274 and IDO1. Joint effects use shared held-out genes.
    shared <- intersect(rownames(discovery_programs), rownames(effects))
    program_shared[[id]] <- bind_rows(
      mutate(program_differences(discovery_programs[shared, , drop = FALSE],
                                 programs, cfg$minimum_program_coverage), cohort = "discovery"),
      mutate(program_differences(effects[shared, , drop = FALSE],
                                 programs, cfg$minimum_program_coverage), cohort = "GSE181919")
    ) %>% mutate(analysis = id)
  }
}

empty_replication <- data.frame(analysis = character(), n_shared = integer(), rho = numeric(),
                                concordant_shared = integer(), concordance_rate = numeric(),
                                n_discovery_hits = integer(), concordant_hits = integer(),
                                concordant_hit_rate = numeric(), external_fdr_hits = integer())
cross_summary <- sensitivity %>%
  select(analysis, pool, discovery_contrast, min_cells, n_donors, status, selected_for_validation) %>%
  left_join(bind_rows(c(list(empty_replication), replication)), by = "analysis")
empty_cross <- data.frame(gene = character(), log2_fc_discovery = numeric(), fdr_discovery = numeric(),
                          log2_fc_external = numeric(), fdr_external = numeric(),
                          annotation_gene = logical(), concordant = logical(), analysis = character())
empty_panel <- data.frame(gene = character(), patient = character(), log2_cpm_difference = numeric(),
                          analysis = character())
write_result(cross_summary, "cross_cohort_summary")
write_result(bind_rows(c(list(empty_cross), cross_cohort)), "cross_cohort_genes")
write_result(bind_rows(c(list(empty_panel), panel_results)), "panel_donor_effects")

empty_scores <- data.frame(program = character(), patient = character(), n_genes = integer(),
                           n_target_genes = integer(), genes = character(), effect = numeric(),
                           analysis = character())
scores <- bind_rows(c(list(empty_scores), program_results))
shared_scores <- bind_rows(c(list(mutate(empty_scores, cohort = character())), program_shared))
write_result(scores, "program_donor_effects")
write_result(shared_scores, "shared_program_donor_effects")

summarize_for_contrasts <- function(effects, analyses, panels) {
  unavailable <- data.frame(program = names(panels), n_genes = NA_integer_,
                            n_target_genes = lengths(panels), effect = NA_real_)
  bind_rows(lapply(analyses, function(id) {
    rows <- effects[effects$analysis == id, , drop = FALSE]
    cbind(analysis = id, summarize_programs(if (nrow(rows)) rows else unavailable))
  }))
}

empty_fry <- data.frame(analysis = character(), program = character(), n_genes_fry = integer(),
                        direction_fry = character(), p_fry = numeric(), fdr_fry = numeric(),
                        p_fry_mixed = numeric(), fdr_fry_mixed = numeric())
fry_results <- bind_rows(c(list(empty_fry), program_fry))
program_summary <- summarize_for_contrasts(scores, sensitivity$analysis, gate_programs) %>%
  left_join(select(sensitivity, analysis, n_eligible = n_donors, status), by = "analysis") %>%
  left_join(fry_results, by = c("analysis", "program")) %>%
  mutate(status = if_else(status == "fitted" & n_donors == 0L, "insufficient_program_coverage", status))
write_result(program_summary, "program_summary")
write_result(select(program_summary, analysis, program, n_eligible, status, n_genes_fry,
                    direction_fry, p_fry, fdr_fry, p_fry_mixed, fdr_fry_mixed), "program_fry")

validation_models <- sensitivity %>%
  filter(analysis %in% combined_ids) %>%
  select(analysis, min_cells, n_external = n_donors, status, selected_for_validation)
combined <- summarize_for_contrasts(shared_scores, combined_ids, programs) %>%
  left_join(validation_models, by = "analysis") %>%
  left_join(discovery_fry, by = "program") %>%
  left_join(select(fry_results, analysis, program, p_fry_external = p_fry,
                   direction_fry_external = direction_fry), by = c("analysis", "program")) %>%
  mutate(n_discovery = ncol(discovery_programs),
         status = if_else(status == "fitted" & n_donors == 0L, "insufficient_shared_program_coverage", status),
         p_fry_meta = if_else(status == "fitted",
                             pchisq(-2 * (log(p_fry_discovery) + log(p_fry_external)),
                                    df = 4L, lower.tail = FALSE), NA_real_)) %>%
  group_by(analysis) %>%
  mutate(fdr_fry_meta = p.adjust(p_fry_meta, "BH")) %>%
  ungroup() %>%
  relocate(min_cells, n_discovery, n_external, selected_for_validation, status, .after = program)
stopifnot(with(subset(combined, status == "fitted"), all(n_donors == n_discovery + n_external)))
write_result(combined, "combined_program_summary")

plot_ids <- intersect(c("primary", hybrid_id), fitted_ids)
if (!length(plot_ids)) plot_ids <- selected_id
if (length(plot_ids)) {
  plot_info <- cross_summary %>%
    filter(analysis %in% plot_ids) %>%
    mutate(heading = case_when(
      analysis == "primary" ~ "Primary myeloid APC pool",
      analysis == hybrid_id ~ "Author migDC and CD163-negative cDC2",
      TRUE ~ "Myeloid APC pool, five-cell minimum"
    ), panel = sprintf("%s (n = %d)\nSpearman rho = %.3f  |  %.1f%% of discovery hits concordant",
                       heading, n_donors, rho, 100 * concordant_hit_rate))
  panel_order <- plot_info$panel[match(plot_ids, plot_info$analysis)]
  plot_data <- bind_rows(cross_cohort[plot_ids]) %>%
    filter(!annotation_gene) %>%
    left_join(select(plot_info, analysis, panel), by = "analysis") %>%
    mutate(panel = factor(panel, levels = panel_order))
  core_genes <- c("CD83", "BIRC3", "PDCD1LG2", "CD200", "SOCS1", "SOCS2", "TNFAIP3", "NFKBIA")
  markers <- filter(plot_data, gene %in% core_genes)
  p <- ggplot(plot_data, aes(log2_fc_discovery, log2_fc_external)) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_point(aes(colour = fdr_discovery < cfg$fdr & fdr_external < cfg$fdr & concordant),
               alpha = 0.54, size = 1.15) +
    ggrepel::geom_text_repel(data = markers, aes(label = gene), size = 3,
                             colour = "#203541", seed = SEED, max.overlaps = Inf,
                             min.segment.length = 0, box.padding = 0.32) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = c("FALSE" = "#9CA9B2", "TRUE" = "#117B75"),
                        labels = c("Other shared genes", "Concordant, FDR < 0.05 in both cohorts")) +
    labs(x = "Discovery log2 fold change", y = "GSE181919 log2 fold change",
         colour = NULL, title = "Cross-cohort validation of migDC vs cDC2",
         subtitle = "Shared genes outside the assignment gene set") +
    theme_pub() + theme(legend.position = "bottom", strip.text = element_text(lineheight = 1.12))
} else {
  gate_counts <- runs$primary$eligibility %>%
    select(sample, migDC = n_positive, cDC2 = n_negative) %>%
    pivot_longer(-sample, names_to = "subtype", values_to = "n_cells") %>%
    mutate(subtype = factor(subtype, levels = c("migDC", "cDC2")))
  thresholds <- filter(validation_models, analysis %in% validation_ids)
  eligibility_labels <- paste0(thresholds$n_external, " pairs at ", thresholds$min_cells, " cells")
  p <- ggplot(gate_counts, aes(n_cells, sample)) +
    geom_vline(xintercept = thresholds$min_cells, colour = "grey55", linetype = "dashed") +
    geom_col(fill = "#117B75", width = 0.63) +
    geom_text(aes(label = n_cells), hjust = -0.35, size = 3) +
    facet_wrap(~ subtype) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = "Cells per tumour library", y = NULL,
         title = "GSE181919 paired validation eligibility",
         subtitle = paste0(paste(eligibility_labels, collapse = " | "), ". ",
                           cfg$min_donors, " matched donors required."),
         caption = "Dashed lines mark the cell thresholds for the two-marker gate.") +
    theme_pub() + theme(panel.grid.major.y = element_blank())
}
save_fig(p, "12_cross_cohort_effects", if (length(plot_ids) > 1L) 11 else 7, 5.5)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(TAB_DIR, "12_sessionInfo.txt"))
log_msg("STEP 12 done")

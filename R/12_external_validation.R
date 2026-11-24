# Stage 12: independent validation in primary oral-cavity tumours from GSE181919.
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

write_result <- function(x, suffix) {
  write.csv(x, file.path(TAB_DIR, paste0("12_", suffix, ".csv")), row.names = FALSE)
}
authors <- read.delim(gzfile(inputs[1]), row.names = NULL)
stopifnot(!anyDuplicated(authors$row.names),
          all(c("patient.id", "sample.id", "tissue.type", "subsite", "cell.type") %in% names(authors)))
oral <- subset(authors, tissue.type == "CA" & subsite == "OC")
author_dc <- subset(oral, cell.type == "Dendritic.cells")
sample_table <- unique(data.frame(sample = oral$sample.id, patient = paste0("Choi_", oral$patient.id),
                                  dataset = "GSE181919", group = "tumour"))
sample_table <- sample_table[order(sample_table$sample), ]
stopifnot(!anyDuplicated(sample_table$patient), !anyDuplicated(sample_table$sample))
inventory <- sample_table %>%
  left_join(count(oral, sample.id, name = "n_cells"), by = c("sample" = "sample.id")) %>%
  left_join(count(author_dc, sample.id, name = "n_author_dc"), by = c("sample" = "sample.id")) %>%
  mutate(n_author_dc = replace_na(n_author_dc, 0L))
write_result(inventory, "cohort_inventory")
write_result(data.frame(file = basename(inputs), bytes = file.size(inputs),
                        md5 = unname(tools::md5sum(inputs)),
                        available = file.exists(inputs)), "inputs")
log_msg("  ", sum(inventory$n_cells), " cells, ", sum(inventory$n_author_dc),
        " author-annotated DCs, ", nrow(inventory), " oral tumours")
if ("--inventory-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no")
if (!file.exists(inputs[2])) stop("Missing GEO UMI counts. Only the metadata inventory was generated.")
md <- sample_table[match(author_dc$sample.id, sample_table$sample), ]
rownames(md) <- author_dc$row.names

# Read only relevant columns in bounded chunks. The source is a wide TSV whose
# decompressed size makes a dense read of all 54,239 cells impractical.
read_dc_columns <- function(path, barcode_ids) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  first_gene <- readLines(con, n = 1)
  fields <- length(strsplit(first_gene, "\t", fixed = TRUE)[[1]])
  if (length(header) == fields) header <- header[-1L]
  header <- gsub('"', "", header, fixed = TRUE)
  at <- match(barcode_ids, header)
  if (length(header) != fields - 1L || anyNA(at) || anyDuplicated(header)) {
    stop("Cell barcodes do not match uniquely between the UMI matrix and author metadata.")
  }
  select_columns <- c(1L, at + 1L)
  parts <- list()
  lines <- c(first_gene, readLines(con, n = 95L))
  repeat {
    if (!length(lines)) break
    chunk <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                               sep = "\t", select = select_columns, showProgress = FALSE,
                               nThread = 1L)
    if (nrow(chunk) != length(lines) || ncol(chunk) != length(select_columns)) {
      stop("Unexpected number of rows or columns while reading GEO counts.")
    }
    parts[[length(parts) + 1L]] <- chunk
    lines <- readLines(con, n = 96L)
  }
  tab <- data.table::rbindlist(parts)
  genes <- tab[[1L]]
  storage <- as.matrix(tab[, -1L])
  dimnames(storage) <- list(genes, barcode_ids)
  if (!is.numeric(storage) || anyNA(genes) || anyDuplicated(genes) ||
      any(!is.finite(storage) | storage < 0 | storage != floor(storage))) {
    stop("Expected unique gene IDs and non-negative integer UMI counts.")
  }
  Matrix::Matrix(storage, sparse = TRUE)
}

log_msg("  Selecting ", nrow(author_dc), " author-annotated DCs in ", nrow(sample_table), " oral tumours")
counts <- read_dc_columns(inputs[2], rownames(md))

fit_validation <- function(labels, min_cells, positive, negative) {
  cells <- labels %in% unique(gates$subtype)
  annotated <- md[cells, , drop = FALSE]
  annotated$dc_subtype <- labels[cells]
  bulk <- if (any(cells)) aggregate_umi(counts[, cells, drop = FALSE], annotated) else {
    list(metadata = transform(annotated, n_cells = integer()))
  }
  eligibility <- paired_cell_counts(bulk, sample_table, positive, negative, min_cells)
  chosen <- subset(eligibility, eligible)
  fit <- if (nrow(chosen) >= cfg$min_donors) {
    fit_paired_bulk(paired_bulk(bulk, chosen, positive, negative))
  } else NULL
  list(eligibility = eligibility, chosen = chosen, fit = fit)
}

main_labels <- classify_dc(counts, gates, n_detected = 2L)
loose_labels <- classify_dc(counts, gates, n_detected = 1L)
specs <- data.frame(analysis = c("primary", "min_5", "min_15", "min_20", "loose_gate",
                                 "migDC_vs_cDC1", "cDC1_vs_cDC2"),
                    gate = c(2L, 2L, 2L, 2L, 1L, 2L, 2L),
                    min_cells = c(10L, 5L, 15L, 20L, 10L, 5L, 5L),
                    positive = c(rep("migDC", 6L), "cDC1"),
                    negative = c(rep("cDC2", 5L), "cDC1", "cDC2"))
runs <- setNames(lapply(seq_len(nrow(specs)), function(i) {
  spec <- specs[i, ]
  labels <- if (spec$gate == 2L) main_labels else loose_labels
  fit_validation(labels, spec$min_cells, spec$positive, spec$negative)
}), specs$analysis)

write_result(bind_rows(
  count(mutate(md, subtype = main_labels), sample, patient, subtype, name = "n") %>%
    mutate(gate_min_detected = 2L),
  count(mutate(md, subtype = loose_labels), sample, patient, subtype, name = "n") %>%
    mutate(gate_min_detected = 1L)
), "classification")
write_result(bind_rows(lapply(names(runs), function(id) mutate(runs[[id]]$eligibility, analysis = id))),
             "donor_eligibility")

summary <- bind_rows(lapply(names(runs), function(id) {
  run <- runs[[id]]
  de <- run$fit$de
  cbind(specs[specs$analysis == id, ], data.frame(
    n_donors = nrow(run$chosen), donors = paste(run$chosen$patient, collapse = " "),
    status = if (is.null(de)) "insufficient_donors" else "fitted",
    n_genes = if (is.null(de)) NA_integer_ else nrow(de),
    n_fdr = if (is.null(de)) NA_integer_ else sum(de$fdr < cfg$fdr)
  ))
}))
write_result(summary, "sensitivity")
print(summary, row.names = FALSE)
if (is.null(runs$primary$fit) && is.null(runs$min_5$fit)) {
  stop("The strict external gates yielded too few donors for the 10- or 5-cell paired comparisons.")
}

discovery <- read.csv(file.path(TAB_DIR, "11_de_migDC_vs_cDC2.csv"))
discovery_programs <- read.csv(file.path(TAB_DIR, "11_program_gene_donor_effects.csv"))
discovery_programs <- as.matrix(xtabs(log2_cpm_difference ~ gene + patient, data = discovery_programs))
discovery_effects <- read.csv(file.path(TAB_DIR, "11_panel_donor_effects.csv"))
annotation_genes <- union(gates$gene, read.csv(file.path("config", "dc_programs.csv")) %>%
                            filter(role == "annotation") %>% pull(gene))
cross_cohort <- list()
replication <- list()
program_results <- list()
program_shared <- list()
panel_results <- list()
for (id in names(runs)) {
  run <- runs[[id]]
  if (is.null(run$fit)) next
  de <- run$fit$de
  if (id %in% c("primary", "min_5", "migDC_vs_cDC1", "cDC1_vs_cDC2")) {
    write_result(de, paste0("de_", id))
  }
  spec <- specs[specs$analysis == id, ]
  effects <- paired_effects(run$fit$y, spec$positive, spec$negative)
  panel_genes <- intersect(unique(discovery_effects$gene), rownames(effects))
  panel_results[[id]] <- mutate(gene_effects_long(effects, panel_genes), analysis = id)
  if (spec$positive != "migDC" || spec$negative != "cDC2") next

  joined <- inner_join(select(discovery, gene, log2_fc, fdr), select(de, gene, log2_fc, fdr),
                       by = "gene", suffix = c("_discovery", "_external")) %>%
    mutate(annotation_gene = gene %in% annotation_genes,
           concordant = sign(log2_fc_discovery) == sign(log2_fc_external))
  compare <- subset(joined, !annotation_gene)
  hits <- subset(compare, fdr_discovery < cfg$fdr)
  replication[[id]] <- data.frame(analysis = id, n_shared = nrow(compare),
                                  rho = cor(compare$log2_fc_discovery, compare$log2_fc_external,
                                            method = "spearman"),
                                  n_discovery_hits = nrow(hits),
                                  concordant_hits = sum(hits$concordant),
                                  external_fdr_hits = sum(hits$concordant & hits$fdr_external < cfg$fdr))
  if (id %in% c("primary", "min_5")) cross_cohort[[id]] <- mutate(joined, analysis = id)

  scores <- program_differences(effects, programs, cfg$minimum_program_coverage)
  program_results[[id]] <- mutate(scores, analysis = id)
  if (id %in% c("primary", "min_5")) {
    shared <- intersect(rownames(discovery_programs), rownames(effects))
    program_shared[[id]] <- bind_rows(
      mutate(program_differences(discovery_programs[shared, , drop = FALSE],
                                 programs, cfg$minimum_program_coverage), cohort = "discovery"),
      mutate(program_differences(effects[shared, , drop = FALSE],
                                 programs, cfg$minimum_program_coverage), cohort = "GSE181919")
    ) %>% mutate(analysis = id)
  }
}
write_result(bind_rows(replication), "cross_cohort_summary")
write_result(bind_rows(cross_cohort), "cross_cohort_genes")
write_result(bind_rows(panel_results), "panel_donor_effects")
scores <- bind_rows(program_results)
write_result(scores, "program_donor_effects")
write_result(bind_rows(program_shared), "shared_program_donor_effects")
write_result(bind_rows(lapply(split(scores, scores$analysis), function(d) {
  cbind(analysis = d$analysis[1], summarize_programs(d))
})), "program_summary")

plot_id <- if (!is.null(runs$primary$fit)) "primary" else "min_5"
selected <- cross_cohort[[plot_id]]
plot_data <- filter(selected, !annotation_gene)
p <- ggplot(plot_data, aes(log2_fc_discovery, log2_fc_external)) +
  geom_hline(yintercept = 0, colour = "grey78", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey78", linetype = "dashed") +
  geom_point(aes(colour = fdr_discovery < cfg$fdr & fdr_external < cfg$fdr & concordant),
             alpha = 0.58, size = 1.7) +
  scale_colour_manual(values = c("FALSE" = "#9CA9B2", "TRUE" = "#117B75"),
                      labels = c("Other shared genes", "Concordant, FDR < 0.05 in both")) +
  labs(x = "Discovery log2 fold change", y = "GSE181919 log2 fold change",
       colour = NULL, title = "Independent oral-cavity tumour cohort",
       subtitle = paste0("migDC vs cDC2, excluding ", length(annotation_genes),
                         " assignment genes; ", runs[[plot_id]]$eligibility$min_cells[1],
                         " cells per donor and subtype minimum")) +
  theme_pub() + theme(legend.position = "bottom")
save_fig(p, "12_cross_cohort_effects", 7, 5.5)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(TAB_DIR, "12_sessionInfo.txt"))
log_msg("STEP 12 done")

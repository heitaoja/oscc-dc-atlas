# =============================================================================
# 01_load_qc.R  --  read per-sample count matrices, map IDs, QC, de-doublet
# -----------------------------------------------------------------------------
# Input : data/raw/<GSE>/<GSM>_seurat.rds      raw UMI counts (integer)
#         data/raw/<GSE>/<GSM>_metaData.rds    depositor QC (cross-check only)
# Output: results/objects/01_qc_filtered.rds   merged Seurat object (counts only)
#         results/tables/01_sample_qc.csv      per-sample cell numbers
#         results/tables/01_cell_numbers.csv   called / QC-pass / doublets
#         results/objects/tmp/<GSM>.rds        per-sample cache (skipped on re-run)
#
# Design notes
#  * the deposited `<GSM>_seurat.rds` objects carry integer UMI counts in the
#    `counts` layer; the separate `<GSM>_mat.rds` files are log1p(CPM) and must
#    NOT be used as counts (see docs/data_provenance.md).
#  * samples are processed in parallel (future multisession); each worker only
#    holds one sample, so peak memory is workers x ~3 GB.
#  * QC combines absolute bounds (PARAMS$qc) with a per-sample adaptive upper
#    bound (median + 3 MAD of log10 nFeature) to handle library-depth variation.
#  * doublets are called per sample with scDblFinder on the raw counts.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr)
  library(future); library(future.apply)
})

log_msg("STEP 01 -- load & QC")
samples <- PARAMS$samples
dir.create(file.path(OBJ_DIR, "tmp"), showWarnings = FALSE)
write.csv(samples, file.path(TAB_DIR, "00_samples.csv"), row.names = FALSE)
write.csv(data.frame(parameter = names(unlist(PARAMS$qc)),
                     value = as.character(unlist(PARAMS$qc))),
          file.path(TAB_DIR, "00_params_qc.csv"), row.names = FALSE)

# ---- 1. common gene space ----------------------------------------------
read_counts <- function(sample) {
  ds <- samples$dataset[samples$sample == sample]
  o  <- readRDS(file.path(RAW_DIR, ds, paste0(sample, "_seurat.rds")))
  as(LayerData(o, assay = "RNA", layer = "counts"), "CsparseMatrix")
}

gene_cache <- file.path(OBJ_DIR, "tmp", "genes.rds")
if (!file.exists(gene_cache)) {
  m1 <- read_counts("GSM5258385"); s1 <- ensembl_to_symbol(rownames(m1))
  m2 <- read_counts("GSM5820668"); s2 <- ensembl_to_symbol(rownames(m2))
  genes <- intersect(s1[!is.na(s1)], s2[!is.na(s2)])
  rm(m1, m2); gc()
  saveRDS(genes, gene_cache)
}
genes <- readRDS(gene_cache)
log_msg("common gene space: ", length(genes), " symbols")

# ---- 2. per-sample processing (parallel) -------------------------------
process_sample <- function(sample) {
  cache <- file.path(OBJ_DIR, "tmp", paste0(sample, ".rds"))
  if (file.exists(cache)) { return(readRDS(cache)@misc$qc_stats) }

  ds <- samples$dataset[samples$sample == sample]
  mat <- read_counts(sample)
  n_called <- ncol(mat)

  sym <- ensembl_to_symbol(rownames(mat))
  mat <- collapse_duplicate_symbols(mat, sym)
  mat <- mat[intersect(genes, rownames(mat)), , drop = FALSE]
  gc()

  md <- readRDS(file.path(RAW_DIR, ds, paste0(sample, "_metaData.rds")))
  md <- md[colnames(mat), , drop = FALSE]
  md$sample  <- sample
  md$dataset <- ds
  md$group   <- samples$group[samples$sample == sample]
  md$patient <- samples$patient[samples$sample == sample]

  obj <- CreateSeuratObject(counts = mat, meta.data = md, project = sample)
  rm(mat); gc()
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[AB]")

  # ---- QC ----
  qc <- PARAMS$qc
  n_genes <- obj$nFeature_RNA; n_counts <- obj$nCount_RNA
  lg <- log10(pmax(n_genes, 1))
  mad_lg <- mad(lg, constant = 1.4826)
  hi <- if (mad_lg > 0) 10^(median(lg) + qc$mad_n * mad_lg) else Inf
  keep <- n_genes >= qc$min_genes_per_cell & n_genes <= qc$max_genes_per_cell &
          n_genes <= hi &
          n_counts >= qc$min_counts_per_cell &
          obj$percent.mt <= qc$max_mito_pct &
          obj$percent.hb <= qc$max_hb_pct
  obj$qc_pass <- keep
  n_after_qc <- sum(keep)
  obj <- obj[, keep]

  # ---- doublets ----
  n_before_dbl <- ncol(obj)
  if (isTRUE(qc$remove_doublets) && requireNamespace("scDblFinder", quietly = TRUE)) {
    res <- scDblFinder::scDblFinder(obj[["RNA"]]$counts, verbose = FALSE)
    cls <- setNames(as.character(res$scDblFinder.class), colnames(res))
    scr <- setNames(as.numeric(res$scDblFinder.score), colnames(res))
    obj$scDblFinder_class <- unname(cls[colnames(obj)])
    obj$scDblFinder_score <- unname(scr[colnames(obj)])
    n_doublet <- sum(obj$scDblFinder_class == "doublet")
    obj <- obj[, obj$scDblFinder_class != "doublet"]
  } else {
    obj$scDblFinder_class <- NA_character_; obj$scDblFinder_score <- NA_real_
    n_doublet <- NA_integer_
  }

  obj@misc$qc_stats <- data.frame(
    sample = sample, dataset = ds, cells_called = n_called,
    pass_qc = n_after_qc,
    remove_pct_qc = round(100 * (1 - n_after_qc / n_called), 1),
    cells_final = ncol(obj),
    doublets_removed = n_doublet,
    remove_pct_doublet = round(100 * n_doublet / n_before_dbl, 1)
  )
  saveRDS(obj, cache)
  obj@misc$qc_stats
}

n_workers <- max(1, min(4, parallel::detectCores() - 2))
plan(multisession, workers = n_workers)
log_msg("  processing ", nrow(samples), " samples with ", n_workers, " workers")
stats_list <- future_lapply(samples$sample, process_sample, future.seed = TRUE)
plan(sequential)
qc_raw <- do.call(rbind, lapply(stats_list, function(d) {
  # tolerate caches written by an earlier version of this script
  if (!"doublets_removed" %in% names(d)) d$doublets_removed <- d$pass_qc - d$cells_final
  if (!"remove_pct_doublet" %in% names(d)) d$remove_pct_doublet <- NA
  d$remove_pct_doublet <- round(100 * d$doublets_removed / d$pass_qc, 1)
  d[, c("sample", "dataset", "cells_called", "pass_qc", "remove_pct_qc",
        "cells_final", "doublets_removed", "remove_pct_doublet")]
}))
write.csv(qc_raw, file.path(TAB_DIR, "01_cell_numbers.csv"), row.names = FALSE)
print(qc_raw)

# ---- 3. merge -----------------------------------------------------------
objs <- lapply(samples$sample, function(s) readRDS(file.path(OBJ_DIR, "tmp", paste0(s, ".rds"))))
names(objs) <- samples$sample
genes_all <- Reduce(intersect, lapply(objs, rownames))
log_msg("genes retained in all samples: ", length(genes_all))
objs <- lapply(objs, function(o) o[genes_all, ])
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "OSCC")
# Seurat v5 keeps one counts layer per sample after merge(); DE and marker
# detection require a single joined layer, so join here (not in every script).
merged <- JoinLayers(merged)
merged$sample  <- factor(merged$sample, levels = samples$sample)
merged$dataset <- factor(merged$dataset, levels = unique(samples$dataset))
merged$group   <- factor(merged$group, levels = c("tumour", "normal"))
rm(objs); gc()
saveRDS(merged, file.path(OBJ_DIR, "01_qc_filtered.rds"))

# ---- 4. summary table ---------------------------------------------------
qc_tab <- merged@meta.data %>%
  group_by(sample, dataset, group, patient) %>%
  summarise(n_cells = n(),
            median_counts = median(nCount_RNA),
            median_genes  = median(nFeature_RNA),
            median_mt     = round(median(percent.mt), 2), .groups = "drop") %>%
  as.data.frame()
write.csv(qc_tab, file.path(TAB_DIR, "01_sample_qc.csv"), row.names = FALSE)
print(qc_tab)
log_msg("STEP 01 done: ", ncol(merged), " cells x ", nrow(merged), " genes")

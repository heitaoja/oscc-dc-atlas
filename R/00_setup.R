# =============================================================================
# 00_setup.R  --  project configuration and shared helpers
# -----------------------------------------------------------------------------
# Sourced by every analysis script.  Nothing here should depend on results of
# earlier steps; it only defines paths, parameters and small utility functions.
# =============================================================================

# ---- paths --------------------------------------------------------------
PROJ      <- normalizePath(".", winslash = "/")
RAW_DIR   <- file.path(PROJ, "data", "raw")
OBJ_DIR   <- file.path(PROJ, "results", "objects")
TAB_DIR   <- file.path(PROJ, "results", "tables")
FIG_DIR   <- file.path(PROJ, "results", "figures")
LOG_DIR   <- file.path(PROJ, "results", "logs")
for (d in c(OBJ_DIR, TAB_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- global options -----------------------------------------------------
options(stringsAsFactors = FALSE, future.globals.maxSize = 8000 * 1024^2)
set.seed(20240601)                       # single seed, used everywhere
SEED <- 20240601

# ---- analysis parameters (single source of truth) -----------------------
# Every number quoted in README / report is defined HERE and written to
# results/tables/00_params.csv by 01_load_qc.R.  No magic numbers in scripts.
PARAMS <- list(
  # sample table: dataset -> samples.  "group" is tumour / normal.
  samples = data.frame(
    sample  = c("GSM5258385", "GSM5258386", "GSM5258387", "GSM5258388",
                "GSM5258389", "GSM5258390", "GSM5820668", "GSM5820669",
                "GSM5820670", "GSM5820671"),
    dataset = c(rep("GSE172577", 6), rep("GSE193766", 4)),
    group   = c(rep("tumour", 6), "tumour", "normal", "tumour", "normal"),
    patient = c(paste0("P", 1:6), "P7", "P7", "P8", "P8"),
    stringsAsFactors = FALSE
  ),

  # ---- QC ----
  qc = list(
    min_genes_per_cell = 200,    # absolute floor
    max_genes_per_cell = 6000,   # absolute ceiling (doublet / ambient guard)
    min_counts_per_cell = 500,
    max_mito_pct = 15,           # %
    max_hb_pct = 3,              # % haemoglobin genes -> erythrocytes
    min_cells_per_gene = 3,      # gene filter
    mad_n = 3,                   # adaptive outlier threshold (median + n*MAD)
    doublet_expected_rate = NULL,  # NULL -> scDblFinder default per sample
    remove_doublets = TRUE
  ),

  # ---- normalisation / feature selection ----
  hvg_n = 2000,
  hvg_method = "vst",
  n_pcs = 50,
  n_pcs_use = 30,

  # ---- integration ----
  integrate_method = "harmony",
  batch_var = "sample",
  harmony_theta = 2,

  # ---- clustering ----
  cluster_method = "leiden",
  resolution_main = 0.5,
  resolution_dc = 0.4,

  # ---- differential expression ----
  de = list(test = "wilcox", min_pct = 0.1, logfc = 0.25,
            padj_cut = 0.05, logfc_cut = 0.5, top_n_label = 15),

  # ---- enrichment ----
  enrich = list(padj_cut = 0.05, q_cut = 0.25, min_gs = 10, max_gs = 500),

  # ---- trajectory ----
  trajectory = list(use_clusters = c("cDC2", "migDC"), root_cluster = "cDC2",
                    n_top_trending = 40)
)

# ---- helper: logging ----------------------------------------------------
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n"); flush.console()
  invisible(msg)
}

# ---- helper: save a ggplot to png + pdf ---------------------------------
save_fig <- function(p, name, width = 6, height = 5, dpi = 300) {
  for (ext in c("png", "pdf")) {
    f <- file.path(FIG_DIR, paste0(name, ".", ext))
    if (ext == "png") {
      ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white")
    } else {
      ggsave(f, p, width = width, height = height, device = cairo_pdf, bg = "white")
    }
  }
  invisible(name)
}

# ---- helper: Ensembl -> SYMBOL -----------------------------------------
# Returns a character vector of symbols, NA where unmapped.
ensembl_to_symbol <- function(ids) {
  ids <- sub("\\..*$", "", ids)                     # drop version suffix
  sym <- rep(NA_character_, length(ids))
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    m <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
      column = "SYMBOL", multiVals = "first"))
    sym <- unname(m[ids])
  }
  sym
}

# ---- helper: collapse duplicated symbols by summing counts --------------
# Uses a sparse aggregation matrix: rownames of the result are unique symbols.
collapse_duplicate_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat  <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) {
    u <- unique(symbols)
    A <- Matrix::sparseMatrix(i = match(symbols, u), j = seq_along(symbols),
                              x = 1, dims = c(length(u), length(symbols)))
    mat <- as(A %*% mat, "CsparseMatrix")
    rownames(mat) <- u
  } else {
    rownames(mat) <- symbols
  }
  mat
}

# ---- helper: figure theme ----------------------------------------------
theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA))
}

# ---- helper: palette for cell types (stable across figures) -------------
celltype_pal <- function(types) {
  base <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
            "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
            "#8CD17D", "#F1CE63", "#A0CBE8", "#FFBE7D", "#D4A6C8", "#79706E")
  setNames(rep(base, length.out = length(types)), types)
}

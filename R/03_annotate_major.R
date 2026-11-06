# =============================================================================
# 03_annotate_major.R  --  reference-based + marker-based annotation
# -----------------------------------------------------------------------------
# Input : results/objects/02_clustered.rds
#         config/annotation_major.csv   (optional, cluster -> cell type)
# Output: results/objects/03_annotated.rds
#         results/tables/03_cluster_markers.csv
#         results/tables/03_singler_labels.csv
#         results/figures/03_*.png
#
# Strategy: three independent references (Blueprint/ENCODE, Monaco immune,
# Human Primary Cell Atlas) plus cluster markers.  Concordance between
# references and marker evidence is reported instead of trusting one label.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(SingleR); library(celldex); library(SummarizedExperiment)
})

log_msg("STEP 03 -- annotation")
obj <- readRDS(file.path(OBJ_DIR, "02_clustered.rds"))

# ---- 1. cluster markers (used as independent evidence) -------------------
mk_file <- file.path(TAB_DIR, "03_cluster_markers.csv")
recompute <- !file.exists(mk_file) || file.size(mk_file) < 1000
m <- NULL
if (!recompute) {
  m <- tryCatch(read.csv(mk_file), error = function(e) NULL)
  if (is.null(m) || !"cluster" %in% colnames(m) || nrow(m) == 0) recompute <- TRUE
}
if (recompute) {
  Idents(obj) <- "seurat_clusters"
  m <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1,
                      logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE)
  if (!is.data.frame(m) || !"cluster" %in% colnames(m) || nrow(m) == 0)
    stop("FindAllMarkers returned nothing -- are RNA layers joined? (JoinLayers)")
  write.csv(m, mk_file, row.names = FALSE)
  log_msg("  markers: ", nrow(m), " rows")
}
top <- m %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()

# ---- 2. SingleR with three references -----------------------------------
# SingleR is run on the *cluster-average* expression profile rather than on
# 85k single cells: it is memory-safe (20k x 23 instead of 20k x 85k), faster,
# and answers the question we actually ask ("which cell type is this cluster?").
singler_file <- file.path(TAB_DIR, "03_singler_labels.csv")
if (!file.exists(singler_file)) {
  avg <- AverageExpression(obj, assays = "RNA", layer = "data",
                           group.by = "seurat_clusters")$RNA
  # AverageExpression returns the mean on the *linear* scale (expm1 of the
  # log-normalised data); SingleR expects log-scale input, so re-log here.
  avg <- log1p(as.matrix(avg))
  colnames(avg) <- sub("^g", "", colnames(avg))
  refs <- list(
    blueprint = celldex::BlueprintEncodeData(),
    monaco    = celldex::MonacoImmuneData(),
    hpca      = celldex::HumanPrimaryCellAtlasData()
  )
  out <- lapply(names(refs), function(rn) {
    ref <- refs[[rn]]
    common <- intersect(rownames(avg), rownames(ref))
    log_msg("  ", rn, ": ", length(common), " common genes")
    lab <- SingleR(test = avg[common, , drop = FALSE],
                   ref = ref[common, , drop = FALSE],
                   labels = ref$label.main)
    data.frame(cluster = colnames(avg), reference = rn,
               label = as.character(lab$labels),
               score = apply(lab$scores, 1, max))
  })
  sl <- bind_rows(out) %>% arrange(cluster, reference)
  write.csv(sl, singler_file, row.names = FALSE)
} else {
  sl <- read.csv(singler_file)
}
print(sl, row.names = FALSE)

# ---- 3. apply curated annotation if available ---------------------------
ann_file <- file.path("config", "annotation_major.csv")
if (file.exists(ann_file)) {
  ann <- read.csv(ann_file, stringsAsFactors = FALSE)
  map <- setNames(ann$cell_type, ann$cluster)
  obj$cell_type <- unname(map[as.character(obj$seurat_clusters)])
  obj$cell_type[is.na(obj$cell_type)] <- "Unassigned"
  obj$cell_type <- factor(obj$cell_type)
  Idents(obj) <- "cell_type"
  pal <- celltype_pal(levels(obj$cell_type))
  p <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE,
               raster = TRUE, cols = pal) + theme_pub() + ggtitle("Major cell types")
  save_fig(p, "03_umap_celltype", 8, 6)
  p2 <- DotPlot(obj, features = unique(unlist(list(
    T = c("CD3D","CD8A","CD4"), NK = c("NKG7","GNLY"), B = c("MS4A1","CD79A"),
    Plasma = c("MZB1","JCHAIN"), Myeloid = c("LYZ","CD14","CD68"),
    DC = c("CD1C","CLEC9A","LILRA4"), Mast = c("TPSAB1","CPA3"),
    Epithelial = c("EPCAM","KRT5"), Fibroblast = c("COL1A1","DCN"),
    Endothelial = c("PECAM1","VWF"), Muscle = c("MYL1","DES"),
    Proliferating = c("MKI67","TOP2A")))),
    group.by = "cell_type") + RotatedAxis() + theme_pub(9)
  save_fig(p2, "03_marker_dotplot_celltype", 14, 6)
  write.csv(table(obj$cell_type, obj$sample) %>% as.data.frame(),
            file.path(TAB_DIR, "03_celltype_by_sample.csv"), row.names = FALSE)
  log_msg("  annotation applied: ", length(unique(obj$cell_type)), " cell types")
} else {
  log_msg("  config/annotation_major.csv not found -- run without final labels")
}

saveRDS(obj, file.path(OBJ_DIR, "03_annotated.rds"))
log_msg("STEP 03 done")

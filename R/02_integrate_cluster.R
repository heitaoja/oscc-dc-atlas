# =============================================================================
# 02_integrate_cluster.R  --  normalise, integrate (harmony), cluster, UMAP
# -----------------------------------------------------------------------------
# Input : results/objects/01_qc_filtered.rds
# Output: results/objects/02_clustered.rds
#         results/figures/02_umap_*.png, results/tables/02_cluster_counts.csv
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(patchwork)
})

log_msg("STEP 02 -- normalise / integrate / cluster")
obj <- readRDS(file.path(OBJ_DIR, "01_qc_filtered.rds"))
# make the script robust to objects saved by an older version of 01 (split layers)
if (length(Layers(obj[["RNA"]])) > 1) obj <- JoinLayers(obj)
log_msg("input: ", ncol(obj), " cells x ", nrow(obj), " genes")

obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = PARAMS$hvg_method,
                            nfeatures = PARAMS$hvg_n, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = PARAMS$n_pcs, verbose = FALSE)

# ---- integration ---------------------------------------------------------
# harmony on the PCA embedding; batch = sample (10 libraries).
obj <- RunHarmony(obj, group.by.vars = PARAMS$batch_var,
                  reduction.use = "pca", dims.use = 1:PARAMS$n_pcs_use,
                  theta = PARAMS$harmony_theta, plot_convergence = FALSE,
                  verbose = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:PARAMS$n_pcs_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = PARAMS$resolution_main,
                    algorithm = 4, random.seed = SEED, verbose = FALSE)   # Leiden
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:PARAMS$n_pcs_use,
               seed.use = SEED, verbose = FALSE)

saveRDS(obj, file.path(OBJ_DIR, "02_clustered.rds"))

# ---- figures -------------------------------------------------------------
pal <- celltype_pal(levels(obj$seurat_clusters))
p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
              raster = TRUE, cols = pal) + ggtitle("Leiden clusters") + theme_pub()
p2 <- DimPlot(obj, group.by = "sample", raster = TRUE) + ggtitle("Sample") + theme_pub()
p3 <- DimPlot(obj, group.by = "dataset", raster = TRUE) + ggtitle("Dataset") + theme_pub()
p4 <- DimPlot(obj, group.by = "group", raster = TRUE) + ggtitle("Tumour vs normal") + theme_pub()
save_fig(p1 + p2 + p3 + p4 + plot_layout(ncol = 2), "02_umap_overview", 12, 10)

p5 <- DimPlot(obj, group.by = "sample", split.by = "dataset", ncol = 5,
              raster = TRUE) + theme_pub(9)
save_fig(p5, "02_umap_by_sample", 15, 6)

p6 <- ElbowPlot(obj, ndims = PARAMS$n_pcs) + theme_pub()
save_fig(p6, "02_elbow", 5, 4)

# ---- cluster composition -------------------------------------------------
tab <- obj@meta.data %>%
  count(seurat_clusters, sample, group) %>%
  group_by(sample) %>% mutate(pct = round(100 * n / sum(n), 2)) %>% ungroup()
write.csv(tab, file.path(TAB_DIR, "02_cluster_composition.csv"), row.names = FALSE)

# ---- quick marker check for major lineages ------------------------------
canonical <- list(
  T        = c("CD3D", "CD3E", "TRAC"),
  NK       = c("NKG7", "GNLY", "KLRD1"),
  B        = c("MS4A1", "CD79A"),
  Plasma   = c("MZB1", "JCHAIN", "DERL3"),
  Myeloid  = c("LYZ", "CD68", "CD14"),
  DC       = c("CD1C", "CLEC9A", "LILRA4", "CD207"),
  Mast     = c("TPSAB1", "CPA3"),
  Epithelial = c("EPCAM", "KRT5", "KRT13"),
  Fibroblast = c("COL1A1", "DCN", "LUM"),
  Endothelial = c("PECAM1", "VWF", "ACKR1"),
  Muscle = c("MYL1", "DES"),
  Proliferating = c("MKI67", "TOP2A")
)
mrk <- unlist(canonical, use.names = FALSE)
mrk <- mrk[mrk %in% rownames(obj)]
p7 <- DotPlot(obj, features = mrk, group.by = "seurat_clusters") +
  RotatedAxis() + theme_pub(9) + ggtitle("Canonical lineage markers")
save_fig(p7, "02_marker_dotplot", 16, 6)

log_msg("STEP 02 done: ", nlevels(obj$seurat_clusters), " clusters")

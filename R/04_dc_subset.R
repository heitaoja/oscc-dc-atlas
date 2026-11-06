# =============================================================================
# 04_dc_subset.R  --  isolate dendritic cells and re-embed them
# -----------------------------------------------------------------------------
# Input : results/objects/03_annotated.rds   (needs obj$cell_type)
# Output: results/objects/04_dc_clustered.rds
#         results/tables/04_dc_cluster_markers.csv
#         results/figures/04_dc_*.png
#
# Rationale: DCs are ~2% of the atlas, so the global embedding is dominated by
# other lineages.  Sub-clustering on DC only, with a lower resolution and the
# sample as batch variable, is required before any subtype claim is made.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(patchwork)
})

log_msg("STEP 04 -- DC sub-clustering")
obj <- readRDS(file.path(OBJ_DIR, "03_annotated.rds"))
stopifnot("cell_type" %in% colnames(obj@meta.data))

dc_labels <- c("DC", "Dendritic cell", "Dendritic cells", "cDC", "pDC", "migDC")
sel <- obj$cell_type %in% dc_labels
log_msg("  DC cells: ", sum(sel), " / ", ncol(obj))
dc <- subset(obj, cells = colnames(obj)[sel])
dc <- DietSeurat(dc, assays = "RNA", dimreducs = NULL, graphs = NULL)
dc <- NormalizeData(dc, verbose = FALSE)
dc <- FindVariableFeatures(dc, nfeatures = 1500, verbose = FALSE)
dc <- ScaleData(dc, features = VariableFeatures(dc), verbose = FALSE)
dc <- RunPCA(dc, npcs = 30, verbose = FALSE)

# harmony only makes sense when several libraries contribute
if (nlevels(droplevels(dc$sample)) > 1) {
  dc <- RunHarmony(dc, group.by.vars = "sample", dims.use = 1:20,
                   theta = 2, verbose = FALSE)
  red <- "harmony"
} else red <- "pca"

dc <- FindNeighbors(dc, reduction = red, dims = 1:20, verbose = FALSE)
dc <- FindClusters(dc, resolution = PARAMS$resolution_dc, algorithm = 4,
                   random.seed = SEED, verbose = FALSE)
dc <- RunUMAP(dc, reduction = red, dims = 1:20, seed.use = SEED, verbose = FALSE)

# ---- markers -------------------------------------------------------------
Idents(dc) <- "seurat_clusters"
mk <- FindAllMarkers(dc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                     test.use = "wilcox")
write.csv(mk, file.path(TAB_DIR, "04_dc_cluster_markers.csv"), row.names = FALSE)

# ---- figures -------------------------------------------------------------
dc_markers <- c("CLEC9A", "XCR1", "BATF3", "IRF8",       # cDC1
                "CD1C", "CLEC10A", "SIRPA", "FCER1A",     # cDC2
                "CCR7", "LAMP3", "CCL19", "CCL17", "CCL22", # migratory / mregDC
                "LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7", # pDC
                "CD207", "CD1A",                          # LC-like
                "CD14", "CD163", "LYZ",                   # monocyte-like
                "MKI67", "TOP2A")                         # proliferating
dc_markers <- dc_markers[dc_markers %in% rownames(dc)]

pal <- celltype_pal(levels(dc$seurat_clusters))
p1 <- DimPlot(dc, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
              cols = pal) + theme_pub() + ggtitle("DC sub-clusters")
p2 <- DimPlot(dc, group.by = "sample") + theme_pub() + ggtitle("Sample")
p3 <- DimPlot(dc, group.by = "group") + theme_pub() + ggtitle("Tumour / normal")
save_fig(p1 + p2 + p3 + plot_layout(ncol = 3), "04_dc_umap", 15, 5)

p4 <- DotPlot(dc, features = dc_markers, group.by = "seurat_clusters") +
  RotatedAxis() + theme_pub(9) + ggtitle("Canonical DC markers by sub-cluster")
save_fig(p4, "04_dc_marker_dotplot", 16, 5.5)

p5 <- FeaturePlot(dc, features = intersect(c("CLEC9A", "CD1C", "CCR7", "LILRA4", "CD207"),
                                           rownames(dc)),
                  ncol = 3, order = TRUE, raster = TRUE) & theme_pub(9)
save_fig(p5, "04_dc_marker_features", 12, 7)

saveRDS(dc, file.path(OBJ_DIR, "04_dc_clustered.rds"))
tab <- dc@meta.data %>% count(seurat_clusters, sample, group) %>% as.data.frame()
write.csv(tab, file.path(TAB_DIR, "04_dc_cluster_by_sample.csv"), row.names = FALSE)
log_msg("STEP 04 done: ", ncol(dc), " DCs in ", nlevels(dc$seurat_clusters), " sub-clusters")

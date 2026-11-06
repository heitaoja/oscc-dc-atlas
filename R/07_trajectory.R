# =============================================================================
# 07_trajectory.R  --  pseudotime along the cDC2 -> migratory axis
# -----------------------------------------------------------------------------
# Input : results/objects/05_dc_annotated.rds
# Output: results/tables/07_*.csv, results/objects/07_trajectory.rds,
#         results/figures/07_*.png
#
# The original analysis used monocle3 with a root chosen by clicking on a
# UMAP plot (not recorded, not reproducible).  Here the root is fixed in code
# (cDC2) and the trajectory is fitted with slingshot on the DC-only embedding,
# which is deterministic given the seed and the cluster labels.
# Pseudotime is a *transcriptional ordering*, not evidence of differentiation;
# see docs/limitations.md.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(slingshot); library(dplyr)
  library(ggplot2); library(patchwork)
})

log_msg("STEP 07 -- trajectory")
dc <- readRDS(file.path(OBJ_DIR, "05_dc_annotated.rds"))
keep <- PARAMS$trajectory$use_clusters
dc <- subset(dc, subset = dc_subtype %in% keep)
dc$dc_subtype <- droplevels(factor(dc$dc_subtype))
log_msg("  cells in trajectory: ", ncol(dc), " (", paste(keep, collapse = ", "), ")")

sce <- as.SingleCellExperiment(dc)
reducedDims(sce)$UMAP <- Embeddings(dc, "umap")[colnames(sce), 1:2]
colData(sce)$cluster <- dc$dc_subtype[colnames(sce)]

set.seed(SEED)
sds <- slingshot(sce, clusterLabels = "cluster", reducedDim = "UMAP",
                 start.clus = PARAMS$trajectory$root_cluster, approx_points = 150)

pt <- slingPseudotime(sds)
pt_mean <- rowMeans(pt, na.rm = TRUE)
pt_mean[!is.finite(pt_mean)] <- NA
pseudo <- data.frame(cell = colnames(sce), pseudotime = pt_mean,
                     cluster = as.character(colData(sce)$cluster))
write.csv(pseudo, file.path(TAB_DIR, "07_pseudotime_by_cell.csv"), row.names = FALSE)
log_msg("  pseudotime range: ", paste(round(range(pt_mean, na.rm = TRUE), 2), collapse = " - "))

# ---- plot trajectory on UMAP -------------------------------------------
emb <- as.data.frame(reducedDims(sce)$UMAP)
colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$pseudotime <- pt_mean
emb$cluster <- pseudo$cluster
p1 <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = pseudotime)) +
  geom_point(size = 0.4) + scale_colour_viridis_c(na.value = "grey85") +
  labs(title = "Pseudotime (slingshot)") + theme_pub()
p2 <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = 0.4) + labs(title = "Subtype") + theme_pub()
save_fig(p1 + p2, "07_pseudotime_umap", 11, 5)

# ---- genes that change along pseudotime --------------------------------
expr <- LayerData(dc, assay = "RNA", layer = "data")[, pseudo$cell]
hv <- VariableFeatures(dc)
hv <- hv[hv %in% rownames(expr)]
ok <- !is.na(pt_mean)
cors <- apply(expr[hv, ok], 1, function(v) {
  if (sd(v) == 0) return(NA_real_)
  suppressWarnings(cor(v, pt_mean[ok], method = "spearman"))
})
cor_df <- data.frame(gene = names(cors), spearman_rho = as.numeric(cors)) %>%
  filter(!is.na(spearman_rho)) %>% arrange(desc(spearman_rho))
write.csv(cor_df, file.path(TAB_DIR, "07_gene_pseudotime_correlation.csv"), row.names = FALSE)

# ---- marker trajectories ------------------------------------------------
mk <- intersect(c("CLEC10A", "CD1C", "CCR7", "LAMP3", "CCL19", "CD274", "IDO1",
                  "CLEC9A", "XCR1", "CD207", "LILRA4"), rownames(expr))
if (length(mk)) {
  long <- lapply(mk, function(g) data.frame(cell = colnames(expr), gene = g,
                                            expr = as.numeric(expr[g, ]),
                                            pseudotime = pt_mean)) %>% bind_rows()
  long <- long %>% filter(!is.na(pseudotime))
  long$bin <- cut(long$pseudotime, breaks = 20, labels = FALSE)
  sm <- long %>% group_by(gene, bin) %>%
    summarise(mean_expr = mean(expr), .groups = "drop") %>%
    group_by(gene) %>% mutate(bin_pos = (bin - 0.5) / max(bin)) %>% ungroup()
  pg <- ggplot(sm, aes(bin_pos, mean_expr, colour = gene)) +
    geom_line(linewidth = 1) + labs(x = "Pseudotime (binned)", y = "Mean expression (log-normalised)",
                                    title = "Marker dynamics along pseudotime") +
    theme_pub()
  save_fig(pg, "07_marker_pseudotime", 8, 5)
}

saveRDS(list(sds = sds, pseudotime = pseudo), file.path(OBJ_DIR, "07_trajectory.rds"))
log_msg("STEP 07 done")

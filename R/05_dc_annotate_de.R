# =============================================================================
# 05_dc_annotate_de.R  --  annotate DC subtypes, then differential expression
# -----------------------------------------------------------------------------
# Input : results/objects/04_dc_clustered.rds
#         config/annotation_dc.csv  (cluster, cell_type, evidence)
# Output: results/objects/05_dc_annotated.rds
#         results/tables/05_dc_markers_<A>_vs_<B>.csv, 05_dc_marker_top.csv
#         results/figures/05_*.png
#
# Every DE test is run on the RNA assay (log-normalised counts), never on the
# batch-corrected embedding, and up- and down-regulated genes are kept apart so
# that enrichment can be direction-aware.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(ggrepel); library(patchwork)
})

log_msg("STEP 05 -- DC annotation & DE")
dc <- readRDS(file.path(OBJ_DIR, "04_dc_clustered.rds"))

# ---- 1. apply annotation table -----------------------------------------
ann_file <- file.path("config", "annotation_dc.csv")
if (!file.exists(ann_file)) {
  stop("config/annotation_dc.csv missing -- inspect results/figures/04_dc_marker_dotplot.png ",
       "and results/tables/04_dc_cluster_markers.csv first, then create the file.")
}
ann <- read.csv(ann_file, stringsAsFactors = FALSE)
dc$dc_subtype <- unname(setNames(ann$cell_type, ann$cluster)[as.character(dc$seurat_clusters)])
dc$dc_subtype[is.na(dc$dc_subtype)] <- "Other"
dc$dc_subtype <- factor(dc$dc_subtype)
# LowQuality clusters (ambient / low-complexity cells) are annotated but not
# carried into differential expression, trajectory or composition analyses.
n_lowq <- sum(dc$dc_subtype == "LowQuality")
dc <- subset(dc, subset = dc_subtype != "LowQuality")
dc$dc_subtype <- droplevels(dc$dc_subtype)
log_msg("  removed ", n_lowq, " low-quality DC-like cells; ", ncol(dc), " retained")
# Recompute the embedding after removing the low-quality cells: their positions
# dominate the UMAP axes and make the remaining subtypes unreadable.
dc <- RunUMAP(dc, reduction = "harmony", dims = 1:20, seed.use = SEED, verbose = FALSE)
Idents(dc) <- "dc_subtype"
log_msg("  subtypes: ", paste(levels(dc$dc_subtype), collapse = ", "))

pal <- celltype_pal(levels(dc$dc_subtype))
p <- DimPlot(dc, group.by = "dc_subtype", label = TRUE, repel = TRUE, cols = pal) +
  theme_pub() + ggtitle("DC subtypes")
save_fig(p, "05_dc_subtype_umap", 7, 5.5)

# ---- 2. subtype markers (all clusters) ----------------------------------
mk <- FindAllMarkers(dc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                     test.use = "wilcox", verbose = FALSE)
mk <- mk %>% filter(p_val_adj < PARAMS$de$padj_cut)
write.csv(mk, file.path(TAB_DIR, "05_dc_marker_all.csv"), row.names = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
write.csv(top, file.path(TAB_DIR, "05_dc_marker_top.csv"), row.names = FALSE)

# ---- 3. focal comparison: migratory vs cDC2 -----------------------------
# This is the axis the original thesis claimed; tested here as a two-group
# contrast on the RNA assay, reported with effect size and detection rate.
contrast <- c(ident.1 = "migDC", ident.2 = "cDC2")
if (all(contrast %in% levels(dc$dc_subtype))) {
  de <- FindMarkers(dc, ident.1 = contrast[1], ident.2 = contrast[2],
                    assay = "RNA", test.use = PARAMS$de$test,
                    min.pct = PARAMS$de$min_pct, logfc.threshold = PARAMS$de$logfc,
                    only.pos = FALSE, verbose = FALSE)
  de$gene <- rownames(de)
  de$direction <- ifelse(de$avg_log2FC > 0, "up_in_migDC", "up_in_cDC2")
  write.csv(de, file.path(TAB_DIR, "05_de_migDC_vs_cDC2.csv"), row.names = FALSE)

  sig <- de %>% filter(p_val_adj < PARAMS$de$padj_cut,
                       abs(avg_log2FC) >= PARAMS$de$logfc_cut)
  sig$score <- ifelse(sig$direction == "up_in_migDC", sig$avg_log2FC, -sig$avg_log2FC)
  lab <- sig %>% group_by(direction) %>%
    slice_max(score, n = PARAMS$de$top_n_label) %>% ungroup()
  de$label <- ifelse(de$gene %in% lab$gene, de$gene, NA)
  pv <- ggplot(de, aes(avg_log2FC, -log10(p_val_adj))) +
    geom_point(aes(colour = direction, size = pct.1), alpha = 0.6) +
    geom_vline(xintercept = c(-PARAMS$de$logfc_cut, PARAMS$de$logfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(PARAMS$de$padj_cut), linetype = "dashed") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    scale_colour_manual(values = c(up_in_migDC = "#D62728", up_in_cDC2 = "#1F77B4")) +
    labs(x = "log2 fold change (migDC / cDC2)", y = "-log10 adjusted P",
         title = "migratory DC vs cDC2", size = "fraction in migDC") +
    theme_pub()
  save_fig(pv, "05_volcano_migDC_vs_cDC2", 7, 6)
  log_msg("  DE genes (padj<0.05): ", sum(de$p_val_adj < PARAMS$de$padj_cut),
          " | |log2FC|>", PARAMS$de$logfc_cut, ": ", nrow(sig))
}

# ---- 4. marker dotplot across subtypes ---------------------------------
dc_markers <- c("CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C", "CLEC10A", "SIRPA",
                "FCER1A", "CCR7", "LAMP3", "CCL19", "CCL17", "CCL22", "CD274",
                "IDO1", "LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7", "CD207",
                "CD1A", "CD14", "CD163", "MKI67")
dc_markers <- dc_markers[dc_markers %in% rownames(dc)]
pd <- DotPlot(dc, features = dc_markers, group.by = "dc_subtype") +
  RotatedAxis() + theme_pub(9) + ggtitle("DC subtype markers")
save_fig(pd, "05_dc_subtype_dotplot", 15, 5.5)

# ---- 5. marker-level sanity: are the defining genes actually different? --
chk <- intersect(c("CLEC9A", "CD1C", "CLEC10A", "CCR7", "LAMP3", "LILRA4", "CD207"),
                 rownames(dc))
pv2 <- VlnPlot(dc, features = chk, group.by = "dc_subtype", pt.size = 0, ncol = 4,
               cols = pal) & theme_pub(9)
save_fig(pv2, "05_dc_subtype_marker_violin", 14, 7)

saveRDS(dc, file.path(OBJ_DIR, "05_dc_annotated.rds"))
tab <- dc@meta.data %>% count(sample, group, dc_subtype) %>% as.data.frame()
write.csv(tab, file.path(TAB_DIR, "05_dc_subtype_by_sample.csv"), row.names = FALSE)
log_msg("STEP 05 done")

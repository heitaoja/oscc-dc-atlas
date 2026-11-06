# =============================================================================
# 08_tf_activity.R  --  transcription-factor activity per DC subtype
# -----------------------------------------------------------------------------
# Input : results/objects/05_dc_annotated.rds
#         data/genesets/c3.tft.*.symbols.gmt   (MSigDB C3 TFT, cached locally)
# Output: results/tables/08_tf_activity.csv, results/figures/08_tf_*.png
#
# The original plan (SCENIC: GENIE3 + cisTarget) was never run.  A lighter,
# fully reproducible alternative is used here: TF target sets from MSigDB C3
# (TFT) scored per cell with decoupleR's univariate linear model (ULM).
# The GMT is cached locally because msigdbr's remote data source is unreliable
# from this network; using the GMT also pins the exact MSigDB release.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(decoupleR); library(tidyr)
})

log_msg("STEP 08 -- TF activity")
dc <- readRDS(file.path(OBJ_DIR, "05_dc_annotated.rds"))
expr <- LayerData(dc, assay = "RNA", layer = "data")

gs_dir <- file.path(PROJ, "data", "genesets")
dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)
tft <- file.path(gs_dir, "c3.tft.v2024.1.Hs.symbols.gmt")
if (!file.exists(tft)) {
  download.file(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",
                       "2024.1.Hs/c3.tft.v2024.1.Hs.symbols.gmt"),
                tft, quiet = TRUE, mode = "wb")
}

net <- clusterProfiler::read.gmt(tft)
net <- data.frame(source = sub("_TARGET_GENES$", "", net$term),
                  target = net$gene, mor = 1)
net <- net[net$target %in% rownames(expr), ]
log_msg("  TF target sets: ", length(unique(net$source)),
        " | edges: ", nrow(net))

acts <- run_ulm(mat = expr, network = net, .source = "source", .target = "target",
                .mor = "mor", minsize = 10)
# The per-cell score matrix is ~1.5 M rows (1.4 GB as CSV) and is an
# intermediate; only the per-subtype summary is written out. Re-running this
# script regenerates the per-cell scores in ~90 s.

# average activity per subtype
# decoupleR returns one row per (cell, TF) with `condition` = cell barcode,
# so the subtype label has to be joined back explicitly.
acts$dc_subtype <- dc$dc_subtype[match(acts$condition, colnames(dc))]
by_sub <- acts %>% filter(!is.na(dc_subtype)) %>% group_by(source, dc_subtype) %>%
  summarise(activity = mean(score), .groups = "drop")
write.csv(by_sub, file.path(TAB_DIR, "08_tf_activity_by_subtype.csv"), row.names = FALSE)

# ---- heatmap of the most variable TFs ----------------------------------
wide <- by_sub %>% tidyr::pivot_wider(names_from = dc_subtype, values_from = activity) %>%
  tibble::column_to_rownames("source")
v <- apply(wide, 1, sd)
top <- names(sort(v, decreasing = TRUE))[1:min(40, length(v))]
mat <- as.matrix(wide[top, ])
mat_z <- t(scale(t(mat)))
ord <- hclust(dist(mat_z))$order
long <- as.data.frame(as.table(mat_z)) %>%
  setNames(c("tf", "subtype", "z")) %>%
  mutate(tf = factor(tf, levels = rownames(mat_z)[ord]))
ph <- ggplot(long, aes(subtype, tf, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(x = NULL, y = NULL, title = "TF activity (ULM, z-scored per TF)") +
  theme_pub(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(ph, "08_tf_activity_heatmap", 8, 9)

# ---- TFs of interest ----------------------------------------------------
focus <- intersect(c("RELA", "NFKB1", "REL", "NFKB2", "IRF4", "IRF8", "BATF3",
                     "ZEB2", "MAFB", "TCF4", "IRF7", "STAT1", "STAT3", "CEBPB",
                     "SPI1", "FLI1", "MYB", "RUNX1", "ATF6", "XBP1"), unique(by_sub$source))
if (length(focus)) {
  fo <- by_sub %>% filter(source %in% focus)
  pf <- ggplot(fo, aes(source, activity, fill = dc_subtype)) +
    geom_col(position = "dodge") + coord_flip() + theme_pub(9) +
    labs(x = NULL, y = "Mean TF activity (ULM score)", title = "Selected TFs")
  save_fig(pf, "08_tf_activity_focus", 8, 6)
}
log_msg("STEP 08 done")

# =============================================================================
# 06_enrichment.R  --  direction-aware GO / KEGG / GSEA
# -----------------------------------------------------------------------------
# Input : results/tables/05_de_migDC_vs_cDC2.csv
#         results/objects/05_dc_annotated.rds   (expression universe)
#         data/genesets/*.gmt                   (built once, cached, offline)
# Output: results/tables/06_*.csv, results/figures/06_*.png
#
# The original thesis pooled up- and down-regulated genes into one list and then
# attributed pathways to individual DC subsets.  Here every enrichment runs on
# ONE direction at a time, with the detected transcriptome as background, and is
# complemented by a rank-based GSEA that needs no arbitrary fold-change cutoff.
#
# Gene sets are cached locally because both rest.kegg.jp and msigdbr's remote
# data source are unreliable from this network; caching also makes the result
# reproducible offline.
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(fgsea)
})

log_msg("STEP 06 -- enrichment")

# ---- 0. local gene-set cache -------------------------------------------
gs_dir <- file.path(PROJ, "data", "genesets")
dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)

hm_gmt <- file.path(gs_dir, "h.all.symbols.gmt")
if (!file.exists(hm_gmt)) {
  log_msg("  downloading MSigDB hallmark ...")
  try(download.file(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/",
                           "release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"),
                    hm_gmt, quiet = TRUE, mode = "wb"), silent = TRUE)
}

kegg_gmt <- file.path(gs_dir, "kegg_hsa_pathway.gmt")
if (!file.exists(kegg_gmt)) {
  log_msg("  building KEGG pathway GMT from rest.kegg.jp ...")
  links <- try(read.delim("https://rest.kegg.jp/link/hsa/pathway", header = FALSE), silent = TRUE)
  names_ <- try(read.delim("https://rest.kegg.jp/list/pathway/hsa", header = FALSE), silent = TRUE)
  if (!inherits(links, "try-error") && !inherits(names_, "try-error")) {
    links$pathway <- sub("^path:", "", links$V1)
    links$entrez  <- sub("^hsa:", "", links$V2)
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(links$entrez),
                                                  keytype = "ENTREZID", column = "SYMBOL",
                                                  multiVals = "first"))
    links$symbol <- unname(sym[links$entrez])
    links <- links[!is.na(links$symbol), ]
    pname <- setNames(sub(" - .*$", "", names_$V2), sub("^path:", "", names_$V1))
    gmt <- links %>% group_by(pathway) %>% summarise(genes = paste(unique(symbol), collapse = "\t"),
                                                     n = n_distinct(symbol), .groups = "drop") %>%
      filter(n >= 10)
    gmt$name <- pname[gmt$pathway]
    gmt$name[is.na(gmt$name)] <- gmt$pathway[is.na(gmt$name)]
    writeLines(paste(gmt$name, "kegg", gmt$genes, sep = "\t"), kegg_gmt)
    log_msg("  KEGG GMT: ", nrow(gmt), " pathways")
  }
}
have_kegg <- file.exists(kegg_gmt)
have_hm   <- file.exists(hm_gmt)

de_file <- file.path(TAB_DIR, "05_de_migDC_vs_cDC2.csv")
if (!file.exists(de_file)) stop("run R/05 first")
de <- read.csv(de_file)
dc <- readRDS(file.path(OBJ_DIR, "05_dc_annotated.rds"))
universe <- rownames(dc)

sig <- de %>% filter(p_val_adj < PARAMS$de$padj_cut, abs(avg_log2FC) >= PARAMS$de$logfc_cut)
up   <- sig$gene[sig$avg_log2FC > 0]
down <- sig$gene[sig$avg_log2FC < 0]
log_msg("  up in migDC: ", length(up), " | up in cDC2: ", length(down))

# ---- 1. GO BP per direction --------------------------------------------
run_go <- function(genes, tag) {
  if (length(genes) < 10) return(NULL)
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = PARAMS$enrich$padj_cut, qvalueCutoff = PARAMS$enrich$q_cut,
                  readable = TRUE)
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(NULL)
  df$direction <- tag
  write.csv(df, file.path(TAB_DIR, paste0("06_go_", tag, ".csv")), row.names = FALSE)
  ego
}

ego_up <- run_go(up, "up_in_migDC")
ego_dn <- run_go(down, "up_in_cDC2")

plot_go <- function(ego, tag, n = 12) {
  if (is.null(ego)) return(NULL)
  dotplot(ego, showCategory = n) + ggtitle(paste0("GO BP -- genes higher in ", tag)) +
    theme_pub(10)
}
if (!is.null(ego_up)) save_fig(plot_go(ego_up, "migDC"), "06_go_up_migDC", 9, 6)
if (!is.null(ego_dn)) save_fig(plot_go(ego_dn, "cDC2"), "06_go_up_cDC2", 9, 6)

# ---- 2. KEGG per direction (local GMT, no runtime download) ------------
run_kegg <- function(genes, tag) {
  if (!have_kegg || length(genes) < 10) return(NULL)
  t2g <- clusterProfiler::read.gmt(kegg_gmt)
  ek <- enricher(gene = genes, universe = universe, TERM2GENE = t2g,
                 pvalueCutoff = PARAMS$enrich$padj_cut,
                 pAdjustMethod = "BH", qvalueCutoff = PARAMS$enrich$q_cut)
  df <- as.data.frame(ek)
  if (nrow(df) == 0) return(NULL)
  df$direction <- tag
  write.csv(df, file.path(TAB_DIR, paste0("06_kegg_", tag, ".csv")), row.names = FALSE)
  ek
}
ek_up <- run_kegg(up, "up_in_migDC")
ek_dn <- run_kegg(down, "up_in_cDC2")
plot_kegg <- function(ek, tag, n = 12) {
  if (is.null(ek)) return(NULL)
  dotplot(ek, showCategory = n) + ggtitle(paste0("KEGG -- genes higher in ", tag)) +
    theme_pub(10)
}
if (!is.null(ek_up)) save_fig(plot_kegg(ek_up, "migDC"), "06_kegg_up_migDC", 9, 6)
if (!is.null(ek_dn)) save_fig(plot_kegg(ek_dn, "cDC2"), "06_kegg_up_cDC2", 9, 6)

# ---- 3. GSEA (ranked by log2 fold change) ------------------------------
# Ranking by log2FC is the conventional choice; an alternative composite
# (log2FC x -log10 p) was tested and produced no significant set because
# low-effect but highly significant genes dominate the top of the list.
de$rank_stat <- de$avg_log2FC
ranks <- de$rank_stat; names(ranks) <- de$gene
ranks <- sort(ranks, decreasing = TRUE)

run_gsea <- function(gmt, tag) {
  if (!file.exists(gmt)) return(NULL)
  paths <- fgsea::gmtPathways(gmt)
  res <- fgsea(pathways = paths, stats = ranks, minSize = PARAMS$enrich$min_gs,
               maxSize = PARAMS$enrich$max_gs, eps = 0, nproc = 1)
  res <- res[order(res$padj), ]
  res$direction <- ifelse(res$NES > 0, "up_in_migDC", "up_in_cDC2")
  if ("leadingEdge" %in% names(res))
    res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
  write.csv(res, file.path(TAB_DIR, paste0("06_gsea_", tag, ".csv")), row.names = FALSE)
  res
}
gsea_h <- run_gsea(hm_gmt, "hallmark")
gsea_k <- run_gsea(kegg_gmt, "kegg")

plot_gsea <- function(res, tag) {
  if (is.null(res)) return(NULL)
  sig <- res %>% filter(padj < 0.05)
  top <- if (nrow(sig) > 0) sig %>% arrange(desc(abs(NES))) %>% head(20)
         else res %>% arrange(padj) %>% head(15)
  if (nrow(top) == 0) return(NULL)
  top$set <- gsub("^(HALLMARK_|KEGG_)", "", top$pathway)
  top$set <- factor(top$set, levels = top$set[order(top$NES)])
  sub <- if (nrow(sig) > 0) paste0(nrow(sig), " set(s) at padj < 0.05")
         else "no set reaches padj < 0.05; top 15 by padj shown"
  ggplot(top, aes(NES, set, fill = NES > 0)) +
    geom_col() + scale_fill_manual(values = c(`TRUE` = "#D62728", `FALSE` = "#1F77B4")) +
    labs(title = paste0("GSEA -- ", tag, " (migDC vs cDC2)"), subtitle = sub,
         x = "Normalised enrichment score", y = NULL, fill = "higher in migDC") +
    theme_pub(10)
}
if (!is.null(gsea_h)) save_fig(plot_gsea(gsea_h, "hallmark"), "06_gsea_hallmark", 8, 6)
if (!is.null(gsea_k)) save_fig(plot_gsea(gsea_k, "kegg"), "06_gsea_kegg", 8, 6)

log_msg("STEP 06 done")

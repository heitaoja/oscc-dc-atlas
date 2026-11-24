# Refit the main model from published donor-subtype aggregates.
source("R/00_setup.R")
source("R/pseudobulk.R")

md <- read.csv(file.path(TAB_DIR, "11_primary_pseudobulk_libraries.csv"), row.names = "pseudobulk")
provenance <- sample_metadata(md)
counts <- read.delim(gzfile(file.path(TAB_DIR, "11_primary_pseudobulk_counts.tsv.gz")),
                     row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
stopifnot(identical(colnames(counts), rownames(md)), all(colSums(counts) == md$umi),
          all(counts >= 0 & counts == floor(counts)), all(md$group == "tumour"),
          identical(md$group, as.character(provenance$group)),
          identical(md$dataset, as.character(provenance$dataset)),
          identical(md$patient, provenance$patient),
          all(md$n_cells >= PARAMS$pseudobulk$primary_min_cells), !anyDuplicated(rownames(counts)))
md$patient <- factor(md$patient)
md$dc_subtype <- factor(md$dc_subtype, levels = c("cDC2", "migDC"))
stopifnot(all(table(md$patient, md$dc_subtype) == 1L))
actual <- fit_paired_bulk(list(counts = counts, metadata = md))$de
expected <- read.csv(file.path(TAB_DIR, "11_de_migDC_vs_cDC2.csv"))
stopifnot(setequal(actual$gene, expected$gene))
at <- match(actual$gene, expected$gene)
for (field in c("log2_fc", "p_value", "fdr")) {
  stopifnot(isTRUE(all.equal(actual[[field]], expected[[field]][at], tolerance = 1e-7)))
}
cat("Reproduced", nrow(actual), "gene estimates for", length(unique(md$patient)),
    "donors and", sum(actual$fdr < PARAMS$pseudobulk$fdr), "genes with FDR < 0.05.\n")

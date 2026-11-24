# Stage 99: audit the reviewed atlas outputs and the manuscript figure exports.
source("R/verify_primary.R")

log_msg("VERIFY ALL -- discovery, external validation, TME and manuscript figures")
cfg <- PARAMS$pseudobulk
require_check <- function(ok, label) {
  if (!isTRUE(ok)) stop("Verification failed: ", label, call. = FALSE)
}
read_table <- function(name) read.csv(file.path(TAB_DIR, paste0(name, ".csv")))
finite_fields <- function(rows, fields) {
  nrow(rows) > 0L && all(fields %in% names(rows)) && all(is.finite(as.matrix(rows[fields])))
}
same_numbers <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-9, check.attributes = FALSE))
nonempty <- function(paths) {
  info <- file.info(paths)
  good <- !is.na(info$size) & !is.na(info$isdir) & !info$isdir & info$size > 0
  require_check(all(good), paste("missing or empty files", paste(paths[!good], collapse = ", ")))
}

table_files <- c(
  "03_celltype_by_sample.csv", "05_dc_subtype_by_sample.csv", "05_de_migDC_vs_cDC2.csv",
  "09_dc_fraction_by_sample.csv", "09_dc_subtype_fraction.csv", "09_dc_subtype_pooled.csv",
  "09_paired_tissue_descriptive.csv", "09_cDC1_cDC2_balance.csv",
  "11_balanced_cell_sensitivity.csv", "11_de_cDC1_vs_cDC2.csv", "11_de_LClike_vs_cDC2.csv",
  "11_de_migDC_vs_cDC1.csv", "11_de_migDC_vs_cDC2.csv", "11_de_migDC_vs_LClike.csv",
  "11_donor_eligibility.csv", "11_gene_detection.csv", "11_input_snapshot.csv",
  "11_label_agreement.csv", "11_leave_one_donor_out.csv", "11_marker_gate_eligibility.csv",
  "11_marker_gate_programs.csv", "11_marker_gate_sensitivity.csv", "11_paired_hallmark_camera.csv",
  "11_panel_donor_effects.csv", "11_panel_model_results.csv", "11_primary_pseudobulk_counts.tsv.gz",
  "11_primary_pseudobulk_libraries.csv", "11_program_by_cohort.csv", "11_program_donor_effects.csv",
  "11_program_fry.csv", "11_program_gene_donor_effects.csv", "11_program_summary.csv",
  "11_sensitivity.csv", "11_sessionInfo.txt",
  "12_classification.csv", "12_cohort_inventory.csv", "12_combined_program_summary.csv",
  "12_cross_cohort_genes.csv", "12_cross_cohort_summary.csv", "12_donor_eligibility.csv",
  "12_inputs.csv", "12_panel_donor_effects.csv", "12_program_donor_effects.csv",
  "12_program_fry.csv", "12_program_summary.csv", "12_sensitivity.csv", "12_sessionInfo.txt",
  "12_shared_program_donor_effects.csv",
  "13_tme_cell_proportions.csv", "13_tme_donor_profiles.csv", "13_tme_ligand_receptor_correlations.csv"
)
manuscript_figures <- c(
  "Fig1_discovery_atlas", "Fig2_pseudobulk_de_and_robustness", "Fig3_programs_and_inflammation_paradox",
  "Fig4_lineage_and_tme_context", "Fig5_external_validation"
)
figure_stems <- c(
  "09_dc_infiltration", "09_dc_subtype_composition", "09_dc_subtype_pooled",
  "09_dc_subtype_tumor_vs_normal", "09_cdc1_cdc2_ratio",
  "11_paired_dc_programs", "11_dc_lineage_comparators",
  "12_cross_cohort_effects", "13_tme_covariation", manuscript_figures
)
figure_paths <- file.path(FIG_DIR, c(paste0(figure_stems, ".png"), paste0(figure_stems, ".pdf")))
nonempty(c(file.path(TAB_DIR, table_files), figure_paths))

png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
pdf_signature <- charToRaw("%PDF-")
for (path in figure_paths) {
  signature <- if (endsWith(path, ".png")) png_signature else pdf_signature
  require_check(identical(readBin(path, "raw", n = length(signature)), signature), basename(path))
}
png_dpi <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  readBin(con, "raw", n = length(png_signature))
  repeat {
    size <- readBin(con, "integer", n = 1L, size = 4L, endian = "big")
    type <- rawToChar(readBin(con, "raw", n = 4L))
    if (!length(size) || type %in% c("IDAT", "IEND") || size < 0L) return(NA_real_)
    if (type == "pHYs" && size == 9L) {
      pixels_per_metre <- readBin(con, "integer", n = 2L, size = 4L, endian = "big")
      unit <- readBin(con, "integer", n = 1L, size = 1L, signed = FALSE)
      return(if (identical(unit, 1L)) pixels_per_metre * 0.0254 else NA_real_)
    }
    seek(con, where = size + 4L, origin = "current")
  }
}
for (name in manuscript_figures) {
  dpi <- png_dpi(file.path(FIG_DIR, paste0(name, ".png")))
  require_check(length(dpi) == 2L && all(round(dpi) == 300), paste(name, "PNG resolution"))
}
log_msg("  Key stage outputs and five 300 dpi PNG/PDF pairs passed.")

atlas <- read_table("03_celltype_by_sample")
subtypes <- read_table("05_dc_subtype_by_sample")
fractions <- read_table("09_dc_fraction_by_sample")
paired_tissues <- read_table("09_paired_tissue_descriptive")
pooled_dc <- read_table("09_dc_subtype_pooled")
provenance <- sample_metadata(fractions)
require_check(sum(atlas$Freq) == 85164L && length(unique(atlas$Var1)) == 24L &&
              sum(subtypes$n) == 1372L && nrow(pooled_dc) == 6L && sum(pooled_dc$n) == 1372L &&
              nrow(fractions) == 10L && sum(fractions$n_atlas) == 85164L &&
              sum(fractions$n_retained_dc) == 1372L, "atlas and DC composition")
require_check(identical(fractions$group, as.character(provenance$group)) &&
              identical(fractions$patient, provenance$patient) &&
              identical(fractions$dataset, as.character(provenance$dataset)) &&
              setequal(paired_tissues$patient, c("P7", "P8")) &&
              all(table(subset(fractions, patient %in% paired_tissues$patient)[c("patient", "group")]) == 1L) &&
              finite_fields(paired_tissues, c("frac_dc_normal", "frac_dc_tumour",
                                             "frac_retained_dc_normal", "frac_retained_dc_tumour")),
              "corrected matched-tissue labels and fractions")
log_msg("  Stage 09 composition and two matched tissue pairs passed.")

program_order <- c("Maturation", "NF-kB feedback", "Regulatory expression", "MHC class II")
positive_programs <- head(program_order, 3L)
primary_programs <- subset(read_table("11_program_summary"), analysis == "primary")
primary_panel <- subset(read_table("11_panel_model_results"),
                        analysis == "primary" & gene %in% c("PDCD1LG2", "CD200"))
sensitivity <- read_table("11_sensitivity")
strict <- subset(sensitivity, analysis == "strict_filter_min10")
require_check(nrow(actual) == 4151L && sum(actual$fdr < cfg$fdr) == 933L &&
              length(unique(md$patient)) == 5L, "primary paired DE totals")
require_check(nrow(primary_programs) == 4L && setequal(primary_programs$program, program_order) &&
              all(primary_programs$n_donors == 5L) &&
              finite_fields(primary_programs, c("mean_effect", "p_paired_t", "p_fry", "fdr_fry")) &&
              subset(primary_programs, program == "Regulatory expression")$n_genes >= 5L,
              "four primary programs and Regulatory expression coverage")
significant_programs <- subset(primary_programs, program %in% positive_programs)
tests <- c("p_paired_t", "fdr_paired_t", "p_fry", "fdr_fry")
require_check(finite_fields(significant_programs, tests) && all(as.matrix(significant_programs[tests]) < cfg$fdr),
              "primary paired t and fry results")
require_check(nrow(primary_panel) == 2L && setequal(primary_panel$gene, c("PDCD1LG2", "CD200")) &&
              finite_fields(primary_panel, c("log2_fc", "fdr")), "restored PDCD1LG2 and CD200")
require_check(nrow(strict) == 1L && finite_fields(strict, "rho_with_primary") &&
              strict$rho_with_primary > 0.99, "strict filter rank concordance")
for (contrast in c("migDC_vs_LClike", "LClike_vs_cDC2")) {
  lc <- read_table(paste0("11_de_", contrast))
  require_check(nrow(lc) > 1000L && !anyDuplicated(lc$gene) &&
                finite_fields(lc, c("log2_fc", "fdr")), paste("LC-like contrast", contrast))
}
log_msg("  Stage 11 model, held-out programs, marker estimates and LC-like contrasts passed.")


# Validate the real GEO inputs and their myeloid APC cache against the audited snapshot.
geo_fingerprints <- c(
  "GSE181919_Barcode_metadata.txt.gz" = "17c1447050efbc0e83ffdf3d27acdf66",
  "GSE181919_UMI_counts.txt.gz" = "6f416ecec9c25f7e00aa481cfb553970"
)
geo_paths <- file.path(RAW_DIR, "GSE181919", names(geo_fingerprints))
cache_path <- file.path(OBJ_DIR, "12_gse181919_oral_myeloid_counts.rds")
nonempty(c(geo_paths, cache_path))
input_record <- read_table("12_inputs")
require_check(nrow(input_record) == 2L && setequal(input_record$file, names(geo_fingerprints)),
              "GEO input ledger schema")
input_record <- input_record[match(names(geo_fingerprints), input_record$file), ]
require_check(all(input_record$available) &&
              identical(unname(tools::md5sum(geo_paths)), unname(geo_fingerprints)) &&
              identical(input_record$md5, unname(geo_fingerprints)) &&
              identical(file.size(geo_paths), c(339614, 127878601)) &&
              identical(as.numeric(input_record$bytes), file.size(geo_paths)), "verified GEO source fingerprints")
metadata <- read.delim(gzfile(geo_paths[1]), row.names = NULL)
oral <- subset(metadata, tissue.type == "CA" & subsite == "OC")
apc <- subset(oral, cell.type %in% c("Dendritic.cells", "Macrophages"))
myeloid <- readRDS(cache_path)
require_check(inherits(myeloid$counts, "dgCMatrix") && identical(dim(myeloid$counts), c(20000L, 2902L)) &&
              identical(myeloid$md5, input_record$md5[2]) && identical(myeloid$bytes, file.size(geo_paths[2])) &&
              identical(myeloid$barcodes, apc$row.names) && identical(colnames(myeloid$counts), apc$row.names) &&
              !anyDuplicated(myeloid$source_barcodes) &&
              identical(sub("\\.([0-9]+)$", "-\\1", myeloid$source_barcodes), apc$row.names) &&
              all(is.finite(myeloid$counts@x) & myeloid$counts@x >= 0 &
                    myeloid$counts@x == floor(myeloid$counts@x)), "integer myeloid cache and barcode provenance")
external_inventory <- read_table("12_cohort_inventory")
require_check(nrow(external_inventory) == 10L && sum(external_inventory$n_cells) == 10397L &&
              sum(external_inventory$n_author_dc) == 688L && sum(external_inventory$n_author_macrophages) == 2214L &&
              sum(external_inventory$n_myeloid_apc) == nrow(apc), "external cohort inventory")
gates <- read.csv(file.path("config", "external_dc_gates.csv"))
annotation <- read.csv(file.path("config", "dc_programs.csv"))
annotation_genes <- union(gates$gene, annotation$gene[annotation$role == "annotation"])
require_check(!any(c("CD274", "IDO1") %in% gates$gene), "extended Regulatory panel gate independence")
classification <- read_table("12_classification")
for (cutoff in c(1L, 2L)) {
  labels <- classify_dc(myeloid$counts, gates, n_detected = cutoff)
  reference <- table(apc$sample.id, apc$cell.type, labels)
  published <- xtabs(n ~ sample + author_cell_type + subtype,
                      subset(classification, gate_min_detected == cutoff))
  require_check(identical(unname(dimnames(reference)), unname(dimnames(published))) &&
                same_numbers(as.vector(reference), as.vector(published)),
                paste("external gate counts at", cutoff, "detected markers"))
}
strict_counts <- xtabs(n ~ author_cell_type + subtype, subset(classification, gate_min_detected == 2L))
require_check(strict_counts["Macrophages", "cDC2"] == 243L && sum(strict_counts[, "cDC2"]) == 245L &&
              strict_counts["Dendritic.cells", "migDC"] == 209L &&
              strict_counts["Dendritic.cells", "cDC1"] == 52L, "cross-annotation DC counts")
log_msg("  GSE181919 source, myeloid cache and author-label audit passed.")

hybrid_id <- "hybrid_author_migDC_vs_cd163neg_cDC2"
external_n <- c(primary = 6L, min_5 = 7L, cd163_neg_min10 = 4L, cd163_neg_min5 = 7L,
                hybrid_author_migDC_vs_cd163neg_cDC2 = 4L, author_dc_only_min10 = 0L,
                author_dc_only_min5 = 0L, author_dc_migDC_vs_cDC1_min3 = 6L,
                migDC_vs_cDC1 = 4L, cDC1_vs_cDC2 = 3L)
external_sensitivity <- read_table("12_sensitivity")
replication <- read_table("12_cross_cohort_summary")
replication <- replication[match(names(external_n), replication$analysis), ]
require_check(nrow(external_sensitivity) == length(external_n) && !anyNA(replication$analysis) &&
              !anyDuplicated(external_sensitivity$analysis) &&
              identical(replication$n_donors, unname(external_n)) &&
              identical(external_sensitivity$n_donors[match(names(external_n), external_sensitivity$analysis)],
                        unname(external_n)),
              "external strict gate donor eligibility")
fitted <- subset(replication, status == "fitted")
not_fitted <- subset(replication, status == "insufficient_donors")
require_check(setequal(fitted$analysis, names(external_n)[external_n >= cfg$min_donors]) &&
              setequal(not_fitted$analysis, names(external_n)[external_n < cfg$min_donors]) &&
              all(is.na(not_fitted$rho)) && all(is.na(not_fitted$concordant_hit_rate)) &&
              finite_fields(fitted, c("rho", "n_shared", "concordant_hit_rate", "external_fdr_hits")),
              "external fitted rows and intentional eligibility gaps")
external_de_files <- paste0("12_de_", fitted$analysis, ".csv")
nonempty(file.path(TAB_DIR, external_de_files))
replication_genes <- read_table("12_cross_cohort_genes")
require_check(!anyDuplicated(replication_genes[c("analysis", "gene")]) &&
              setequal(replication_genes$analysis, fitted$analysis), "external shared gene keys")
for (i in seq_len(nrow(fitted))) {
  model <- fitted[i, ]
  genes <- subset(replication_genes, analysis == model$analysis)
  require_check(identical(genes$annotation_gene, genes$gene %in% annotation_genes),
                paste("assignment flag", model$analysis))
  genes <- subset(genes, !annotation_gene)
  hit <- genes$fdr_discovery < cfg$fdr
  concordant <- sign(genes$log2_fc_discovery) == sign(genes$log2_fc_external)
  require_check(model$n_shared == nrow(genes) && same_numbers(model$concordance_rate, mean(concordant)) &&
                same_numbers(model$rho, cor(genes$log2_fc_discovery, genes$log2_fc_external, method = "spearman")) &&
                model$n_discovery_hits == sum(hit) && model$concordant_hits == sum(hit & concordant) &&
                same_numbers(model$concordant_hit_rate, mean(concordant[hit])) &&
                model$external_fdr_hits == sum(hit & concordant & genes$fdr_external < cfg$fdr),
                paste("cross-cohort metrics", model$analysis))
}
external_primary <- subset(fitted, analysis == "primary")
hybrid <- subset(fitted, analysis == hybrid_id)
require_check(round(external_primary$rho, 3) == 0.671 && round(100 * external_primary$concordant_hit_rate, 1) == 93.4 &&
              round(hybrid$rho, 3) == 0.751 && hybrid$concordant_hits == 860L &&
              hybrid$n_discovery_hits == 888L && hybrid$external_fdr_hits == 618L,
              "primary and hybrid replication estimates")
external_programs <- read_table("12_program_summary")
external_fry <- read_table("12_program_fry")
for (id in c("primary", hybrid_id)) {
  external <- subset(external_programs, analysis == id)
  require_check(nrow(external) == 5L && all(external$status == "fitted") &&
                "Regulatory expression (incl. CD274/IDO1)" %in% external$program &&
                finite_fields(external, c("mean_effect", "p_paired_t", "p_fry", "fdr_fry")) &&
                same_numbers(external$p_fry, subset(external_fry, analysis == id)$p_fry),
                paste(id, "external held-out and extended programs"))
}
log_msg("  All seven fitted external comparisons and their replication estimates passed.")

combined <- read_table("12_combined_program_summary")
shared <- read_table("12_shared_program_donor_effects")
joint_n <- external_n[c("primary", "min_5", "cd163_neg_min10", "cd163_neg_min5", hybrid_id)] + 5L
joint_fields <- c("mean_effect", "lower_95_t", "upper_95_t", "p_sign_flip", "fdr",
                  "p_paired_t", "fdr_paired_t", "p_fry_meta", "fdr_fry_meta")
require_check(nrow(combined) == 4L * length(joint_n) && setequal(combined$analysis, names(joint_n)) &&
              setequal(shared$analysis, names(joint_n)) && all(combined$status == "fitted") &&
              finite_fields(combined, joint_fields), "joint external model estimates")
for (id in names(joint_n)) {
  output <- subset(combined, analysis == id)
  effects <- subset(shared, analysis == id)
  require_check(nrow(output) == 4L && setequal(output$program, program_order) &&
                all(output$n_discovery == 5L) && all(output$n_external == external_n[[id]]) &&
                all(output$n_donors == joint_n[[id]]) &&
                !anyDuplicated(effects[c("cohort", "patient", "program")]) &&
                all(is.finite(effects$effect)), paste(id, "joint cohort size and shared effects"))
  recomputed <- summarize_programs(effects)
  recomputed <- recomputed[match(output$program, recomputed$program), ]
  core <- c("mean_effect", "lower_95_t", "upper_95_t", "positive_donors", "p_paired_t", "p_sign_flip")
  require_check(same_numbers(output[core], recomputed[core]), paste(id, "pooled effects and tests"))
  fisher_p <- with(output, pchisq(-2 * (log(p_fry_discovery) + log(p_fry_external)),
                                    df = 4L, lower.tail = FALSE))
  require_check(same_numbers(output$p_fry_meta, fisher_p) &&
                same_numbers(output$fdr_fry_meta, p.adjust(fisher_p, "BH")),
                paste(id, "Fisher meta-analysis"))
  positive <- subset(output, program %in% positive_programs)
  require_check(all(positive$positive_donors == joint_n[[id]]) &&
                all(positive$p_sign_flip == 2^(1 - joint_n[[id]])) &&
                all(positive$fdr_fry_meta < cfg$fdr), paste(id, "unanimous positive program effects"))
}
log_msg("  Stage 12 joint estimates passed, including 11/11 and 12/12 concordant program effects.")

tme <- read_table("13_tme_ligand_receptor_correlations")
profiles <- read_table("13_tme_donor_profiles")
tme_cells <- read_table("13_tme_cell_proportions")
tested <- subset(tme, status == "tested")
undercovered <- subset(tme, status == "insufficient_program_coverage")
require_check(nrow(profiles) == 8L && !anyDuplicated(profiles$patient) &&
              nrow(tme_cells) == 10L && sum(tme_cells$n_total) == 85164L &&
              nrow(tme) == 46L && nrow(tested) == 42L && nrow(undercovered) == 4L &&
              all(undercovered$predictor == "min_3_regulatory_expression") &&
              all(is.na(undercovered[c("rho", "p_value", "fdr")])) &&
              finite_fields(tested, c("n_donors", "rho", "p_value", "fdr")) &&
              !anyDuplicated(tme[c("predictor", "response")]), "TME cohort and genuine coverage gaps")
for (i in seq_len(nrow(tested))) {
  row <- tested[i, ]
  x <- profiles[[row$predictor]]
  y <- profiles[[row$response]]
  keep <- is.finite(x) & is.finite(y)
  result <- cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  require_check(row$n_donors == sum(keep) && same_numbers(row$rho, unname(result$estimate)) &&
                same_numbers(row$p_value, result$p.value), paste(row$predictor, row$response, "Spearman test"))
}
ccl17 <- subset(tested, predictor == "dc_CCL17" & response == "treg_CCR4")
ccl19 <- subset(tested, predictor == "dc_CCL19" & response == "tcell_CCR7")
require_check(same_numbers(tested$fdr, p.adjust(tested$p_value, "BH")) &&
              nrow(ccl17) == 1L && ccl17$n_donors == 8L && round(ccl17$rho, 2) == 0.93 &&
              round(ccl17$fdr, 3) == 0.036 && nrow(ccl19) == 1L && round(ccl19$rho, 2) == 0.88,
              "TME chemokine/receptor axes and multiple testing")
log_msg("  Stage 13 verified 42 TME tests and four documented program coverage gaps.")

documents <- c("report/manuscript.md", "report/report.md", "docs/methods.md", "README.md")
nonempty(documents)
legacy_patterns <- c(
  "external\\s+UMI\\s+file\\s+(?:was|is)\\s+unavailable",
  "replication\\s+remains\\s+pending",
  "\\b1,?556\\s+(?:(?:genes\\s+)?tested|genes)\\b",
  "\\b457\\s+(?:significant\\s+genes|(?:genes\\s+)?(?:met|passed|reached))\\b",
  "regulatory\\s+(?:expression\\s+)?(?:program|panel)\\s+had\\s+insufficient\\s+expression\\s+coverage",
  "1,?556\\s*\u4e2a\u57fa\u56e0",
  "457\\s*\u4e2a(?:\u663e\u8457\u57fa\u56e0|\u8fbe\u5230)",
  "\u8868\u8fbe\u590d\u73b0\u7b49\u5f85|\u4e0b\u8f7d\u672a\u5b8c\u6210"
)

markdown_targets <- function(lines, prefix) {
  destination <- "(?:<([^>]+)>|([^\\s)<]+))"
  pattern <- paste0(prefix, destination)
  matches <- unlist(regmatches(lines, gregexpr(pattern, lines, perl = TRUE)), use.names = FALSE)
  sub(pattern, "\\1\\2", matches, perl = TRUE)
}

joint_sizes <- unname(joint_n[c("primary", "min_5")])
core_values <- c(formatC(nrow(actual), format = "d", big.mark = ","),
                 as.character(sum(actual$fdr < cfg$fdr)),
                 sprintf("%.4f", ccl17$fdr),
                 sprintf("%.3f", c(external_primary$rho, hybrid$rho)),
                 paste(joint_sizes, joint_sizes, sep = "/"))
links_checked <- 0L
for (document in documents) {
  lines <- readLines(document, encoding = "UTF-8", warn = FALSE)
  prose <- paste(lines, collapse = " ")
  stale <- legacy_patterns[vapply(legacy_patterns, grepl, logical(1), x = prose,
                                  ignore.case = TRUE, perl = TRUE)]
  require_check(!length(stale), paste(document, "contains obsolete text", paste(stale, collapse = ", ")))
  missing_values <- core_values[!vapply(core_values, grepl, logical(1), x = prose, fixed = TRUE)]
  require_check(!length(missing_values),
                paste(document, "missing current summary values", paste(missing_values, collapse = ", ")))

  targets <- unique(c(markdown_targets(lines, "\\]\\(\\s*"),
                      markdown_targets(lines, "^ {0,3}\\[[^]]+\\]:\\s*")))
  local <- targets[!grepl("^(?:[[:alpha:]][[:alnum:].+-]*:|//|#)", targets, perl = TRUE)]
  local <- unique(utils::URLdecode(sub("[?#].*$", "", local)))
  absent <- local[!file.exists(file.path(dirname(document), local))]
  require_check(!length(absent), paste(document, "broken relative links", paste(absent, collapse = ", ")))
  links_checked <- links_checked + length(local)
}
log_msg("  Documentation passed: ", length(documents), " files, ", links_checked,
        " local links and shared summary values verified.")
cat("ALL VERIFICATION CHECKS PASSED\n")

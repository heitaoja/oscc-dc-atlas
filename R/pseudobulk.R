# Shared donor-paired expression analysis on original UMI counts.

aggregate_umi <- function(counts, md) {
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      !identical(colnames(counts), rownames(md))) {
    stop("Expected unique gene names and count columns aligned with cell metadata.")
  }
  fields <- c("sample", "patient", "dataset", "group", "dc_subtype")
  if (!all(fields %in% names(md)) || anyNA(md[fields])) stop("Incomplete cell metadata.")
  values <- if (inherits(counts, "sparseMatrix")) counts@x else counts
  if (any(!is.finite(values) | values < 0 | values != floor(values))) {
    stop("Pseudobulk requires non-negative integer UMI counts.")
  }
  key <- paste(md$sample, md$dc_subtype, sep = "/")
  keys <- unique(key)
  at <- match(key, keys)
  bulk_info <- md[match(keys, key), fields, drop = FALSE]
  rownames(bulk_info) <- keys
  bulk_info$n_cells <- tabulate(at, length(keys))
  cell_to_bulk <- Matrix::sparseMatrix(i = seq_along(at), j = at, x = 1,
                                     dims = c(length(at), length(keys)),
                                     dimnames = list(rownames(md), keys))
  bulk_counts <- counts %*% cell_to_bulk
  list(counts = bulk_counts, metadata = bulk_info)
}

paired_cell_counts <- function(bulk, samples, positive, negative, min_cells) {
  samples <- samples[, c("sample", "patient", "dataset", "group")]
  if (anyDuplicated(samples$sample) || anyDuplicated(samples[c("patient", "group")])) {
    stop("Select one library per donor per tissue before forming paired contrasts.")
  }
  md <- bulk$metadata
  n <- setNames(md$n_cells, paste(md$sample, md$dc_subtype, sep = "/"))
  count <- function(label) {
    x <- unname(n[paste(samples$sample, label, sep = "/")])
    replace(x, is.na(x), 0L)
  }
  samples$positive <- positive
  samples$negative <- negative
  samples$n_positive <- count(positive)
  samples$n_negative <- count(negative)
  samples$min_cells <- min_cells
  samples$eligible <- pmin(samples$n_positive, samples$n_negative) >= min_cells
  samples
}

paired_bulk <- function(bulk, eligible, positive, negative) {
  md <- bulk$metadata
  keep <- md$sample %in% eligible$sample & md$dc_subtype %in% c(negative, positive)
  md <- md[keep, , drop = FALSE]
  md$patient <- factor(md$patient, levels = eligible$patient)
  md$dc_subtype <- factor(md$dc_subtype, levels = c(negative, positive))
  md <- md[order(md$patient, md$dc_subtype), , drop = FALSE]
  if (nrow(md) != 2L * nrow(eligible) ||
      !all(table(md$patient, md$dc_subtype) == 1L)) {
    stop("Expected one pseudobulk per donor per subtype.")
  }
  list(counts = bulk$counts[, rownames(md), drop = FALSE], metadata = md)
}

fit_paired_bulk <- function(pair, norm = "TMM") {
  md <- pair$metadata
  design <- model.matrix(~ patient + dc_subtype, md)
  stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design))
  y <- edgeR::DGEList(as.matrix(pair$counts), samples = md)
  keep <- edgeR::filterByExpr(y, group = md$dc_subtype)
  if (!any(keep)) stop("No genes pass the pseudobulk expression filter.")
  y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE], method = norm)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = ncol(design))$table
  de <- data.frame(gene = rownames(test), log2_fc = test$logFC, log_cpm = test$logCPM,
                   ql_f = test$F, p_value = test$PValue, fdr = p.adjust(test$PValue, "BH"))
  list(de = de[order(de$p_value, de$gene), ], y = y)
}

paired_effects <- function(y, positive, negative) {
  md <- y$samples
  values <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  hi <- which(md$dc_subtype == positive)
  lo <- which(md$dc_subtype == negative)
  stopifnot(identical(md$patient[hi], md$patient[lo]))
  differences <- values[, hi, drop = FALSE] - values[, lo, drop = FALSE]
  colnames(differences) <- as.character(md$patient[hi])
  differences
}

read_dc_programs <- function(path = file.path("config", "dc_programs.csv")) {
  definition <- read.csv(path)
  stopifnot(!anyDuplicated(definition[c("role", "program", "gene")]),
            !anyNA(definition),
            !any(definition$gene[definition$role == "held_out"] %in%
                 definition$gene[definition$role == "annotation"]))
  split(definition$gene[definition$role == "held_out"],
        definition$program[definition$role == "held_out"])
}

classify_dc <- function(counts, gates, n_detected) {
  gate_sets <- split(gates$gene, gates$subtype)
  if (!all(gates$gene %in% rownames(counts))) {
    stop("Missing gate genes: ", paste(setdiff(gates$gene, rownames(counts)), collapse = ", "))
  }
  pass <- do.call(cbind, lapply(gate_sets, function(genes) {
    Matrix::colSums(counts[genes, , drop = FALSE] > 0) >= n_detected
  }))
  n_pass <- rowSums(pass)
  picked <- colnames(pass)[max.col(pass, ties.method = "first")]
  picked[n_pass == 0L] <- "unassigned"
  picked[n_pass > 1L] <- "ambiguous"
  picked
}

gene_effects_long <- function(effects, genes) {
  out <- as.data.frame(as.table(effects[intersect(genes, rownames(effects)), , drop = FALSE]),
                       stringsAsFactors = FALSE, responseName = "log2_cpm_difference")
  names(out)[1:2] <- c("gene", "patient")
  out
}

program_differences <- function(differences, programs, min_coverage) {
  do.call(rbind, lapply(names(programs), function(program) {
    genes <- programs[[program]]
    used <- intersect(genes, rownames(differences))
    adequate <- length(used) >= ceiling(min_coverage * length(genes))
    data.frame(program = program, patient = colnames(differences),
               n_genes = length(used), n_target_genes = length(genes),
               genes = paste(used, collapse = " "),
               effect = if (adequate) colMeans(differences[used, , drop = FALSE]) else NA_real_)
  }))
}

exact_sign_flip <- function(x) {
  observed <- abs(sum(x))
  distribution <- 0
  for (effect in x) distribution <- c(distribution + effect, distribution - effect)
  mean(abs(distribution) >= observed - 1e-10 * max(1, observed))
}

summarize_programs <- function(effects) {
  out <- do.call(rbind, lapply(split(effects, effects$program), function(program) {
    x <- program$effect[is.finite(program$effect)]
    n <- length(x)
    mean_effect <- if (n) mean(x) else NA_real_
    half_width <- if (n >= 2L) qt(0.975, n - 1L) * sd(x) / sqrt(n) else NA_real_
    data.frame(program = program$program[1], n_donors = n,
               n_genes = program$n_genes[1], n_target_genes = program$n_target_genes[1],
               mean_effect = mean_effect, lower_95_t = mean_effect - half_width,
               upper_95_t = mean_effect + half_width, positive_donors = sum(x > 0),
               p_sign_flip = if (n >= 2L && n <= 20L) exact_sign_flip(x) else NA_real_)
  }))
  out$fdr <- p.adjust(out$p_sign_flip, "BH", n = nrow(out))
  out
}

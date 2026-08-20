#!/usr/bin/env Rscript

# Apply the frozen Wilkerson 2012 506-gene LUAD nearest-centroid classifier.
# This script deliberately separates TCGA and DRA001846 preprocessing and
# never combines expression intensities across cohorts.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Usage: classify_luad_states.R <state_audit_root> <figure1_root> <figure2_root> [--tcga-only]")
}
state_root <- normalizePath(args[[1]], mustWork = TRUE)
figure1_root <- normalizePath(args[[2]], mustWork = TRUE)
figure2_root <- normalizePath(args[[3]], mustWork = TRUE)
tcga_only <- length(args) == 4L && identical(args[[4]], "--tcga-only")
if (length(args) == 4L && !tcga_only) stop("Only optional flag is --tcga-only")

out_dir <- file.path(state_root, "results", "state_classifier")
audit_dir <- file.path(state_root, "audit", "state_classifier")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  answer <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(answer)) return(NA_character_)
  strsplit(answer[[1]], "[[:space:]]+")[[1]][[1]]
}

centroid_path <- file.path(
  state_root, "reference", "lungCancerSubtypes", "lung_adenocarcinoma_subtypes",
  "wilkerson.2012.LAD.predictor.centroids.csv"
)
if (!file.exists(centroid_path)) stop("Missing frozen centroid file: ", centroid_path)
centroids <- fread(centroid_path, check.names = FALSE)
setnames(centroids, 1L, "gene")
expected_centroids <- c("bronchioid", "magnoid", "squamoid")
if (nrow(centroids) != 506L || !all(expected_centroids %in% names(centroids))) {
  stop("Frozen Wilkerson centroid file is not 506 genes x 3 centroids")
}
if (anyDuplicated(centroids$gene)) stop("Duplicated predictor genes in centroid file")
centroid_matrix <- as.matrix(centroids[, ..expected_centroids])
storage.mode(centroid_matrix) <- "double"
rownames(centroid_matrix) <- centroids$gene
subtype_map <- c(bronchioid = "TRU", magnoid = "PP", squamoid = "PI")
lineage_genes <- c("NKX2-1", "NAPSA", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC")

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  cor(x[keep], y[keep], method = "pearson")
}

classify_expression <- function(log2_expression, cohort, input_scale) {
  if (is.null(rownames(log2_expression)) || is.null(colnames(log2_expression))) {
    stop("Expression matrix requires gene row names and sample column names")
  }
  common <- intersect(rownames(centroid_matrix), rownames(log2_expression))
  if (length(common) < 405L) {
    stop(cohort, " has only ", length(common), "/506 predictor genes; frozen minimum is 405")
  }
  x <- log2_expression[common, , drop = FALSE]
  storage.mode(x) <- "double"
  gene_medians <- apply(x, 1L, median, na.rm = TRUE)
  centered <- sweep(x, 1L, gene_medians, "-")
  correlations <- matrix(
    NA_real_, nrow = ncol(centered), ncol = length(expected_centroids),
    dimnames = list(colnames(centered), expected_centroids)
  )
  available <- integer(ncol(centered))
  for (sample_index in seq_len(ncol(centered))) {
    values <- centered[, sample_index]
    available[[sample_index]] <- sum(is.finite(values))
    for (centroid_name in expected_centroids) {
      correlations[sample_index, centroid_name] <- safe_cor(
        values, centroid_matrix[common, centroid_name]
      )
    }
  }

  result_rows <- lapply(seq_len(nrow(correlations)), function(i) {
    scores <- correlations[i, ]
    ordered <- names(sort(scores, decreasing = TRUE, na.last = TRUE))
    classifiable <- available[[i]] >= 405L && is.finite(scores[[ordered[[1]]]])
    data.table(
      sample_id = rownames(correlations)[[i]],
      cohort = cohort,
      input_scale = input_scale,
      available_predictor_genes = available[[i]],
      coverage_fraction = available[[i]] / 506,
      correlation_bronchioid = scores[["bronchioid"]],
      correlation_magnoid = scores[["magnoid"]],
      correlation_squamoid = scores[["squamoid"]],
      top_centroid = if (classifiable) ordered[[1]] else "UNCLASSIFIABLE",
      second_centroid = if (classifiable) ordered[[2]] else "UNCLASSIFIABLE",
      correlation_margin = if (classifiable) scores[[ordered[[1]]]] - scores[[ordered[[2]]]] else NA_real_,
      subtype = if (classifiable) unname(subtype_map[[ordered[[1]]]]) else "UNCLASSIFIABLE"
    )
  })
  result <- rbindlist(result_rows)

  # Orthogonal NKX2-1/alveolar-lineage audit uses uncentered expression within
  # each cohort. Missing genes remain explicit and never receive imputed values.
  present_lineage <- intersect(lineage_genes, rownames(log2_expression))
  lineage_z <- matrix(
    NA_real_, nrow = length(present_lineage), ncol = ncol(log2_expression),
    dimnames = list(present_lineage, colnames(log2_expression))
  )
  for (gene in present_lineage) {
    values <- as.numeric(log2_expression[gene, ])
    gene_sd <- sd(values, na.rm = TRUE)
    if (is.finite(gene_sd) && gene_sd > 0) {
      lineage_z[gene, ] <- (values - mean(values, na.rm = TRUE)) / gene_sd
    } else {
      lineage_z[gene, is.finite(values)] <- 0
    }
  }
  lineage_score <- if (length(present_lineage) >= 4L) {
    colMeans(lineage_z, na.rm = TRUE)
  } else {
    rep(NA_real_, ncol(log2_expression))
  }
  nkx_values <- if ("NKX2-1" %in% rownames(log2_expression)) {
    as.numeric(log2_expression["NKX2-1", ])
  } else {
    rep(NA_real_, ncol(log2_expression))
  }
  nkx_median <- if (any(is.finite(nkx_values))) median(nkx_values, na.rm = TRUE) else NA_real_
  lineage_by_sample <- data.table(
    sample_id = colnames(log2_expression),
    lineage_genes_available = length(present_lineage),
    nkx2_1_log2_expression = nkx_values,
    nkx2_1_cohort_median = nkx_median,
    nkx2_1_above_median = ifelse(is.finite(nkx_values), nkx_values > nkx_median, NA),
    alveolar_lineage_z_mean = as.numeric(lineage_score)
  )
  result <- merge(result, lineage_by_sample, by = "sample_id", sort = FALSE)
  result[, nkx2_1_high_tru_concordant :=
    subtype == "TRU" & nkx2_1_above_median %in% TRUE & alveolar_lineage_z_mean > 0]
  setcolorder(result, c(
    "sample_id", "cohort", "input_scale", "available_predictor_genes",
    "coverage_fraction", "correlation_bronchioid", "correlation_magnoid",
    "correlation_squamoid", "top_centroid", "second_centroid",
    "correlation_margin", "subtype", "lineage_genes_available",
    "nkx2_1_log2_expression", "nkx2_1_cohort_median", "nkx2_1_above_median",
    "alveolar_lineage_z_mean", "nkx2_1_high_tru_concordant"
  ))
  list(result = result, common = common, centered = centered)
}

entropy_bits <- function(labels) {
  labels <- labels[!is.na(labels)]
  if (!length(labels)) return(NA_real_)
  probabilities <- as.numeric(table(labels)) / length(labels)
  -sum(probabilities * log2(probabilities))
}

wilson <- function(successes, total, z = 1.95996398454005) {
  if (!total) return(c(lower = NA_real_, upper = NA_real_))
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}

summarize_labels <- function(result, cohort) {
  levels <- c("TRU", "PP", "PI", "UNCLASSIFIABLE")
  counts <- table(factor(result$subtype, levels = levels))
  rows <- lapply(levels, function(label) {
    interval <- wilson(unname(counts[[label]]), nrow(result))
    data.table(
      cohort = cohort, subtype = label, n = unname(counts[[label]]),
      denominator = nrow(result), fraction = unname(counts[[label]]) / nrow(result),
      wilson_95_lower = interval[["lower"]], wilson_95_upper = interval[["upper"]],
      classification_entropy_bits = entropy_bits(result$subtype)
    )
  })
  rbindlist(rows)
}

# TCGA expression written by Figure 1 is on a scale of 2^Xena_log2 - 0.001.
# The official downloaded STAR TPM matrix has zero encoded as 0 and is
# log2(TPM+1), so log2(processed + 0.001) exactly restores that log scale even
# after duplicate-symbol averaging. This also avoids propagating the historic
# +0.999 TPM offset in the Figure 1 preparation comment.
tcga_expression_path <- file.path(figure1_root, "data", "processed", "tcga_aracne_expression.tsv")
tcga_manifest_path <- file.path(figure1_root, "data", "processed", "tcga_manifest.tsv")
tcga_atac_manifest_path <- file.path(
  figure2_root, "audit", "manifests", "tcga", "tcga_luad_atac_samples.tsv"
)
for (path in c(tcga_expression_path, tcga_manifest_path, tcga_atac_manifest_path)) {
  if (!file.exists(path)) stop("Missing TCGA input: ", path)
}
tcga_tab <- fread(tcga_expression_path, check.names = FALSE)
setnames(tcga_tab, 1L, "gene")
tcga_manifest <- fread(tcga_manifest_path)
luad_samples <- tcga_manifest[group == "LUAD", sample_id]
if (!length(luad_samples) || !all(luad_samples %in% names(tcga_tab))) {
  stop("TCGA LUAD manifest and expression columns do not match")
}
tcga_values <- as.matrix(tcga_tab[, ..luad_samples])
storage.mode(tcga_values) <- "double"
rownames(tcga_values) <- tcga_tab$gene
tcga_log2 <- log2(tcga_values + 0.001)
tcga_classification <- classify_expression(
  tcga_log2, "TCGA_LUAD", "restored_Xena_STAR_log2_TPM_plus_1"
)$result
fwrite(tcga_classification, file.path(out_dir, "tcga_luad_wilkerson506_classification.tsv"), sep = "\t")
fwrite(
  summarize_labels(tcga_classification, "TCGA_LUAD"),
  file.path(out_dir, "tcga_luad_wilkerson506_summary.tsv"), sep = "\t"
)

atac_manifest <- fread(tcga_atac_manifest_path)
atac_manifest[, patient_id := substr(sample_id, 1L, 12L)]
tcga_classification[, patient_id := substr(sample_id, 1L, 12L)]
matched <- merge(
  atac_manifest[, .(patient_id, atac_sample_id = sample_id, has_rna_seq_published)],
  tcga_classification,
  by = "patient_id", all.x = TRUE, sort = FALSE
)
fwrite(matched, file.path(out_dir, "tcga_luad_atac_rna_matched_state_labels.tsv"), sep = "\t")

cell_result <- NULL
cell_quant_root <- file.path(state_root, "data", "processed", "dra001846_quant")
cell_manifest_path <- file.path(
  state_root, "audit", "manifests", "dra_rna", "ddbj_dra001846_baseline_libraries.tsv"
)
if (!tcga_only) {
  if (!file.exists(cell_manifest_path)) stop("Missing frozen DRA001846 manifest")
  cell_manifest <- fread(cell_manifest_path)
  if (nrow(cell_manifest) != 19L || uniqueN(cell_manifest$cell_line) != 19L) {
    stop("Frozen DRA001846 manifest must contain 19 unique cell lines")
  }
  quant_tables <- lapply(seq_len(nrow(cell_manifest)), function(i) {
    model <- cell_manifest$cell_line[[i]]
    slug <- gsub("/", "_", model, fixed = TRUE)
    path <- file.path(cell_quant_root, slug, "salmon", "quant.sf")
    if (!file.exists(path)) stop("Missing completed Salmon quantification for ", model, ": ", path)
    tab <- fread(path, select = c("Name", "TPM"))
    fields <- tstrsplit(tab$Name, "|", fixed = TRUE, keep = 6L)
    tab[, gene := fields[[1]]]
    if (anyNA(tab$gene) || any(tab$gene == "")) stop("Cannot parse GENCODE gene symbol for ", model)
    aggregated <- tab[, .(TPM = sum(TPM)), by = gene]
    setnames(aggregated, "TPM", model)
    aggregated
  })
  cell_tpm <- Reduce(function(left, right) merge(left, right, by = "gene", all = TRUE), quant_tables)
  for (column in setdiff(names(cell_tpm), "gene")) set(cell_tpm, which(is.na(cell_tpm[[column]])), column, 0)
  setorder(cell_tpm, gene)
  fwrite(cell_tpm, file.path(out_dir, "dra001846_19_cellline_gene_tpm.tsv.gz"), sep = "\t")
  cell_columns <- setdiff(names(cell_tpm), "gene")
  cell_values <- as.matrix(cell_tpm[, ..cell_columns])
  storage.mode(cell_values) <- "double"
  rownames(cell_values) <- cell_tpm$gene
  cell_result <- classify_expression(
    log2(cell_values + 1), "DRA001846_19_LUAD_cell_lines", "Salmon_gene_TPM_log2_plus_1"
  )$result
  fwrite(cell_result, file.path(out_dir, "dra001846_19_cellline_wilkerson506_classification.tsv"), sep = "\t")
  fwrite(
    summarize_labels(cell_result, "DRA001846_19_LUAD_cell_lines"),
    file.path(out_dir, "dra001846_19_cellline_wilkerson506_summary.tsv"), sep = "\t"
  )
}

receipt <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  tcga_only = tcga_only,
  classifier = list(
    name = "Wilkerson_2012_LAD_nearest_centroid",
    predictor_genes = nrow(centroids),
    minimum_available_genes = 405L,
    centroids_sha256 = sha256(centroid_path),
    subtype_mapping = as.list(subtype_map),
    transformation = "within-cohort gene median centering; Pearson sample-to-centroid correlation; maximum correlation label",
    confidence_exclusion = FALSE
  ),
  tcga = list(
    samples = nrow(tcga_classification),
    expression_sha256 = sha256(tcga_expression_path),
    manifest_sha256 = sha256(tcga_manifest_path),
    atac_manifest_sha256 = sha256(tcga_atac_manifest_path),
    matched_atac_with_state_label = sum(!is.na(matched$subtype)),
    restored_scale_formula = "log2(figure1_processed_value + 0.001) = official Xena STAR log2(TPM+1)"
  ),
  cell_lines = if (is.null(cell_result)) NULL else list(
    samples = nrow(cell_result),
    manifest_sha256 = sha256(cell_manifest_path),
    quantification = "Salmon transcript TPM summed by GENCODE v36 gene symbol; log2(TPM+1)"
  )
)
writeLines(
  toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, null = "null"),
  file.path(audit_dir, if (tcga_only) "tcga_state_classifier_receipt.json" else "state_classifier_receipt.json")
)
cat(toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")

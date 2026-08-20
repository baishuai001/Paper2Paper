#!/usr/bin/env Rscript

# Transfer the validated TCGA promoter-ATAC proxy to the frozen 13 LUAD PDX.
# Requires count_pdx_atac_proxy_intervals.py to have completed without errors.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: classify_pdx_atac_states.R <state_audit_root> <figure2_root>")
}
state_root <- normalizePath(args[[1]], mustWork = TRUE)
figure2_root <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- file.path(state_root, "results", "state_classifier")
audit_dir <- file.path(state_root, "audit", "state_classifier")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  answer <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  strsplit(answer[[1]], "[[:space:]]+")[[1]][[1]]
}

centroid_path <- file.path(
  state_root, "reference", "lungCancerSubtypes", "lung_adenocarcinoma_subtypes",
  "wilkerson.2012.LAD.predictor.centroids.csv"
)
validation_path <- file.path(
  state_root, "audit", "atac_state_proxy", "tcga_atac_state_proxy_validation_receipt.json"
)
count_receipt_path <- file.path(
  state_root, "audit", "pdx_atac_proxy_counts", "pdx_atac_proxy_counting_receipt.json"
)
bed_path <- file.path(
  state_root, "reference", "atac_proxy", "tcga_luad_139135_consensus_intervals.bed"
)
tss_path <- file.path(
  figure2_root, "reference", "figure2_features", "gencode.v47.protein_coding.gene_tss.bed"
)
qc_path <- file.path(figure2_root, "audit", "raw_atac_qc", "figure2_raw_atac_qc.tsv")
for (path in c(centroid_path, validation_path, count_receipt_path, bed_path, tss_path, qc_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}
validation <- fromJSON(validation_path, simplifyVector = FALSE)
if (!isTRUE(validation$tru_vs_rest$pdx_transfer_authorized)) {
  stop("Frozen TCGA validation did not authorize PDX transfer")
}
count_receipt <- fromJSON(count_receipt_path, simplifyVector = FALSE)
if (length(count_receipt$errors) || length(count_receipt$completed_models) != 13L) {
  stop("PDX fixed-interval counting receipt is incomplete")
}

centroids <- fread(centroid_path, check.names = FALSE)
setnames(centroids, 1L, "gene")
centroid_names <- c("bronchioid", "magnoid", "squamoid")
centroid_matrix <- as.matrix(centroids[, ..centroid_names])
storage.mode(centroid_matrix) <- "double"
rownames(centroid_matrix) <- centroids$gene
subtype_map <- c(bronchioid = "TRU", magnoid = "PP", squamoid = "PI")

intervals <- fread(
  bed_path, header = FALSE,
  col.names = c("chromosome", "start", "end", "peak_id")
)
if (nrow(intervals) != 139135L || anyDuplicated(intervals$peak_id)) {
  stop("Frozen coordinate BED is not 139,135 unique intervals")
}
tss <- fread(
  tss_path, header = FALSE,
  col.names = c("chromosome", "tss_start", "tss_end", "gene", "score", "strand")
)
tss <- tss[gene %in% centroids$gene]
tss[, promoter_start := pmax(0L, as.integer(tss_start) - 2000L)]
tss[, promoter_end := as.integer(tss_end) + 2000L]
promoters <- unique(tss[, .(gene, chromosome, start = promoter_start, end = promoter_end)])
setkey(promoters, chromosome, start, end)
interval_ranges <- copy(intervals)
setkey(interval_ranges, chromosome, start, end)
overlaps <- foverlaps(
  interval_ranges, promoters,
  by.x = c("chromosome", "start", "end"),
  by.y = c("chromosome", "start", "end"),
  type = "any", nomatch = NULL
)
gene_peak <- unique(overlaps[, .(gene, peak_id)])
mapped_genes <- intersect(centroids$gene, unique(promoters$gene))
if (length(mapped_genes) < 405L) stop("Fewer than 405 predictor genes map to TSS")

count_dir <- file.path(state_root, "data", "processed", "pdx_atac_proxy_counts")
count_files <- sort(list.files(count_dir, pattern = "\\.counts\\.tsv\\.gz$", full.names = TRUE))
if (length(count_files) != 13L) stop("Expected 13 PDX count tables, found ", length(count_files))
models <- sub("\\.tcga_luad_intervals\\.counts\\.tsv\\.gz$", "", basename(count_files))
count_matrix <- matrix(
  0, nrow = nrow(intervals), ncol = length(models),
  dimnames = list(intervals$peak_id, models)
)
for (i in seq_along(count_files)) {
  tab <- fread(count_files[[i]])
  if (nrow(tab) != nrow(intervals) || !identical(tab$peak_id, intervals$peak_id)) {
    stop("Coordinate/order mismatch in ", count_files[[i]])
  }
  count_matrix[, i] <- tab$count
}
library_sizes <- colSums(count_matrix)
if (any(!is.finite(library_sizes) | library_sizes <= 0)) stop("Invalid fixed-interval library size")
cpm <- sweep(count_matrix, 2L, library_sizes / 1e6, "/")

gene_cpm <- matrix(
  0, nrow = length(mapped_genes), ncol = ncol(cpm),
  dimnames = list(mapped_genes, colnames(cpm))
)
peak_row <- setNames(seq_len(nrow(count_matrix)), rownames(count_matrix))
for (gene_name in mapped_genes) {
  peak_ids <- gene_peak[gene == gene_name, peak_id]
  indices <- unname(peak_row[peak_ids])
  indices <- indices[!is.na(indices)]
  if (length(indices)) gene_cpm[gene_name, ] <- colSums(cpm[indices, , drop = FALSE])
}
gene_log2 <- log2(gene_cpm + 1)
gene_medians <- apply(gene_log2, 1L, median, na.rm = TRUE)
centered <- sweep(gene_log2, 1L, gene_medians, "-")

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  cor(x[keep], y[keep], method = "pearson")
}

rows <- lapply(seq_len(ncol(centered)), function(i) {
  scores <- vapply(
    centroid_names,
    function(name) safe_cor(centered[, i], centroid_matrix[mapped_genes, name]),
    numeric(1)
  )
  ordered <- names(sort(scores, decreasing = TRUE, na.last = TRUE))
  data.table(
    sample_id = colnames(centered)[[i]],
    cohort = "GSE269746_13_LUAD_PDX_ATAC_proxy",
    available_predictor_genes = length(mapped_genes),
    predictor_genes_with_nonzero_accessibility = sum(gene_cpm[, i] > 0),
    fixed_interval_library_size = library_sizes[[i]],
    correlation_bronchioid = scores[["bronchioid"]],
    correlation_magnoid = scores[["magnoid"]],
    correlation_squamoid = scores[["squamoid"]],
    top_centroid = ordered[[1]],
    second_centroid = ordered[[2]],
    correlation_margin = scores[[ordered[[1]]]] - scores[[ordered[[2]]]],
    subtype = unname(subtype_map[[ordered[[1]]]])
  )
})
result <- rbindlist(rows)

qc <- fread(qc_path)[system == "PDX"]
qc_columns <- c(
  "sample_id", "filtered_human_reads", "filtered_human_units", "peak_count",
  "canonical_peak_count", "frip", "peak_call_status", "worker_diagnostic_qc_status"
)
result <- merge(result, qc[, ..qc_columns], by = "sample_id", all.x = TRUE, sort = FALSE)
if (anyNA(result$filtered_human_reads)) stop("PDX state labels did not fully match raw ATAC QC")

entropy_bits <- function(labels) {
  probabilities <- as.numeric(table(labels)) / length(labels)
  -sum(probabilities * log2(probabilities))
}
wilson <- function(successes, total, z = 1.95996398454005) {
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}
summary_rows <- lapply(c("TRU", "PP", "PI", "UNCLASSIFIABLE"), function(label) {
  n <- sum(result$subtype == label)
  interval <- wilson(n, nrow(result))
  data.table(
    cohort = "GSE269746_13_LUAD_PDX_ATAC_proxy", subtype = label,
    n = n, denominator = nrow(result), fraction = n / nrow(result),
    wilson_95_lower = interval[["lower"]], wilson_95_upper = interval[["upper"]],
    classification_entropy_bits = entropy_bits(result$subtype)
  )
})
summary <- rbindlist(summary_rows)

fwrite(result, file.path(out_dir, "gse269746_13_pdx_wilkerson506_atac_proxy_classification.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "gse269746_13_pdx_wilkerson506_atac_proxy_summary.tsv"), sep = "\t")
fwrite(
  data.table(gene = rownames(gene_log2), gene_log2),
  file.path(out_dir, "gse269746_13_pdx_promoter_atac_wilkerson506_log2_cpm.tsv.gz"), sep = "\t"
)

receipt <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  transfer_authorization = list(
    validation_receipt_sha256 = sha256(validation_path),
    balanced_accuracy = validation$tru_vs_rest$balanced_accuracy,
    one_sided_permutation_p = validation$tru_vs_rest$one_sided_permutation_p,
    passed = validation$tru_vs_rest$passed
  ),
  pdx_models = nrow(result),
  predictor_genes_with_gencode_v47_tss = length(mapped_genes),
  fixed_coordinate_intervals = nrow(intervals),
  interval_bed_sha256 = sha256(bed_path),
  cpm_library_size_definition = "sum of alignment counts over the same 139,135 TCGA-LUAD intervals",
  gene_feature_definition = "sum interval CPM within union of all TSS +/-2 kb; log2(sum(CPM)+1); within-PDX-cohort gene median centering",
  labels = as.list(table(result$subtype)),
  all_same_state = uniqueN(result$subtype) == 1L,
  classification_entropy_bits = entropy_bits(result$subtype),
  no_qc_based_exclusions = TRUE,
  caution = "Labels are validated ATAC-proxy assignments, not matched-RNA assignments; raw ATAC QC is retained per model."
)
writeLines(
  toJSON(receipt, auto_unbox = TRUE, pretty = TRUE),
  file.path(audit_dir, "pdx_atac_state_classifier_receipt.json")
)
cat(toJSON(receipt, auto_unbox = TRUE, pretty = TRUE), "\n")

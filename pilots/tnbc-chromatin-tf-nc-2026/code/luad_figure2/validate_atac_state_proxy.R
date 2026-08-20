#!/usr/bin/env Rscript

# Validate the prespecified promoter-ATAC nearest-centroid proxy on TCGA-LUAD
# tumors with independently assigned RNA subtypes. The proxy is not transferred
# to PDX unless both frozen pass criteria are met.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: validate_atac_state_proxy.R <state_audit_root> <figure2_root>")
}
state_root <- normalizePath(args[[1]], mustWork = TRUE)
figure2_root <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- file.path(state_root, "results", "atac_state_proxy")
audit_dir <- file.path(state_root, "audit", "atac_state_proxy")
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
rna_label_path <- file.path(
  state_root, "results", "state_classifier", "tcga_luad_atac_rna_matched_state_labels.tsv"
)
atac_rds_path <- file.path(
  figure2_root, "data", "processed", "tcga_luad", "tcga_luad_accessibility_matrices.rds"
)
tss_path <- file.path(
  figure2_root, "reference", "figure2_features", "gencode.v47.protein_coding.gene_tss.bed"
)
for (path in c(centroid_path, rna_label_path, atac_rds_path, tss_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

centroids <- fread(centroid_path, check.names = FALSE)
setnames(centroids, 1L, "gene")
centroid_names <- c("bronchioid", "magnoid", "squamoid")
if (nrow(centroids) != 506L || !all(centroid_names %in% names(centroids))) {
  stop("Unexpected Wilkerson centroid file")
}
centroid_matrix <- as.matrix(centroids[, ..centroid_names])
storage.mode(centroid_matrix) <- "double"
rownames(centroid_matrix) <- centroids$gene
subtype_map <- c(bronchioid = "TRU", magnoid = "PP", squamoid = "PI")

atac <- readRDS(atac_rds_path)
required_objects <- c("coordinates", "log2_merged_cpm")
if (!all(required_objects %in% names(atac))) stop("Unexpected TCGA ATAC RDS schema")
coordinates <- as.data.table(atac$coordinates)
log2_cpm <- atac$log2_merged_cpm
if (nrow(coordinates) != 139135L || nrow(log2_cpm) != nrow(coordinates) || ncol(log2_cpm) != 22L) {
  stop("Frozen TCGA coordinate universe must be 139,135 intervals x 22 tumors")
}
cpm <- pmax(2^log2_cpm - 1, 0)

tss <- fread(
  tss_path, header = FALSE,
  col.names = c("chromosome", "tss_start", "tss_end", "gene", "score", "strand")
)
tss <- tss[gene %in% centroids$gene]
tss[, promoter_start := pmax(0L, as.integer(tss_start) - 2000L)]
tss[, promoter_end := as.integer(tss_end) + 2000L]
promoters <- unique(tss[, .(
  gene, chromosome, start = promoter_start, end = promoter_end
)])

peaks <- data.table(
  peak_index = seq_len(nrow(coordinates)),
  chromosome = as.character(coordinates$seqnames),
  start = as.integer(coordinates$start),
  end = as.integer(coordinates$end)
)
setkey(promoters, chromosome, start, end)
setkey(peaks, chromosome, start, end)
overlaps <- foverlaps(
  peaks, promoters,
  by.x = c("chromosome", "start", "end"),
  by.y = c("chromosome", "start", "end"),
  type = "any", nomatch = NULL
)
gene_peak <- unique(overlaps[, .(gene, peak_index)])
mapped_genes <- intersect(centroids$gene, unique(promoters$gene))
if (length(mapped_genes) < 405L) {
  stop("Only ", length(mapped_genes), "/506 predictor genes map to frozen GENCODE TSS")
}

gene_cpm <- matrix(
  0, nrow = length(mapped_genes), ncol = ncol(cpm),
  dimnames = list(mapped_genes, colnames(cpm))
)
for (gene_name in mapped_genes) {
  indices <- gene_peak[gene == gene_name, peak_index]
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

prediction_rows <- lapply(seq_len(ncol(centered)), function(i) {
  scores <- vapply(
    centroid_names,
    function(name) safe_cor(centered[, i], centroid_matrix[mapped_genes, name]),
    numeric(1)
  )
  ordered <- names(sort(scores, decreasing = TRUE, na.last = TRUE))
  data.table(
    atac_sample_id = colnames(centered)[[i]],
    predictor_genes_with_tss = length(mapped_genes),
    predictor_genes_with_overlapping_interval = sum(gene_cpm[, i] > 0),
    correlation_bronchioid = scores[["bronchioid"]],
    correlation_magnoid = scores[["magnoid"]],
    correlation_squamoid = scores[["squamoid"]],
    atac_top_centroid = ordered[[1]],
    atac_second_centroid = ordered[[2]],
    atac_correlation_margin = scores[[ordered[[1]]]] - scores[[ordered[[2]]]],
    atac_proxy_subtype = unname(subtype_map[[ordered[[1]]]])
  )
})
predictions <- rbindlist(prediction_rows)
rna_labels <- fread(rna_label_path)
validation <- merge(
  predictions,
  rna_labels[!is.na(subtype) & nzchar(subtype), .(
    patient_id, atac_sample_id, rna_subtype = subtype,
    rna_correlation_margin = correlation_margin
  )],
  by = "atac_sample_id", all.x = TRUE, sort = FALSE
)
evaluation <- validation[!is.na(rna_subtype) & nzchar(rna_subtype)]
if (nrow(evaluation) != 21L) stop("Frozen validation set must contain 21 RNA-ATAC matched tumors")

truth <- evaluation$rna_subtype == "TRU"
prediction <- evaluation$atac_proxy_subtype == "TRU"
if (!any(truth) || !any(!truth)) stop("TRU-vs-rest truth labels are not evaluable")
sensitivity <- mean(prediction[truth])
specificity <- mean(!prediction[!truth])
balanced_accuracy <- (sensitivity + specificity) / 2
set.seed(5062012L)
n_permutations <- 10000L
permuted_ba <- replicate(n_permutations, {
  permuted_truth <- sample(truth, replace = FALSE)
  permuted_sensitivity <- mean(prediction[permuted_truth])
  permuted_specificity <- mean(!prediction[!permuted_truth])
  (permuted_sensitivity + permuted_specificity) / 2
})
permutation_p <- (1 + sum(permuted_ba >= balanced_accuracy)) / (n_permutations + 1)
passed <- balanced_accuracy >= 0.70 && permutation_p < 0.05

fwrite(validation, file.path(out_dir, "tcga_atac_proxy_predictions.tsv"), sep = "\t")
fwrite(
  data.table(
    gene = mapped_genes,
    tss_records = vapply(mapped_genes, function(g) nrow(promoters[gene == g]), integer(1)),
    distinct_consensus_intervals = vapply(mapped_genes, function(g) uniqueN(gene_peak[gene == g, peak_index]), integer(1)),
    variable_across_tcga = apply(gene_log2, 1L, function(x) sd(x) > 0)
  ),
  file.path(out_dir, "atac_proxy_gene_mapping.tsv"), sep = "\t"
)
fwrite(
  data.table(gene = rownames(gene_log2), gene_log2),
  file.path(out_dir, "tcga_promoter_atac_wilkerson506_log2_cpm.tsv.gz"), sep = "\t"
)

receipt <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  frozen_coordinate_intervals = nrow(coordinates),
  tcga_atac_tumors = ncol(log2_cpm),
  matched_validation_tumors = nrow(evaluation),
  rna_truth_counts = as.list(table(evaluation$rna_subtype)),
  predictor_genes_total = 506L,
  predictor_genes_with_gencode_v47_tss = length(mapped_genes),
  predictor_genes_with_any_consensus_interval = sum(rowSums(gene_cpm) > 0),
  promoter_definition = "union of all GENCODE v47 protein-coding TSS +/-2 kb per gene",
  feature_definition = "sum merged-replicate CPM over distinct overlapping TCGA-LUAD consensus intervals; log2(sum(CPM)+1); within-cohort gene median centering",
  classifier = "Pearson correlation to frozen Wilkerson 2012 centroids; maximum correlation label",
  tru_vs_rest = list(
    true_tru = sum(truth),
    true_rest = sum(!truth),
    predicted_tru = sum(prediction),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_accuracy,
    permutations = n_permutations,
    permutation_seed = 5062012L,
    one_sided_permutation_p = permutation_p,
    pass_rule = "balanced_accuracy >= 0.70 AND one-sided permutation p < 0.05",
    passed = passed,
    pdx_transfer_authorized = passed
  ),
  inputs = list(
    centroid_sha256 = sha256(centroid_path),
    rna_labels_sha256 = sha256(rna_label_path),
    tcga_atac_rds_sha256 = sha256(atac_rds_path),
    tss_bed_sha256 = sha256(tss_path)
  )
)
writeLines(
  toJSON(receipt, auto_unbox = TRUE, pretty = TRUE),
  file.path(audit_dir, "tcga_atac_state_proxy_validation_receipt.json")
)
cat(toJSON(receipt, auto_unbox = TRUE, pretty = TRUE), "\n")

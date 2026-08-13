#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(viper)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: viper_msviper.R BULK_TPM_RDS CONSOLIDATED_NETWORK",
    "ATLAS_LOG2CPM PATIENT_METADATA OUTPUT_DIR THREADS PERMUTATIONS"
  ))
}
bulk_path <- args[[1]]
network_path <- args[[2]]
atlas_path <- args[[3]]
metadata_path <- args[[4]]
output_dir <- normalizePath(args[[5]], mustWork = FALSE)
threads <- as.integer(args[[6]])
permutations <- as.integer(args[[7]])
if (!is.finite(permutations) || permutations < 50) {
  stop("PERMUTATIONS must be >=50: viper::aecdf requires a sufficiently diverse empirical null")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_expression <- function(path) {
  tab <- fread(path, data.table = FALSE, check.names = FALSE)
  genes <- as.character(tab[[1]])
  values <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- genes
  if (anyDuplicated(rownames(values)) || anyDuplicated(colnames(values))) {
    stop(sprintf("Duplicate rows or columns in %s", path))
  }
  values
}

network <- fread(network_path, data.table = FALSE)
required_network <- c("regulator.values", "target.values", "mi.values")
if (!all(required_network %in% colnames(network))) {
  stop(sprintf("ARACNe3 author-style subnetwork missing columns: %s", paste(setdiff(required_network, colnames(network)), collapse = ", ")))
}
if (nrow(network) == 0 || anyDuplicated(network[, required_network])) {
  stop("ARACNe3 author-style subnetwork is empty or has duplicate edge records")
}
fwrite(network, file.path(output_dir, "aracne3_author_subnetwork_edges.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

bulk_tpm <- readRDS(bulk_path)
if (anyDuplicated(rownames(bulk_tpm)) || anyDuplicated(colnames(bulk_tpm))) {
  stop("Bulk TPM matrix contains duplicate genes or participants")
}
network_3col <- tempfile(pattern = "aracne3_author_subnetwork_", fileext = ".tsv")
fwrite(
  network[, required_network],
  network_3col, sep = "\t", quote = FALSE, col.names = FALSE
)
regulon <- aracne2regulon(
  network_3col,
  bulk_tpm,
  format = "3col",
  verbose = TRUE
)
unlink(network_3col)
saveRDS(regulon, file.path(output_dir, "crc_aracne3_regulon.rds"), compress = "xz")

regulon_edges <- rbindlist(lapply(names(regulon), function(tf) {
  data.table(
    TF = tf,
    target = names(regulon[[tf]]$tfmode),
    tfmode = as.numeric(regulon[[tf]]$tfmode),
    likelihood = as.numeric(regulon[[tf]]$likelihood)
  )
}))
fwrite(regulon_edges, file.path(output_dir, "crc_aracne3_regulon_edges.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

atlas <- read_expression(atlas_path)
metadata <- fread(metadata_path, data.table = FALSE)
if (!all(c("donor_id", "dataset", "group", "analysis_cells_aggregated") %in% colnames(metadata))) {
  stop("Patient metadata lacks required columns")
}
position <- match(colnames(atlas), metadata$donor_id)
if (anyNA(position) || anyDuplicated(position)) {
  stop("Atlas expression columns do not map one-to-one to metadata")
}
metadata <- metadata[position, , drop = FALSE]
if (!identical(as.character(metadata$donor_id), colnames(atlas))) {
  stop("Atlas expression/metadata order mismatch")
}

measured_targets <- vapply(regulon, function(x) sum(names(x$tfmode) %in% rownames(atlas)), numeric(1))
eligible_regulons <- names(measured_targets)[measured_targets >= 1]
if (length(eligible_regulons) < 400) {
  stop(sprintf("Only %d regulons have >=1 measured atlas target", length(eligible_regulons)))
}

activity <- viper(
  atlas,
  regulon,
  minsize = 1,
  nes = TRUE,
  eset.filter = TRUE,
  cores = threads,
  verbose = TRUE
)
if (is.null(dim(activity))) {
  stop("VIPER returned a non-matrix activity object")
}
activity_tab <- data.table(TF = rownames(activity), activity, keep.rownames = FALSE)
fwrite(activity_tab, file.path(output_dir, "viper_activity.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

primary <- metadata$analysis_cells_aggregated >= 50
primary_metadata <- metadata[primary, , drop = FALSE]
primary_expression <- atlas[, primary, drop = FALSE]
groups <- primary_metadata$group == "case"
datasets <- as.character(primary_metadata$dataset)
dataset_names <- sort(unique(datasets))
informative <- dataset_names[vapply(dataset_names, function(dataset) {
  index <- datasets == dataset
  sum(groups[index]) >= 3 && sum(!groups[index]) >= 3
}, logical(1))]
if (length(informative) < 3) {
  stop(sprintf("Only %d informative independent studies at the frozen 50-cell threshold", length(informative)))
}

welch_t <- function(matrix, label) {
  case <- matrix[, label, drop = FALSE]
  control <- matrix[, !label, drop = FALSE]
  if (ncol(case) < 2 || ncol(control) < 2) {
    stop("Welch t statistic requires at least two patients per group")
  }
  mean_difference <- rowMeans(case) - rowMeans(control)
  row_variance <- function(x) {
    n <- ncol(x)
    means <- rowMeans(x)
    pmax((rowSums(x * x) - n * means * means) / (n - 1), 0)
  }
  denominator <- sqrt(
    row_variance(case) / ncol(case) +
      row_variance(control) / ncol(control)
  )
  statistic <- mean_difference / denominator
  statistic[!is.finite(statistic)] <- 0
  statistic
}

meta_signature <- function(label) {
  numer <- rep(0, nrow(primary_expression))
  denominator <- 0
  for (dataset in informative) {
    index <- datasets == dataset
    local <- label[index]
    weight <- sqrt(sum(local) * sum(!local) / length(local))
    numer <- numer + weight * welch_t(primary_expression[, index, drop = FALSE], local)
    denominator <- denominator + weight^2
  }
  result <- numer / sqrt(denominator)
  names(result) <- rownames(primary_expression)
  result
}

observed_signature <- meta_signature(groups)
set.seed(1729)
null_signatures <- matrix(
  0,
  nrow = length(observed_signature),
  ncol = permutations,
  dimnames = list(names(observed_signature), sprintf("perm_%04d", seq_len(permutations)))
)
for (permutation in seq_len(permutations)) {
  permuted <- groups
  for (dataset in informative) {
    index <- which(datasets == dataset)
    permuted[index] <- sample(permuted[index], length(index), replace = FALSE)
  }
  null_signatures[, permutation] <- meta_signature(permuted)
}
saveRDS(
  list(observed = observed_signature, null = null_signatures, informative_datasets = informative),
  file.path(output_dir, "msviper_gene_signatures.rds"),
  compress = "xz"
)

ms <- msviper(
  observed_signature,
  regulon,
  nullmodel = null_signatures,
  minsize = 1,
  adaptive.size = FALSE,
  ges.filter = TRUE,
  cores = threads,
  verbose = TRUE
)
ms_table <- data.frame(
  TF = names(ms$es$nes),
  NES = as.numeric(ms$es$nes),
  size = as.numeric(ms$es$size),
  p.value = as.numeric(ms$es$p.value),
  stringsAsFactors = FALSE
)
ms_table$FDR <- p.adjust(ms_table$p.value, method = "BH")
ms_table <- ms_table[order(ms_table$FDR, -abs(ms_table$NES)), , drop = FALSE]
fwrite(ms_table, file.path(output_dir, "msviper_results.tsv"), sep = "\t", quote = FALSE)

receipt <- list(
  status = "passed",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  network = list(
    author_subnetwork_edges = nrow(network),
    author_subnetwork_regulators = length(unique(network$regulator.values)),
    ARACNe3_subnetwork_FDR_alpha = 0.05,
    second_consensus_BH = FALSE,
    regulons_after_conversion = length(regulon),
    regulons_ge1_measured_target = length(eligible_regulons)
  ),
  viper = list(
    patients_scored = ncol(activity),
    TFs_scored = nrow(activity),
    method = "none (viper default; anchor code omits method)",
    minsize = 1,
    threads = threads
  ),
  msviper = list(
    primary_patients = sum(primary),
    case_patients = sum(groups),
    control_patients = sum(!groups),
    informative_datasets = informative,
    validation_unit = "study_id (stored in the generic dataset field)",
    permutations = permutations,
    seed = 1729,
    TFs_tested = nrow(ms_table),
    TFs_FDR_0_01 = sum(ms_table$FDR <= 0.01)
  ),
  package_versions = list(
    R = as.character(getRversion()),
    viper = as.character(packageVersion("viper")),
    data.table = as.character(packageVersion("data.table"))
  )
)
write_json(receipt, file.path(output_dir, "viper_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

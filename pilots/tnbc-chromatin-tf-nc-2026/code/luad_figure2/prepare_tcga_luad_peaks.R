#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: prepare_tcga_luad_peaks.R <figure2_run_root>")
run_root <- normalizePath(args[[1]], mustWork = TRUE)
raw_dir <- file.path(run_root, "data", "raw", "tcga_atac_gdc")
manifest_path <- file.path(run_root, "audit", "manifests", "tcga", "tcga_luad_atac_samples.tsv")
out_dir <- file.path(run_root, "data", "processed", "tcga_luad")
peak_dir <- file.path(out_dir, "sample_peaks", "cpm1_both_technical_replicates")
peak_dir_05 <- file.path(out_dir, "sample_peaks", "cpm0.5_both_technical_replicates")
peak_dir_2 <- file.path(out_dir, "sample_peaks", "cpm2_both_technical_replicates")
peak_dir_merged <- file.path(out_dir, "sample_peaks", "cpm1_merged_technical_replicates")
audit_dir <- file.path(run_root, "audit", "tcga_luad_accessibility")
for (path in c(peak_dir, peak_dir_05, peak_dir_2, peak_dir_merged)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(raw_dir, "TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip")
if (!file.exists(zip_path) || !file.exists(manifest_path)) stop("Missing count ZIP or frozen manifest")
manifest <- fread(manifest_path)
if (nrow(manifest) != 22L || !all(manifest$figure2_primary_include)) {
  stop("TCGA manifest must contain 22 included independent tumors")
}

message("Streaming LUAD_raw_counts.txt from the official GDC ZIP...")
counts <- fread(
  cmd = sprintf("unzip -p %s LUAD_raw_counts.txt", shQuote(zip_path)),
  sep = "\t",
  header = TRUE,
  showProgress = TRUE
)
metadata_columns <- c("seqnames", "start", "end", "name", "score")
if (!all(metadata_columns %in% names(counts))) stop("Unexpected GDC LUAD matrix schema")
library_columns <- setdiff(names(counts), metadata_columns)
raw_matrix <- as.matrix(counts[, ..library_columns])
storage.mode(raw_matrix) <- "double"
library_sizes <- colSums(raw_matrix)
if (any(!is.finite(library_sizes) | library_sizes <= 0)) stop("Invalid ATAC library size")
cpm <- sweep(raw_matrix, 2L, library_sizes / 1e6, "/")

binary_primary <- matrix(FALSE, nrow = nrow(counts), ncol = nrow(manifest))
binary_05 <- matrix(FALSE, nrow = nrow(counts), ncol = nrow(manifest))
binary_2 <- matrix(FALSE, nrow = nrow(counts), ncol = nrow(manifest))
binary_merged <- matrix(FALSE, nrow = nrow(counts), ncol = nrow(manifest))
log2_merged_cpm <- matrix(NA_real_, nrow = nrow(counts), ncol = nrow(manifest))
colnames(binary_primary) <- manifest$sample_id
colnames(binary_05) <- manifest$sample_id
colnames(binary_2) <- manifest$sample_id
colnames(binary_merged) <- manifest$sample_id
colnames(log2_merged_cpm) <- manifest$sample_id
peak_count_rows <- vector("list", nrow(manifest))
peak_manifest_rows <- list()

write_peak_bed <- function(mask, directory, sample_id) {
  bed <- counts[mask, .(seqnames, start, end, name)]
  bed_path <- file.path(directory, paste0(sample_id, ".bed.gz"))
  con <- gzfile(bed_path, open = "wt")
  write.table(bed, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)
  normalizePath(bed_path, mustWork = TRUE)
}

for (i in seq_len(nrow(manifest))) {
  libs <- strsplit(manifest$count_columns[[i]], ";", fixed = TRUE)[[1]]
  if (length(libs) != 2L || !all(libs %in% library_columns)) {
    stop("Expected two count columns for ", manifest$sample_id[[i]])
  }
  ix <- match(libs, library_columns)
  both_05 <- cpm[, ix[[1]]] >= 0.5 & cpm[, ix[[2]]] >= 0.5
  both_1 <- cpm[, ix[[1]]] >= 1.0 & cpm[, ix[[2]]] >= 1.0
  both_2 <- cpm[, ix[[1]]] >= 2.0 & cpm[, ix[[2]]] >= 2.0
  merged_cpm <- rowSums(raw_matrix[, ix, drop = FALSE]) / sum(library_sizes[ix]) * 1e6
  merged_1 <- merged_cpm >= 1.0
  binary_05[, i] <- both_05
  binary_primary[, i] <- both_1
  binary_2[, i] <- both_2
  binary_merged[, i] <- merged_1
  log2_merged_cpm[, i] <- log2(merged_cpm + 1)

  definitions <- list(
    cpm0.5_both_reps = list(mask = both_05, directory = peak_dir_05),
    cpm1_both_reps = list(mask = both_1, directory = peak_dir),
    cpm2_both_reps = list(mask = both_2, directory = peak_dir_2),
    cpm1_merged_reps = list(mask = merged_1, directory = peak_dir_merged)
  )
  for (definition in names(definitions)) {
    item <- definitions[[definition]]
    peak_manifest_rows[[length(peak_manifest_rows) + 1L]] <- data.frame(
      system = "patient",
      sample_id = manifest$sample_id[[i]],
      definition = definition,
      peak_count = sum(item$mask),
      peak_path = write_peak_bed(item$mask, item$directory, manifest$sample_id[[i]]),
      stringsAsFactors = FALSE
    )
  }

  peak_count_rows[[i]] <- data.frame(
    sample_id = manifest$sample_id[[i]],
    cpm0.5_both_reps = sum(both_05),
    cpm1_both_reps = sum(both_1),
    cpm2_both_reps = sum(both_2),
    cpm1_merged_reps = sum(merged_1),
    tech1_library_size = library_sizes[ix[[1]]],
    tech2_library_size = library_sizes[ix[[2]]],
    stringsAsFactors = FALSE
  )
}

peak_counts <- rbindlist(peak_count_rows)
fwrite(peak_counts, file.path(audit_dir, "tcga_luad_peak_counts_threshold_sensitivity.tsv"), sep = "\t")
fwrite(
  rbindlist(peak_manifest_rows),
  file.path(audit_dir, "tcga_luad_peak_file_manifest.tsv"),
  sep = "\t"
)
fwrite(
  data.table(peak_id = counts$name, binary_primary * 1L),
  file.path(out_dir, "tcga_luad_binary_accessibility_cpm1_both_reps.tsv.gz"),
  sep = "\t"
)
saveRDS(
  list(
    coordinates = counts[, ..metadata_columns],
    binary_05 = binary_05,
    binary_primary = binary_primary,
    binary_2 = binary_2,
    binary_merged = binary_merged,
    log2_merged_cpm = log2_merged_cpm,
    library_sizes = library_sizes,
    thresholds = list(
      sensitivity_low = "CPM >= 0.5 in both technical replicates",
      primary = "CPM >= 1 in both technical replicates",
      sensitivity_high = "CPM >= 2 in both technical replicates",
      sensitivity_merged = "merged-replicate CPM >= 1"
    )
  ),
  file.path(out_dir, "tcga_luad_accessibility_matrices.rds"),
  compress = "xz"
)

correlation <- cor(log2_merged_cpm, method = "pearson", use = "pairwise.complete.obs")
fwrite(
  data.table(sample_id = rownames(correlation), correlation),
  file.path(out_dir, "tcga_luad_sample_pearson_correlation.tsv.gz"),
  sep = "\t"
)

receipt <- list(
  independent_tumors = nrow(manifest),
  technical_replicates = length(library_columns),
  candidate_fixed_width_peaks = nrow(counts),
  primary_definition = "CPM >= 1 independently in both technical replicates",
  sensitivity_definitions = c(
    "CPM >= 0.5 independently in both technical replicates",
    "CPM >= 2 independently in both technical replicates",
    "merged-replicate CPM >= 1"
  ),
  primary_peak_count_min = min(peak_counts$cpm1_both_reps),
  primary_peak_count_median = median(peak_counts$cpm1_both_reps),
  primary_peak_count_max = max(peak_counts$cpm1_both_reps),
  coordinate_interpretation = "GDC hg38 fixed-width intervals retained exactly as supplied",
  exact_fastq_reproduction = FALSE,
  reason = "TCGA raw/aligned ATAC is controlled; open technical-replicate counts and published QC were used"
)
json <- jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE)
writeLines(json, file.path(audit_dir, "tcga_luad_accessibility_receipt.json"))
cat(json, "\n")

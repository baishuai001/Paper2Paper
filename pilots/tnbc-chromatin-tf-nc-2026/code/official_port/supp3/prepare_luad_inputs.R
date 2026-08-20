#!/usr/bin/env Rscript

# Pure LUAD -> official-schema adapter for Supplementary Figure 3.
# This file contains no plotting constructor, theme, colour, geometry or scale.

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: prepare_luad_inputs.R LUAD_FIGURE2_ROOT ADAPTER_OUTDIR")
}

figure2_root <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("missing input: ", path)
  read.delim(path, check.names = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

count_text_rows <- function(path) {
  if (!file.exists(path)) stop("missing peak file: ", path)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    block <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(block)) break
    n <- n + length(block)
  }
  n
}

expected_n <- c(patient = 22L, PDX = 13L, cell_line = 19L)
cohort_for_system <- c(patient = "TCGA", PDX = "PDX", cell_line = "CellLines")

manifest_path <- file.path(
  figure2_root, "audit", "manifests", "figure2_atomic", "figure2_atomic_sample_manifest.tsv"
)
raw_qc_path <- file.path(figure2_root, "audit", "raw_atac_qc", "figure2_raw_atac_qc.tsv")
supp3_dir <- file.path(figure2_root, "results", "supplementary_figure3")
manifest <- read_tsv(manifest_path)
raw_qc <- read_tsv(raw_qc_path)

for (system_name in names(expected_n)) {
  observed <- length(unique(manifest$sample_id[manifest$system == system_name]))
  if (observed != expected_n[[system_name]]) {
    stop(system_name, " sample count mismatch: ", observed, " != ", expected_n[[system_name]])
  }
}

# A-C: counts in the exact two-column schema consumed by official 02J2.
peak_receipts <- list()
for (system_name in names(expected_n)) {
  cohort <- cohort_for_system[[system_name]]
  samples <- manifest[manifest$system == system_name, , drop = FALSE]
  if (system_name == "patient") {
    counts <- vapply(samples$peak_path, count_text_rows, integer(1))
    semantic_status <- "PARTIAL_FIXED_GDC_OPEN_COUNT_MATRIX_NOT_NEW_MACS2_IDR"
    idr_status <- "NOT_RERUN_PUBLIC_FIXED_MATRIX"
  } else {
    idx <- match(samples$sample_id, raw_qc$sample_id)
    if (anyNA(idx)) stop("raw-QC peak counts missing for ", system_name)
    counts <- as.integer(raw_qc$canonical_peak_count[idx])
    if (anyNA(counts)) stop("non-numeric raw-QC peak counts for ", system_name)
    semantic_status <- "MACS2_Q0.01_CANONICAL_PEAK_COUNT"
    idr_status <- "NOT_APPLICABLE_NO_TRUE_REPLICATE_SET"
  }
  number_peaks <- data.frame(Number_Peaks = counts, Sample = samples$sample_id)
  write_tsv(number_peaks, file.path(outdir, "number_peaks", paste0(cohort, "_Number_Peaks.tsv")))
  peak_receipts[[length(peak_receipts) + 1L]] <- data.frame(
    system = system_name,
    cohort = cohort,
    n_samples = nrow(number_peaks),
    min_peaks = min(counts),
    median_peaks = median(counts),
    max_peaks = max(counts),
    semantic_status = semantic_status,
    idr_status = idr_status
  )
}
write_tsv(do.call(rbind, peak_receipts), file.path(outdir, "audit", "peak_count_adapter_receipt.tsv"))

# D-F: translate the frozen 1000-permutation objects into the exact 02C schema.
permutations <- read_tsv(file.path(supp3_dir, "saturation_1000_permutations.tsv"))
summary_df <- read_tsv(file.path(supp3_dir, "saturation_summary.tsv"))
saturation_receipts <- list()
for (system_name in names(expected_n)) {
  cohort <- cohort_for_system[[system_name]]
  p <- permutations[permutations$system == system_name, , drop = FALSE]
  s <- summary_df[summary_df$system == system_name, , drop = FALSE]
  if (length(unique(p$permutation)) != 1000L) stop("expected 1000 permutations for ", system_name)
  if (max(p$sample_number) != expected_n[[system_name]]) stop("saturation sample count mismatch")
  s <- s[order(s$sample_number), , drop = FALSE]
  fit <- nls(
    mean_cumulative_consensus_peaks ~ SSasymp(sample_number, Asym, R0, lrc),
    data = s,
    control = nls.control(maxiter = 10000)
  )
  co <- coef(fit)
  asym <- unname(co[["Asym"]])
  r0 <- unname(co[["R0"]])
  lrc <- unname(co[["lrc"]])
  predicted <- function(n) asym + (r0 - asym) * exp(-exp(lrc) * n)
  max_samples <- 200L
  sample_ests <- data.frame(
    N_Samples = seq_len(max_samples),
    Est_N_Peaks = predicted(seq_len(max_samples)),
    Frac_Saturation = predicted(seq_len(max_samples)) / asym,
    Model = "All_Samples"
  )
  fractions <- c(0.5, 0.9, 0.95, 0.99)
  n_required <- -log((fractions - 1) * asym / (r0 - asym)) / exp(lrc)
  saturation_ests <- data.frame(
    Frac_Saturation = fractions,
    N_Peaks = asym * fractions,
    N_Samples = n_required,
    Model = "All_Samples"
  )
  bootstraps_long <- data.frame(
    N_Samples = p$sample_number,
    Iteration = p$permutation,
    N_Peaks = p$cumulative_consensus_peaks,
    Model = "All_Samples"
  )
  identified <- s$mean_cumulative_consensus_peaks[s$sample_number == max(s$sample_number)][[1]]
  residual_sum_squares <- sum(residuals(fit)^2)
  total_sum_squares <- sum((s$mean_cumulative_consensus_peaks - mean(s$mean_cumulative_consensus_peaks))^2)
  saturation_model <- list(
    bootstraps = bootstraps_long,
    bootstrap_summary = data.frame(
      N_Samples = s$sample_number,
      N_Peaks_Mean = s$mean_cumulative_consensus_peaks,
      N_Peaks_SEM = s$sem_cumulative_consensus_peaks
    ),
    model = list(
      Asym = asym, R0 = r0, lrc = lrc,
      RSS = residual_sum_squares,
      TSS = total_sum_squares,
      R_Squared = 1 - residual_sum_squares / total_sum_squares
    ),
    peaks = list(identified = identified, saturation = identified / asym),
    estimates = list(saturation = saturation_ests, samples = sample_ests),
    Frac_Saturation = identified / asym
  )
  saturation_out <- file.path(outdir, "saturation")
  dir.create(saturation_out, recursive = TRUE, showWarnings = FALSE)
  write_tsv(bootstraps_long, file.path(saturation_out, paste0(cohort, "_peak-saturation.bootstraps.tsv")))
  write_tsv(sample_ests, file.path(saturation_out, paste0(cohort, "_peak-saturation.model-estimates.tsv")))
  write_tsv(saturation_ests, file.path(saturation_out, paste0(cohort, "_peak-saturation.saturation-estimates.tsv")))
  saveRDS(saturation_model, file.path(saturation_out, paste0(cohort, "_saturation_model.rds")))
  saturation_receipts[[length(saturation_receipts) + 1L]] <- data.frame(
    system = system_name,
    cohort = cohort,
    n_samples = expected_n[[system_name]],
    permutations = length(unique(p$permutation)),
    asymptote = asym,
    identified = identified,
    fraction_saturation = identified / asym,
    model_r_squared = saturation_model$model$R_Squared
  )
}
write_tsv(do.call(rbind, saturation_receipts), file.path(outdir, "audit", "saturation_adapter_receipt.tsv"))

# G-I: materialise symmetric Pearson matrices consumed by official 03D.
cor_long <- read_tsv(file.path(supp3_dir, "sample_accessibility_pearson_correlations.tsv"))
cor_receipts <- list()
for (system_name in names(expected_n)) {
  cohort <- cohort_for_system[[system_name]]
  z <- cor_long[cor_long$system == system_name, , drop = FALSE]
  samples <- sort(unique(c(z$sample_1, z$sample_2)))
  if (length(samples) != expected_n[[system_name]]) stop("correlation sample count mismatch")
  mat <- matrix(NA_real_, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
  for (i in seq_len(nrow(z))) mat[z$sample_1[[i]], z$sample_2[[i]]] <- z$pearson_r[[i]]
  missing_reverse <- which(is.na(mat) & !is.na(t(mat)), arr.ind = TRUE)
  if (nrow(missing_reverse)) mat[missing_reverse] <- t(mat)[missing_reverse]
  diag(mat) <- 1
  if (anyNA(mat)) stop("incomplete correlation matrix for ", system_name)
  if (max(abs(mat - t(mat))) > 1e-8) stop("asymmetric correlation matrix for ", system_name)
  dir.create(file.path(outdir, "correlation"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(mat, file.path(outdir, "correlation", paste0(cohort, "_pearson_cor.rds")))
  cor_receipts[[length(cor_receipts) + 1L]] <- data.frame(
    system = system_name,
    cohort = cohort,
    n_samples = nrow(mat),
    min_r = min(mat),
    median_off_diagonal_r = median(mat[row(mat) != col(mat)]),
    max_r = max(mat)
  )
}
write_tsv(do.call(rbind, cor_receipts), file.path(outdir, "audit", "correlation_adapter_receipt.tsv"))

cell_qc <- raw_qc[raw_qc$system == "cell_line", , drop = FALSE]
cell_meta <- data.frame(
  sample_name = cell_qc$sample_id,
  Cohort = ifelse(toupper(cell_qc$layout) == "SINGLE", "Nergiz_se", "Other")
)
cell_meta <- cell_meta[cell_meta$sample_name %in% manifest$sample_id[manifest$system == "cell_line"], ]
if (nrow(cell_meta) != expected_n[["cell_line"]]) stop("cell-line sequencing metadata incomplete")
write_tsv(cell_meta, file.path(outdir, "correlation", "CellLines_metadata.tsv"))

# J-L: expand the existing per-sample count table into the official row schema.
annotation <- read_tsv(file.path(supp3_dir, "peak_genomic_annotation.tsv"))
annotation_map <- c(promoter = "Promoter", exonic = "Exonic", intronic = "Intronic", distal = "Distal")
annotation$official_annotation <- unname(annotation_map[tolower(annotation$category)])
if (anyNA(annotation$official_annotation)) stop("unknown genomic annotation category")
genomic_receipts <- list()
for (system_name in names(expected_n)) {
  cohort <- cohort_for_system[[system_name]]
  z <- annotation[annotation$system == system_name, , drop = FALSE]
  if (length(unique(z$sample_id)) != expected_n[[system_name]]) stop("genomic annotation sample count mismatch")
  target_dir <- file.path(outdir, "consensus", paste0(cohort, "_results"))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(target_dir, "PeakCalls_unionPeakSet.bed.gz")
  con <- gzfile(target, "wt", compression = 6)
  writeLines("GroupReplicate\tpeakType", con)
  for (i in seq_len(nrow(z))) {
    line <- paste(z$sample_id[[i]], z$official_annotation[[i]], sep = "\t")
    remaining <- as.integer(z$peak_count[[i]])
    while (remaining > 0L) {
      n_block <- min(remaining, 100000L)
      writeLines(rep.int(line, n_block), con)
      remaining <- remaining - n_block
    }
  }
  close(con)
  genomic_receipts[[length(genomic_receipts) + 1L]] <- data.frame(
    system = system_name,
    cohort = cohort,
    n_samples = length(unique(z$sample_id)),
    n_nonzero_samples = length(unique(z$sample_id[z$peak_count > 0])),
    expanded_rows = sum(z$peak_count),
    expanded_file = target,
    expanded_bytes = file.info(target)$size
  )
}
write_tsv(do.call(rbind, genomic_receipts), file.path(outdir, "audit", "genomic_annotation_adapter_receipt.tsv"))

baseline <- file.path(supp3_dir, "genome_annotation_baseline.tsv")
if (!file.exists(baseline)) stop("real hg38 500-bp baseline is missing")
header <- strsplit(readLines(baseline, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
if (!all(c("seqnames", "start", "end", "width", "strand", "type") %in% header)) {
  stop("500-bp baseline schema mismatch")
}
baseline_bins <- count_text_rows(baseline) - 1L
if (baseline_bins != 6062095L) stop("unexpected hg38 500-bp baseline bin count: ", baseline_bins)
writeLines(normalizePath(baseline), file.path(outdir, "GENOME_BASELINE_PATH.txt"))
write_tsv(
  data.frame(
    path = normalizePath(baseline), bytes = file.info(baseline)$size,
    bin_width_bp = 500L, bins = baseline_bins
  ),
  file.path(outdir, "audit", "genome_baseline_receipt.tsv")
)

write_tsv(
  data.frame(
    system = names(expected_n),
    cohort = unname(cohort_for_system[names(expected_n)]),
    n_samples = unname(expected_n),
    adapter_role = "schema_only_no_plotting"
  ),
  file.path(outdir, "audit", "cohort_receipt.tsv")
)

cat("Prepared strict official Supplementary Figure 3 inputs for 22 patients, 13 PDX and 19 cell lines.\n")

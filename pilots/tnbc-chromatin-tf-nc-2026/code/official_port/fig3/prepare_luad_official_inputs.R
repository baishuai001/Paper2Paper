#!/usr/bin/env Rscript

# Data-only adapter for the unmodified TNBC Figure 3 data flow.
#
# This file intentionally contains no plotting code. It converts frozen LUAD
# results to the column names, row-name conventions and subtype manifests read
# by the official CodeOcean Figure3 scripts.

suppressPackageStartupMessages({
  library(data.table)
})

parse_named_args <- function(values) {
  if (length(values) %% 2L != 0L || any(!grepl("^--", values[seq(1L, length(values), 2L)]))) {
    stop(
      paste(
        "Usage: prepare_luad_official_inputs.R",
        "--fig1-root PATH --fig2-root PATH --fig3-root PATH",
        "--publication-root PATH --output-root PATH"
      ),
      call. = FALSE
    )
  }
  keys <- sub("^--", "", values[seq(1L, length(values), 2L)])
  setNames(as.list(values[seq(2L, length(values), 2L)]), keys)
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
required <- c("fig1-root", "fig2-root", "fig3-root", "publication-root", "output-root")
missing_args <- setdiff(required, names(args))
if (length(missing_args)) {
  stop(sprintf("Missing argument(s): %s", paste(missing_args, collapse = ", ")), call. = FALSE)
}

fig1_root <- normalizePath(args[["fig1-root"]], mustWork = TRUE)
fig2_root <- normalizePath(args[["fig2-root"]], mustWork = TRUE)
fig3_root <- normalizePath(args[["fig3-root"]], mustWork = TRUE)
publication_root <- normalizePath(args[["publication-root"]], mustWork = TRUE)
output_root <- args[["output-root"]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)

input_dir <- file.path(output_root, "inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(input_dir, "MKI67"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(input_dir, "ESTIMATE"), recursive = TRUE, showWarnings = FALSE)

assert_columns <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) {
    stop(sprintf("%s lacks: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  }
}

read_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  if (ncol(x) < 2L) stop(sprintf("Matrix has fewer than two columns: %s", path), call. = FALSE)
  ids <- as.character(x[[1L]])
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

write_official_matrix <- function(x, path) {
  # Official scripts 03 and 06 call read.table(..., header=TRUE) without a
  # row.names argument, and then subset TFs by character row names. R only
  # infers the first field as row names when data rows have one more field than
  # the header. Therefore the header must contain sample names only. Contrary
  # to the usual write.table round-trip recipe, col.names=NA adds a leading
  # empty header field here and leaves numeric row names on readback.
  write.table(
    as.data.frame(x, check.names = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
  )
  probe_n <- min(2L, nrow(x))
  probe <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    nrows = probe_n
  )
  if (!identical(rownames(probe), head(rownames(x), probe_n)) ||
      !identical(colnames(probe), colnames(x))) {
    stop(sprintf(
      "Official read.table row-name round-trip failed for %s.",
      path
    ), call. = FALSE)
  }
}

first_existing <- function(paths, label) {
  selected <- paths[file.exists(paths)]
  if (!length(selected)) {
    stop(sprintf("No %s file found. Tried: %s", label, paste(paths, collapse = "; ")), call. = FALSE)
  }
  selected[[1L]]
}

# -------------------------------------------------------------------------
# HC-TF list and the official Figure 3A / 3B / 3D / 3E regulon schema.
# -------------------------------------------------------------------------

edge_path <- file.path(
  fig3_root, "results", "tables", "Figure3_HC_TF_signed_regulon_edges.tsv.gz"
)
edges <- fread(edge_path)
assert_columns(edges, c("TF", "target", "tfmode", "direction", "target_class"), "signed regulon")

hc_tfs <- sort(unique(edges$TF))
if (length(hc_tfs) != 31L) {
  stop(sprintf("Expected the frozen 31 LUAD HC-TFs; found %d.", length(hc_tfs)), call. = FALSE)
}
writeLines(hc_tfs, file.path(input_dir, "High_RNA_NES.tsv"), useBytes = TRUE)

official_viper <- unique(edges[, .(
  TR_Source = TF,
  Target_Gene = target,
  VIPER_Weight = as.numeric(tfmode)
)])
if (anyNA(official_viper$VIPER_Weight)) {
  stop("The frozen regulon contains non-numeric VIPER/tfmode weights.", call. = FALSE)
}
fwrite(
  official_viper,
  file.path(input_dir, "TCGA_Regulon_TF_Target_Weights_Viper.tsv"),
  sep = "\t"
)

# The official Figure 3A motif tile marks an HC-TF if the motif passed in at
# least one system. This is deliberately not the stricter three-system subset.
motif_path <- file.path(
  fig2_root, "results", "motif_primary", "anchor_consensus_author",
  "anchor_consensus_system_tf_motif_summary.tsv"
)
motif <- fread(motif_path)
assert_columns(
  motif,
  c("system", "TF", "TF_role", "support_fixed_original_denominator"),
  "motif summary"
)
motif_supported <- motif[
  TF %in% hc_tfs & support_fixed_original_denominator == TRUE,
  .(hc_group = "Motif"),
  by = TF
]
motif_supported <- unique(motif_supported, by = "TF")
motif_table <- data.frame(hc_group = motif_supported$hc_group, row.names = motif_supported$TF)
write.table(
  motif_table,
  file.path(input_dir, "JASPAR_CISBP_HC_TF_motif.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

# -------------------------------------------------------------------------
# Activity matrices and subtype files in the exact row-name conventions used
# by official scripts 03 and 06.
# -------------------------------------------------------------------------

activity_specs <- list(
  TCGA = list(
    matrix = file.path(fig1_root, "results", "tcga", "TCGA_viper_activity.tsv.gz"),
    manifest = file.path(fig1_root, "data", "processed", "tcga_manifest.tsv")
  ),
  PDX = list(
    matrix = file.path(fig1_root, "results", "pdmr", "PDMR_viper_activity.tsv.gz"),
    manifest = file.path(fig1_root, "data", "processed", "pdmr_manifest.tsv")
  ),
  CellLines = list(
    matrix = file.path(fig1_root, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
    manifest = file.path(fig1_root, "data", "processed", "depmap_manifest.tsv")
  )
)

matrix_receipts <- list()
for (cohort in names(activity_specs)) {
  spec <- activity_specs[[cohort]]
  matrix <- read_matrix(spec$matrix)
  if (length(setdiff(hc_tfs, rownames(matrix)))) {
    stop(sprintf("%s activity is missing HC-TFs: %s", cohort,
                 paste(setdiff(hc_tfs, rownames(matrix)), collapse = ", ")), call. = FALSE)
  }
  manifest <- fread(spec$manifest, check.names = FALSE)
  assert_columns(manifest, c("sample_id", "group"), paste(cohort, "manifest"))
  manifest <- manifest[sample_id %in% colnames(matrix)]
  if (!nrow(manifest) || !all(manifest$sample_id %in% colnames(matrix))) {
    stop(sprintf("%s manifest/activity alignment failed.", cohort), call. = FALSE)
  }
  write_official_matrix(matrix, file.path(input_dir, paste0(cohort, "_TR_activities_viper.tsv")))

  if (cohort == "TCGA") {
    subtype <- manifest[group == "LUAD", .(Sample.ID = sample_id)]
    fwrite(subtype, file.path(input_dir, "TCGA_LUAD.tsv"), sep = "\t")
  } else {
    subtype <- manifest[, .(
      Sample = sample_id,
      PAM50 = ifelse(group == "LUAD", "Basal", "Other"),
      SCMOD2 = ifelse(group == "LUAD", "Basal", "Other")
    )]
    fwrite(subtype, file.path(input_dir, paste0(cohort, "_molecular_subtype.tsv")), sep = "\t")
  }

  matrix_receipts[[cohort]] <- data.table(
    cohort = cohort,
    tf_rows = nrow(matrix),
    samples = ncol(matrix),
    luad_samples = sum(manifest$group == "LUAD"),
    non_luad_samples = sum(manifest$group != "LUAD")
  )
}
fwrite(rbindlist(matrix_receipts), file.path(output_root, "activity_adapter_receipt.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# Linked Supplementary Figure 4D. The correlations are frozen upstream
# results; they are converted to the official 07B/07C plotting schema. Four
# systems mirror the anchor (two patient cohorts, PDX, cell lines). GSE81089 is
# retained in the source table but not substituted for the independent cohort.
# -------------------------------------------------------------------------

mki67_source <- fread(file.path(
  publication_root, "results", "tables", "SupplementaryFigure4D_MKI67_correlations.tsv"
))
assert_columns(
  mki67_source,
  c("cohort", "TF", "pearson_r", "p_value", "FDR"),
  "MKI67 correlations"
)
mki67_map <- c(
  TCGA_LUAD = "TCGA",
  GSE41271_LUAD = "GSE41271",
  PDMR_LUAD = "PDX",
  DEPMAP_LUAD = "CellLines"
)
mki67_source <- mki67_source[cohort %in% names(mki67_map) & TF %in% hc_tfs]
mki67_source[, cohort_official := unname(mki67_map[cohort])]
for (cohort_name in unname(mki67_map)) {
  x <- mki67_source[cohort_official == cohort_name, .(
    TR = TF,
    r = as.numeric(pearson_r),
    p = as.numeric(p_value),
    fdr = as.numeric(FDR),
    sig = as.numeric(FDR) < 0.05 & abs(as.numeric(pearson_r)) > 0.4,
    label = fifelse(as.numeric(FDR) < 0.05 & abs(as.numeric(pearson_r)) > 0.4, TF, NA_character_),
    cohort = cohort_name
  )]
  if (nrow(x) != length(hc_tfs)) {
    stop(sprintf("MKI67 %s expected %d HC-TFs; found %d.", cohort_name, length(hc_tfs), nrow(x)), call. = FALSE)
  }
  fwrite(
    x,
    file.path(input_dir, "MKI67", paste0(cohort_name, "_MKi67_Correlation_HC-TR_results.tsv")),
    sep = "\t"
  )
}

# -------------------------------------------------------------------------
# Linked Supplementary Figure 5B. Convert frozen Pearson/FDR results to the
# official 08D/08E schema; official plot constructors then reapply thresholds.
# -------------------------------------------------------------------------

estimate_source <- fread(file.path(
  publication_root, "results", "tables", "SupplementaryFigure5B_ESTIMATE_correlations.tsv"
))
assert_columns(
  estimate_source,
  c("cohort", "score", "TF", "pearson_r", "p_value", "FDR"),
  "ESTIMATE correlations"
)
estimate_map <- c(TCGA_LUAD = "TCGA", GSE41271_LUAD = "GSE41271")
estimate_source <- estimate_source[cohort %in% names(estimate_map) & TF %in% hc_tfs]
estimate_source[, cohort_official := unname(estimate_map[cohort])]
for (cohort_name in unname(estimate_map)) {
  x <- estimate_source[cohort_official == cohort_name, .(
    TR = TF,
    marker = score,
    r = as.numeric(pearson_r),
    p = as.numeric(p_value),
    fdr = as.numeric(FDR),
    sig = as.numeric(FDR) < 0.05 & abs(as.numeric(pearson_r)) > 0.4,
    label = fifelse(as.numeric(FDR) < 0.05 & abs(as.numeric(pearson_r)) > 0.4, TF, NA_character_),
    cohort = cohort_name
  )]
  expected <- length(hc_tfs) * 3L
  if (nrow(x) != expected) {
    stop(sprintf("ESTIMATE %s expected %d TF-score rows; found %d.", cohort_name, expected, nrow(x)), call. = FALSE)
  }
  fwrite(
    x,
    file.path(input_dir, "ESTIMATE", paste0(cohort_name, "_ESTIMATE_all_results.tsv")),
    sep = "\t"
  )
}

adapter_manifest <- data.table(
  artifact = c(
    "High_RNA_NES.tsv",
    "TCGA_Regulon_TF_Target_Weights_Viper.tsv",
    "JASPAR_CISBP_HC_TF_motif.tsv",
    "TCGA_TR_activities_viper.tsv",
    "PDX_TR_activities_viper.tsv",
    "CellLines_TR_activities_viper.tsv",
    "MKI67/{TCGA,GSE41271,PDX,CellLines}_*.tsv",
    "ESTIMATE/{TCGA,GSE41271}_*.tsv"
  ),
  role = c(
    "official HC-TR list",
    "official Figure3A regulon input",
    "official Figure3A motif annotation",
    "official Figure3C/F/G patient activity",
    "official Figure3G/Supp5A PDX activity",
    "official Figure3G/Supp5A cell-line activity",
    "official SupplementaryFigure4D result schema",
    "official SupplementaryFigure5B result schema"
  ),
  biological_filter = c(
    "31 frozen LUAD HC-TFs",
    "all frozen signed ARACNe3/VIPER regulon edges",
    "motif passes original denominator rule in >=1 system",
    "LUAD versus non-LUAD samples from frozen TCGA manifest",
    "LUAD versus non-LUAD samples from frozen PDMR manifest",
    "LUAD versus non-LUAD samples from frozen DepMap manifest",
    "TCGA, GSE41271, PDMR, DepMap; |r|>0.4 and BH FDR<0.05 reapplied by official code",
    "TCGA and GSE41271; |r|>0.4 and BH FDR<0.05 reapplied by official code"
  )
)
fwrite(adapter_manifest, file.path(output_root, "adapter_manifest.tsv"), sep = "\t")

message(sprintf(
  "Prepared strict official Figure 3 inputs: %d HC-TFs, %d regulon edges, %d any-system motif TFs.",
  length(hc_tfs), nrow(official_viper), nrow(motif_supported)
))

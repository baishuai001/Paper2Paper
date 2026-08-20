#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: prepare_gse41271_microarray.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
raw_dir <- file.path(root, "data", "raw", "gse41271")
processed <- file.path(root, "data", "processed")
audit <- file.path(root, "audit", "gse41271")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

soft_path <- file.path(raw_dir, "GSE41271_family.soft.gz")
matrix_path <- file.path(raw_dir, "GSE41271_series_matrix.txt.gz")
if (!file.exists(soft_path) || !file.exists(matrix_path)) {
  stop("GSE41271 SOFT or Series Matrix is missing")
}

soft <- readLines(gzfile(soft_path), warn = FALSE)

sample_starts <- grep("^\\^SAMPLE = ", soft)
sample_ends <- c(sample_starts[-1L] - 1L, length(soft))
parse_sample <- function(index) {
  block <- soft[sample_starts[index]:sample_ends[index]]
  first_value <- function(prefix) {
    hit <- grep(prefix, block)
    if (!length(hit)) return(NA_character_)
    sub(prefix, "", block[hit[[1L]]])
  }
  chars <- sub(
    "^!Sample_characteristics_ch1 = ", "",
    block[grep("^!Sample_characteristics_ch1 = ", block)]
  )
  pairs <- strsplit(chars, ": ", fixed = TRUE)
  fields <- setNames(
    vapply(pairs, function(x) paste(x[-1L], collapse = ": "), character(1)),
    vapply(pairs, `[[`, character(1), 1L)
  )
  get_field <- function(name) {
    value <- unname(fields[name])
    if (!length(value)) NA_character_ else value
  }
  data.table(
    geo_accession = first_value("^!Sample_geo_accession = "),
    sample_title = first_value("^!Sample_title = "),
    histology = get_field("histology"),
    sex = get_field("gender"),
    race = get_field("race"),
    smoking = get_field("tobacco history"),
    stage = get_field("final patient stage"),
    date_birth = get_field("date of birth"),
    date_surgery = get_field("date of surgery"),
    vital_status = get_field("vital statistics"),
    date_last_survival = get_field("last follow-up survival"),
    recurrence = get_field("recurrence"),
    date_last_recurrence = get_field("last follow-up recurrence")
  )
}
manifest_all <- rbindlist(lapply(seq_along(sample_starts), parse_sample), fill = TRUE)
# All downstream VIPER helpers define `sample_id` as the expression-matrix
# column key.  GSE41271 matrix columns are GEO accessions, while the human
# specimen name is retained separately as `sample_title`.
manifest_all[, sample_id := geo_accession]
# Freeze the two principal histologies by exact label.  Substring matching is
# unsafe here because the 12-study-specified "other" tumours include two
# adenosquamous and one sarcomatoid-squamous specimen.
manifest_all[, histology_normalized := toupper(trimws(histology))]
manifest_all[, group := fcase(
  histology_normalized == "ADENOCARCINOMA", "LUAD",
  histology_normalized == "SQUAMOUS", "LUSC",
  default = "OTHER"
)]
if (nrow(manifest_all) != 275L ||
    manifest_all[group == "LUAD", .N] != 183L ||
    manifest_all[group == "LUSC", .N] != 80L ||
    manifest_all[group == "OTHER", .N] != 12L) {
  stop("GSE41271 histology counts do not match the official 183/80/12 design")
}
if (anyDuplicated(manifest_all$geo_accession) || anyDuplicated(manifest_all$sample_id)) {
  stop("Duplicate GSE41271 sample identifiers")
}

safe_date <- function(x) suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
birth <- safe_date(manifest_all$date_birth)
surgery <- safe_date(manifest_all$date_surgery)
last_survival <- safe_date(manifest_all$date_last_survival)
last_recurrence <- safe_date(manifest_all$date_last_recurrence)
manifest_all[, age_at_surgery := as.numeric(surgery - birth) / 365.25]
manifest_all[, os_days := as.numeric(last_survival - surgery)]
manifest_all[, os_event := fifelse(toupper(vital_status) == "D", 1L,
                                    fifelse(toupper(vital_status) == "A", 0L, NA_integer_))]
manifest_all[, rfs_days := as.numeric(last_recurrence - surgery)]
manifest_all[, rfs_event := fifelse(toupper(recurrence) == "Y", 1L,
                                     fifelse(toupper(recurrence) == "N", 0L, NA_integer_))]
manifest_all[, `:=`(system = "patient_microarray_validation", source = "GSE41271")]
manifest <- manifest_all[group %in% c("LUAD", "LUSC")]

matrix_lines <- readLines(gzfile(matrix_path), warn = FALSE)
matrix_start <- grep("^!series_matrix_table_begin", matrix_lines)
matrix_end <- grep("^!series_matrix_table_end", matrix_lines)
if (length(matrix_start) != 1L || length(matrix_end) != 1L || matrix_end <= matrix_start) {
  stop("Cannot locate GSE41271 Series Matrix table")
}
probe_expression <- fread(
  text = paste(matrix_lines[(matrix_start + 1L):(matrix_end - 1L)], collapse = "\n"),
  sep = "\t", header = TRUE, check.names = FALSE, na.strings = c("null", "NA")
)
setnames(probe_expression, 1L, "probe_id")
target_samples <- manifest$geo_accession
missing_samples <- setdiff(target_samples, names(probe_expression))
if (length(missing_samples)) {
  stop(sprintf("GSE41271 target samples missing from Series Matrix: %s",
               paste(missing_samples, collapse = ", ")))
}
probe_ids <- probe_expression$probe_id
expression <- as.matrix(probe_expression[, ..target_samples])
storage.mode(expression) <- "double"
rownames(expression) <- probe_ids
if (any(!is.finite(expression))) stop("Non-finite GSE41271 expression values")

platform_start <- grep("^!platform_table_begin", soft)
platform_end <- grep("^!platform_table_end", soft)
if (!length(platform_start) || !length(platform_end)) stop("Cannot locate GPL6884 platform table")
platform_end <- platform_end[platform_end > platform_start[[1L]]][[1L]]
platform <- fread(
  text = paste(soft[(platform_start[[1L]] + 1L):(platform_end - 1L)], collapse = "\n"),
  sep = "\t", header = TRUE, quote = "", fill = TRUE, check.names = FALSE
)
if (!all(c("ID", "Symbol") %in% names(platform))) stop("GPL6884 ID/Symbol columns missing")
platform <- platform[ID %in% rownames(expression), .(probe_id = ID, gene = trimws(Symbol))]
platform <- platform[
  !is.na(gene) & gene != "" & gene != "." &
    grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", gene)
]
platform <- unique(platform, by = "probe_id")

probe_mad <- apply(expression[platform$probe_id, , drop = FALSE], 1L, mad, na.rm = TRUE)
platform[, MAD := as.numeric(probe_mad[probe_id])]
setorder(platform, gene, -MAD, probe_id)
selected_probes <- platform[, .SD[1L], by = gene]
gene_expression <- expression[selected_probes$probe_id, , drop = FALSE]
rownames(gene_expression) <- selected_probes$gene
if (anyDuplicated(rownames(gene_expression))) stop("Duplicate genes after probe selection")
if (!identical(colnames(gene_expression), manifest$geo_accession)) {
  stop("GSE41271 expression columns do not match manifest order")
}

output <- as.data.table(gene_expression)
output[, gene := rownames(gene_expression)]
setcolorder(output, "gene")
fwrite(output, file.path(processed, "gse41271_aracne_expression.tsv"), sep = "\t")
fwrite(manifest, file.path(processed, "gse41271_manifest.tsv"), sep = "\t")
fwrite(manifest_all, file.path(processed, "gse41271_all_samples.tsv"), sep = "\t")
fwrite(selected_probes, file.path(processed, "gse41271_selected_probes.tsv"), sep = "\t")

pango_path <- file.path(root, "reference", "pango_regulators.txt")
if (!file.exists(pango_path)) stop("PAN-GO regulator list missing")
pango <- scan(pango_path, what = "character", quiet = TRUE)
regulators <- sort(intersect(pango, rownames(gene_expression)))
writeLines(regulators, file.path(processed, "gse41271_aracne_regulators.txt"))

receipt <- list(
  status = "passed",
  accession = "GSE41271",
  platform = "GPL6884 Illumina HumanWG-6 v3.0",
  all_tumors = nrow(manifest_all),
  group_counts_all = as.list(table(manifest_all$group)),
  analysis_samples = nrow(manifest),
  group_counts = as.list(table(manifest$group)),
  matrix_probes = nrow(expression),
  single_symbol_mapped_probes = nrow(platform),
  selected_genes = nrow(gene_expression),
  pango_regulators = length(pango),
  regulators_in_expression = length(regulators),
  expression_scale = "GEO processed log2-normalized microarray intensity",
  duplicate_probe_rule = "highest MAD across the frozen 263-sample LUAD/LUSC cohort; probe ID breaks exact ties",
  ambiguous_symbol_rule = "exclude blank and multi-symbol probe annotations",
  package_versions = list(R = as.character(getRversion()), data.table = as.character(packageVersion("data.table")))
)
write_json(receipt, file.path(audit, "gse41271_preparation_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

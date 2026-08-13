#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: prepare_tcga_crc.R OUTPUT_DIR")
}
output_dir <- normalizePath(args[[1]], mustWork = FALSE)
raw_dir <- file.path(output_dir, "gdc_raw")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

query <- GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
manifest <- as.data.frame(getResults(query))
fwrite(manifest, file.path(output_dir, "gdc_query_manifest.tsv"), sep = "\t", quote = FALSE)
if (nrow(manifest) < 600) {
  stop(sprintf("Only %d GDC primary-tumor files returned", nrow(manifest)))
}

download_marker <- file.path(output_dir, ".gdc_download_complete")
if (!file.exists(download_marker)) {
  GDCdownload(query, directory = raw_dir, method = "api", files.per.chunk = 20)
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), download_marker)
}

se <- GDCprepare(query, directory = raw_dir, summarizedExperiment = TRUE)
assay_names <- assayNames(se)
tpm_candidates <- grep("^tpm", assay_names, ignore.case = TRUE, value = TRUE)
if (length(tpm_candidates) != 1) {
  stop(sprintf("Expected exactly one TPM assay, found: %s", paste(tpm_candidates, collapse = ", ")))
}
tpm <- as.matrix(assay(se, tpm_candidates[[1]]))
if (any(!is.finite(tpm)) || any(tpm < 0)) {
  stop("GDC TPM assay contains non-finite or negative values")
}

annotation <- as.data.frame(rowData(se))
symbol_field <- intersect(c("gene_name", "gene_symbol", "external_gene_name"), colnames(annotation))
if (length(symbol_field) == 0) {
  stop(sprintf("No gene symbol field in rowData: %s", paste(colnames(annotation), collapse = ", ")))
}
symbols <- trimws(as.character(annotation[[symbol_field[[1]]]]))
valid <- !is.na(symbols) & nzchar(symbols) & symbols != "NA"
tpm <- tpm[valid, , drop = FALSE]
symbols <- symbols[valid]

# GDC annotations occasionally contain multiple Ensembl rows for one symbol.
# Sum them before any filtering so every output row has one unique gene symbol.
tpm_by_symbol <- rowsum(tpm, group = symbols, reorder = TRUE, na.rm = TRUE)
sample_barcodes <- colnames(tpm_by_symbol)
if (is.null(sample_barcodes) || any(nchar(sample_barcodes) < 12)) {
  stop("Invalid or missing TCGA sample barcodes")
}
participants <- substr(sample_barcodes, 1, 12)
participant_sums <- t(rowsum(t(tpm_by_symbol), group = participants, reorder = TRUE, na.rm = TRUE))
participant_n <- table(participants)
participant_tpm <- sweep(participant_sums, 2, as.numeric(participant_n[colnames(participant_sums)]), "/")

prevalent <- rowMeans(participant_tpm >= 1) >= 0.10
variable <- apply(participant_tpm, 1, var) > 0
keep <- prevalent & variable & is.finite(rowMeans(participant_tpm))
filtered_tpm <- participant_tpm[keep, , drop = FALSE]
if (ncol(filtered_tpm) < 550 || nrow(filtered_tpm) < 10000) {
  stop(sprintf("Frozen data minimum failed: %d participants, %d genes", ncol(filtered_tpm), nrow(filtered_tpm)))
}
if (anyDuplicated(rownames(filtered_tpm)) || anyDuplicated(colnames(filtered_tpm))) {
  stop("Duplicate gene or participant identifiers after collapse")
}

sample_metadata <- data.frame(
  sample_barcode = sample_barcodes,
  participant = participants,
  project = manifest$project[match(sample_barcodes, manifest$cases)],
  stringsAsFactors = FALSE
)
fwrite(sample_metadata, file.path(output_dir, "tcga_sample_to_participant.tsv"), sep = "\t", quote = FALSE)
fwrite(
  data.frame(gene = rownames(participant_tpm), retained = keep, stringsAsFactors = FALSE),
  file.path(output_dir, "tcga_gene_filter.tsv"), sep = "\t", quote = FALSE
)
saveRDS(filtered_tpm, file.path(output_dir, "tcga_crc_tpm.rds"), compress = "xz")
fwrite(
  data.table(gene = rownames(filtered_tpm), filtered_tpm, keep.rownames = FALSE),
  file.path(output_dir, "tcga_crc_tpm.tsv"), sep = "\t", quote = FALSE
)

receipt <- list(
  status = "passed",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  projects = as.list(table(manifest$project)),
  gdc_files = nrow(manifest),
  tpm_assay = tpm_candidates[[1]],
  input_rows = nrow(tpm),
  unique_gene_symbols_before_filter = nrow(participant_tpm),
  retained_genes = nrow(filtered_tpm),
  aliquots = length(sample_barcodes),
  participants = ncol(filtered_tpm),
  participants_with_multiple_aliquots = sum(participant_n > 1),
  filter = list(tpm_ge_1_fraction = 0.10, variance_gt = 0),
  package_versions = list(
    R = as.character(getRversion()),
    TCGAbiolinks = as.character(packageVersion("TCGAbiolinks")),
    SummarizedExperiment = as.character(packageVersion("SummarizedExperiment")),
    data.table = as.character(packageVersion("data.table"))
  )
)
write_json(receipt, file.path(output_dir, "tcga_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

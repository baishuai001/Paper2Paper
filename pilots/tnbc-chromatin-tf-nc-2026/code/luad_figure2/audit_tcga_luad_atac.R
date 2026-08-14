#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readxl))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: audit_tcga_luad_atac.R <figure2_run_root> <figure1_run_root>")
}
run_root <- normalizePath(args[[1]], mustWork = TRUE)
figure1_root <- normalizePath(args[[2]], mustWork = TRUE)
raw_dir <- file.path(run_root, "data", "raw", "tcga_atac_gdc")
audit_dir <- file.path(run_root, "audit", "manifests", "tcga")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "TCGA_identifier_mapping.txt",
  "TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx",
  "TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip"
)
missing <- required[!file.exists(file.path(raw_dir, required))]
if (length(missing)) stop("Missing TCGA inputs: ", paste(missing, collapse = ", "))

mapping <- read.delim(
  file.path(raw_dir, "TCGA_identifier_mapping.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
mapping <- mapping[!is.na(mapping$bam_prefix) & grepl("^LUAD-", mapping$bam_prefix), , drop = FALSE]

workbook <- file.path(raw_dir, "TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx")
donors <- as.data.frame(
  read_excel(workbook, sheet = "Donor Characteristics", skip = 23),
  stringsAsFactors = FALSE
)
stats <- as.data.frame(
  read_excel(workbook, sheet = "Sequencing Statistics", skip = 37),
  stringsAsFactors = FALSE
)
donors <- donors[!is.na(donors$cohort) & donors$cohort == "LUAD", , drop = FALSE]
stats <- stats[!is.na(stats$cohort) & stats$cohort == "LUAD", , drop = FALSE]

zip_path <- file.path(raw_dir, "TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip")
# Base R's unz connection fails on this large ZIP64 archive on some builds;
# Info-ZIP is already required to inspect the official GDC bundle and streams
# only the single header line here.
con <- pipe(sprintf("unzip -p %s LUAD_raw_counts.txt", shQuote(zip_path)), open = "rt")
header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
close(con)
library_columns <- header[-seq_len(5L)]

stanford_uuid <- vapply(
  library_columns,
  function(x) {
    pieces <- strsplit(sub("^LUAD_", "", x), "_", fixed = TRUE)[[1]]
    if (length(pieces) < 5L) stop("Cannot parse LUAD count column: ", x)
    paste(pieces[1:5], collapse = "-")
  },
  character(1)
)

map_by_uuid <- split(mapping, toupper(mapping$stanfordUUID))
rows <- lapply(seq_along(library_columns), function(i) {
  uuid <- toupper(stanford_uuid[[i]])
  hit <- map_by_uuid[[uuid]]
  if (is.null(hit) || nrow(hit) < 1L) stop("No identifier mapping for ", uuid)
  hit <- hit[grepl(paste0("-T", sub(".*_T([0-9]+)_.*", "\\1", library_columns[[i]]), "-"), hit$bam_prefix), , drop = FALSE]
  if (nrow(hit) != 1L) {
    # Pool suffixes can differ between the count header and mapping.  The
    # Stanford UUID and technical-replicate token still identify one row.
    hit <- map_by_uuid[[uuid]]
    tech <- sub(".*_(T[0-9]+)_.*", "\\1", library_columns[[i]])
    hit <- hit[grepl(paste0("-", tech, "-"), hit$bam_prefix), , drop = FALSE]
  }
  if (nrow(hit) != 1L) stop("Ambiguous mapping for ", library_columns[[i]])
  st <- stats[toupper(stats$Stanford_UUID) == uuid &
                stats$Technical_Replicate_Number == sub(".*_(T[0-9]+)_.*", "\\1", library_columns[[i]]), , drop = FALSE]
  if (nrow(st) != 1L) stop("Missing sequencing-statistics row for ", library_columns[[i]])
  donor <- donors[!is.na(donors$submitter_id) & donors$submitter_id == st$submitter_id, , drop = FALSE]
  if (nrow(donor) != 1L) stop("Missing donor row for ", st$submitter_id)
  data.frame(
    count_column = library_columns[[i]],
    count_column_index = i + 5L,
    stanford_uuid = stanford_uuid[[i]],
    technical_replicate = st$Technical_Replicate_Number,
    donor_id = st$submitter_id,
    case_id = st$case_id,
    tissue_barcode = st$Tissue_Barcode,
    tissue_sample_type = st$Tissue_Sample_Type,
    published_qc_pass = as.logical(donor[["Passed_QC?"]]),
    has_rna_seq = as.logical(donor[["Has RNA-seq?"]]),
    raw_reads = as.numeric(st$Raw_reads),
    final_dedup_reads = as.numeric(st$Final_DeDup_reads),
    tss_enrichment_score = as.numeric(st$TSS_Enrichment_Score),
    frip = as.numeric(st$FRIP),
    aliquot_id = hit$aliquot_id,
    stringsAsFactors = FALSE
  )
})
technical <- do.call(rbind, rows)
technical <- technical[order(technical$donor_id, technical$technical_replicate), ]

activity_path <- file.path(
  figure1_root, "results", "tables", "Figure1C_TCGA_activity_matrix.tsv.gz"
)
if (!file.exists(activity_path)) stop("Missing frozen Figure 1 TCGA activity matrix")
activity_header <- strsplit(readLines(gzfile(activity_path), n = 1L), "\t", fixed = TRUE)[[1]][-1L]
activity_donor <- substr(activity_header, 1L, 12L)

by_donor <- split(technical, technical$donor_id)
sample_rows <- lapply(by_donor, function(d) {
  data.frame(
    sample_id = d$donor_id[[1]],
    case_id = d$case_id[[1]],
    stanford_uuid = d$stanford_uuid[[1]],
    tissue_barcode = d$tissue_barcode[[1]],
    tissue_sample_type = d$tissue_sample_type[[1]],
    independent_biological_sample = TRUE,
    technical_replicate_count = nrow(d),
    count_columns = paste(d$count_column, collapse = ";"),
    published_qc_pass = all(d$published_qc_pass),
    has_rna_seq_published = all(d$has_rna_seq),
    matched_frozen_figure1_activity = d$donor_id[[1]] %in% activity_donor,
    median_tss_enrichment_score = median(d$tss_enrichment_score),
    median_frip = median(d$frip),
    total_final_dedup_reads = sum(d$final_dedup_reads),
    figure2_primary_include = (
      all(d$published_qc_pass) &&
      all(d$tissue_sample_type == "Primary Solid Tumor") &&
      nrow(d) == 2L
    ),
    stringsAsFactors = FALSE
  )
})
samples <- do.call(rbind, sample_rows)
samples <- samples[order(samples$sample_id), ]

write.table(
  technical,
  file.path(audit_dir, "tcga_luad_atac_technical_replicates.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)
write.table(
  samples,
  file.path(audit_dir, "tcga_luad_atac_samples.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)

if (nrow(technical) != 44L) stop("Expected 44 LUAD technical replicates; observed ", nrow(technical))
if (nrow(samples) != 22L) stop("Expected 22 independent LUAD samples; observed ", nrow(samples))
if (!all(samples$figure2_primary_include)) stop("Not all 22 LUAD samples passed the frozen data gate")

receipt <- list(
  technical_replicates = nrow(technical),
  independent_primary_tumors = nrow(samples),
  technical_replicates_per_tumor = as.integer(table(samples$technical_replicate_count)),
  published_qc_pass = sum(samples$published_qc_pass),
  published_rna_seq_available = sum(samples$has_rna_seq_published),
  matched_frozen_figure1_activity = sum(samples$matched_frozen_figure1_activity),
  median_tss_enrichment = median(samples$median_tss_enrichment_score),
  median_frip = median(samples$median_frip),
  count_matrix = "TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip::LUAD_raw_counts.txt",
  count_matrix_genome = "hg38",
  independence_rule = "one Stanford tissue UUID / TCGA donor; two ATAC reactions are technical replicates",
  primary_inclusion_rule = "published QC pass; Primary Solid Tumor; exactly two technical replicates"
)
json <- jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE)
writeLines(json, file.path(audit_dir, "tcga_luad_atac_manifest_receipt.json"))
cat(json, "\n")

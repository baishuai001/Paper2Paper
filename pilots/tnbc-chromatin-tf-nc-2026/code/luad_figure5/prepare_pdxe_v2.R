#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_pdxe_v2.R PDXE_TSV_DIR OUTPUT_ROOT")
source_dir <- normalizePath(args[[1]], mustWork = TRUE)
out <- normalizePath(args[[2]], mustWork = TRUE)
processed <- file.path(out, "data", "processed")
audit <- file.path(out, "audit")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

rna_path <- file.path(source_dir, "molecular_profiles", "RNASeq.tsv")
model_path <- file.path(source_dir, "model.tsv")
if (!file.exists(rna_path) || !file.exists(model_path)) stop("PDXE v2 TSV inputs missing")

header <- names(fread(rna_path, nrows = 0L, check.names = FALSE))
if (length(header) < 2L || header[[1]] != "feature_rowname") stop("Unexpected RNASeq.tsv schema")
samples <- header[-1]
model <- fread(model_path, na.strings = c("", "NA"))
required <- c("model.id", "patient.id", "tissue", "tissue.name", "drug")
if (!all(required %in% names(model))) stop("Unexpected model.tsv schema")

patient_meta <- unique(model[, .(patient_id = as.character(patient.id), tissue, tissue_name = tissue.name)])
patient_meta <- patient_meta[, .(
  tissue = paste(sort(unique(na.omit(tissue))), collapse = ";"),
  tissue_name = paste(sort(unique(na.omit(tissue_name))), collapse = ";")
), by = patient_id]
manifest <- data.table(sample_id = samples, patient_id = samples)
manifest <- merge(manifest, patient_meta, by = "patient_id", all.x = TRUE, sort = FALSE)
setcolorder(manifest, c("sample_id", "patient_id", "tissue", "tissue_name"))
manifest[, group := fifelse(tissue == "NSCLC", "NSCLC", fifelse(is.na(tissue) | tissue == "", "unknown", tissue))]
manifest[, system := "PDXE_v2_PDX"]
manifest[, source := "Xeva PDXE v2 / Zenodo 10.5281/zenodo.19354439"]
manifest[, sample_order := match(sample_id, samples)]
setorder(manifest, sample_order)
manifest[, sample_order := NULL]

if (!identical(manifest$sample_id, samples)) stop("Manifest order differs from RNA matrix")
fwrite(manifest, file.path(processed, "pdxe_v2_manifest.tsv"), sep = "\t")

receipt <- list(
  status = "passed",
  expression_file = normalizePath(rna_path),
  genes = nrow(fread(rna_path, select = 1L)),
  expression_samples = length(samples),
  metadata_matched = sum(!is.na(manifest$tissue)),
  NSCLC_samples = sum(manifest$tissue == "NSCLC", na.rm = TRUE),
  tissue_counts = as.list(table(manifest$group, useNA = "ifany")),
  response_models = uniqueN(model$model.id),
  response_patients = uniqueN(model$patient.id)
)
write_json(receipt, file.path(audit, "pdxe_v2_preparation_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

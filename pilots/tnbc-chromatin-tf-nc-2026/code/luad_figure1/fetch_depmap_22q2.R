#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: fetch_depmap_22q2.R OUTPUT_DIR AUDIT_DIR")
}
output_dir <- args[[1]]
audit_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

resources <- list(
  metadata = list(id = "EH7558", fetch_id = "7608", object = "metadata_22Q2"),
  tpm = list(id = "EH7556", fetch_id = "7606", object = "TPM_22Q2")
)

# ExperimentHub stores these fixed 22Q2 resources as xz-compressed R
# workspaces (RDX3), not as single-object RDS files.  Loading the fixed fetch
# endpoints directly avoids a moving ExperimentHub snapshot while preserving
# the exact EH7558/EH7556 resource versions used by the TNBC capsule.
fetch_workspace <- function(resource) {
  compressed <- file.path(output_dir, sprintf("%s_22Q2.RData.xz", resource$id))
  workspace <- tempfile(fileext = ".RData")
  on.exit(unlink(workspace), add = TRUE)
  url <- sprintf("https://experimenthub.bioconductor.org/fetch/%s", resource$fetch_id)
  if (!file.exists(compressed) || file.info(compressed)$size == 0) {
    download.file(url, compressed, mode = "wb", quiet = FALSE)
  }
  status <- system2("xz", c("-dc", shQuote(compressed)), stdout = workspace)
  if (status != 0) stop(sprintf("xz failed for %s", resource$id))
  environment <- new.env(parent = emptyenv())
  loaded <- load(workspace, envir = environment)
  if (!identical(loaded, resource$object)) {
    stop(sprintf("Expected object %s in %s; found %s", resource$object, resource$id, paste(loaded, collapse = ", ")))
  }
  get(resource$object, envir = environment, inherits = FALSE)
}

metadata <- fetch_workspace(resources$metadata)
tpm <- fetch_workspace(resources$tpm)

metadata_source_rows <- nrow(metadata)
tpm_source_rows <- nrow(tpm)
strict_luad <- !is.na(metadata$lineage) & !is.na(metadata$subtype_disease) &
  metadata$lineage == "lung" & metadata$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma"
strict_lusc <- !is.na(metadata$lineage) & !is.na(metadata$subtype_disease) &
  metadata$lineage == "lung" & metadata$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma"
metadata <- metadata[strict_luad | strict_lusc, , drop = FALSE]
tpm <- tpm[tpm$depmap_id %in% metadata$depmap_id, , drop = FALSE]

fwrite(as.data.table(metadata), file.path(output_dir, "depmap_metadata_22Q2.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")
fwrite(as.data.table(tpm), file.path(output_dir, "depmap_TPM_22Q2.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

receipt <- list(
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  release = "22Q2",
  metadata_resource = "EH7558",
  tpm_resource = "EH7556",
  metadata_fetch_url = "https://experimenthub.bioconductor.org/fetch/7608",
  tpm_fetch_url = "https://experimenthub.bioconductor.org/fetch/7606",
  metadata_source_rows = metadata_source_rows,
  tpm_source_rows = tpm_source_rows,
  metadata_rows = nrow(metadata),
  tpm_rows = nrow(tpm),
  metadata_columns = names(metadata),
  tpm_columns = names(tpm),
  package_versions = list(
    R = as.character(getRversion())
  )
)
write_json(receipt, file.path(audit_dir, "depmap_receipt.json"), pretty = TRUE, auto_unbox = TRUE)

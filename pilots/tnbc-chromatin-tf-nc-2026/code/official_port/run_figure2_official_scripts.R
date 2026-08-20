#!/usr/bin/env Rscript

# Invoke the materialised CodeOcean Figure 2 plotting scripts directly.
# This wrapper contains no plotting constructor.

suppressPackageStartupMessages(library(data.table))

usage <- "Usage: run_figure2_official_scripts.R RUNTIME_ROOT OFFICIAL_PORT_ROOT"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(usage, call. = FALSE)

runtime_root <- normalizePath(args[[1L]], mustWork = TRUE)
port_root <- normalizePath(args[[2L]], mustWork = TRUE)
official_aesthetics <- file.path(
  port_root, "upstream_edf5314", "code", "Static_Scripts", "plotting_aesthetics.R"
)
Sys.setenv(PAPER2PAPER_OFFICIAL_AESTHETICS = official_aesthetics)

rscript <- file.path(R.home("bin"), "Rscript")
runtime_code <- file.path(runtime_root, "code")
log_dir <- file.path(runtime_root, "logs")
audit_dir <- file.path(runtime_root, "audit")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

jobs <- list(
  list(panel = "Figure2B", source = "Figure2/07B-UpSet_Plot_Promoter_Accessibility.R", args = character()),
  list(panel = "Figure2C_CAPSULE_OUTPUT", source = "Figure2/06C-PromoterAccessibility_MotifEnrichment.R", args = character()),
  list(panel = "Figure2D", source = "Figure2/07A-UpSet_Plot_Motif_Enrichment.R", args = c("--category", "JASPAR_CISBP")),
  list(panel = "Figure2E", source = "Figure2/08-Plot_TR_Examples.R", args = character()),
  list(panel = "SupplementaryFigure4B", source = "Figure3/02-Determine_HC_TR.R", args = character()),
  list(panel = "SupplementaryFigure4C", source = "Figure3/01-Compare_JASPAR_CISBP_Motifs.R", args = character())
)

receipt <- list()
for (job in jobs) {
  script <- file.path(runtime_code, job$source)
  if (!file.exists(script)) stop("Materialised source missing: ", script)
  log_path <- file.path(log_dir, paste0(job$panel, ".log"))
  started <- Sys.time()
  status <- system2(
    rscript,
    args = c(script, job$args),
    stdout = log_path,
    stderr = log_path,
    env = paste0("PAPER2PAPER_OFFICIAL_AESTHETICS=", official_aesthetics)
  )
  receipt[[length(receipt) + 1L]] <- data.table(
    panel = job$panel,
    official_source = job$source,
    status = if (identical(status, 0L)) "GENERATED_OFFICIAL_CONSTRUCTOR" else "FAILED",
    exit_status = status,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    log = log_path,
    publication_boundary = if (job$panel == "Figure2C_CAPSULE_OUTPUT") {
      "CAPSULE_PRINT_MISMATCH;PAIRED_TCGA_ATAC_RNA_N21_OF_N22"
    } else if (job$panel == "SupplementaryFigure4C") {
      "ANCHOR_EMPTY_SET_BUGFIX"
    } else {
      "NONE"
    }
  )
  if (!identical(status, 0L)) {
    fwrite(rbindlist(receipt), file.path(audit_dir, "figure2_official_script_receipt.tsv"), sep = "\t")
    stop(job$panel, " failed; inspect ", log_path, call. = FALSE)
  }
}

fwrite(rbindlist(receipt), file.path(audit_dir, "figure2_official_script_receipt.tsv"), sep = "\t")
cat("Figure 2 official scripts completed.\n")

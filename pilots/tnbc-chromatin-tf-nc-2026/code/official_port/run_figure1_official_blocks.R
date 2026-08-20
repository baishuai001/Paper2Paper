#!/usr/bin/env Rscript

# Execute only the official Figure 1 plotting blocks from CodeOcean edf5314.
# Plot expressions are read from the materialised official source; none are
# reimplemented in this wrapper.  Dendrogram export side calculations between
# plot blocks are intentionally not evaluated because they are not constructors.

suppressPackageStartupMessages(library(data.table))

usage <- paste(
  "Usage: run_figure1_official_blocks.R RUNTIME_ROOT OFFICIAL_PORT_ROOT",
  "[VALIDATION_ROLE] [VALIDATION_DATASET] [ANCHOR_BUG_STATUS]"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 5L) stop(usage, call. = FALSE)

runtime_root <- normalizePath(args[[1L]], mustWork = TRUE)
port_root <- normalizePath(args[[2L]], mustWork = TRUE)
validation_role <- if (length(args) >= 3L) args[[3L]] else "PRIMARY_METABRIC_ANALOGUE"
validation_dataset <- if (length(args) >= 4L) args[[4L]] else "GSE41271"
anchor_bug_status <- if (length(args) >= 5L) args[[5L]] else "CAPSULE_BUG_RETAINED"
official_aesthetics <- file.path(
  port_root, "upstream_edf5314", "code", "Static_Scripts", "plotting_aesthetics.R"
)
if (!file.exists(official_aesthetics)) stop("Frozen official aesthetics file is missing.")
Sys.setenv(PAPER2PAPER_OFFICIAL_AESTHETICS = official_aesthetics)

runtime_f1 <- file.path(runtime_root, "code", "Figure1")
source09 <- file.path(runtime_f1, "09-Compare_TCGA_METABRIC.R")
source10 <- file.path(runtime_f1, "10-Plot_Figure1.R")
if (!all(file.exists(c(source09, source10)))) stop("Materialised Figure 1 sources are missing.")

receipt <- list()
run_full_source <- function(panel, path) {
  started <- Sys.time()
  sys.source(path, envir = .GlobalEnv)
  receipt[[length(receipt) + 1L]] <<- data.table(
    panel = panel,
    official_source = basename(path),
    official_line_block = "full source",
    status = "GENERATED_OFFICIAL_CONSTRUCTOR",
    validation_role = validation_role,
    validation_dataset = validation_dataset,
    anchor_bug_status = anchor_bug_status,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}

lines <- readLines(source10, warn = FALSE)
run_block <- function(panel, start, end) {
  started <- Sys.time()
  code <- paste(lines[start:end], collapse = "\n")
  eval(parse(text = code, keep.source = TRUE), envir = .GlobalEnv)
  receipt[[length(receipt) + 1L]] <<- data.table(
    panel = panel,
    official_source = basename(source10),
    official_line_block = paste0(start, "-", end),
    status = "GENERATED_OFFICIAL_CONSTRUCTOR",
    validation_role = validation_role,
    validation_dataset = validation_dataset,
    anchor_bug_status = anchor_bug_status,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}

# Supplementary Figure 2B: official eulerr constructor.
run_full_source("SupplementaryFigure2B", source09)

# Package imports, official aesthetics source, and output-directory creation.
run_block("Figure1_prelude", 14L, 26L)

# Exact plot blocks in the official Figure1/10 source.
run_block("Figure1B", 28L, 52L)
run_block("Figure1C", 54L, 129L)
run_block("SupplementaryFigure2C", 150L, 175L)
run_block("Figure1D", 177L, 217L)
run_block("SupplementaryFigure2D", 238L, 262L)
run_block("Figure1E", 264L, 315L)
run_block("SupplementaryFigure2E", 336L, 353L)
run_block("SupplementaryFigure2A", 355L, 398L)
run_block("SupplementaryFigure2F", 400L, 607L)

audit_dir <- file.path(runtime_root, "audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(receipt), file.path(audit_dir, "figure1_official_block_receipt.tsv"), sep = "\t")
cat("Figure 1 official blocks completed.\n")

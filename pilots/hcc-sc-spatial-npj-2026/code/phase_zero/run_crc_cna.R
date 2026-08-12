#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) stop("arguments must be --key value pairs")
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

required_arg <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || !nzchar(value)) stop(sprintf("missing --%s", key))
  value
}

read_input <- function(input_dir) {
  required <- file.path(input_dir, c(
    "counts_genes_by_cells.mtx.gz", "genes.tsv", "cells.tsv", "cell_metadata.tsv"
  ))
  missing <- required[!file.exists(required) | file.info(required)$size == 0]
  if (length(missing)) stop(sprintf("missing or empty CNA input: %s", paste(missing, collapse = ", ")))
  counts <- readMM(gzfile(required[[1L]]))
  genes <- read.delim(required[[2L]], header = FALSE, stringsAsFactors = FALSE)[[1L]]
  cells <- read.delim(required[[3L]], header = FALSE, stringsAsFactors = FALSE)[[1L]]
  metadata <- read.delim(required[[4L]], check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(counts) != length(genes) || ncol(counts) != length(cells)) stop("matrix dimensions disagree with genes/cells")
  if (anyDuplicated(genes)) stop("gene symbols are not unique")
  if (anyDuplicated(cells)) stop("cell IDs are not unique")
  if (!identical(as.character(metadata$cell_id), as.character(cells))) stop("metadata cell order disagrees with matrix")
  if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-6)) stop("CNA matrix is not nonnegative integer counts")
  rownames(counts) <- genes
  colnames(counts) <- cells
  list(counts = counts, metadata = metadata)
}

write_input_receipt <- function(input, output_dir, method, sample_name) {
  receipt <- data.frame(
    method = method,
    sample_name = sample_name,
    genes = nrow(input$counts),
    cells = ncol(input$counts),
    known_normal_reference_cells = sum(input$metadata$is_known_normal_reference %in% c(TRUE, "TRUE", "True", "true")),
    author_cancer_cells = sum(input$metadata$cna_input_role == "author_cancer"),
    other_epithelial_cells = sum(input$metadata$cna_input_role == "other_epithelial"),
    stringsAsFactors = FALSE
  )
  write.table(receipt, file.path(output_dir, "P0_cna_runner_input_receipt.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

validate_method_runtime <- function(method) {
  if (method == "copykat") {
    if (!requireNamespace("copykat", quietly = TRUE)) stop("copykat package is not installed")
    suppressPackageStartupMessages(library(copykat))
    if (as.character(utils::packageVersion("copykat")) != "1.1.0") {
      stop(sprintf("copykat version must be 1.1.0, got %s", utils::packageVersion("copykat")))
    }
    required_data <- c("full.anno", "DNA.hg20", "cyclegenes")
    unavailable <- required_data[!vapply(required_data, exists, logical(1), inherits = TRUE)]
    if (length(unavailable)) {
      stop(sprintf("copykat LazyData are unavailable after library(copykat): %s", paste(unavailable, collapse = ", ")))
    }
  }
  if (method == "scevan") {
    if (!requireNamespace("SCEVAN", quietly = TRUE)) stop("SCEVAN package is not installed")
    if (as.character(utils::packageVersion("SCEVAN")) != "1.0.3") {
      stop(sprintf("SCEVAN version must be 1.0.3, got %s", utils::packageVersion("SCEVAN")))
    }
    actual <- names(formals(SCEVAN::pipelineCNA))
    required <- c("count_mtx", "norm_cell", "FIXED_NORMAL_CELLS", "SUBCLONES", "ClonalCN", "output_dir")
    missing <- setdiff(required, actual)
    if (length(missing)) {
      stop(sprintf("SCEVAN pipelineCNA API changed; missing formals: %s", paste(missing, collapse = ", ")))
    }
  }
}

run_copykat <- function(input, output_dir, sample_name, cores, seed) {
  known_normal <- input$metadata$cell_id[
    input$metadata$is_known_normal_reference %in% c(TRUE, "TRUE", "True", "true")
  ]
  if (length(known_normal) < 20L) stop("fewer than 20 known-normal reference cells")
  set.seed(seed)
  old <- setwd(output_dir)
  on.exit(setwd(old), add = TRUE)
  result <- copykat(
    rawmat = as.matrix(input$counts),
    id.type = "S",
    cell.line = "no",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = sample_name,
    distance = "euclidean",
    norm.cell.names = known_normal,
    output.seg = FALSE,
    plot.genes = FALSE,
    genome = "hg20",
    n.cores = cores
  )
  saveRDS(result, file.path(output_dir, "copykat_result.rds"))
  prediction <- as.data.frame(result$prediction, stringsAsFactors = FALSE)
  write.table(prediction, file.path(output_dir, "copykat_prediction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(result$CNAmat)) {
    saveRDS(result$CNAmat, file.path(output_dir, "copykat_cna_matrix.rds"))
  }
}

run_scevan <- function(input, output_dir, sample_name, cores, seed) {
  if (!requireNamespace("SCEVAN", quietly = TRUE)) stop("SCEVAN package is not installed")
  known_normal <- input$metadata$cell_id[
    input$metadata$is_known_normal_reference %in% c(TRUE, "TRUE", "True", "true")
  ]
  if (length(known_normal) < 20L) stop("fewer than 20 known-normal reference cells")
  # Matrix::readMM() returns a dgTMatrix.  SCEVAN 1.0.3 only recognises
  # dgCMatrix in annotateGenes(); other sparse classes fall through to a
  # data.frame cbind and fail.  Convert the sparse representation without
  # changing values, cell order, or gene order.  SCEVAN then performs its own
  # bounded dense conversion after filtering genes.
  scevan_counts <- as(input$counts, "CsparseMatrix")
  if (!inherits(scevan_counts, "dgCMatrix")) {
    stop(sprintf("SCEVAN requires dgCMatrix after sparse conversion, got %s", class(scevan_counts)[[1L]]))
  }
  if (!identical(dim(scevan_counts), dim(input$counts)) ||
      !identical(rownames(scevan_counts), rownames(input$counts)) ||
      !identical(colnames(scevan_counts), colnames(input$counts))) {
    stop("SCEVAN sparse-class conversion changed matrix dimensions or dimnames")
  }
  set.seed(seed)
  old <- setwd(output_dir)
  on.exit(setwd(old), add = TRUE)
  # SCEVAN 1.0.3 passes output_dir to most writers, but plotCNclonal() calls
  # getScevanCNV() without forwarding that path and therefore reads from the
  # hard-coded relative directory ./output.  Expose the already-selected
  # output directory under that legacy name for the duration of the official
  # pipeline.  This changes path resolution only; no SCEVAN function or result
  # is patched.  Refuse to reuse a pre-existing path so stale files cannot be
  # consumed silently.
  legacy_output_alias <- file.path(output_dir, "output")
  alias_target <- Sys.readlink(legacy_output_alias)
  alias_is_link <- !is.na(alias_target) && nzchar(alias_target)
  if (file.exists(legacy_output_alias) || alias_is_link) {
    stop("SCEVAN legacy ./output alias already exists")
  }
  if (!isTRUE(file.symlink(output_dir, legacy_output_alias))) {
    stop("failed to create SCEVAN legacy ./output path alias")
  }
  on.exit(unlink(legacy_output_alias), add = TRUE)
  result <- SCEVAN::pipelineCNA(
    scevan_counts,
    sample = sample_name,
    par_cores = cores,
    norm_cell = known_normal,
    FIXED_NORMAL_CELLS = FALSE,
    SUBCLONES = FALSE,
    ClonalCN = FALSE,
    plotTree = FALSE,
    organism = "human",
    output_dir = output_dir
  )
  if (!is.data.frame(result) || is.null(rownames(result)) || anyDuplicated(rownames(result))) {
    stop("SCEVAN pipelineCNA did not return a uniquely indexed classification data.frame")
  }
  saveRDS(result, file.path(output_dir, "scevan_result.rds"))
  prediction <- data.frame(cell.names = rownames(result), result, check.names = FALSE, stringsAsFactors = FALSE)
  write.table(prediction, file.path(output_dir, "scevan_prediction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(capture.output(str(result)), file.path(output_dir, "scevan_result_structure.txt"))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  input_dir <- normalizePath(required_arg(args, "input-dir"), mustWork = TRUE)
  output_dir <- required_arg(args, "output-dir")
  method <- required_arg(args, "method")
  sample_name <- required_arg(args, "sample-name")
  cores <- as.integer(if (is.null(args$cores)) "4" else args$cores)
  seed <- as.integer(if (is.null(args$seed)) "20260810" else args$seed)
  validate_only <- identical(if (is.null(args$`validate-only`)) "false" else tolower(args$`validate-only`), "true")
  if (!method %in% c("copykat", "scevan")) stop("--method must be copykat or scevan")
  if (is.na(cores) || cores < 1L) stop("--cores must be positive")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  input <- read_input(input_dir)
  validate_method_runtime(method)
  write_input_receipt(input, output_dir, method, sample_name)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  if (!validate_only) {
    if (method == "copykat") run_copykat(input, output_dir, sample_name, cores, seed)
    if (method == "scevan") run_scevan(input, output_dir, sample_name, cores, seed)
  }
  cat(sprintf("status=%s method=%s genes=%d cells=%d\n", if (validate_only) "validated" else "completed", method, nrow(input$counts), ncol(input$counts)))
}

main()

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(viper)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste(
    "Usage: run_viper_anchor.R NETWORK_TSV NETWORK_EXPRESSION MANIFEST",
    "COHORT OUTPUT_DIR MODE THREADS PERMUTATIONS"
  ))
}
network_path <- args[[1]]
expression_path <- args[[2]]
manifest_path <- args[[3]]
cohort <- args[[4]]
output_dir <- args[[5]]
mode <- args[[6]]
threads <- as.integer(args[[7]])
permutations <- as.integer(args[[8]])
if (!mode %in% c("discover", "project")) stop("MODE must be discover or project")
if (mode == "discover" && permutations != 1000L) stop("Anchor discovery requires exactly 1000 permutations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2025)

read_expression <- function(path) {
  tab <- fread(path, check.names = FALSE)
  genes <- tab[[1]]
  tab[[1]] <- NULL
  matrix <- as.matrix(tab)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  if (anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop(sprintf("Duplicate genes/samples in %s", path))
  }
  matrix
}

# viper 1.38.0 has a dimension-drop bug in aREA's weighted loop when a
# regulon has exactly one retained target and the permutation matrix contains
# missing values.  It applies over columns, obtains a vector, and then calls
# colSums() on that vector.  Restoring the known 1 x n matrix shape changes no
# values and preserves the anchor's minsize = 1 definition.
patch_area_single_target_dimension <- function() {
  if (as.character(packageVersion("viper")) != "1.38.0") {
    stop("The aREA compatibility patch is validated only for viper 1.38.0")
  }
  namespace <- asNamespace("viper")
  original <- get("aREA", envir = namespace)
  source <- paste(deparse(original, width.cutoff = 500L), collapse = "\n")
  target <- "sqrt(colSums(apply(tmp, 2, function(x) x/max(x))^2))"
  replacement <- paste0(
    "sqrt(colSums(matrix(apply(tmp, 2, function(x) x/max(x))^2, ",
    "nrow = nrow(tmp), ncol = ncol(tmp))))"
  )
  if (!grepl(target, source, fixed = TRUE)) {
    stop("Cannot locate the viper 1.38.0 aREA dimension-drop expression")
  }
  patched_source <- sub(target, replacement, source, fixed = TRUE)
  patched <- eval(parse(text = patched_source), envir = namespace)
  assignInNamespace("aREA", patched, ns = "viper")
  TRUE
}

area_single_target_compat_patch <- patch_area_single_target_dimension()

network <- fread(network_path)
required <- c("regulator.values", "target.values", "mi.values")
if (!all(required %in% names(network))) stop("ARACNe3 subnetwork columns are not anchor-compatible")
# The published Figure1/06 script reads the header and then removes the first
# data row while commenting "Remove the header". Preserve that executable-code
# behavior exactly, even though it discards one valid edge.
network_input_edges <- nrow(network)
network <- network[-1]
expression <- read_expression(expression_path)
regulon_path <- file.path(output_dir, sprintf("%s_regulon.rds", cohort))
if (file.exists(regulon_path) && file.info(regulon_path)$size > 0) {
  message("Loading checkpointed regulon: ", regulon_path)
  regulon <- readRDS(regulon_path)
} else {
  network_file <- tempfile(fileext = ".tsv")
  fwrite(network[, ..required], network_file, sep = "\t", col.names = FALSE)
  regulon <- aracne2regulon(network_file, expression, format = "3col", verbose = TRUE)
  unlink(network_file)
  saveRDS(regulon, regulon_path, compress = "xz")
}
single_target_regulons <- sum(vapply(regulon, function(x) length(x$tfmode), integer(1)) == 1L)

activity_path <- file.path(output_dir, sprintf("%s_viper_activity.tsv.gz", cohort))
if (file.exists(activity_path) && file.info(activity_path)$size > 0) {
  message("Preserving checkpointed VIPER activity: ", activity_path)
  activity_TF_count <- nrow(fread(activity_path, select = 1L))
} else {
  activity <- viper(
    expression,
    regulon,
    verbose = TRUE,
    minsize = 1,
    cores = threads
  )
  activity_TFs <- rownames(activity)
  activity_TF_count <- nrow(activity)
  activity_output <- as.data.table(activity)
  activity_output[, TF := activity_TFs]
  setcolorder(activity_output, "TF")
  fwrite(activity_output, activity_path, sep = "\t")
}

manifest <- fread(manifest_path)
manifest <- manifest[match(colnames(expression), sample_id)]
if (anyNA(manifest$sample_id) || !identical(manifest$sample_id, colnames(expression))) {
  stop("Expression columns do not map one-to-one to manifest")
}

ms_table <- NULL
if (mode == "discover") {
  msviper_rds_path <- file.path(output_dir, sprintf("%s_msviper.rds", cohort))
  msviper_table_path <- file.path(output_dir, sprintf("%s_msviper.tsv", cohort))
  if (file.exists(msviper_rds_path) && file.info(msviper_rds_path)$size > 0 &&
      file.exists(msviper_table_path) && file.info(msviper_table_path)$size > 0) {
    message("Preserving checkpointed msVIPER result: ", msviper_table_path)
    ms_table <- fread(msviper_table_path)
  } else {
    luad <- expression[, manifest$group == "LUAD", drop = FALSE]
    lusc <- expression[, manifest$group == "LUSC", drop = FALSE]
    if (min(ncol(luad), ncol(lusc)) < 10) stop("Discovery group has fewer than 10 independent samples")
    signature_raw <- rowTtest(luad, lusc, "N")
    signature <- (qnorm(signature_raw$p.value / 2, lower.tail = FALSE) * sign(signature_raw$statistic))[, 1]
    # The author's scientific definition is preserved (1,000 with-replacement
    # permutations, seed 1). The official cores argument distributes independent
    # permutations because this 1,017-patient transfer is much larger than the
    # anchor cohort. L'Ecuyer-CMRG makes the parallel streams reproducible.
    RNGkind("L'Ecuyer-CMRG")
    nullmodel <- ttestNull(
      luad, lusc, per = permutations, repos = TRUE, seed = 1,
      cores = threads, verbose = TRUE
    )
    mrs <- msviper(signature, regulon, nullmodel, verbose = TRUE, minsize = 1, cores = threads)
    saveRDS(mrs, msviper_rds_path, compress = "xz")
    ms_table <- data.table(
      TF = names(mrs$es$nes),
      NES = as.numeric(mrs$es$nes),
      p_value = as.numeric(mrs$es$p.value),
      FDR = p.adjust(as.numeric(mrs$es$p.value), method = "BH"),
      regulon_size = as.numeric(mrs$es$size)
    )
    fwrite(ms_table, msviper_table_path, sep = "\t")
  }
}

receipt <- list(
  status = "passed",
  cohort = cohort,
  mode = mode,
  samples = ncol(expression),
  genes = nrow(expression),
  group_counts = as.list(table(manifest$group)),
  subnetwork_input_edges = network_input_edges,
  subnetwork_edges_after_anchor_first_row_drop = nrow(network),
  anchor_first_row_drop_reproduced = TRUE,
  regulons = length(regulon),
  viper_TFs = activity_TF_count,
  permutations = if (mode == "discover") permutations else 0,
  msviper_TFs = if (is.null(ms_table)) 0 else nrow(ms_table),
  viper_script_seed = 2025,
  ttestNull_seed = if (mode == "discover") 1 else 0,
  ttestNull_cores = if (mode == "discover") threads else 0,
  ttestNull_rng_kind = if (mode == "discover") "L'Ecuyer-CMRG" else "not_applicable",
  minsize = 1,
  single_target_regulons = single_target_regulons,
  area_single_target_dimension_patch = area_single_target_compat_patch,
  area_patch_effect = "shape restoration only; minsize and numeric definition unchanged",
  package_versions = list(R = as.character(getRversion()), viper = as.character(packageVersion("viper")))
)
write_json(receipt, file.path(output_dir, sprintf("%s_viper_receipt.json", cohort)), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

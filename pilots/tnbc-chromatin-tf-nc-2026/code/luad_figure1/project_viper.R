#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(viper)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: project_viper.R REGULON_RDS EXPRESSION_TSV MANIFEST COHORT OUTPUT_DIR THREADS")
}
regulon_path <- args[[1]]
expression_path <- args[[2]]
manifest_path <- args[[3]]
cohort <- args[[4]]
output_dir <- args[[5]]
threads <- as.integer(args[[6]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2025)

tab <- fread(expression_path, check.names = FALSE)
genes <- tab[[1]]
tab[[1]] <- NULL
expression <- as.matrix(tab)
storage.mode(expression) <- "double"
rownames(expression) <- genes
manifest <- fread(manifest_path)
manifest <- manifest[match(colnames(expression), sample_id)]
if (anyNA(manifest$sample_id)) stop("Expression/manifest mismatch")
regulon <- readRDS(regulon_path)
activity <- viper(expression, regulon, verbose = TRUE, minsize = 1, cores = threads)
activity_TFs <- rownames(activity)
activity_output <- as.data.table(activity)
activity_output[, TF := activity_TFs]
setcolorder(activity_output, "TF")
fwrite(activity_output, file.path(output_dir, sprintf("%s_viper_activity.tsv.gz", cohort)), sep = "\t")
receipt <- list(
  status = "passed", cohort = cohort, samples = ncol(expression), genes = nrow(expression),
  group_counts = as.list(table(manifest$group)), regulons = length(regulon), TFs = nrow(activity),
  seed = 2025, minsize = 1,
  package_versions = list(R = as.character(getRversion()), viper = as.character(packageVersion("viper")))
)
write_json(receipt, file.path(output_dir, sprintf("%s_projection_receipt.json", cohort)), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

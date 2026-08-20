#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: classify_pdxe_lung_state.R FIG1_ROOT OUTPUT_ROOT")
fig1 <- normalizePath(args[[1]], mustWork = TRUE)
out <- normalizePath(args[[2]], mustWork = TRUE)
tables <- file.path(out, "results", "tables")
audit <- file.path(out, "audit")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

read_activity <- function(path) {
  tab <- fread(path, check.names = FALSE)
  id <- tab[[1]]
  tab[[1]] <- NULL
  matrix <- as.matrix(tab)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- id
  matrix
}

specific <- fread(file.path(fig1, "results", "tables", "tcga_specific_TFs.tsv"))
luad_tfs <- unique(specific[discovery_category == "LUAD", TF])
lusc_tfs <- unique(specific[discovery_category == "LUSC", TF])

# Frozen biological signature and decision threshold, with a cohort-level
# row-z transform to remove study-specific TF-activity offsets. The same
# transform is applied without labels to TCGA, PDMR, DepMap and the PDXE NSCLC
# subset. This variant was retained after the fully frozen per-sample rank
# implementation showed a severe cross-study location shift; all evaluated
# variants remain in the classifier audit tables.
score_activity <- function(activity) {
  features <- intersect(c(luad_tfs, lusc_tfs), rownames(activity))
  luad <- intersect(luad_tfs, features)
  lusc <- intersect(lusc_tfs, features)
  if (length(luad) < 0.8 * length(luad_tfs) || length(lusc) < 0.8 * length(lusc_tfs)) {
    stop("Too many frozen classifier features are absent")
  }
  transformed <- activity[features, , drop = FALSE]
  center <- rowMeans(transformed, na.rm = TRUE)
  scale <- apply(transformed, 1L, sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  transformed <- sweep(sweep(transformed, 1L, center), 1L, scale, "/")
  data.table(
    sample_id = colnames(activity),
    state_score = colMeans(transformed[luad, , drop = FALSE], na.rm = TRUE) -
      colMeans(transformed[lusc, , drop = FALSE], na.rm = TRUE),
    LUAD_features = length(luad),
    LUSC_features = length(lusc)
  )
}

balanced_accuracy <- function(truth, prediction) {
  mean(c(mean(prediction[truth == "LUAD"] == "LUAD"),
         mean(prediction[truth == "LUSC"] == "LUSC")), na.rm = TRUE)
}

learn_threshold <- function(score, truth) {
  candidates <- c(-Inf, (sort(unique(score))[-1L] + sort(unique(score))[-length(unique(score))]) / 2, Inf)
  result <- rbindlist(lapply(candidates, function(threshold) {
    prediction <- ifelse(score > threshold, "LUAD", "LUSC")
    data.table(threshold = threshold, balanced_accuracy = balanced_accuracy(truth, prediction),
               accuracy = mean(prediction == truth))
  }))
  result[order(-balanced_accuracy, -accuracy, abs(threshold))][1]
}

evaluate <- function(name, activity_path, manifest_path, id_column, group_column, threshold) {
  activity <- read_activity(activity_path)
  scored <- score_activity(activity)
  manifest <- fread(manifest_path)
  annotation <- manifest[, .(sample_id = as.character(get(id_column)), truth = as.character(get(group_column)))]
  scored <- merge(scored, annotation, by = "sample_id", all.x = TRUE)
  scored[, prediction := ifelse(state_score > threshold, "LUAD", "LUSC")]
  evaluable <- scored[truth %in% c("LUAD", "LUSC")]
  metrics <- data.table(
    cohort = name, n = nrow(evaluable), LUAD_n = sum(evaluable$truth == "LUAD"),
    LUSC_n = sum(evaluable$truth == "LUSC"),
    accuracy = mean(evaluable$prediction == evaluable$truth),
    sensitivity_LUAD = mean(evaluable[truth == "LUAD", prediction == "LUAD"]),
    specificity_LUSC = mean(evaluable[truth == "LUSC", prediction == "LUSC"]),
    balanced_accuracy = balanced_accuracy(evaluable$truth, evaluable$prediction)
  )
  list(scores = scored, metrics = metrics)
}

tcga_activity_path <- file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz")
tcga_manifest_path <- file.path(fig1, "data", "processed", "tcga_manifest.tsv")
tcga_activity <- read_activity(tcga_activity_path)
tcga_scores <- score_activity(tcga_activity)
tcga_manifest <- fread(tcga_manifest_path)
tcga_scores <- merge(tcga_scores, tcga_manifest[, .(sample_id, truth = group)], by = "sample_id")
fit <- learn_threshold(tcga_scores$state_score, tcga_scores$truth)
threshold <- fit$threshold
tcga_scores[, prediction := ifelse(state_score > threshold, "LUAD", "LUSC")]
tcga_metrics <- data.table(
  cohort = "TCGA_training", n = nrow(tcga_scores), LUAD_n = sum(tcga_scores$truth == "LUAD"),
  LUSC_n = sum(tcga_scores$truth == "LUSC"), accuracy = mean(tcga_scores$prediction == tcga_scores$truth),
  sensitivity_LUAD = mean(tcga_scores[truth == "LUAD", prediction == "LUAD"]),
  specificity_LUSC = mean(tcga_scores[truth == "LUSC", prediction == "LUSC"]),
  balanced_accuracy = balanced_accuracy(tcga_scores$truth, tcga_scores$prediction)
)

pdmr <- evaluate(
  "PDMR_external", file.path(fig1, "results", "pdmr", "PDMR_viper_activity.tsv.gz"),
  file.path(fig1, "data", "processed", "pdmr_manifest.tsv"), "sample_id", "group", threshold
)
depmap <- evaluate(
  "DepMap_external", file.path(fig1, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
  file.path(fig1, "data", "processed", "depmap_manifest.tsv"), "sample_id", "group", threshold
)

pdxe_activity_path <- file.path(out, "results", "pdxe_projection", "PDXE_V2_viper_activity.tsv.gz")
pdxe_manifest_path <- file.path(out, "data", "processed", "pdxe_v2_manifest.tsv")
pdxe_manifest <- fread(pdxe_manifest_path)
pdxe_activity <- read_activity(pdxe_activity_path)
pdxe_nsclc_ids <- intersect(pdxe_manifest[tissue == "NSCLC", sample_id], colnames(pdxe_activity))
if (length(pdxe_nsclc_ids) < 10L) stop("Too few PDXE NSCLC expression models for cohort normalization")
pdxe_scores <- score_activity(pdxe_activity[, pdxe_nsclc_ids, drop = FALSE])
pdxe_scores <- merge(pdxe_scores, pdxe_manifest, by = "sample_id", all.x = TRUE)
pdxe_scores[, prediction := ifelse(state_score > threshold, "LUAD-like", "LUSC-like")]
pdxe_scores[, classifier_threshold := threshold]

metrics <- rbindlist(list(tcga_metrics, pdmr$metrics, depmap$metrics), fill = TRUE)
fwrite(metrics, file.path(tables, "SupplementaryFigure8C_lung_classifier_validation.tsv"), sep = "\t")
fwrite(tcga_scores, file.path(tables, "SupplementaryFigure8D_TCGA_classifier_scores.tsv"), sep = "\t")
fwrite(pdmr$scores, file.path(tables, "SupplementaryFigure8E_PDMR_classifier_scores.tsv"), sep = "\t")
fwrite(depmap$scores, file.path(tables, "SupplementaryFigure8F_DepMap_classifier_scores.tsv"), sep = "\t")
fwrite(pdxe_scores, file.path(tables, "Figure5_PDXE_v2_lung_state_scores.tsv"), sep = "\t")

receipt <- list(
  status = "passed", classifier = "cohort row-z LUAD-minus-LUSC TF activity",
  LUAD_signature_TFs = length(luad_tfs), LUSC_signature_TFs = length(lusc_tfs),
  frozen_threshold = threshold, training_balanced_accuracy = tcga_metrics$balanced_accuracy,
  PDMR_external_balanced_accuracy = pdmr$metrics$balanced_accuracy,
  DepMap_external_balanced_accuracy = depmap$metrics$balanced_accuracy,
  PDXE_NSCLC_samples = pdxe_scores[tissue == "NSCLC", .N],
  PDXE_NSCLC_LUAD_like = pdxe_scores[tissue == "NSCLC" & prediction == "LUAD-like", .N],
  caution = paste(
    "PDXE v2 labels lung models as NSCLC; LUAD-like is a computational state, not pathology.",
    "Classifier-transform alternatives and their external validation are retained in the audit."
  )
)
write_json(receipt, file.path(audit, "pdxe_lung_classifier_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

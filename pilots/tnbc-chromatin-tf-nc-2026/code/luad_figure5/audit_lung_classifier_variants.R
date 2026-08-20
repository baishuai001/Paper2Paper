#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: audit_lung_classifier_variants.R FIG1_ROOT OUTPUT_ROOT")
fig1 <- normalizePath(args[[1]], mustWork = TRUE)
out <- normalizePath(args[[2]], mustWork = TRUE)

read_activity <- function(path) {
  x <- fread(path, check.names = FALSE); ids <- x[[1]]; x[[1]] <- NULL
  m <- as.matrix(x); storage.mode(m) <- "double"; rownames(m) <- ids; m
}
row_z <- function(m) {
  center <- rowMeans(m, na.rm = TRUE)
  scale <- apply(m, 1L, sd, na.rm = TRUE); scale[!is.finite(scale) | scale == 0] <- 1
  sweep(sweep(m, 1L, center), 1L, scale, "/")
}
row_robust_z <- function(m) {
  center <- apply(m, 1L, median, na.rm = TRUE)
  scale <- apply(m, 1L, mad, na.rm = TRUE); scale[!is.finite(scale) | scale == 0] <- 1
  sweep(sweep(m, 1L, center), 1L, scale, "/")
}
rank_rows_in_sample <- function(m) apply(m, 2L, function(x) rank(x, ties.method = "average") / length(x))
auc <- function(score, truth) {
  y <- truth == "LUAD"; (sum(rank(score)[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
balacc <- function(truth, pred) mean(c(mean(pred[truth == "LUAD"] == "LUAD"), mean(pred[truth == "LUSC"] == "LUSC")))
learn_threshold <- function(score, truth) {
  u <- sort(unique(score)); cuts <- c(-Inf, (u[-1L] + u[-length(u)]) / 2, Inf)
  z <- rbindlist(lapply(cuts, function(t) data.table(threshold=t, balacc=balacc(truth, ifelse(score > t,"LUAD","LUSC")))))
  z[order(-balacc, abs(threshold))][1, threshold]
}

spec <- fread(file.path(fig1, "results", "tables", "tcga_specific_TFs.tsv"))
luad <- spec[discovery_category == "LUAD", unique(TF)]
lusc <- spec[discovery_category == "LUSC", unique(TF)]
features <- c(luad, lusc)
hc <- fread(file.path(out, "..", "luad-figure2", "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
hc <- hc[HC_TF_promoter_activity_definition == TRUE, TF]

paths <- list(
  TCGA = list(activity=file.path(fig1,"results/tcga/TCGA_viper_activity.tsv.gz"), manifest=file.path(fig1,"data/processed/tcga_manifest.tsv"), id="sample_id"),
  PDMR = list(activity=file.path(fig1,"results/pdmr/PDMR_viper_activity.tsv.gz"), manifest=file.path(fig1,"data/processed/pdmr_manifest.tsv"), id="sample_id"),
  DepMap = list(activity=file.path(fig1,"results/depmap/DEPMAP_viper_activity.tsv.gz"), manifest=file.path(fig1,"data/processed/depmap_manifest.tsv"), id="sample_id")
)
data <- lapply(paths, function(z) {
  m <- read_activity(z$activity); a <- fread(z$manifest)[, .(sample_id=as.character(get(z$id)), truth=group)]
  list(m=m, truth=a$truth[match(colnames(m),a$sample_id)])
})
features <- Reduce(intersect, c(list(features), lapply(data, function(z) rownames(z$m))))
luad <- intersect(luad, features)
lusc <- intersect(lusc, features)
hc <- Reduce(intersect, c(list(hc), lapply(data, function(z) rownames(z$m))))
if (length(luad) < 0.8 * sum(spec$discovery_category == "LUAD") ||
    length(lusc) < 0.8 * sum(spec$discovery_category == "LUSC") || length(hc) < 20L) {
  stop("Too many classifier features are missing from the common projection space")
}

tcga_z <- row_z(data$TCGA$m[features,,drop=FALSE])
centroids <- sapply(c("LUAD","LUSC"), function(g) rowMeans(tcga_z[,data$TCGA$truth==g,drop=FALSE]))

make_scores <- function(z, method) {
  m <- z$m[features,,drop=FALSE]
  if (method == "within_sample_rank") {
    q <- rank_rows_in_sample(m); colMeans(q[luad,,drop=FALSE])-colMeans(q[lusc,,drop=FALSE])
  } else if (method == "cohort_row_z_signature") {
    q <- row_z(m); colMeans(q[luad,,drop=FALSE])-colMeans(q[lusc,,drop=FALSE])
  } else if (method == "cohort_robust_z_signature") {
    q <- row_robust_z(m); colMeans(q[luad,,drop=FALSE])-colMeans(q[lusc,,drop=FALSE])
  } else if (method == "cohort_row_z_centroid") {
    q <- row_z(m)
    apply(q,2L,function(x) cor(x,centroids[,"LUAD"],use="pairwise")-cor(x,centroids[,"LUSC"],use="pairwise"))
  } else if (method == "HC_cohort_row_z") {
    q <- row_z(z$m[hc,,drop=FALSE]); colMeans(q)
  } else stop("Unknown method")
}

methods <- c("within_sample_rank","cohort_row_z_signature","cohort_robust_z_signature","cohort_row_z_centroid","HC_cohort_row_z")
all_scores <- rbindlist(lapply(methods, function(method) rbindlist(lapply(names(data), function(cohort) {
  score <- make_scores(data[[cohort]], method)
  data.table(method, cohort, sample_id=colnames(data[[cohort]]$m), truth=data[[cohort]]$truth, score)
}))))
thresholds <- all_scores[cohort=="TCGA", .(threshold=learn_threshold(score,truth)), by=method]
all_scores <- merge(all_scores, thresholds, by="method")
all_scores[, prediction := ifelse(score > threshold,"LUAD","LUSC")]
metrics <- all_scores[truth %in% c("LUAD","LUSC"), .(
  n=.N, LUAD_n=sum(truth=="LUAD"), LUSC_n=sum(truth=="LUSC"),
  AUC=auc(score,truth), balanced_accuracy=balacc(truth,prediction), accuracy=mean(truth==prediction),
  score_LUAD_median=median(score[truth=="LUAD"]), score_LUSC_median=median(score[truth=="LUSC"]),
  threshold=unique(threshold)
), by=.(method,cohort)]
fwrite(metrics, file.path(out,"results/tables/SupplementaryFigure8_classifier_variant_audit.tsv"), sep="\t")
fwrite(all_scores, file.path(out,"results/tables/SupplementaryFigure8_classifier_variant_scores.tsv.gz"), sep="\t")
print(metrics)

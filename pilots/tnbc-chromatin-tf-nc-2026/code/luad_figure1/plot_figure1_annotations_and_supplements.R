#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: plot_figure1_annotations_and_supplements.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
processed <- file.path(root, "data", "processed")
raw <- file.path(root, "data", "raw")
results <- file.path(root, "results")
tables <- file.path(results, "tables")
figures <- file.path(results, "figures")
annotations <- file.path(processed, "annotations")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  tab <- fread(path, check.names = FALSE)
  ids <- tab[[1]]
  tab[[1]] <- NULL
  answer <- as.matrix(tab)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

unavailable <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x %in% c("NA", "N/A", "Not Reported", "not reported", "Unknown", "-")] <- "Unavailable"
  x
}

site_group <- function(x) {
  x <- tolower(unavailable(x))
  fifelse(x == "unavailable", "Unavailable",
    fifelse(grepl("lymph", x), "Lymph node",
      fifelse(grepl("pleur", x), "Pleura",
        fifelse(grepl("lung|bronch", x), "Lung",
          fifelse(grepl("brain", x), "Brain", "Other")))))
}

growth_group <- function(x) {
  x <- tolower(unavailable(x))
  fifelse(x == "unavailable", "Unavailable",
    fifelse(grepl("adherent", x), "2D adherent",
      fifelse(grepl("suspension", x), "Suspension", "Other")))
}

make_colors <- function(values, preferred = NULL) {
  levels <- sort(unique(as.character(values)))
  base <- setNames(grDevices::hcl.colors(max(3, length(levels)), palette = "Dark 3")[seq_along(levels)], levels)
  if (!is.null(preferred)) {
    hit <- intersect(names(preferred), names(base))
    base[hit] <- preferred[hit]
  }
  if ("Unavailable" %in% names(base)) base[["Unavailable"]] <- "#D9D9D9"
  base
}

preferred_colors <- list(
  Histology = c(LUAD = "#B2182B", LUSC = "#2166AC"),
  Published_subtype = c(
    `LUAD: TRU` = "#1B9E77", `LUAD: PI` = "#D95F02", `LUAD: PP` = "#7570B3",
    `LUSC: Basal` = "#66C2A5", `LUSC: Classical` = "#FC8D62",
    `LUSC: Primitive` = "#8DA0CB", `LUSC: Secretory` = "#E78AC3"
  ),
  Sex = c(Female = "#CC79A7", Male = "#0072B2"),
  Age = c(`<=60` = "#D9F0A3", `61-70` = "#78C679", `>70` = "#238443"),
  Stage = c(I = "#D9F0D3", II = "#A6DBA0", III = "#5AAE61", IV = "#1B7837"),
  Pack_years = c(`0` = "#F7FCF5", `>0-20` = "#C7E9C0", `>20-40` = "#74C476", `>40` = "#238B45"),
  Smoking = c(Current = "#B2182B", `Former >1 year` = "#EF8A62", Never = "#67A9CF", Yes = "#B2182B", No = "#67A9CF"),
  Primary_or_metastasis = c(Primary = "#1B9E77", Metastasis = "#D95F02"),
  Metastatic_disease = c(Yes = "#D95F02", No = "#1B9E77", `Not Reported` = "#D9D9D9"),
  Molecular_record = c(Reported = "#4DAF4A", Unavailable = "#D9D9D9")
)

annotation_for <- function(cohort, sample_ids) {
  if (cohort == "TCGA") {
    x <- fread(file.path(annotations, "tcga_figure1_annotations.tsv"))[match(sample_ids, sample_id)]
    answer <- data.frame(
      Histology = unavailable(x$group),
      Published_subtype = unavailable(x$published_expression_subtype),
      Sex = unavailable(x$sex),
      Age = unavailable(x$age_group),
      Stage = unavailable(x$stage),
      Pack_years = unavailable(x$pack_year_group),
      Oncogene_mutation = unavailable(x$oncogene_mutation),
      Suppressor_mutation = unavailable(x$suppressor_pathway_mutation),
      check.names = FALSE
    )
  } else if (cohort == "PDMR_PDX") {
    x <- fread(file.path(annotations, "pdmr_figure1_annotations.tsv"))[match(sample_ids, sample_id)]
    answer <- data.frame(
      Histology = unavailable(x$group),
      Sex = unavailable(x$sex),
      Age = unavailable(x$age_group),
      Metastatic_disease = unavailable(x$known_metastatic_disease),
      Smoking = unavailable(x$has_smoked_100_cigarettes),
      Biopsy_site = site_group(x$biopsy_site),
      Tissue_type = unavailable(x$tissue_type),
      Molecular_record = unavailable(x$molecular_record),
      check.names = FALSE
    )
  } else if (cohort == "DepMap_22Q2") {
    x <- fread(file.path(annotations, "depmap_figure1_annotations.tsv"))[match(sample_ids, sample_id)]
    answer <- data.frame(
      Histology = unavailable(x$group),
      Sex = unavailable(x$sex),
      Age = unavailable(x$age_group),
      Primary_or_metastasis = unavailable(x$primary_or_metastasis),
      Collection_site = site_group(x$sample_collection_site),
      Growth = growth_group(x$growth_pattern),
      Oncogene_mutation = unavailable(x$oncogene_mutation),
      Suppressor_mutation = unavailable(x$suppressor_pathway_mutation),
      check.names = FALSE
    )
  } else stop(sprintf("Unsupported cohort: %s", cohort))
  rownames(answer) <- sample_ids
  answer
}

selected <- fread(file.path(tables, "tcga_specific_TFs.tsv"))
activity_colors <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B"))
specificity_colors <- c(LUAD = "#B2182B", LUSC = "#2166AC")

plot_annotated_heatmap <- function(matrix_path, cohort, panel, width) {
  matrix <- read_matrix(matrix_path)
  matrix[!is.finite(matrix)] <- NA_real_
  anno_df <- annotation_for(cohort, colnames(matrix))
  anno_colors <- lapply(names(anno_df), function(name) make_colors(anno_df[[name]], preferred_colors[[name]]))
  names(anno_colors) <- names(anno_df)
  top <- HeatmapAnnotation(
    df = anno_df,
    col = anno_colors,
    na_col = "#D9D9D9",
    simple_anno_size = unit(2.5, "mm"),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 7),
    annotation_legend_param = lapply(names(anno_df), function(x) list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
  )
  row_group <- selected$discovery_category[match(rownames(matrix), selected$TF)]
  left <- rowAnnotation(
    Specificity = factor(row_group, levels = c("LUAD", "LUSC")),
    col = list(Specificity = specificity_colors),
    simple_anno_size = unit(3, "mm")
  )
  heatmap <- Heatmap(
    matrix,
    name = "TF activity",
    col = activity_colors,
    cluster_rows = panel != "C",
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = ncol(matrix) > 200,
    raster_quality = 2,
    top_annotation = top,
    left_annotation = left,
    column_title = sprintf("Figure 1%s annotated — %s (%d samples; frozen activity matrix)", panel, cohort, ncol(matrix)),
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  base_name <- sprintf("Figure1%s_%s_TF_activity_annotated", panel, cohort)
  cairo_pdf(file.path(figures, paste0(base_name, ".pdf")), width = width, height = 8.2)
  draw(heatmap, merge_legends = TRUE, padding = unit(c(4, 10, 4, 10), "mm"))
  dev.off()
  png(file.path(figures, paste0(base_name, ".png")), width = round(width * 220), height = round(8.2 * 220), res = 220)
  draw(heatmap, merge_legends = TRUE, padding = unit(c(4, 10, 4, 10), "mm"))
  dev.off()
}

plot_annotated_heatmap(file.path(tables, "Figure1C_TCGA_activity_matrix.tsv.gz"), "TCGA", "C", 12.5)
plot_annotated_heatmap(file.path(tables, "Figure1D_PDMR_PDX_activity_matrix.tsv.gz"), "PDMR_PDX", "D", 10.5)
plot_annotated_heatmap(file.path(tables, "Figure1E_DepMap_22Q2_activity_matrix.tsv.gz"), "DepMap_22Q2", "E", 12.0)

# Supplementary Figure 1: auditable inclusion flow from downloaded expression matrices.
matrix_samples <- function(path) names(fread(path, nrows = 0, check.names = FALSE))[-1]
flow_rows <- list()
for (histology in c("LUAD", "LUSC")) {
  counts_samples <- matrix_samples(file.path(raw, "tcga", sprintf("TCGA-%s.star_counts.tsv.gz", histology)))
  tpm_samples <- matrix_samples(file.path(raw, "tcga", sprintf("TCGA-%s.star_tpm.tsv.gz", histology)))
  common <- intersect(counts_samples, tpm_samples)
  primary <- common[substr(common, 14, 15) == "01"]
  final <- fread(file.path(processed, "tcga_manifest.tsv"))[group == histology, sample_id]
  flow_rows[[histology]] <- data.table(
    cohort = paste0("TCGA-", histology),
    downloaded_count_columns = length(counts_samples),
    downloaded_tpm_columns = length(tpm_samples),
    count_tpm_intersection = length(common),
    primary_tumor_columns = length(primary),
    unique_primary_patients = uniqueN(substr(primary, 1, 12)),
    final_analysis_samples = length(final)
  )
}
gse_all <- matrix_samples(file.path(raw, "gse81089", "GSE81089_readcounts_featurecounts.tsv.gz"))
gse_final <- fread(file.path(processed, "gse81089_manifest.tsv"))
flow_rows[["GSE81089"]] <- data.table(
  cohort = "GSE81089",
  downloaded_count_columns = length(gse_all),
  downloaded_tpm_columns = NA_integer_,
  count_tpm_intersection = NA_integer_,
  primary_tumor_columns = nrow(gse_final),
  unique_primary_patients = uniqueN(gse_final$sample_id),
  final_analysis_samples = nrow(gse_final)
)
flow <- rbindlist(flow_rows, fill = TRUE)
fwrite(flow, file.path(tables, "SupplementaryFigure1_inclusion_counts.tsv"), sep = "\t")

draw_flow_box <- function(x, y, label, fill, width = 0.22, height = 0.15) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"), width = unit(width, "npc"), height = unit(height, "npc"),
    r = unit(0.012, "npc"), gp = gpar(fill = fill, col = "#555555", lwd = 1))
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"), gp = gpar(fontsize = 8.5))
}
draw_arrow <- function(x0, x1, y) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y, y), "npc"), arrow = arrow(length = unit(0.09, "inches")), gp = gpar(col = "#666666"))
}
draw_inclusion_flow <- function(path, png_device = FALSE) {
  if (png_device) png(path, width = 2400, height = 1150, res = 220) else cairo_pdf(path, width = 11, height = 5.2)
  grid.newpage()
  grid.text("Supplementary Figure 1 — cohort inclusion and exclusion audit", x = unit(0.5, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
  for (i in seq_len(nrow(flow))) {
    y <- c(0.72, 0.47, 0.22)[i]
    row <- flow[i]
    draw_flow_box(0.12, y, sprintf("%s\ndownloaded RNA matrices\ncounts n=%d%s", row$cohort, row$downloaded_count_columns,
      ifelse(is.na(row$downloaded_tpm_columns), "", sprintf("; TPM n=%d", row$downloaded_tpm_columns))), "#E8F1FA")
    draw_flow_box(0.39, y, if (row$cohort == "GSE81089") sprintf("Tumor histology metadata\nLUAD/LUSC retained n=%d", row$primary_tumor_columns)
      else sprintf("Count/TPM intersection\nprimary tumor n=%d\nunique patients n=%d", row$primary_tumor_columns, row$unique_primary_patients), "#FFF1D6")
    draw_flow_box(0.66, y, sprintf("Frozen Figure 1 cohort\nn=%d", row$final_analysis_samples), "#E5F5E0")
    excluded <- row$downloaded_count_columns - row$final_analysis_samples
    draw_flow_box(0.89, y, sprintf("Not analyzed\nn=%d\n(normal/non-target/duplicate)", excluded), "#FEE2E2", width = 0.18)
    draw_arrow(0.235, 0.275, y); draw_arrow(0.505, 0.545, y); draw_arrow(0.775, 0.795, y)
  }
  dev.off()
}
draw_inclusion_flow(file.path(figures, "SupplementaryFigure1_inclusion_flow.pdf"), FALSE)
draw_inclusion_flow(file.path(figures, "SupplementaryFigure1_inclusion_flow.png"), TRUE)

# Supplementary Figure 2A: independent-patient expression/activity agreement.
gse <- fread(file.path(tables, "gse81089_limma_msviper.tsv"))
gse_colors <- c(LUAD = "#B2182B", LUSC = "#2166AC", Significant = "#7F7F7F", `Non-Significant` = "#D9D9D9",
                `LogFC Significant` = "#F4A582", `VIPER Significant` = "#92C5DE")
gse_cor <- cor.test(gse$NES, gse$logFC, method = "pearson")
p_a <- ggplot(gse, aes(NES, logFC, color = category)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey75") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey75") +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.45) +
  scale_color_manual(values = gse_colors, na.value = "#D9D9D9") +
  labs(title = "Supplementary Figure 2A — GSE81089", x = "msVIPER NES (LUAD vs LUSC)", y = "log2 fold change",
       subtitle = sprintf("Pearson r = %.2f; p = %.2g", unname(gse_cor$estimate), gse_cor$p.value)) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figures, "SupplementaryFigure2A_GSE81089_logFC_msviper.pdf"), p_a, width = 6.4, height = 5.4, useDingbats = FALSE)
ggsave(file.path(figures, "SupplementaryFigure2A_GSE81089_logFC_msviper.png"), p_a, width = 6.4, height = 5.4, dpi = 220)

# Supplementary Figure 2B: discovery/validation overlap, preserving direction.
tcga <- fread(file.path(tables, "figure1b_tcga_limma_msviper.tsv"))
overlap <- rbindlist(lapply(c("LUAD", "LUSC"), function(direction) {
  a <- tcga[category == direction, TF]
  b <- gse[category == direction, TF]
  data.table(
    direction = direction,
    component = c("TCGA only", "Overlap", "GSE81089 only"),
    n = c(length(setdiff(a, b)), length(intersect(a, b)), length(setdiff(b, a)))
  )
}))
fwrite(overlap, file.path(tables, "SupplementaryFigure2B_TF_overlap_counts.tsv"), sep = "\t")
p_b <- ggplot(overlap, aes(component, n, fill = component)) +
  geom_col(width = 0.72) + geom_text(aes(label = n), vjust = -0.3, size = 3.4) +
  facet_wrap(~direction, scales = "free_y") +
  scale_fill_manual(values = c(`TCGA only` = "#BDBDBD", Overlap = "#4DAF4A", `GSE81089 only` = "#80B1D3")) +
  labs(title = "Supplementary Figure 2B — direction-concordant TF overlap", x = NULL, y = "TF count") +
  theme_classic(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(figures, "SupplementaryFigure2B_TF_overlap.pdf"), p_b, width = 7.0, height = 4.4, useDingbats = FALSE)
ggsave(file.path(figures, "SupplementaryFigure2B_TF_overlap.png"), p_b, width = 7.0, height = 4.4, dpi = 220)

plot_sample_correlations <- function(matrix_path, cohort, panel, width = 8.2) {
  matrix <- read_matrix(matrix_path)
  corr <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
  corr[!is.finite(corr)] <- 0
  distance_matrix <- 1 - corr
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  distance <- as.dist(distance_matrix)
  tree <- hclust(distance, method = "average")
  anno <- annotation_for(cohort, colnames(matrix))[, "Histology", drop = FALSE]
  hist_colors <- list(Histology = preferred_colors$Histology)
  top <- HeatmapAnnotation(df = anno, col = hist_colors, simple_anno_size = unit(3, "mm"))
  left <- rowAnnotation(df = anno, col = hist_colors, simple_anno_size = unit(3, "mm"))
  heatmap <- Heatmap(
    corr, name = "Pearson r", col = colorRamp2(c(-0.25, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = tree, cluster_columns = tree, show_row_names = FALSE, show_column_names = FALSE,
    use_raster = ncol(corr) > 150, raster_quality = 2, top_annotation = top, left_annotation = left,
    column_title = sprintf("Supplementary Figure 2%s — %s sample correlations", panel, cohort)
  )
  name <- sprintf("SupplementaryFigure2%s_%s_sample_correlations", panel, cohort)
  cairo_pdf(file.path(figures, paste0(name, ".pdf")), width = width, height = 7.6)
  draw(heatmap, merge_legends = TRUE)
  dev.off()
  png(file.path(figures, paste0(name, ".png")), width = round(width * 200), height = 1520, res = 200)
  draw(heatmap, merge_legends = TRUE)
  dev.off()
}
plot_sample_correlations(file.path(tables, "Figure1C_TCGA_activity_matrix.tsv.gz"), "TCGA", "C", 9.2)
plot_sample_correlations(file.path(tables, "Figure1D_PDMR_PDX_activity_matrix.tsv.gz"), "PDMR_PDX", "D", 7.8)
plot_sample_correlations(file.path(tables, "Figure1E_DepMap_22Q2_activity_matrix.tsv.gz"), "DepMap_22Q2", "E", 8.2)

# Supplementary Figure 2F: effect reproducibility of the 97 frozen LUAD TFs.
replicated <- fread(file.path(tables, "externally_replicated_TFs.tsv"))[discovery_category == "LUAD", TF]
effects <- fread(file.path(tables, "cross_system_TF_effects.tsv.gz"))[TF %in% replicated]
required_effect_cohorts <- c("TCGA", "GSE81089", "PDMR_PDX", "DepMap_22Q2")
missing_effect_cohorts <- setdiff(required_effect_cohorts, unique(effects$cohort))
if (length(missing_effect_cohorts)) {
  stop(sprintf("Missing Figure 2F effect cohorts: %s", paste(missing_effect_cohorts, collapse = ", ")))
}
effects[, cohort := factor(cohort, levels = required_effect_cohorts)]
effect_matrix <- data.table::dcast(
  effects,
  TF ~ cohort,
  value.var = "mean_difference",
  fun.aggregate = mean
)
effect_ids <- effect_matrix$TF
effect_matrix[, TF := NULL]
effect_raw <- as.matrix(effect_matrix[, ..required_effect_cohorts])
rownames(effect_raw) <- effect_ids
effect_z <- vapply(
  seq_len(ncol(effect_raw)),
  function(j) as.numeric(scale(effect_raw[, j])),
  numeric(nrow(effect_raw))
)
rownames(effect_z) <- effect_ids
colnames(effect_z) <- required_effect_cohorts
if (!identical(dim(effect_z), c(length(replicated), length(required_effect_cohorts)))) {
  stop(sprintf(
    "Figure 2F matrix dimension mismatch: observed %s, expected %d x %d",
    paste(dim(effect_z), collapse = " x "), length(replicated), length(required_effect_cohorts)
  ))
}
effect_z[!is.finite(effect_z)] <- 0
effect_z[effect_z > 2.5] <- 2.5
effect_z[effect_z < -2.5] <- -2.5
effect_z_table <- cbind(data.table(TF = rownames(effect_z)), as.data.table(effect_z))
if (!identical(names(effect_z_table), c("TF", required_effect_cohorts))) {
  stop(sprintf("Figure 2F output columns are invalid: %s", paste(names(effect_z_table), collapse = ", ")))
}
fwrite(effect_z_table, file.path(tables, "SupplementaryFigure2F_replicated_LUAD_TF_effect_zscores.tsv"), sep = "\t")
heatmap_f <- Heatmap(
  effect_z, name = "cohort-scaled\neffect", col = activity_colors,
  cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE,
  column_title = sprintf("Supplementary Figure 2F — %d externally replicated LUAD TFs", length(replicated))
)
cairo_pdf(file.path(figures, "SupplementaryFigure2F_replicated_LUAD_TF_effects.pdf"), width = 6.2, height = 7.5)
draw(heatmap_f)
dev.off()
png(file.path(figures, "SupplementaryFigure2F_replicated_LUAD_TF_effects.png"), width = 1360, height = 1650, res = 220)
draw(heatmap_f)
dev.off()

core_index <- data.table(
  panel = c("S1", "S2A", "S2B", "S2C", "S2D", "S2E", "S2F"),
  content = c(
    "Cohort inclusion/exclusion flow", "Independent GSE81089 expression/activity concordance",
    "Direction-concordant TF overlap", "TCGA sample activity correlations",
    "PDMR PDX sample activity correlations", "DepMap sample activity correlations",
    "Cross-system effects of the frozen 97 LUAD TFs"
  ),
  status = "generated"
)
fwrite(core_index, file.path(tables, "SupplementaryFigures1_2_core_index.tsv"), sep = "\t")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: select_and_plot_figure1.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
processed <- file.path(root, "data", "processed")
results <- file.path(root, "results")
figures <- file.path(results, "figures")
tables <- file.path(results, "tables")
audit <- file.path(root, "audit")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  tab <- fread(path, check.names = FALSE)
  genes <- tab[[1]]
  tab[[1]] <- NULL
  matrix <- as.matrix(tab)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  matrix
}

run_limma_voom <- function(counts, manifest) {
  manifest <- manifest[match(colnames(counts), sample_id)]
  group <- factor(manifest$group, levels = c("LUSC", "LUAD"))
  design <- model.matrix(~ group)
  # The anchor calls voom(counts) without passing the group design, then uses
  # the design only in lmFit. Preserve that executable-code behavior.
  fit <- eBayes(lmFit(voom(round(counts), plot = FALSE), design))
  tab <- topTable(fit, coef = "groupLUAD", number = Inf, adjust.method = "fdr", sort.by = "none")
  tab$TF <- rownames(tab)
  as.data.table(tab)
}

author_discovery_category <- function(tab) {
  tab[, category := ifelse(p_value > 0.01 & logFC < 0.5 & logFC > -0.5, "Non-Significant", "Significant")]
  tab[category == "Significant" & logFC < 0.5 & logFC > -0.5, category := "LogFC Non-Significant"]
  nonsig <- tab[category == "Non-Significant", NES]
  if (length(nonsig) == 0) stop("Anchor category rule has no non-significant NES reference interval")
  tab[category == "Significant" & NES < max(nonsig) & NES > min(nonsig), category := "VIPER Significant"]
  tab[category == "Significant" & NES > 0, category := "LUAD"]
  tab[category == "Significant" & NES < 0, category := "LUSC"]
  tab
}

validation_category <- function(tab) {
  tab[, category := ifelse(p_value <= 0.05 & adj.P.Val <= 0.05, "Significant", "Non-Significant")]
  tab[category == "Non-Significant" & p_value > 0.05 & adj.P.Val <= 0.05, category := "LogFC Significant"]
  tab[category == "Non-Significant" & p_value <= 0.05 & adj.P.Val > 0.05, category := "VIPER Significant"]
  tab[category == "Significant" & NES > 0 & logFC > 0, category := "LUAD"]
  tab[category == "Significant" & NES < 0 & logFC < 0, category := "LUSC"]
  tab
}

tcga_counts <- read_matrix(file.path(processed, "tcga_counts_symbols.tsv.gz"))
tcga_manifest <- fread(file.path(processed, "tcga_manifest.tsv"))
tcga_limma <- run_limma_voom(tcga_counts, tcga_manifest)
tcga_ms <- fread(file.path(results, "tcga", "TCGA_msviper.tsv"))
tcga <- merge(tcga_ms, tcga_limma, by = "TF")
tcga <- author_discovery_category(tcga)
fwrite(tcga, file.path(tables, "figure1b_tcga_limma_msviper.tsv"), sep = "\t")

gse_counts <- read_matrix(file.path(processed, "gse81089_counts_symbols.tsv.gz"))
gse_manifest <- fread(file.path(processed, "gse81089_manifest.tsv"))
gse_limma <- run_limma_voom(gse_counts, gse_manifest)
gse_ms <- fread(file.path(results, "gse81089", "GSE81089_msviper.tsv"))
gse <- validation_category(merge(gse_ms, gse_limma, by = "TF"))
fwrite(gse, file.path(tables, "gse81089_limma_msviper.tsv"), sep = "\t")

tcga_luad <- tcga[category == "LUAD", TF]
tcga_lusc <- tcga[category == "LUSC", TF]
gse_luad <- gse[category == "LUAD", TF]
gse_lusc <- gse[category == "LUSC", TF]
replicated_luad <- intersect(tcga_luad, gse_luad)
replicated_lusc <- intersect(tcga_lusc, gse_lusc)
selected <- rbind(
  data.table(TF = tcga_luad, discovery_category = "LUAD", externally_replicated = tcga_luad %in% replicated_luad),
  data.table(TF = tcga_lusc, discovery_category = "LUSC", externally_replicated = tcga_lusc %in% replicated_lusc)
)
fwrite(selected, file.path(tables, "tcga_specific_TFs.tsv"), sep = "\t")
fwrite(
  selected[externally_replicated == TRUE],
  file.path(tables, "externally_replicated_TFs.tsv"), sep = "\t"
)

plot_workflow <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    cairo_pdf(path, width = 11.2, height = 4.2)
  } else {
    png(path, width = 2460, height = 920, res = 220)
  }
  grid.newpage()
  box <- function(x, y, width, height, label, fill, fontsize = 9.5) {
    grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
      width = unit(width, "npc"), height = unit(height, "npc"),
      r = unit(0.018, "npc"), gp = gpar(fill = fill, col = "#4D4D4D", lwd = 1.1))
    grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
      gp = gpar(fontsize = fontsize, col = "#202020"))
  }
  arrow <- function(x0, y0, x1, y1) {
    grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
      arrow = grid::arrow(length = unit(0.12, "inches")), gp = gpar(col = "#555555", lwd = 1.3))
  }
  grid.text("A", x = unit(0.015, "npc"), y = unit(0.96, "npc"), just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold"))
  grid.text("LUAD histology-specific transcriptional-regulator discovery and validation",
    x = unit(0.50, "npc"), y = unit(0.91, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
  box(0.11, 0.61, 0.18, 0.24,
    sprintf("TCGA discovery\n%d LUAD / %d LUSC\nprimary tumors; TPM",
      sum(tcga_manifest$group == "LUAD"), sum(tcga_manifest$group == "LUSC")), "#E8F1FA")
  box(0.34, 0.61, 0.18, 0.24, "ARACNe3\nPAN-GO TF regulons\none anchor-default subnet", "#F3E8FA")
  box(0.57, 0.61, 0.18, 0.24, "VIPER / msVIPER\n1,000-permutation null\nlimma expression effect", "#FFF1D6")
  box(0.80, 0.61, 0.18, 0.24,
    sprintf("TCGA-specific TFs\n%d LUAD / %d LUSC",
      length(tcga_luad), length(tcga_lusc)), "#FBE2E2")
  arrow(0.20, 0.61, 0.25, 0.61); arrow(0.43, 0.61, 0.48, 0.61); arrow(0.66, 0.61, 0.71, 0.61)
  box(0.22, 0.23, 0.25, 0.22,
    sprintf("Independent patient validation\nGSE81089: %d LUAD / %d LUSC\nown ARACNe3 network",
      sum(gse_manifest$group == "LUAD"), sum(gse_manifest$group == "LUSC")), "#E5F5E0")
  pdmr_manifest_workflow <- fread(file.path(processed, "pdmr_manifest.tsv"))
  depmap_manifest_workflow <- fread(file.path(processed, "depmap_manifest.tsv"))
  box(0.52, 0.23, 0.25, 0.22,
    sprintf("PDX projection\nPDMR: %d LUAD / %d LUSC\nTCGA regulon fixed",
      sum(pdmr_manifest_workflow$group == "LUAD"), sum(pdmr_manifest_workflow$group == "LUSC")), "#E0F2F1")
  box(0.82, 0.23, 0.25, 0.22,
    sprintf("Cell-line projection\nDepMap 22Q2: %d LUAD / %d LUSC\nTCGA regulon fixed",
      sum(depmap_manifest_workflow$group == "LUAD"), sum(depmap_manifest_workflow$group == "LUSC")), "#FDE9D9")
  arrow(0.77, 0.49, 0.22, 0.35); arrow(0.80, 0.49, 0.52, 0.35); arrow(0.83, 0.49, 0.82, 0.35)
  dev.off()
}

plot_workflow(file.path(figures, "Figure1A_workflow.pdf"), "pdf")
plot_workflow(file.path(figures, "Figure1A_workflow.png"), "png")

scatter_colors <- c(
  LUAD = "#B2182B", LUSC = "#2166AC", `Non-Significant` = "#D9D9D9",
  `LogFC Non-Significant` = "#F4A582", `VIPER Significant` = "#92C5DE"
)
label_rows <- tcga[category %in% c("LUAD", "LUSC")][order(-abs(NES))][1:min(.N, 40)]
tcga_correlation <- cor.test(tcga$NES, tcga$logFC, method = "pearson")
correlation_p_label <- if (tcga_correlation$p.value < .Machine$double.eps) {
  "p < 2.2e-16"
} else {
  sprintf("p = %.2g", tcga_correlation$p.value)
}
scatter <- ggplot(tcga, aes(NES, logFC, color = category)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey75") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey75") +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.45) +
  geom_text_repel(data = label_rows, aes(label = TF), size = 2.4, max.overlaps = 40) +
  scale_color_manual(values = scatter_colors, drop = FALSE) +
  labs(
    x = "msVIPER NES (LUAD vs LUSC)", y = "log2 fold change (LUAD vs LUSC)",
    subtitle = sprintf("Pearson r = %.2f; %s", unname(tcga_correlation$estimate), correlation_p_label)
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
ggsave(file.path(figures, "Figure1B_TCGA_logFC_msviper.pdf"), scatter, width = 6.4, height = 5.4, useDingbats = FALSE)
ggsave(file.path(figures, "Figure1B_TCGA_logFC_msviper.png"), scatter, width = 6.4, height = 5.4, dpi = 220)

activity_colors <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B"))
group_colors <- c(LUAD = "#B2182B", LUSC = "#2166AC")
row_colors <- c(LUAD = "#B2182B", LUSC = "#2166AC")
cluster_statistics <- list()

plot_heatmap <- function(activity_path, manifest_path, cohort, panel) {
  activity <- read_matrix(activity_path)
  manifest <- fread(manifest_path)
  manifest <- manifest[match(colnames(activity), sample_id)]
  present <- selected$TF[selected$TF %in% rownames(activity)]
  matrix <- activity[present, , drop = FALSE]
  keep <- apply(matrix, 1, sd, na.rm = TRUE) > 0
  matrix <- matrix[keep, , drop = FALSE]
  # Anchor Figure 1C plots TCGA VIPER activity directly with fixed row order;
  # Figure 1D/E row-standardize model activity and cluster rows.
  if (panel != "C") matrix <- t(scale(t(matrix)))
  row_group <- selected$discovery_category[match(rownames(matrix), selected$TF)]
  column_group <- factor(manifest$group, levels = c("LUAD", "LUSC"))
  column_tree <- hclust(dist(t(matrix)))
  two_clusters <- factor(cutree(column_tree, k = 2))
  contingency <- table(two_clusters, column_group)
  fisher <- fisher.test(contingency)
  orientations <- c(
    contingency[1, "LUAD"] + contingency[2, "LUSC"],
    contingency[1, "LUSC"] + contingency[2, "LUAD"]
  )
  cluster_statistics[[cohort]] <<- data.table(
    cohort = cohort, samples = ncol(matrix), TFs = nrow(matrix),
    fisher_p_value = fisher$p.value,
    best_two_cluster_accuracy = max(orientations) / sum(contingency),
    contingency = paste(capture.output(print(contingency)), collapse = " | ")
  )
  annotation <- HeatmapAnnotation(
    Histology = column_group,
    col = list(Histology = group_colors),
    annotation_name_side = "left"
  )
  row_annotation <- rowAnnotation(
    Specificity = factor(row_group, levels = c("LUAD", "LUSC")),
    col = list(Specificity = row_colors)
  )
  heatmap <- Heatmap(
    matrix,
    name = "TF activity",
    col = activity_colors,
    cluster_rows = panel != "C",
    cluster_columns = column_tree,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = annotation,
    left_annotation = row_annotation,
    column_title = sprintf("%s: %d LUAD / %d LUSC", cohort, sum(column_group == "LUAD"), sum(column_group == "LUSC"))
  )
  pdf_path <- file.path(figures, sprintf("Figure1%s_%s_TF_activity.pdf", panel, cohort))
  png_path <- file.path(figures, sprintf("Figure1%s_%s_TF_activity.png", panel, cohort))
  device_width <- if (panel == "C") 8.6 else if (panel == "E") 12.0 else 10.0
  horizontal_padding <- if (panel == "E") 18 else 12
  cairo_pdf(pdf_path, width = device_width, height = 6.6)
  draw(heatmap, padding = unit(c(5, horizontal_padding, 5, horizontal_padding), "mm"))
  dev.off()
  png(png_path, width = round(device_width * 220), height = 1450, res = 220)
  draw(heatmap, padding = unit(c(5, horizontal_padding, 5, horizontal_padding), "mm"))
  dev.off()
  matrix_TFs <- rownames(matrix)
  output_matrix <- as.data.table(matrix)
  output_matrix[, TF := matrix_TFs]
  setcolorder(output_matrix, "TF")
  fwrite(output_matrix, file.path(tables, sprintf("Figure1%s_%s_activity_matrix.tsv.gz", panel, cohort)), sep = "\t")
}

plot_heatmap(
  file.path(results, "tcga", "TCGA_viper_activity.tsv.gz"),
  file.path(processed, "tcga_manifest.tsv"), "TCGA", "C"
)
plot_heatmap(
  file.path(results, "pdmr", "PDMR_viper_activity.tsv.gz"),
  file.path(processed, "pdmr_manifest.tsv"), "PDMR_PDX", "D"
)
plot_heatmap(
  file.path(results, "depmap", "DEPMAP_viper_activity.tsv.gz"),
  file.path(processed, "depmap_manifest.tsv"), "DepMap_22Q2", "E"
)
cluster_statistics <- rbindlist(cluster_statistics)
cluster_statistics[, FDR := p.adjust(fisher_p_value, method = "BH")]
fwrite(cluster_statistics, file.path(tables, "unsupervised_cluster_histology_statistics.tsv"), sep = "\t")

effect_by_system <- function(activity_path, manifest_path, cohort) {
  activity <- read_matrix(activity_path)
  manifest <- fread(manifest_path)[match(colnames(activity), sample_id)]
  common <- intersect(selected$TF, rownames(activity))
  rbindlist(lapply(common, function(tf) {
    a <- activity[tf, manifest$group == "LUAD"]
    b <- activity[tf, manifest$group == "LUSC"]
    test <- tryCatch(t.test(a, b), error = function(error) NULL)
    data.table(
      cohort = cohort, TF = tf, mean_LUAD = mean(a), mean_LUSC = mean(b),
      mean_difference = mean(a) - mean(b), p_value = if (is.null(test)) 1 else test$p.value,
      expected_sign = ifelse(selected$discovery_category[match(tf, selected$TF)] == "LUAD", 1, -1),
      sign_concordant = sign(mean(a) - mean(b)) == ifelse(selected$discovery_category[match(tf, selected$TF)] == "LUAD", 1, -1)
    )
  }))
}

effects <- rbind(
  effect_by_system(file.path(results, "tcga", "TCGA_viper_activity.tsv.gz"), file.path(processed, "tcga_manifest.tsv"), "TCGA"),
  effect_by_system(file.path(results, "gse81089", "GSE81089_viper_activity.tsv.gz"), file.path(processed, "gse81089_manifest.tsv"), "GSE81089"),
  effect_by_system(file.path(results, "pdmr", "PDMR_viper_activity.tsv.gz"), file.path(processed, "pdmr_manifest.tsv"), "PDMR_PDX"),
  effect_by_system(file.path(results, "depmap", "DEPMAP_viper_activity.tsv.gz"), file.path(processed, "depmap_manifest.tsv"), "DepMap_22Q2")
)
effects[, FDR := p.adjust(p_value, method = "BH"), by = cohort]
fwrite(effects, file.path(tables, "cross_system_TF_effects.tsv.gz"), sep = "\t")

program_by_system <- function(activity_path, manifest_path, cohort) {
  activity <- read_matrix(activity_path)
  manifest <- fread(manifest_path)[match(colnames(activity), sample_id)]
  program_tfs <- intersect(replicated_luad, rownames(activity))
  if (length(program_tfs) == 0) {
    return(list(summary = data.table(
      cohort = cohort, TFs = 0L, n_LUAD = sum(manifest$group == "LUAD"),
      n_LUSC = sum(manifest$group == "LUSC"), mean_LUAD = NA_real_, mean_LUSC = NA_real_,
      mean_difference = NA_real_, hedges_g = NA_real_, p_value = 1
    ), scores = data.table()))
  }
  matrix <- activity[program_tfs, , drop = FALSE]
  matrix <- t(scale(t(matrix)))
  matrix[!is.finite(matrix)] <- 0
  score <- colMeans(matrix)
  luad_score <- score[manifest$group == "LUAD"]
  lusc_score <- score[manifest$group == "LUSC"]
  test <- t.test(luad_score, lusc_score)
  pooled_sd <- sqrt(((length(luad_score) - 1) * var(luad_score) +
    (length(lusc_score) - 1) * var(lusc_score)) /
    (length(luad_score) + length(lusc_score) - 2))
  cohen_d <- (mean(luad_score) - mean(lusc_score)) / pooled_sd
  correction <- 1 - 3 / (4 * (length(luad_score) + length(lusc_score)) - 9)
  score_table <- data.table(
    cohort = cohort, sample_id = manifest$sample_id, group = manifest$group,
    LUAD_program_score = as.numeric(score)
  )
  summary <- data.table(
    cohort = cohort, TFs = length(program_tfs), n_LUAD = length(luad_score), n_LUSC = length(lusc_score),
    mean_LUAD = mean(luad_score), mean_LUSC = mean(lusc_score),
    mean_difference = mean(luad_score) - mean(lusc_score),
    hedges_g = cohen_d * correction, p_value = test$p.value
  )
  list(summary = summary, scores = score_table)
}

program_results <- list(
  program_by_system(file.path(results, "gse81089", "GSE81089_viper_activity.tsv.gz"), file.path(processed, "gse81089_manifest.tsv"), "GSE81089"),
  program_by_system(file.path(results, "pdmr", "PDMR_viper_activity.tsv.gz"), file.path(processed, "pdmr_manifest.tsv"), "PDMR_PDX"),
  program_by_system(file.path(results, "depmap", "DEPMAP_viper_activity.tsv.gz"), file.path(processed, "depmap_manifest.tsv"), "DepMap_22Q2")
)
program_summary <- rbindlist(lapply(program_results, `[[`, "summary"))
program_summary[, FDR := p.adjust(p_value, method = "BH")]
program_scores <- rbindlist(lapply(program_results, `[[`, "scores"), fill = TRUE)
fwrite(program_summary, file.path(tables, "LUAD_program_cross_system_statistics.tsv"), sep = "\t")
fwrite(program_scores, file.path(tables, "LUAD_program_sample_scores.tsv.gz"), sep = "\t")

replicated <- selected[externally_replicated == TRUE]
cross_summary <- effects[TF %in% replicated$TF, .(
  TFs = uniqueN(TF),
  sign_concordant = sum(sign_concordant),
  FDR_0_05_and_sign = sum(FDR <= 0.05 & sign_concordant),
  sign_concordance_fraction = mean(sign_concordant)
), by = cohort]
fwrite(cross_summary, file.path(tables, "cross_system_summary.tsv"), sep = "\t")

external_models <- program_summary[cohort %in% c("PDMR_PDX", "DepMap_22Q2")]
biological_gate <- length(replicated_luad) >= 3 && nrow(external_models) == 2 &&
  all(external_models$mean_difference > 0) && all(external_models$FDR <= 0.05) &&
  all(external_models$hedges_g >= 0.5)
receipt <- list(
  status = if (biological_gate) "passed" else "failed",
  comparison = "LUAD_vs_LUSC",
  discovery = list(
    LUAD_specific_TFs = length(tcga_luad),
    LUSC_specific_TFs = length(tcga_lusc),
    NES_expression_Pearson_r = unname(tcga_correlation$estimate),
    NES_expression_Pearson_p = tcga_correlation$p.value,
    author_rule = "verbatim logical structure from Figure1/08-Prepare_Annotation_ForPlotting.R"
  ),
  independent_patient_replication = list(
    cohort = "GSE81089",
    LUAD_specific_replicated = length(replicated_luad),
    LUSC_specific_replicated = length(replicated_lusc),
    rule = "msVIPER raw p<=0.05; limma BH<=0.05; direction concordant"
  ),
  cross_system = cross_summary,
  unsupervised_clustering = cluster_statistics,
  LUAD_program_statistics = program_summary,
  stop_rule = list(
    minimum_independent_patient_replicated_LUAD_TFs = 3,
    external_models_required = c("PDMR_PDX", "DepMap_22Q2"),
    model_requirements = "positive mean difference, BH FDR<=0.05, Hedges g>=0.5 in both models",
    biological_gate_passed = biological_gate
  ),
  method_execution = list(
    ttestNull = "viper::ttestNull",
    permutations = 1000,
    replacement = TRUE,
    seed = 1,
    cores = 32,
    rng_kind = "L'Ecuyer-CMRG",
    adaptation_timing = "frozen after serial runtime measurement at 2 percent and before any msVIPER/TF result",
    scientific_definition_changed = FALSE,
    viper_1_38_aREA_single_target_dimension_patch = TRUE,
    compatibility_patch_effect = "shape restoration only; minsize=1 and numeric definition unchanged"
  ),
  output_figures = list.files(figures, full.names = FALSE),
  package_versions = list(
    R = as.character(getRversion()), limma = as.character(packageVersion("limma")),
    ComplexHeatmap = as.character(packageVersion("ComplexHeatmap"))
  )
)
write_json(receipt, file.path(audit, "figure1_result_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

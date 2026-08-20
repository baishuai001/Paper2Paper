#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(limma)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: analyze_gse41271_microarray.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
processed <- file.path(root, "data", "processed")
results <- file.path(root, "results")
tables <- file.path(results, "tables")
figures <- file.path(results, "figures")
audit <- file.path(root, "audit", "gse41271")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  tab <- fread(path, check.names = FALSE)
  genes <- tab[[1L]]
  tab[[1L]] <- NULL
  matrix <- as.matrix(tab)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  matrix
}

validation_category <- function(tab) {
  tab[, category := ifelse(p_value <= 0.05 & adj.P.Val <= 0.05,
                            "Significant", "Non-Significant")]
  tab[category == "Non-Significant" & p_value > 0.05 & adj.P.Val <= 0.05,
      category := "LogFC Significant"]
  tab[category == "Non-Significant" & p_value <= 0.05 & adj.P.Val > 0.05,
      category := "VIPER Significant"]
  tab[category == "Significant" & NES > 0 & logFC > 0, category := "LUAD"]
  tab[category == "Significant" & NES < 0 & logFC < 0, category := "LUSC"]
  tab
}

expression <- read_matrix(file.path(processed, "gse41271_aracne_expression.tsv"))
manifest <- fread(file.path(processed, "gse41271_manifest.tsv"))
manifest <- manifest[match(colnames(expression), geo_accession)]
if (anyNA(manifest$geo_accession)) stop("GSE41271 manifest/expression mismatch")
group <- factor(manifest$group, levels = c("LUSC", "LUAD"))
design <- model.matrix(~ group)
fit <- eBayes(lmFit(expression, design))
limma_table <- topTable(fit, coef = "groupLUAD", number = Inf,
                        adjust.method = "fdr", sort.by = "none")
limma_table$TF <- rownames(limma_table)
limma_table <- as.data.table(limma_table)

ms <- fread(file.path(results, "gse41271", "GSE41271_msviper.tsv"))
gse41271 <- validation_category(merge(ms, limma_table, by = "TF"))
fwrite(gse41271, file.path(tables, "gse41271_limma_msviper.tsv"), sep = "\t")

tcga <- fread(file.path(tables, "figure1b_tcga_limma_msviper.tsv"))
gse81089 <- fread(file.path(tables, "gse81089_limma_msviper.tsv"))
specific_sets <- list(
  TCGA_LUAD = tcga[category == "LUAD", TF],
  TCGA_LUSC = tcga[category == "LUSC", TF],
  GSE81089_LUAD = gse81089[category == "LUAD", TF],
  GSE81089_LUSC = gse81089[category == "LUSC", TF],
  GSE41271_LUAD = gse41271[category == "LUAD", TF],
  GSE41271_LUSC = gse41271[category == "LUSC", TF]
)
replicated <- rbind(
  data.table(
    TF = specific_sets$TCGA_LUAD,
    discovery_category = "LUAD",
    replicated_GSE81089 = specific_sets$TCGA_LUAD %in% specific_sets$GSE81089_LUAD,
    replicated_GSE41271 = specific_sets$TCGA_LUAD %in% specific_sets$GSE41271_LUAD
  ),
  data.table(
    TF = specific_sets$TCGA_LUSC,
    discovery_category = "LUSC",
    replicated_GSE81089 = specific_sets$TCGA_LUSC %in% specific_sets$GSE81089_LUSC,
    replicated_GSE41271 = specific_sets$TCGA_LUSC %in% specific_sets$GSE41271_LUSC
  )
)
replicated[, replicated_both := replicated_GSE81089 & replicated_GSE41271]
fwrite(replicated, file.path(tables, "figure1_cross_platform_replication.tsv"), sep = "\t")
fwrite(replicated[replicated_GSE41271 == TRUE],
       file.path(tables, "gse41271_replicated_TFs.tsv"), sep = "\t")

combination_table <- rbindlist(lapply(c("LUAD", "LUSC"), function(label) {
  universe <- sort(unique(c(
    specific_sets[[paste0("TCGA_", label)]],
    specific_sets[[paste0("GSE81089_", label)]],
    specific_sets[[paste0("GSE41271_", label)]]
  )))
  data.table(
    specificity = label,
    TF = universe,
    TCGA = universe %in% specific_sets[[paste0("TCGA_", label)]],
    GSE81089 = universe %in% specific_sets[[paste0("GSE81089_", label)]],
    GSE41271 = universe %in% specific_sets[[paste0("GSE41271_", label)]]
  )
}))
combination_table[, combination := paste0(
  ifelse(TCGA, "TCGA", ""),
  ifelse(TCGA & GSE81089, "+", ""),
  ifelse(GSE81089, "GSE81089", ""),
  ifelse((TCGA | GSE81089) & GSE41271, "+", ""),
  ifelse(GSE41271, "GSE41271", "")
)]
combination_counts <- combination_table[, .(TFs = .N), by = .(specificity, combination)]
fwrite(combination_table, file.path(tables, "figure1_patient_cohort_TF_membership.tsv"), sep = "\t")
fwrite(combination_counts, file.path(tables, "figure1_patient_cohort_TF_overlap_counts.tsv"), sep = "\t")

scatter_colors <- c(
  LUAD = "#B2182B", LUSC = "#2166AC", `Non-Significant` = "#D9D9D9",
  `LogFC Significant` = "#F4A582", `VIPER Significant` = "#92C5DE",
  Significant = "#969696"
)
correlation <- cor.test(gse41271$NES, gse41271$logFC, method = "pearson")
format_p <- function(value) {
  if (!is.finite(value) || value < .Machine$double.eps) return("p < 2.2e-16")
  sprintf("p = %.3g", value)
}
scatter <- ggplot(gse41271, aes(NES, logFC, color = category)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey75") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey75") +
  geom_point(size = 1.05, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.45) +
  scale_color_manual(values = scatter_colors, drop = FALSE) +
  labs(
    title = "GSE41271 microarray validation",
    subtitle = sprintf("Pearson r = %.2f; %s", unname(correlation$estimate), format_p(correlation$p.value)),
    x = "msVIPER NES (LUAD vs LUSC)",
    y = "log2 expression fold change (LUAD vs LUSC)",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures, "SupplementaryFigure2G_GSE41271_logFC_msviper.pdf"),
       scatter, width = 6.4, height = 5.4, useDingbats = FALSE)
ggsave(file.path(figures, "SupplementaryFigure2G_GSE41271_logFC_msviper.png"),
       scatter, width = 6.4, height = 5.4, dpi = 220)

combination_counts[, combination := factor(
  combination,
  levels = c("TCGA+GSE81089+GSE41271", "TCGA+GSE41271", "TCGA+GSE81089",
             "GSE81089+GSE41271", "TCGA", "GSE81089", "GSE41271")
)]
overlap_plot <- ggplot(combination_counts[!is.na(combination)],
                       aes(combination, TFs, fill = specificity)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = TFs), position = position_dodge(width = 0.78),
            vjust = -0.25, size = 3) +
  scale_fill_manual(values = c(LUAD = "#B2182B", LUSC = "#2166AC")) +
  labs(title = "Direction-concordant TF replication across patient cohorts",
       x = NULL, y = "TF count", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "top")
ggsave(file.path(figures, "SupplementaryFigure2H_patient_cohort_TF_overlap.pdf"),
       overlap_plot, width = 8.2, height = 5.0, useDingbats = FALSE)
ggsave(file.path(figures, "SupplementaryFigure2H_patient_cohort_TF_overlap.png"),
       overlap_plot, width = 8.2, height = 5.0, dpi = 220)

summary <- data.table(
  specificity = c("LUAD", "LUSC"),
  TCGA_specific = c(length(specific_sets$TCGA_LUAD), length(specific_sets$TCGA_LUSC)),
  GSE81089_specific = c(length(specific_sets$GSE81089_LUAD), length(specific_sets$GSE81089_LUSC)),
  GSE41271_specific = c(length(specific_sets$GSE41271_LUAD), length(specific_sets$GSE41271_LUSC)),
  TCGA_replicated_GSE81089 = c(
    sum(replicated$discovery_category == "LUAD" & replicated$replicated_GSE81089),
    sum(replicated$discovery_category == "LUSC" & replicated$replicated_GSE81089)
  ),
  TCGA_replicated_GSE41271 = c(
    sum(replicated$discovery_category == "LUAD" & replicated$replicated_GSE41271),
    sum(replicated$discovery_category == "LUSC" & replicated$replicated_GSE41271)
  ),
  TCGA_replicated_both = c(
    sum(replicated$discovery_category == "LUAD" & replicated$replicated_both),
    sum(replicated$discovery_category == "LUSC" & replicated$replicated_both)
  )
)
fwrite(summary, file.path(tables, "figure1_cross_platform_replication_summary.tsv"), sep = "\t")

receipt <- list(
  status = "passed",
  accession = "GSE41271",
  samples = nrow(manifest),
  group_counts = as.list(table(manifest$group)),
  genes = nrow(expression),
  msviper_TFs = nrow(ms),
  expression_activity_Pearson_r = unname(correlation$estimate),
  expression_activity_Pearson_p = correlation$p.value,
  replication_summary = summary,
  validation_rule = "msVIPER raw p<=0.05; limma BH<=0.05; NES and logFC direction concordant",
  downstream_figure2_input_changed = FALSE,
  downstream_reason = "Figure 2 remains defined by all 158 TCGA-discovered LUAD TFs",
  package_versions = list(R = as.character(getRversion()), limma = as.character(packageVersion("limma")))
)
write_json(receipt, file.path(audit, "gse41271_analysis_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

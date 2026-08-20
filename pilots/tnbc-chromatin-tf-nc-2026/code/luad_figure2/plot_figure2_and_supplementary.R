#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: plot_figure2_and_supplementary.R <figure2_run_root> <project_root>")
run_root <- normalizePath(args[[1]], mustWork = TRUE)
project_root <- normalizePath(args[[2]], mustWork = TRUE)
figure_dir <- file.path(run_root, "results", "figures")
table_dir <- file.path(run_root, "results", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

primary_definition <- "PRIMARY_anchor_promoter_tcga_cpm1_both"
system_levels <- c("patient", "PDX", "cell_line")
system_labels <- c(patient = "Primary LUAD", PDX = "LUAD PDX", cell_line = "LUAD cell lines")
system_colors <- c(patient = "#7B2C83", PDX = "#718B32", cell_line = "#0B8585")

save_figure <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
}

save_heatmap <- function(ht, stem, width, height) {
  grDevices::cairo_pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  png(file.path(figure_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

heatmap_grob <- function(ht) {
  grid.grabExpr(draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"))
}

make_upset <- function(combination_dt, set_columns, colors = NULL, title = "", subtitle = "") {
  dt <- copy(combination_dt)
  dt[, combination := apply(.SD, 1L, function(values) paste(as.integer(as.logical(values)), collapse = "")), .SDcols = set_columns]
  if (is.null(colors)) {
    counts <- dt[, .(count = .N), by = combination]
    counts[, fill_group := "TFs"]
  } else {
    counts <- dt[, .(count = .N), by = c("combination", colors)]
    setnames(counts, colors, "fill_group")
  }
  totals <- counts[, .(total = sum(count)), by = combination][order(-total, combination)]
  totals[, order := .I]
  counts <- merge(counts, totals, by = "combination")
  counts[, combination := factor(combination, levels = totals$combination)]
  bar <- ggplot(counts, aes(combination, count, fill = fill_group)) +
    geom_col(width = 0.72) +
    geom_text(data = totals, aes(factor(combination, levels = totals$combination), total, label = total),
              inherit.aes = FALSE, vjust = -0.25, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "TF count", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = if (is.null(colors)) "none" else "right")

  dot_rows <- rbindlist(lapply(seq_len(nrow(totals)), function(i) {
    bits <- strsplit(totals$combination[[i]], "", fixed = TRUE)[[1]]
    data.table(combination = totals$combination[[i]], set = set_columns, present = bits == "1")
  }))
  dot_rows[, combination := factor(combination, levels = totals$combination)]
  dot_rows[, set := factor(set, levels = rev(set_columns))]
  # Use one numeric coordinate system for both dots and connector segments.
  # Recent ggplot2 releases reject mixing a discrete y aesthetic with numeric
  # ymin/ymax values in the same panel.
  dot_rows[, set_index := as.integer(set)]
  connectors <- dot_rows[present == TRUE, .(ymin = min(set_index), ymax = max(set_index)), by = combination]
  dots <- ggplot(dot_rows, aes(combination, set_index)) +
    geom_segment(data = connectors, aes(x = combination, xend = combination, y = ymin, yend = ymax),
                 inherit.aes = FALSE, linewidth = 0.6, color = "#333333") +
    geom_point(aes(color = present), size = 2.8) +
    scale_y_continuous(breaks = seq_along(set_columns), labels = rev(set_columns)) +
    scale_color_manual(values = c(`TRUE` = "#222222", `FALSE` = "#D9D9D9"), guide = "none") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
  bar / dots + plot_layout(heights = c(2.1, 1))
}

promoter_tf <- fread(file.path(run_root, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
promoter_tf <- promoter_tf[analysis_definition == primary_definition]
promoter_system <- fread(file.path(run_root, "results", "promoter_gate", "system_tf_promoter_activity_summary.tsv"))
promoter_system <- promoter_system[analysis_definition == primary_definition]
promoter_sample <- fread(file.path(run_root, "results", "promoter_gate", "sample_tf_promoter_accessibility.tsv"))
promoter_sample <- promoter_sample[analysis_definition == primary_definition]
motif_inventory <- fread(file.path(run_root, "results", "motif_gate", "hc_tf_motif_inventory.tsv"))
motif_primary_dir <- file.path(run_root, "results", "motif_primary", "anchor_consensus_author")
motif_sample <- fread(file.path(motif_primary_dir, "anchor_consensus_sample_tf_best_motif.tsv"))
motif_system_raw <- fread(file.path(motif_primary_dir, "anchor_consensus_system_tf_motif_summary.tsv"))
motif_tf_raw <- fread(file.path(motif_primary_dir, "anchor_consensus_hc_tf_triple_system_summary.tsv"))
motif_lor_stats <- motif_sample[motif_test_status == "TESTED", .(
  mean_log2_odds_ratio = mean(log2_odds_ratio, na.rm = TRUE),
  sd_log2_odds_ratio = sd(log2_odds_ratio, na.rm = TRUE)
), by = .(system, TF)]
motif_lor_stats[!is.finite(sd_log2_odds_ratio), sd_log2_odds_ratio := 0]
motif_system <- merge(
  motif_system_raw[, .(
    system, TF, n_all_samples, n_motif_enriched,
    fraction_enriched_all_samples,
    system_motif_enriched = support_fixed_original_denominator
  )],
  motif_lor_stats,
  by = c("system", "TF"), all.x = TRUE
)
motif_tf <- motif_tf_raw[, .(
  TF,
  patient_motif_enriched = patient_fixed_support,
  PDX_motif_enriched = PDX_fixed_support,
  cell_line_motif_enriched = cell_line_fixed_support,
  triple_system_motif_enriched = triple_system_fixed_original_denominator
)]
motif_tf[, any_system_motif_enriched := fifelse(
  patient_motif_enriched == "TRUE" | PDX_motif_enriched == "TRUE" | cell_line_motif_enriched == "TRUE",
  "TRUE", "FALSE"
)]
atomic_manifest <- fread(file.path(run_root, "audit", "manifests", "figure2_atomic", "figure2_atomic_sample_manifest.tsv"))
atomic_counts <- atomic_manifest[, .N, by = system]
system_n <- setNames(atomic_counts$N, atomic_counts$system)

# Figure 2A: the same evidence order as the anchor, with LUAD counts and frozen TF input.
p2a <- ggplot() +
  annotate("rect", xmin = 0.1, xmax = 3.0, ymin = 2.9, ymax = 4.2, fill = "#F5F5F5", color = "#333333") +
  annotate("text", x = 1.55, y = 3.85, label = "ATAC-seq (LUAD)", fontface = "bold", size = 4.2) +
  annotate("text", x = 0.35, y = c(3.55, 3.25, 2.98), hjust = 0,
           label = c(sprintf("%d primary tumors", system_n[["patient"]]),
                     sprintf("%d LUAD PDX models", system_n[["PDX"]]),
                     sprintf("%d strict LUAD cell lines", system_n[["cell_line"]])),
           color = c(system_colors[["patient"]], system_colors[["PDX"]], system_colors[["cell_line"]]), size = 3.6) +
  annotate("segment", x = 1.55, xend = 1.55, y = 2.85, yend = 2.35, arrow = arrow(length = unit(0.12, "in"))) +
  annotate("rect", xmin = 0.35, xmax = 2.75, ymin = 1.55, ymax = 2.3, fill = "#E7F1F5", color = "#333333") +
  annotate("text", x = 1.55, y = 1.93, label = "Accessible chromatin profiling", size = 3.9) +
  annotate("text", x = 1.55, y = 1.15, label = "158 Figure-1 LUAD-specific TFs", fontface = "bold", size = 3.8) +
  annotate("segment", x = 1.1, xend = 0.75, y = 0.95, yend = 0.45, arrow = arrow(length = unit(0.11, "in"))) +
  annotate("segment", x = 2.0, xend = 2.35, y = 0.95, yend = 0.45, arrow = arrow(length = unit(0.11, "in"))) +
  annotate("text", x = 0.7, y = 0.22, label = "Promoter\naccessibility", size = 3.5) +
  annotate("text", x = 2.4, y = 0.22, label = "HOMER genome-wide\nmotif enrichment", size = 3.5) +
  coord_cartesian(xlim = c(0, 3.1), ylim = c(0, 4.3), clip = "off") +
  labs(title = "A") + theme_void(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 16))

# Figure 2B: promoter accessibility combinations among all 158 Figure 1 TFs.
promoter_combo <- promoter_tf[, .(
  TF,
  patient = patient_promoter_accessible == "TRUE",
  PDX = PDX_promoter_accessible == "TRUE",
  cell_line = cell_line_promoter_accessible == "TRUE"
)]
p2b <- make_upset(promoter_combo, c("PDX", "patient", "cell_line"),
                  title = "B", subtitle = sprintf("Promoter-accessibility patterns across %d input TFs", nrow(promoter_combo)))

# Figure 2C: anchor-style activity-category plus sample-level promoter matrix.
activity_file <- file.path(project_root, "pilots", "tnbc-chromatin-tf-nc-2026", "execution", "luad-figure1", "results", "tables", "cross_system_TF_effects.tsv.gz")
activity <- fread(activity_file)
activity <- activity[TF %in% promoter_tf$TF & cohort %in% c("TCGA", "GSE81089", "PDMR_PDX", "DepMap_22Q2")]
cohort_names <- c(TCGA = "TCGA", GSE81089 = "GSE81089", PDMR_PDX = "PDMR PDX", DepMap_22Q2 = "DepMap")
mean_activity <- dcast(activity, cohort ~ TF, value.var = "mean_LUAD")
activity_mat <- as.matrix(mean_activity[, -"cohort"])
rownames(activity_mat) <- cohort_names[mean_activity$cohort]
for (missing_tf in setdiff(promoter_tf$TF, colnames(activity_mat))) {
  activity_mat <- cbind(activity_mat, setNames(rep(NA_real_, nrow(activity_mat)), missing_tf))
}
activity_category <- matrix("Not present", nrow(activity_mat), ncol(activity_mat), dimnames = dimnames(activity_mat))
activity_category[is.finite(activity_mat) & activity_mat < 0] <- "<0"
activity_category[is.finite(activity_mat) & activity_mat >= 0 & activity_mat < 1] <- "0–1"
activity_category[is.finite(activity_mat) & activity_mat >= 1 & activity_mat < 3] <- "1–3"
activity_category[is.finite(activity_mat) & activity_mat >= 3] <- ">3"

promoter_wide <- dcast(promoter_sample, system + sample_id ~ TF, value.var = "promoter_accessible")
promoter_matrix <- as.matrix(promoter_wide[, -c("system", "sample_id")])
rownames(promoter_matrix) <- promoter_wide$sample_id
promoter_display <- matrix("Closed", nrow(promoter_matrix), ncol(promoter_matrix), dimnames = dimnames(promoter_matrix))
for (i in seq_len(nrow(promoter_matrix))) {
  promoter_display[i, promoter_matrix[i, ] == "TRUE"] <- promoter_wide$system[[i]]
}
tf_mean <- activity[, .(mean_across_activity_cohorts = mean(mean_LUAD, na.rm = TRUE)), by = TF]
tf_order_dt <- merge(promoter_tf[, .(TF, HC = HC_TF_promoter_activity_definition == "TRUE")], tf_mean, by = "TF", all.x = TRUE)
setorder(tf_order_dt, -HC, -mean_across_activity_cohorts, TF)
tf_order <- tf_order_dt$TF
activity_category <- activity_category[, tf_order, drop = FALSE]
promoter_display <- promoter_display[, tf_order, drop = FALSE]
combined_matrix <- rbind(activity_category, promoter_display)
row_group <- factor(c(rep("TF activity (mean NES)", nrow(activity_category)), promoter_wide$system),
                    levels = c("TF activity (mean NES)", "patient", "PDX", "cell_line"))
hc_annotation <- HeatmapAnnotation(
  `TF group` = ifelse(tf_order_dt$HC, "HC-TF", "Not HC-TF"),
  col = list(`TF group` = c(`HC-TF` = "#112D3C", `Not HC-TF` = "#C8C8C8")),
  show_annotation_name = TRUE
)
ht2c <- Heatmap(
  combined_matrix,
  name = "Evidence",
  col = c(
    `Not present` = "white", Closed = "white", `<0` = "#C1C7D0", `0–1` = "#C77A9B",
    `1–3` = "#A43870", `>3` = "#6C003B", patient = system_colors[["patient"]],
    PDX = system_colors[["PDX"]], cell_line = system_colors[["cell_line"]]
  ),
  cluster_columns = FALSE, cluster_rows = FALSE, row_split = row_group,
  row_gap = unit(c(2, 1, 1), "mm"), column_names_rot = 90,
  column_names_gp = gpar(fontsize = 4.2), row_names_gp = gpar(fontsize = 4.5),
  top_annotation = hc_annotation, show_heatmap_legend = TRUE,
  column_title = "C — TF activity categories and sample-level promoter accessibility"
)
save_heatmap(ht2c, "Figure2C_activity_promoter_matrix", 17, 10)
g2c <- heatmap_grob(ht2c)

# Figure 2D: B0 author-style motif combinations among motif-testable HC-TFs.
motif_combo <- merge(motif_tf[any_system_motif_enriched == "TRUE"],
                     motif_inventory[, .(TF, motif_database_category)], by = "TF")
motif_combo[, `:=`(
  patient = patient_motif_enriched == "TRUE",
  PDX = PDX_motif_enriched == "TRUE",
  cell_line = cell_line_motif_enriched == "TRUE"
)]
p2d <- make_upset(motif_combo, c("PDX", "patient", "cell_line"), colors = "motif_database_category",
                  title = "D", subtitle = sprintf("HC-TFs with >=1-system motif support (n = %d/%d testable)",
                                                   nrow(motif_combo), motif_inventory[motif_testable == "TRUE", .N])) +
  plot_annotation(theme = theme(legend.position = "right"))

# Figure 2E: B0 system LOR/SD/frequency and top motif distributions.
supported_tfs <- motif_tf[any_system_motif_enriched == "TRUE", TF]
motif_system <- motif_system[TF %in% supported_tfs]
motif_sample <- motif_sample[TF %in% supported_tfs]
motif_system[, system := factor(system, levels = system_levels)]
motif_system[, overall_mean := mean(mean_log2_odds_ratio), by = TF]
motif_system[, TF_order := factor(TF, levels = unique(TF[order(overall_mean)]))]
p2e_left <- ggplot(motif_system, aes(mean_log2_odds_ratio, TF_order, color = system,
                                     size = pmax(sd_log2_odds_ratio, 0.04), fill = system_motif_enriched)) +
  geom_point(shape = 21, stroke = 0.9) +
  scale_color_manual(values = system_colors, labels = system_labels) +
  scale_fill_manual(values = c(`TRUE` = alpha("black", 0.15), `FALSE` = "white"), guide = "none") +
  scale_size_continuous(range = c(1.4, 5), name = "SD of LOR") +
  labs(title = "E", x = "Mean log2 odds ratio\n(target vs background)", y = "TF with motif", color = NULL) +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")

top_tf_table <- motif_system[, .(overall = mean(mean_log2_odds_ratio)), by = TF][order(-overall)]
top_tfs <- head(top_tf_table$TF, 6L)
top_values <- motif_sample[TF %in% top_tfs & motif_test_status == "TESTED" & is.finite(log2_odds_ratio)]
top_values[, TF := factor(TF, levels = top_tfs)]
p2e_right <- ggplot(top_values, aes(system, log2_odds_ratio, color = system)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.45) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.75) +
  facet_wrap(~TF, ncol = 1, scales = "free_y") +
  scale_color_manual(values = system_colors, labels = system_labels) +
  labs(x = NULL, y = "Log2 odds ratio", color = NULL) +
  theme_bw(base_size = 7) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
p2e <- p2e_left + p2e_right + plot_layout(widths = c(1.3, 1))
save_figure(p2a, "Figure2A_design", 5.2, 4.8)
save_figure(p2b, "Figure2B_promoter_upset", 5.5, 5.5)
save_figure(p2d, "Figure2D_motif_upset", 5.5, 5.5)
save_figure(p2e, "Figure2E_motif_enrichment", 11, 8.5)

top_row <- p2a + p2b + p2d + p2e + plot_layout(widths = c(1, 1, 1, 2.2))
figure2 <- top_row / wrap_elements(full = g2c) + plot_layout(heights = c(1, 1.7)) +
  plot_annotation(
    title = "Figure 2 — Chromatin support for the frozen LUAD TF program",
    subtitle = "Primary motif analysis: B0 shared lung-accessibility background; HOMER q ≤ 1×10⁻⁵; ≥50% of each system"
  )
save_figure(figure2, "Figure2_complete", 24, 17)

# Pairwise system tests for motif LOR, reported independently of plotting.
comparisons <- list(c("patient", "PDX"), c("patient", "cell_line"), c("PDX", "cell_line"))
wilcox_rows <- rbindlist(lapply(unique(motif_sample$TF), function(tf) {
  rbindlist(lapply(comparisons, function(pair) {
    x <- motif_sample[TF == tf & system == pair[[1]], log2_odds_ratio]
    y <- motif_sample[TF == tf & system == pair[[2]], log2_odds_ratio]
    p <- if (length(x) > 0L && length(y) > 0L) wilcox.test(x, y, exact = FALSE)$p.value else NA_real_
    data.table(TF = tf, system_1 = pair[[1]], system_2 = pair[[2]], p_value = p,
               median_1 = median(x), median_2 = median(y))
  }))
}))
wilcox_rows[, FDR := p.adjust(p_value, method = "BH")]
fwrite(wilcox_rows, file.path(table_dir, "Figure2E_pairwise_system_wilcoxon.tsv"), sep = "\t")

# Supplementary Figure 3 A-C: sample peak counts and honest IDR status.
patient_peaks <- fread(file.path(run_root, "audit", "tcga_luad_accessibility", "tcga_luad_peak_counts_threshold_sensitivity.tsv"))
patient_qc <- patient_peaks[, .(system = "patient", sample_id, peak_count = cpm1_both_reps,
                               idr_status = "Public fixed-peak matrix; IDR not rerun")]
raw_qc <- fread(file.path(run_root, "audit", "raw_atac_qc", "figure2_raw_atac_qc.tsv"))
raw_qc <- raw_qc[analysis_included == "TRUE", .(system, sample_id, peak_count,
                                                 idr_status = "Single public library; IDR N/A")]
peak_qc <- rbind(patient_qc, raw_qc, fill = TRUE)
peak_qc[, system := factor(system, levels = system_levels)]
make_peak_plot <- function(system_name, letter) {
  dt <- peak_qc[system == system_name][order(-peak_count)]
  dt[, sample_id := factor(sample_id, levels = sample_id)]
  ggplot(dt, aes(sample_id, peak_count / 1000, fill = system)) +
    geom_col(width = 0.8) + scale_fill_manual(values = system_colors, guide = "none") +
    labs(title = paste0(letter, " — ", system_labels[[system_name]], " ATAC-seq"),
         subtitle = unique(dt$idr_status), x = NULL, y = "Number of peaks (×1000)") +
    theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), panel.grid.minor = element_blank())
}
p3a <- make_peak_plot("patient", "A")
p3b <- make_peak_plot("PDX", "B")
p3c <- make_peak_plot("cell_line", "C")

# Supplementary Figure 3 D-F: nls(SSasymp) saturation fits.
saturation <- fread(file.path(run_root, "results", "supplementary_figure3", "saturation_summary.tsv"))
saturation_raw <- fread(file.path(run_root, "results", "supplementary_figure3", "saturation_1000_permutations.tsv"))
fit_rows <- list()
sat_plots <- list()
letters_sat <- c(patient = "D", PDX = "E", cell_line = "F")
for (system_name in system_levels) {
  dt <- saturation[system == system_name]
  fit <- try(nls(mean_cumulative_consensus_peaks ~ SSasymp(sample_number, Asym, R0, lrc), data = dt), silent = TRUE)
  if (inherits(fit, "try-error")) stop("SSasymp failed for ", system_name)
  pars <- coef(fit)
  grid_x <- seq(0, max(dt$sample_number), length.out = 300)
  grid_y <- predict(fit, newdata = data.frame(sample_number = grid_x))
  observed <- unique(dt$observed_full_consensus_peaks)
  saturation_fraction <- observed / unname(pars[["Asym"]])
  fit_rows[[system_name]] <- data.table(system = system_name, asymptote = unname(pars[["Asym"]]),
                                         observed_consensus_peaks = observed, saturation_fraction,
                                         R0 = unname(pars[["R0"]]), lrc = unname(pars[["lrc"]]))
  plot_dt <- saturation_raw[system == system_name]
  sat_plots[[system_name]] <- ggplot(plot_dt, aes(sample_number, cumulative_consensus_peaks)) +
    geom_boxplot(aes(group = sample_number), outlier.shape = NA, width = 0.55, fill = alpha(system_colors[[system_name]], 0.35), linewidth = 0.25) +
    geom_line(data = data.table(sample_number = grid_x, fitted = grid_y), aes(sample_number, fitted), inherit.aes = FALSE, linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
             label = sprintf("Observed: %s\nEstimated: %s\nSaturation: %.1f%%", comma(observed), comma(round(pars[["Asym"]])), 100 * saturation_fraction), size = 2.6) +
    labs(title = paste0(letters_sat[[system_name]], " — ", system_labels[[system_name]], " saturation"), x = "Number of samples", y = "Unique consensus peaks") +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())
}
fwrite(rbindlist(fit_rows), file.path(table_dir, "SupplementaryFigure3_saturation_nls_SSasymp.tsv"), sep = "\t")

# Supplementary Figure 3 G-I: Pearson correlation matrices.
corr <- fread(file.path(run_root, "results", "supplementary_figure3", "sample_accessibility_pearson_correlations.tsv"))
corr_grobs <- list()
letters_corr <- c(patient = "G", PDX = "H", cell_line = "I")
for (system_name in system_levels) {
  dt <- corr[system == system_name]
  mat <- dcast(dt, sample_1 ~ sample_2, value.var = "pearson_r")
  rn <- mat$sample_1
  mat <- as.matrix(mat[, -"sample_1"])
  rownames(mat) <- rn
  ht <- Heatmap(mat, name = "Pearson r", col = colorRamp2(c(-0.2, 0.4, 0.7, 1), c("#24205C", "#8A3F8D", "#F08A78", "#FFF8BE")),
                cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE,
                column_title = paste0(letters_corr[[system_name]], " — ", system_labels[[system_name]], " Pearson correlation"))
  save_heatmap(ht, paste0("SupplementaryFigure3", letters_corr[[system_name]], "_correlation_", system_name), 5, 4.8)
  corr_grobs[[system_name]] <- heatmap_grob(ht)
}

# Supplementary Figure 3 J-L: hierarchical genomic peak annotation.
annotation <- fread(file.path(run_root, "results", "supplementary_figure3", "peak_genomic_annotation.tsv"))
annotation[, category := factor(category, levels = c("distal", "intronic", "exonic", "promoter"))]
# A zero-peak sample has undefined category fractions. Keep it in the source
# table and in the frozen 13-PDX denominator, but omit undefined proportions
# from this descriptive panel and its pairwise fraction tests.
annotation_valid <- annotation[total_peaks > 0 & is.finite(fraction)]
feature_colors <- c(distal = "#7C8B6B", intronic = "#C98265", exonic = "#6E9CB6", promoter = "#C84532")
letters_ann <- c(patient = "J", PDX = "K", cell_line = "L")
ann_plots <- list()
for (system_name in system_levels) {
  ann_plots[[system_name]] <- ggplot(annotation_valid[system == system_name], aes(category, fraction, fill = category)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.12, size = 0.8) +
    scale_fill_manual(values = feature_colors, guide = "none") +
    labs(title = paste0(letters_ann[[system_name]], " — ", system_labels[[system_name]], " peak annotation"), x = NULL, y = "Proportion of peaks") +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())
}
ann_tests <- annotation_valid[, {
  pairs <- combn(levels(category), 2, simplify = FALSE)
  rbindlist(lapply(pairs, function(pair) {
    x <- fraction[category == pair[[1]]]
    y <- fraction[category == pair[[2]]]
    data.table(category_1 = pair[[1]], category_2 = pair[[2]], p_value = wilcox.test(x, y, exact = FALSE)$p.value)
  }))
}, by = system]
ann_tests[, FDR := p.adjust(p_value, method = "BH")]
fwrite(ann_tests, file.path(table_dir, "SupplementaryFigure3_peak_annotation_wilcoxon.tsv"), sep = "\t")

supp3 <- wrap_plots(
  p3a, p3b, p3c,
  sat_plots[["patient"]], sat_plots[["PDX"]], sat_plots[["cell_line"]],
  wrap_elements(full = corr_grobs[["patient"]]), wrap_elements(full = corr_grobs[["PDX"]]), wrap_elements(full = corr_grobs[["cell_line"]]),
  ann_plots[["patient"]], ann_plots[["PDX"]], ann_plots[["cell_line"]],
  ncol = 3
) + plot_annotation(title = "Supplementary Figure 3 — Genome-wide chromatin accessibility in LUAD systems")
save_figure(supp3, "SupplementaryFigure3_complete", 18, 20)

# Supplementary Figure 4A: deterministic filtering flow.
n_frozen <- nrow(promoter_tf)
n_patient_open <- promoter_tf[patient_promoter_accessible == "TRUE", .N]
n_triple_open <- promoter_tf[triple_system_promoter_accessible == "TRUE", .N]
n_hc <- promoter_tf[HC_TF_promoter_activity_definition == "TRUE", .N]
n_testable <- motif_inventory[motif_testable == "TRUE", .N]
n_triple_motif <- motif_tf[triple_system_motif_enriched == "TRUE", .N]
flow <- data.table(
  x = seq_len(6),
  label = c(sprintf("Frozen LUAD TFs\n(n=%d)", n_frozen), sprintf("Patient promoter-open\n(n=%d)", n_patient_open),
            sprintf("Three-system promoter-open\n(n=%d)", n_triple_open), sprintf("HC-TFs after NES filter\n(n=%d)", n_hc),
            sprintf("Motif-testable HC-TFs\n(n=%d)", n_testable), sprintf("Three-system motif\n(n=%d)", n_triple_motif))
)
p4a <- ggplot(flow, aes(x, 1)) +
  geom_segment(data = flow[x < max(x)], aes(x = x + 0.38, xend = x + 0.62, y = 1, yend = 1),
               arrow = arrow(length = unit(0.09, "in")), linewidth = 0.55) +
  geom_label(aes(label = label, fill = x == 4), label.size = 0.35, size = 3.1, label.padding = unit(0.22, "lines")) +
  scale_fill_manual(values = c(`TRUE` = "#37526B", `FALSE` = "white"), guide = "none") +
  coord_cartesian(xlim = c(0.55, 6.45), ylim = c(0.5, 1.5), clip = "off") +
  labs(title = "A — Frozen chromatin/TF filtering flow") + theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold"))

# Supplementary Figure 4B: promoter combinations with the anchor's sole
# activity exclusion (mean NES <0 in all three activity cohorts) highlighted.
wide_promoter <- promoter_tf[, .(
  TF,
  patient_promoter = patient_promoter_accessible == "TRUE",
  PDX_promoter = PDX_promoter_accessible == "TRUE",
  cell_promoter = cell_line_promoter_accessible == "TRUE",
  activity_filter_group = ifelse(
    all_three_system_mean_NES_negative == "TRUE" & triple_system_promoter_accessible == "TRUE",
    "Excluded: NES<0 in all three", ifelse(HC_TF_promoter_activity_definition == "TRUE", "HC-TF", "Other promoter pattern")
  )
)]
logical_cols <- c("PDX_promoter", "patient_promoter", "cell_promoter")
p4b <- make_upset(wide_promoter, logical_cols, colors = "activity_filter_group", title = "B",
                  subtitle = "Promoter combinations; all-three-negative NES exclusion highlighted")

# Supplementary Figure 4C: JASPAR/CIS-BP motif availability Venn.
cat_counts <- motif_inventory[, .N, by = motif_database_category]
get_count <- function(category) {
  value <- cat_counts[motif_database_category == category, N]
  if (length(value) == 0L) 0L else value[[1]]
}
circle <- rbind(
  data.table(group = "JASPAR", t = seq(0, 2 * pi, length.out = 300), cx = -0.55),
  data.table(group = "CIS-BP", t = seq(0, 2 * pi, length.out = 300), cx = 0.55)
)
circle[, `:=`(x = cx + cos(t), y = sin(t))]
p4c <- ggplot(circle, aes(x, y, group = group, fill = group)) +
  geom_polygon(alpha = 0.28, color = "#333333") +
  annotate("text", x = -1.05, y = 1.18, label = "JASPAR", size = 4) +
  annotate("text", x = 1.05, y = 1.18, label = "CIS-BP", size = 4) +
  annotate("text", x = -0.75, y = 0, label = get_count("JASPAR_only"), size = 5) +
  annotate("text", x = 0, y = 0, label = get_count("both"), size = 5) +
  annotate("text", x = 0.75, y = 0, label = get_count("CIS-BP_only"), size = 5) +
  annotate("text", x = 0, y = -1.35, label = sprintf("No assigned motif: %d", get_count("none")), size = 3.5) +
  scale_fill_manual(values = c(JASPAR = "#84C4C4", `CIS-BP` = "#D98EA9"), guide = "none") +
  coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.55, 1.45), clip = "off") +
  labs(title = "C — HC-TF DNA motif availability") + theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold"))

p4b_nested <- wrap_elements(full = p4b)
supp4 <- p4a / (p4b_nested + p4c + plot_layout(widths = c(1.5, 1))) +
  plot_layout(heights = c(0.65, 1.5)) +
  plot_annotation(title = "Supplementary Figure 4A–C — LUAD HC-TF evidence accounting")
save_figure(p4a, "SupplementaryFigure4A_filtering_flow", 14, 3.3)
save_figure(p4b, "SupplementaryFigure4B_promoter_activity_upset", 8, 6.5)
save_figure(p4c, "SupplementaryFigure4C_motif_availability", 6, 5.5)
save_figure(supp4, "SupplementaryFigure4A-C_complete", 16, 11)

# Supplementary Figure 5: heterogeneity, shared-background robustness, and the
# user-requested q<=0.05 / >=50%-of-samples sensitivity analysis.  None of
# these panels replaces the immutable B0 primary result in Figure 2.
state_root <- file.path(dirname(run_root), "luad-state-audit")
state_membership <- fread(file.path(
  state_root, "results", "tru_state_diagnostic", "state_membership.tsv"
))
state_membership[, subtype_plot := fifelse(
  as.character(classifiable) == "TRUE", subtype, "Unclassified"
)]
state_counts <- state_membership[, .N, by = .(system, subtype_plot)]
state_counts[, system := factor(system, levels = system_levels)]
state_counts[, subtype_plot := factor(subtype_plot, levels = c("TRU", "PP", "PI", "Unclassified"))]
state_palette <- c(TRU = "#2A9D8F", PP = "#E9C46A", PI = "#E76F51", Unclassified = "#B8B8B8")
p5a <- ggplot(state_counts, aes(system, N, fill = subtype_plot)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.25) +
  geom_text(aes(label = N), position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = state_palette, drop = FALSE) +
  scale_x_discrete(labels = system_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A — LUAD intrinsic-state heterogeneity", x = NULL, y = "Fraction of models", fill = "State") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 15, hjust = 1))

background_stability <- fread(file.path(
  run_root, "results", "motif_background_sensitivity", "hc_tf_background_stability.tsv"
))
background_long <- melt(
  background_stability,
  id.vars = c("TF", "backgrounds_with_triple_support", "primary_B0_supported"),
  measure.vars = paste0("B", 0:3, "_triple_support"),
  variable.name = "background", value.name = "supported"
)
background_long[, background := sub("_triple_support$", "", background)]
background_long[, supported := as.character(supported) == "TRUE"]
tf_background_order <- background_stability[
  order(backgrounds_with_triple_support, as.character(primary_B0_supported) == "TRUE", TF), TF
]
background_long[, TF := factor(TF, levels = tf_background_order)]
background_long[, background := factor(background, levels = paste0("B", 0:3))]
p5b <- ggplot(background_long, aes(background, TF, fill = supported)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(supported, "●", "")), color = "white", size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#315B6D", `FALSE` = "#EBEEF0"), guide = "none") +
  labs(title = "B — HC-TF support across shared backgrounds", x = "Background definition", y = NULL) +
  theme_minimal(base_size = 8) + theme(panel.grid = element_blank())

background_summary <- fread(file.path(
  run_root, "results", "motif_background_sensitivity", "background_variant_summary.tsv"
))
background_summary[, variant := factor(variant, levels = paste0("B", 0:3))]
p5c <- ggplot(background_summary, aes(variant, triple_system_hc_tf_count)) +
  geom_col(aes(fill = variant == "B0"), width = 0.68, color = "#333333", linewidth = 0.25) +
  geom_text(aes(label = triple_system_hc_tf_count), vjust = -0.35, fontface = "bold", size = 3.5) +
  geom_text(aes(y = -0.28, label = sprintf("%s peaks", comma(background_peak_count))),
            vjust = 1, size = 2.7, color = "#555555") +
  scale_fill_manual(values = c(`TRUE` = "#A7446B", `FALSE` = "#8DA9B5"), guide = "none") +
  scale_y_continuous(breaks = 0:6, limits = c(-0.65, 6), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "C — Background sensitivity", subtitle = "B0 is primary; B1–B3 are sensitivity analyses",
       x = NULL, y = "Three-system motif-supported HC-TFs") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

threshold_scenarios <- fread(file.path(
  run_root, "results", "motif_threshold_sensitivity", "motif_threshold_sensitivity_scenarios.tsv"
))
threshold_tfs <- fread(file.path(
  run_root, "results", "motif_threshold_sensitivity", "motif_threshold_sensitivity_triple_system_tfs.tsv"
))
scenario_levels <- threshold_scenarios$scenario
threshold_union <- sort(unique(threshold_tfs$TF))
threshold_grid <- CJ(scenario = scenario_levels, TF = threshold_union, unique = TRUE)
threshold_grid <- merge(
  threshold_grid,
  unique(threshold_tfs[, .(scenario, TF, supported = TRUE)]),
  by = c("scenario", "TF"), all.x = TRUE
)
threshold_grid[is.na(supported), supported := FALSE]
scenario_labels <- setNames(
  c("Primary: q≤1×10⁻⁵; ≥50%", "Sensitivity: q≤0.05; ≥50%"),
  scenario_levels
)
threshold_grid[, scenario := factor(scenario, levels = scenario_levels, labels = scenario_labels)]
threshold_grid[, TF := factor(TF, levels = threshold_union)]
p5d <- ggplot(threshold_grid, aes(TF, scenario, fill = supported)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.92, height = 0.78) +
  geom_text(aes(label = ifelse(supported, "●", "")), color = "white", size = 4) +
  scale_fill_manual(values = c(`TRUE` = "#A7446B", `FALSE` = "#E7E7E7"), guide = "none") +
  labs(title = "D — Motif-threshold sensitivity",
       subtitle = "The relaxed result is supplementary and does not replace B0",
       x = "Three-system HC-TF", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))

supp5 <- (p5a + p5b) / (p5c + p5d) +
  plot_layout(widths = c(1, 1.25), heights = c(1, 0.9)) +
  plot_annotation(title = "Supplementary Figure 5 — LUAD heterogeneity and robustness analyses")
save_figure(p5a, "SupplementaryFigure5A_state_heterogeneity", 7, 5)
save_figure(p5b, "SupplementaryFigure5B_background_stability", 6.5, 7)
save_figure(p5c, "SupplementaryFigure5C_background_counts", 6.5, 5)
save_figure(p5d, "SupplementaryFigure5D_threshold_sensitivity", 8, 4.5)
save_figure(supp5, "SupplementaryFigure5_complete", 15, 12)

index <- data.table(
  artifact = c("Figure2_complete", "SupplementaryFigure3_complete", "SupplementaryFigure4A-C_complete", "SupplementaryFigure5_complete"),
  pdf = file.path(figure_dir, c("Figure2_complete.pdf", "SupplementaryFigure3_complete.pdf", "SupplementaryFigure4A-C_complete.pdf", "SupplementaryFigure5_complete.pdf")),
  png = file.path(figure_dir, c("Figure2_complete.png", "SupplementaryFigure3_complete.png", "SupplementaryFigure4A-C_complete.png", "SupplementaryFigure5_complete.png")),
  atomic_manifest = file.path(run_root, "audit", "manifests", "figure2_atomic", "figure2_atomic_sample_manifest.tsv"),
  primary_definition = primary_definition
)
fwrite(index, file.path(table_dir, "Figure2_atomic_artifact_index.tsv"), sep = "\t")
cat("FIGURE2_AND_REQUIRED_SUPPLEMENTS_COMPLETE\n")

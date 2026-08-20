#!/usr/bin/env Rscript

# Publication-grade second pass for the LUAD transfer of the TNBC framework.
#
# This script intentionally does not re-run ARACNe, VIPER, ATAC processing, or
# HOMER.  It consumes the frozen analysis receipts produced on the cloud and
# rebuilds the missing anchor-equivalent panels, canonical supplementary
# figures, and tabular evidence packages.  It is therefore deterministic and
# safe to resume after a disconnected desktop session.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(ComplexUpset)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(estimate)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "Usage: run_publication_rebuild.R FIG1_ROOT FIG2_ROOT FIG3_ROOT",
    "FIG4_ROOT FIG5_ROOT OUTPUT_ROOT"
  ))
}
fig1 <- normalizePath(args[[1]], mustWork = TRUE)
fig2 <- normalizePath(args[[2]], mustWork = TRUE)
fig3 <- normalizePath(args[[3]], mustWork = TRUE)
fig4 <- normalizePath(args[[4]], mustWork = TRUE)
fig5 <- normalizePath(args[[5]], mustWork = TRUE)
out <- normalizePath(args[[6]], mustWork = FALSE)
figures <- file.path(out, "results", "atomic_figures")
tables <- file.path(out, "results", "tables")
supp_root <- file.path(out, "results", "supplementary_data")
estimate_root <- file.path(out, "work", "estimate")
audit_root <- file.path(out, "audit")
for (d in c(figures, tables, supp_root, estimate_root, audit_root)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
set.seed(20260820)

pal <- c(
  LUAD = "#2878B5", LUSC = "#C45A38", PDX = "#799349",
  cell_line = "#178F8A", adverse = "#B2182B", favorable = "#2166AC",
  nonsig = "#B7BDC3", strict = "#6A3D9A", sensitivity = "#D98C10"
)
theme_pub <- function(base_size = 8.2) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(color = "#40464D", size = rel(0.86)),
      axis.title = element_text(face = "plain"),
      strip.background = element_rect(fill = "#F1F3F5", color = "#D3D7DB"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}
save_plot <- function(name, plot, width, height, dpi = 240) {
  ggsave(file.path(figures, paste0(name, ".pdf")), plot,
         width = width, height = height, useDingbats = FALSE, limitsize = FALSE)
  ggsave(file.path(figures, paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi, limitsize = FALSE)
}
read_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  ids <- x[[1L]]
  x[[1L]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "double"
  rownames(m) <- ids
  m
}
safe_copy <- function(src, group, description, dest_name = basename(src)) {
  if (!file.exists(src)) return(NULL)
  dest_dir <- file.path(supp_root, sprintf("Data_%02d", group))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dest_dir, dest_name)
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) stop(sprintf("Failed to copy %s", src))
  info <- file.info(dest)
  data.table(
    data_group = group,
    file_name = dest_name,
    description = description,
    source_path = src,
    output_path = dest,
    bytes = info$size,
    md5 = unname(tools::md5sum(dest))
  )
}
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  path
}
copy_rows <- list()
add_copy <- function(...) {
  z <- safe_copy(...)
  if (!is.null(z)) copy_rows[[length(copy_rows) + 1L]] <<- z
}

# ---------------------------------------------------------------------------
# Frozen TF sets and motif evidence
# ---------------------------------------------------------------------------

promoter_tf <- fread(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
primary_def <- "PRIMARY_anchor_promoter_tcga_cpm1_both"
hc <- sort(promoter_tf[
  analysis_definition == primary_def & HC_TF_promoter_activity_definition == TRUE,
  unique(TF)
])
if (length(hc) != 31L) stop(sprintf("Expected 31 frozen HC-TFs; found %d", length(hc)))

threshold_file <- file.path(
  fig2, "results", "motif_threshold_sensitivity",
  "motif_threshold_sensitivity_triple_system_tfs.tsv"
)
motif_threshold <- fread(threshold_file)
strict_tf <- sort(unique(motif_threshold[scenario == "primary_anchor", TF]))
sensitivity_tf <- sort(unique(motif_threshold[
  q_threshold == 0.05 & prevalence_threshold == 0.5, TF
]))
if (!identical(strict_tf, c("FOXA3", "NFATC4", "XBP1"))) {
  stop(sprintf("Frozen strict motif set changed: %s", paste(strict_tf, collapse = ",")))
}
if (!identical(sensitivity_tf, c("ETV1", "FOXA3", "NFATC4", "XBP1", "ZNF75D"))) {
  stop(sprintf("Frozen q<=0.05/50%% motif set changed: %s", paste(sensitivity_tf, collapse = ",")))
}

system_summary <- fread(file.path(
  fig2, "results", "promoter_gate", "system_tf_promoter_activity_summary.tsv"
))[analysis_definition == primary_def & TF %in% hc]
system_levels <- c("patient", "PDX", "cell_line")
system_labels <- c(patient = "Patients", PDX = "PDX", cell_line = "Cell lines")
system_summary[, system := factor(system, levels = system_levels)]

# ---------------------------------------------------------------------------
# Main Figure 2C: compact aligned activity, promoter, and motif evidence.
# ---------------------------------------------------------------------------

activity_wide <- dcast(system_summary, TF ~ system, value.var = "mean_LUAD_NES")
promoter_wide <- dcast(system_summary, TF ~ system, value.var = "fraction_promoter_accessible")
activity_mat <- as.matrix(activity_wide[, ..system_levels])
promoter_mat <- as.matrix(promoter_wide[, ..system_levels])
rownames(activity_mat) <- activity_wide$TF
rownames(promoter_mat) <- promoter_wide$TF
activity_mat <- activity_mat[hc, , drop = FALSE]
promoter_mat <- promoter_mat[hc, , drop = FALSE]
colnames(activity_mat) <- unname(system_labels[system_levels])
colnames(promoter_mat) <- unname(system_labels[system_levels])

motif_primary_summary <- fread(file.path(
  fig2, "results", "motif_primary", "anchor_consensus_author",
  "anchor_consensus_hc_tf_triple_system_summary.tsv"
))
motif_mat <- matrix(0, nrow = length(hc), ncol = 3L,
                    dimnames = list(hc, unname(system_labels[system_levels])))
for (i in seq_len(nrow(motif_primary_summary))) {
  tf <- motif_primary_summary$TF[[i]]
  if (!tf %in% hc) next
  motif_mat[tf, ] <- as.integer(c(
    motif_primary_summary$patient_fixed_support[[i]],
    motif_primary_summary$PDX_fixed_support[[i]],
    motif_primary_summary$cell_line_fixed_support[[i]]
  ))
}
row_score <- rowMeans(scale(activity_mat), na.rm = TRUE) +
  2 * as.integer(hc %in% strict_tf) + as.integer(hc %in% setdiff(sensitivity_tf, strict_tf))
tf_order <- names(sort(row_score, decreasing = TRUE))
activity_mat <- activity_mat[tf_order, , drop = FALSE]
promoter_mat <- promoter_mat[tf_order, , drop = FALSE]
motif_mat <- motif_mat[tf_order, , drop = FALSE]

motif_class <- factor(
  fifelse(tf_order %in% strict_tf, "Strict three-system",
          fifelse(tf_order %in% sensitivity_tf, "q<=0.05 sensitivity", "No three-system support")),
  levels = c("Strict three-system", "q<=0.05 sensitivity", "No three-system support")
)
row_ha <- rowAnnotation(
  `Motif evidence` = motif_class,
  col = list(`Motif evidence` = c(
    `Strict three-system` = pal[["strict"]],
    `q<=0.05 sensitivity` = pal[["sensitivity"]],
    `No three-system support` = "#E2E5E8"
  )),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3.2, "mm")
)
nes_limit <- max(3, quantile(abs(activity_mat), 0.95, na.rm = TRUE))
activity_display <- activity_mat
activity_display[] <- pmax(-nes_limit, pmin(nes_limit, as.numeric(activity_mat)))
hm_activity <- Heatmap(
  activity_display, name = "Mean NES",
  col = colorRamp2(c(-nes_limit, 0, nes_limit), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6.2), column_names_gp = gpar(fontsize = 7),
  column_title = "TF activity", column_title_gp = gpar(fontface = "bold", fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.35), width = unit(34, "mm")
)
hm_promoter <- Heatmap(
  promoter_mat, name = "Open fraction",
  col = colorRamp2(c(0.5, 0.75, 1), c("#F2F2F2", "#9ECAE1", "#08519C")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 7),
  column_title = "Promoter accessibility", column_title_gp = gpar(fontface = "bold", fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.35), width = unit(34, "mm")
)
hm_motif <- Heatmap(
  motif_mat, name = "Strict motif",
  col = c(`0` = "#F1F3F5", `1` = pal[["strict"]]),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 7),
  column_title = "Distal motif support", column_title_gp = gpar(fontface = "bold", fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.35), width = unit(34, "mm")
)
cairo_pdf(file.path(figures, "Figure2C_compact_HC_TF_evidence.pdf"), width = 7.6, height = 7.8)
draw(row_ha + hm_activity + hm_promoter + hm_motif,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     column_title = "C  Thirty-one chromatin-prioritized LUAD TFs",
     column_title_gp = gpar(fontface = "bold", fontsize = 11))
grid.text("Promoters: open in >= half of samples per system; motif: HOMER q<1e-5 in >= half per system",
          x = unit(0.5, "npc"), y = unit(0.012, "npc"), gp = gpar(fontsize = 7, col = "#4C535A"))
dev.off()
png(file.path(figures, "Figure2C_compact_HC_TF_evidence.png"), width = 1824, height = 1872, res = 240)
draw(row_ha + hm_activity + hm_promoter + hm_motif,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     column_title = "C  Thirty-one chromatin-prioritized LUAD TFs",
     column_title_gp = gpar(fontface = "bold", fontsize = 11))
grid.text("Promoters: open in >= half of samples per system; motif: HOMER q<1e-5 in >= half per system",
          x = unit(0.5, "npc"), y = unit(0.012, "npc"), gp = gpar(fontsize = 7, col = "#4C535A"))
dev.off()

motif_fraction <- rbindlist(lapply(c("primary_anchor", "relaxed_q05_half_samples"), function(scenario_id) {
  z <- motif_threshold[scenario == scenario_id]
  if (!nrow(z)) return(NULL)
  melt(z,
       id.vars = c("scenario", "q_threshold", "prevalence_threshold", "TF"),
       measure.vars = patterns("_n_enriched$", "_n_all_samples$"),
       value.name = c("n_enriched", "n_total"), variable.name = "system_index")
}))
if (nrow(motif_fraction)) {
  motif_fraction[, system := system_levels[system_index]]
  motif_fraction[, fraction := n_enriched / n_total]
  motif_fraction[, rule := fifelse(q_threshold <= 1e-5, "Primary: q<1e-5", "Sensitivity: q<=0.05")]
  motif_fraction[, TF := factor(TF, levels = rev(sensitivity_tf))]
  motif_fraction[, system := factor(system, levels = system_levels, labels = unname(system_labels))]
  panel_2e <- ggplot(motif_fraction, aes(system, TF, size = fraction, fill = rule)) +
    geom_point(shape = 21, color = "#333333", stroke = 0.35) +
    geom_text(aes(label = sprintf("%d/%d", n_enriched, n_total)), size = 2.15, color = "white") +
    scale_size(range = c(7, 12), limits = c(0.5, 1)) +
    scale_fill_manual(values = c(`Primary: q<1e-5` = pal[["strict"]],
                                 `Sensitivity: q<=0.05` = pal[["sensitivity"]])) +
    facet_wrap(~rule, ncol = 1, scales = "free_y") +
    labs(title = "E  Cross-system motif prevalence",
         subtitle = "JASPAR priority; CIS-BP only when JASPAR is unavailable",
         x = NULL, y = NULL, size = "Fraction", fill = NULL) +
    theme_pub(8.5) +
    theme(panel.grid.major = element_line(color = "#E7E9EB", linewidth = 0.3),
          panel.grid.minor = element_blank())
  save_plot("Figure2E_cross_system_motif_prevalence", panel_2e, 5.8, 5.4)
}

# ---------------------------------------------------------------------------
# Figure 3B: target-gene community network, matching the anchor's logic.
# ---------------------------------------------------------------------------

reg_edges <- fread(file.path(fig3, "results", "tables", "Figure3_HC_TF_signed_regulon_edges.tsv.gz"))
activated <- unique(reg_edges[direction == "Activated", .(TF, target)])
target_tf <- activated[, .(TFs = list(sort(unique(TF))), regulator_n = uniqueN(TF)), by = target]
target_tf <- target_tf[regulator_n >= 3L]
target_pairs <- if (nrow(target_tf) >= 2L) {
  cmb <- t(combn(target_tf$target, 2L))
  z <- data.table(from = cmb[, 1], to = cmb[, 2])
  tf_lookup <- setNames(target_tf$TFs, target_tf$target)
  z[, shared_TFs := lengths(Map(intersect, tf_lookup[from], tf_lookup[to]))]
  z[shared_TFs >= 3L]
} else data.table(from = character(), to = character(), shared_TFs = integer())

target_graph <- graph_from_data_frame(target_pairs, directed = FALSE,
                                      vertices = target_tf[, .(name = target, regulator_n)])
if (vcount(target_graph) && ecount(target_graph)) {
  comm <- cluster_louvain(target_graph, weights = E(target_graph)$shared_TFs)
  target_nodes <- data.table(
    target = V(target_graph)$name,
    community = as.integer(membership(comm)[V(target_graph)$name]),
    regulator_n = V(target_graph)$regulator_n,
    degree = degree(target_graph)
  )
} else {
  target_nodes <- target_tf[, .(target, community = seq_len(.N), regulator_n, degree = 0L)]
}
connected_communities <- sort(unique(target_nodes[degree > 0, community]))
target_nodes[, display_group := "Isolated target"]
if (length(connected_communities)) {
  for (i in seq_along(connected_communities)) {
    target_nodes[community == connected_communities[[i]], display_group := paste0("Connected component ", i)]
  }
}
target_pairs <- merge(target_pairs, target_nodes[, .(from = target, from_community = community)],
                      by = "from", all.x = TRUE)
target_pairs <- merge(target_pairs, target_nodes[, .(to = target, to_community = community)],
                      by = "to", all.x = TRUE)
write_table(target_nodes, file.path(tables, "Figure3B_target_gene_network_nodes.tsv"))
write_table(target_pairs, file.path(tables, "Figure3B_target_gene_network_edges.tsv"))

universe <- unique(reg_edges$target)
target_go_path <- file.path(tables, "Figure3B_target_gene_community_GO.tsv")
cached_target_go <- if (file.exists(target_go_path)) fread(target_go_path) else data.table()
if (nrow(cached_target_go) && "community_gene_n" %in% names(cached_target_go)) {
  target_go <- cached_target_go
} else {
  target_go <- rbindlist(lapply(sort(unique(target_nodes$community)), function(k) {
    genes <- target_nodes[community == k, target]
    if (length(genes) < 3L) {
      return(data.table(
        community = k, community_gene_n = length(genes), ID = NA_character_,
        Description = "Not tested: community contains fewer than three genes",
        p.adjust = NA_real_, Count = length(genes)
      ))
    }
    ego <- tryCatch(suppressMessages(enrichGO(
      gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 5, maxGSSize = 500
    )), error = function(e) NULL)
    if (is.null(ego) || !nrow(as.data.frame(ego))) {
      return(data.table(community = k, community_gene_n = length(genes), ID = NA_character_,
                        Description = "No GO term returned", p.adjust = NA_real_, Count = 0L))
    }
    z <- as.data.table(ego@result)
    z[, `:=`(community = k, community_gene_n = length(genes))]
    z
  }), fill = TRUE)
}
setorder(target_go, community, p.adjust, -Count)
write_table(target_go, target_go_path)
community_label <- target_go[!is.na(ID), head(.SD, 1L), by = community]
community_label[, short := ifelse(nchar(Description) > 38,
                                  paste0(substr(Description, 1, 35), "..."), Description)]

if (vcount(target_graph)) {
  display_levels <- c(paste0("Connected component ", seq_along(connected_communities)), "Isolated target")
  display_colors <- c(grDevices::hcl.colors(max(1, length(connected_communities)), "Dark 3"), "#BFC5CA")
  display_colors <- setNames(display_colors[seq_along(display_levels)], display_levels)
  graph_tbl <- as_tbl_graph(target_graph) |>
    activate(nodes) |>
    mutate(
      display_group = factor(target_nodes$display_group[match(name, target_nodes$target)], levels = display_levels),
      regulator_n = target_nodes$regulator_n[match(name, target_nodes$target)]
    )
  panel_3b <- ggraph(graph_tbl, layout = "fr") +
    geom_edge_link(aes(width = shared_TFs), color = "#71808E", alpha = 0.36,
                   show.legend = TRUE) +
    geom_node_point(aes(size = regulator_n, fill = display_group), shape = 21,
                    color = "white", stroke = 0.55) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.4, max.overlaps = Inf) +
    scale_edge_width(range = c(0.35, 2.2)) +
    scale_size(range = c(4.2, 9.0)) +
    scale_fill_manual(values = display_colors) +
    labs(title = "B  Shared-target communities",
         subtitle = sprintf("%d targets meet the node rule; only %d pairs share >=3 TFs",
                            nrow(target_nodes), nrow(target_pairs)),
         size = "Regulating TFs", fill = "Topology", edge_width = "Shared TFs") +
    theme_void(base_size = 8.3) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "#4C535A"), legend.position = "bottom")
} else {
  panel_3b <- ggplot() + annotate("text", 0, 0, label = "No target-target edge met the frozen >=3 shared-TF rule") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + ggtitle("B  Shared-target communities")
}
save_plot("Figure3B_target_gene_community_network", panel_3b, 7.0, 6.0)

# ---------------------------------------------------------------------------
# Supplementary Figure 4D: MKI67 correlations in four anchor-equivalent systems.
# ---------------------------------------------------------------------------

cohort_specs <- list(
  TCGA_LUAD = list(
    activity = file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz"),
    expression = file.path(fig1, "data", "processed", "tcga_aracne_expression.tsv"),
    manifest = file.path(fig1, "data", "processed", "tcga_manifest.tsv"),
    label = "TCGA patients", system = "Patient"
  ),
  GSE41271_LUAD = list(
    activity = file.path(fig4, "results", "gse41271_tcga_projection", "GSE41271_TCGA_PROJECTED_viper_activity.tsv.gz"),
    expression = file.path(fig1, "data", "processed", "gse41271_aracne_expression.tsv"),
    manifest = file.path(fig1, "data", "processed", "gse41271_manifest.tsv"),
    label = "GSE41271 patients", system = "Independent patient"
  ),
  PDMR_LUAD = list(
    activity = file.path(fig1, "results", "pdmr", "PDMR_viper_activity.tsv.gz"),
    expression = file.path(fig1, "data", "processed", "pdmr_viper_expression.tsv.gz"),
    manifest = file.path(fig1, "data", "processed", "pdmr_manifest.tsv"),
    label = "PDMR PDX", system = "PDX"
  ),
  DEPMAP_LUAD = list(
    activity = file.path(fig1, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
    expression = file.path(fig1, "data", "processed", "depmap_viper_expression.tsv.gz"),
    manifest = file.path(fig1, "data", "processed", "depmap_manifest.tsv"),
    label = "DepMap cell lines", system = "Cell line"
  ),
  GSE81089_LUAD = list(
    activity = file.path(fig1, "results", "gse81089", "GSE81089_viper_activity.tsv.gz"),
    expression = file.path(fig1, "data", "processed", "gse81089_aracne_expression.tsv"),
    manifest = file.path(fig1, "data", "processed", "gse81089_manifest.tsv"),
    label = "GSE81089 patients", system = "Additional patient"
  )
)

mki67_rows <- rbindlist(lapply(names(cohort_specs), function(cohort_id) {
  s <- cohort_specs[[cohort_id]]
  a <- read_matrix(s$activity)
  e <- read_matrix(s$expression)
  m <- fread(s$manifest)
  ids <- m[group == "LUAD", sample_id]
  ids <- Reduce(intersect, list(ids, colnames(a), colnames(e)))
  if (!"MKI67" %in% rownames(e)) stop(sprintf("MKI67 missing in %s", cohort_id))
  rbindlist(lapply(hc, function(tf) {
    if (!tf %in% rownames(a)) {
      return(data.table(
        cohort = cohort_id, cohort_label = s$label, system = s$system, TF = tf,
        n = 0L, pearson_r = NA_real_, p_value = NA_real_
      ))
    }
    keep <- is.finite(a[tf, ids]) & is.finite(e["MKI67", ids])
    test <- if (sum(keep) >= 4L) cor.test(a[tf, ids][keep], e["MKI67", ids][keep], method = "pearson") else NULL
    data.table(
      cohort = cohort_id, cohort_label = s$label, system = s$system, TF = tf,
      n = sum(keep), pearson_r = if (is.null(test)) NA_real_ else unname(test$estimate),
      p_value = if (is.null(test)) NA_real_ else test$p.value
    )
  }))
}), fill = TRUE)
mki67_rows[, FDR := p.adjust(p_value, method = "BH"), by = cohort]
mki67_rows[, anchor_display_support := abs(pearson_r) > 0.4 & FDR <= 0.05]
write_table(mki67_rows, file.path(tables, "SupplementaryFigure4D_MKI67_correlations.tsv"))

mki67_main <- mki67_rows[cohort != "GSE81089_LUAD"]
mki67_main[, cohort_label := factor(cohort_label,
  levels = c("TCGA patients", "GSE41271 patients", "PDMR PDX", "DepMap cell lines"))]
mki67_main[, TF := factor(TF, levels = rev(hc))]
panel_s4d <- ggplot(mki67_main, aes(cohort_label, TF, fill = pearson_r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_point(data = mki67_main[anchor_display_support == TRUE], shape = 8, size = 1.6, color = "#1A1A1A") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(title = "D  HC-TF activity versus MKI67 expression",
       subtitle = "Star: |Pearson r|>0.4 and within-cohort BH FDR<=0.05",
       x = NULL, y = NULL, fill = "Pearson r") +
  theme_pub(7.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
save_plot("SupplementaryFigure4D_MKI67_four_system", panel_s4d, 6.6, 7.1)

# ---------------------------------------------------------------------------
# Supplementary Figure 5A: cross-system activity heterogeneity.
# ---------------------------------------------------------------------------

skew_one <- function(cohort_id, s) {
  a <- read_matrix(s$activity)
  m <- fread(s$manifest)
  ids <- intersect(m[group == "LUAD", sample_id], colnames(a))
  rbindlist(lapply(hc, function(tf) {
    if (!tf %in% rownames(a)) {
      return(data.table(cohort = cohort_id, cohort_label = s$label, system = s$system,
                        TF = tf, n = 0L, skewness = NA_real_, z = NA_real_, p_value = NA_real_))
    }
    x <- as.numeric(a[tf, ids])
    x <- x[is.finite(x)]
    n <- length(x)
    g1 <- if (n >= 8L) mean((x - mean(x))^3) / sd(x)^3 else NA_real_
    se <- if (n >= 8L) sqrt(6 / n) else NA_real_
    p <- if (is.finite(g1) && is.finite(se) && se > 0) 2 * pnorm(-abs(g1 / se)) else NA_real_
    data.table(cohort = cohort_id, cohort_label = s$label, system = s$system,
               TF = tf, n = n, skewness = g1, z = g1 / se, p_value = p)
  }))
}
skew_rows <- rbindlist(lapply(names(cohort_specs)[1:4], function(k) skew_one(k, cohort_specs[[k]])))
skew_rows[, FDR := p.adjust(p_value, method = "BH"), by = cohort]
write_table(skew_rows, file.path(tables, "SupplementaryFigure5A_cross_system_skewness.tsv"))
skew_rows[, cohort_label := factor(cohort_label,
  levels = c("TCGA patients", "GSE41271 patients", "PDMR PDX", "DepMap cell lines"))]
skew_rows[, TF := factor(TF, levels = rev(hc))]
panel_s5a <- ggplot(skew_rows, aes(cohort_label, TF, fill = skewness)) +
  geom_tile(color = "white", linewidth = 0.18) +
  geom_point(data = skew_rows[FDR <= 0.05], shape = 8, color = "#1F1F1F", size = 1.45) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish) +
  labs(title = "A  Cross-system HC-TF activity skewness",
       subtitle = "Star: normal-approximation skewness test, within-cohort BH FDR<=0.05",
       x = NULL, y = NULL, fill = "Skewness") +
  theme_pub(7.5) + theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
save_plot("SupplementaryFigure5A_cross_system_skewness", panel_s5a, 6.6, 7.1)

# ---------------------------------------------------------------------------
# Supplementary Figure 5B: ESTIMATE scores in two patient cohorts.
# ---------------------------------------------------------------------------

run_estimate <- function(cohort_id, s) {
  a <- read_matrix(s$activity)
  e <- read_matrix(s$expression)
  m <- fread(s$manifest)
  ids <- Reduce(intersect, list(m[group == "LUAD", sample_id], colnames(a), colnames(e)))
  e <- e[, ids, drop = FALSE]
  if (quantile(e, 0.99, na.rm = TRUE) > 100) e <- log2(e + 1)
  input <- file.path(estimate_root, paste0(cohort_id, "_expression.gct"))
  filtered <- file.path(estimate_root, paste0(cohort_id, "_common_genes.gct"))
  scored <- file.path(estimate_root, paste0(cohort_id, "_estimate_scores.gct"))
  estimate_input <- as.data.table(e, keep.rownames = "NAME")
  fwrite(estimate_input, input, sep = "\t")
  if (!file.exists(scored)) {
    suppressMessages(filterCommonGenes(input.f = input, output.f = filtered, id = "GeneSymbol"))
    suppressMessages(estimateScore(filtered, scored, platform = "illumina"))
  }
  score_tab <- fread(scored, skip = 2, check.names = FALSE)
  score_name <- score_tab[[1L]]
  drop_cols <- intersect(c("NAME", "Name", "Description"), names(score_tab))
  score_values <- as.matrix(score_tab[, setdiff(names(score_tab), drop_cols), with = FALSE])
  storage.mode(score_values) <- "double"
  rownames(score_values) <- score_name
  # ESTIMATE reads tables with check.names=TRUE and therefore converts TCGA
  # hyphens to dots.  Map those syntactic names back to the frozen sample IDs
  # rather than silently losing the complete TCGA cohort.
  estimate_ids <- make.names(ids, unique = TRUE)
  if (!all(estimate_ids %in% colnames(score_values))) {
    stop(sprintf("ESTIMATE sample-name mapping failed for %s: %d/%d found",
                 cohort_id, sum(estimate_ids %in% colnames(score_values)), length(ids)))
  }
  score_values <- score_values[, estimate_ids, drop = FALSE]
  colnames(score_values) <- ids
  if (!"TumorPurity" %in% rownames(score_values) && "ESTIMATEScore" %in% rownames(score_values)) {
    purity <- cos(0.6049872018 + 0.0001467884 * score_values["ESTIMATEScore", ])
    score_values <- rbind(score_values, TumorPurity = purity)
  }
  score_names <- intersect(c("StromalScore", "ImmuneScore", "TumorPurity"), rownames(score_values))
  rbindlist(lapply(score_names, function(score_name) {
    rbindlist(lapply(hc, function(tf) {
      keep <- is.finite(a[tf, ids]) & is.finite(score_values[score_name, ids])
      test <- if (sum(keep) >= 4L) cor.test(a[tf, ids][keep], score_values[score_name, ids][keep]) else NULL
      data.table(
        cohort = cohort_id, cohort_label = s$label, score = score_name, TF = tf,
        n = sum(keep), pearson_r = if (is.null(test)) NA_real_ else unname(test$estimate),
        p_value = if (is.null(test)) NA_real_ else test$p.value
      )
    }))
  }))
}

estimate_rows <- rbindlist(lapply(c("TCGA_LUAD", "GSE41271_LUAD"), function(k) {
  tryCatch(run_estimate(k, cohort_specs[[k]]), error = function(e) {
    warning(sprintf("ESTIMATE failed for %s: %s", k, conditionMessage(e)))
    data.table()
  })
}), fill = TRUE)
if (nrow(estimate_rows)) {
  estimate_rows[, FDR := p.adjust(p_value, method = "BH"), by = .(cohort, score)]
  estimate_rows[, anchor_display_support := abs(pearson_r) > 0.4 & FDR <= 0.05]
  write_table(estimate_rows, file.path(tables, "SupplementaryFigure5B_ESTIMATE_correlations.tsv"))
  estimate_rows[, evidence := fifelse(FDR <= 0.05 & abs(pearson_r) > 0.4, "FDR<=0.05 & |r|>0.4",
                                      fifelse(FDR <= 0.05, "FDR<=0.05", "Not significant"))]
  estimate_rows[, score_label := factor(score,
    levels = c("StromalScore", "ImmuneScore", "TumorPurity"),
    labels = c("Stromal score", "Immune score", "Tumor purity"))]
  panel_s5b <- ggplot(estimate_rows, aes(pearson_r, -log10(pmax(FDR, 1e-300)), color = evidence)) +
    geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "#9EA4AA", linewidth = 0.35) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#6E747A", linewidth = 0.35) +
    geom_point(size = 1.7, alpha = 0.85) +
    geom_text_repel(data = estimate_rows[FDR <= 0.05 & abs(pearson_r) > 0.4],
                    aes(label = TF), size = 2.1, max.overlaps = 12, min.segment.length = 0) +
    facet_grid(cohort_label ~ score_label) +
    scale_color_manual(values = c(`FDR<=0.05 & |r|>0.4` = "#B2182B",
                                  `FDR<=0.05` = "#E69F00", `Not significant` = "#B8BDC2")) +
    labs(title = "B  HC-TF activity versus tumor-microenvironment scores",
         subtitle = "Pearson correlation; BH correction within cohort and score",
         x = "Pearson r", y = expression(-log[10]("BH FDR")), color = NULL) +
    theme_pub(7.2)
  save_plot("SupplementaryFigure5B_ESTIMATE_correlations", panel_s5b, 10.5, 6.4)
}

# ---------------------------------------------------------------------------
# Figure 4 compact main forests and canonical Supplementary Figure 6.
# ---------------------------------------------------------------------------

cox <- fread(file.path(fig4, "results", "tables", "Figure4_all_Cox_models.tsv"))
cox[, panel_id := paste(cohort, endpoint, model, sep = " | ")]
cox[, evidence := fifelse(FDR <= 0.05, "BH FDR<=0.05",
                          fifelse(p_value <= 0.05, "Nominal p<=0.05", "Not significant"))]
cox[, evidence := factor(evidence, levels = c("BH FDR<=0.05", "Nominal p<=0.05", "Not significant"))]

forest_defs <- data.table(
  letter = c("A", "B", "C", "D"),
  cohort = c("GSE41271_LUAD", "GSE41271_LUAD", "TCGA_LUAD", "TCGA_LUAD"),
  endpoint = c("OS", "RFS", "OS", "RFS")
)
for (i in seq_len(nrow(forest_defs))) {
  d <- forest_defs[i]
  z_all <- cox[cohort == d$cohort & endpoint == d$endpoint & model == "multivariable"]
  z <- z_all[p_value <= 0.05]
  if (!nrow(z)) z <- head(z_all[order(p_value)], 5L)
  z[, TF := factor(TF, levels = rev(TF[order(HR)]))]
  p <- ggplot(z, aes(HR, TF, color = evidence)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#777777", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.18, linewidth = 0.55) +
    geom_point(size = 2.1) +
    scale_x_log10() +
    scale_color_manual(values = c(`BH FDR<=0.05` = "#7A0177", `Nominal p<=0.05` = "#D95F0E",
                                  `Not significant` = "#8F969C"), drop = FALSE) +
    labs(title = sprintf("%s  %s %s", d$letter,
                         fifelse(d$cohort == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD"), d$endpoint),
         subtitle = sprintf("Multivariable Cox; %d/31 nominal, %d/31 FDR-significant",
                            sum(z_all$p_value <= 0.05), sum(z_all$FDR <= 0.05)),
         x = "Hazard ratio (log scale)", y = NULL, color = NULL) +
    theme_pub(8.0)
  save_plot(sprintf("Figure4%s_compact_forest", d$letter), p, 4.0, 3.2)
}

volcano_panel <- function(z) {
  z <- copy(z)
  z[, neglogp := -log10(pmax(p_value, 1e-300))]
  lab <- z[FDR <= 0.05]
  if (!nrow(lab)) lab <- head(z[order(p_value)], 2L)
  ggplot(z, aes(log2_HR, neglogp, color = evidence)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#757B80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#B8BDC2", linewidth = 0.3) +
    geom_point(size = 1.45, alpha = 0.82) +
    geom_text_repel(data = lab, aes(label = TF), size = 1.9, max.overlaps = 4, min.segment.length = 0) +
    scale_color_manual(values = c(`BH FDR<=0.05` = "#7A0177", `Nominal p<=0.05` = "#D95F0E",
                                  `Not significant` = "#B8BDC2"), drop = FALSE) +
    labs(title = unique(z$panel_id), x = expression(log[2](HR)), y = expression(-log[10](italic(P)))) +
    theme_pub(6.7) + theme(legend.position = "none")
}
volcano_list <- lapply(unique(cox$panel_id), function(k) volcano_panel(cox[panel_id == k]))
supp6_volcano <- wrap_plots(volcano_list, ncol = 4) +
  plot_annotation(title = "Supplementary Figure 6A-H | Continuous-activity Cox screens",
                  subtitle = "Purple: BH FDR<=0.05; orange: nominal p<=0.05; all 31 HC-TFs are displayed")
save_plot("SupplementaryFigure6A-H_Cox_volcanoes", supp6_volcano, 14.2, 6.7)

perm_counts <- fread(file.path(fig4, "results", "tables", "SupplementaryFigure7_permutation_null_counts.tsv.gz"))
perm_summary <- fread(file.path(fig4, "results", "tables", "SupplementaryFigure7_permutation_summary.tsv"))
perm_counts[, panel_id := paste(cohort, endpoint, model, sep = " | ")]
perm_summary[, panel_id := paste(cohort, endpoint, model, sep = " | ")]
perm_panel <- ggplot(perm_counts, aes(significant_TFs)) +
  geom_histogram(binwidth = 1, boundary = -0.5, fill = "#C9D7E5", color = "white", linewidth = 0.25) +
  geom_vline(aes(xintercept = observed_significant_TFs), color = "#B2182B", linewidth = 0.7) +
  geom_text(data = perm_summary,
            aes(x = Inf, y = Inf, label = sprintf("observed=%d\nempirical p=%.3f",
                                                 observed_significant_TFs, empirical_p)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 2.0) +
  facet_wrap(~panel_id, ncol = 4, scales = "free_y") +
  labs(title = "Supplementary Figure 6I-P | Five-thousand-permutation null distributions",
       subtitle = "Red line: observed number of nominally significant HC-TFs",
       x = "TFs with nominal p<=0.05", y = "Permutations") +
  theme_pub(6.7) + theme(legend.position = "none")
save_plot("SupplementaryFigure6I-P_permutation_nulls", perm_panel, 14.2, 6.7)

# ---------------------------------------------------------------------------
# Figure 5 and canonical Supplementary Figure 7.
# ---------------------------------------------------------------------------

cell_inputs <- readRDS(file.path(fig5, "data", "processed", "figure5_cellline_inputs.rds"))
cell_activity <- cell_inputs$activity
cell_all <- fread(file.path(fig5, "results", "tables", "Figure5ABC_cellline_concordance_all.tsv.gz"))
drug_map <- fread(file.path(fig5, "results", "tables", "Figure5_drug_identity_map.tsv"))
rep_pairs <- fread(file.path(fig5, "results", "tables", "Figure5D_replicated_drug_TF_pairs.tsv"))[replicated == TRUE]
if (nrow(rep_pairs) != 5L) warning(sprintf("Expected 5 replicated drug-TF pairs; found %d", nrow(rep_pairs)))
drug_annotation <- drug_map[, .(
  pathway_or_moa = {
    values <- sort(unique(na.omit(as.character(pathway_or_moa))))
    if (length(values)) paste(values, collapse = "; ") else "Unclassified"
  }
), by = canonical_drug]

membership <- unique(cell_all[significant == TRUE, .(pair = paste(canonical_drug, TF, sep = "::"), dataset)])
membership_wide <- dcast(membership, pair ~ dataset, fun.aggregate = length, value.var = "dataset")
for (d in c("GDSC2", "CTRPv2", "PRISM")) {
  if (!d %in% names(membership_wide)) membership_wide[, (d) := 0L]
  membership_wide[, (d) := get(d) > 0]
}
write_table(membership_wide, file.path(tables, "Figure5D_significant_pair_membership_binary.tsv"))
panel_5d <- upset(
  as.data.frame(membership_wide), intersect = c("GDSC2", "CTRPv2", "PRISM"),
  name = "Dataset", min_size = 1,
  base_annotations = list(
    `Significant drug-TF pairs` = intersection_size(text = list(size = 2.8))
  ),
  width_ratio = 0.22, set_sizes = upset_set_size()
) +
  ggtitle("D  Replication of significant drug-TF associations") +
  theme(plot.title = element_text(face = "bold", size = 10))
save_plot("Figure5D_true_upset", panel_5d, 6.0, 4.2)

rep_detail <- merge(
  cell_all[significant == TRUE],
  rep_pairs[, .(canonical_drug, TF, replicated_direction = direction)],
  by = c("canonical_drug", "TF")
)
rep_detail <- merge(rep_detail, drug_annotation,
                    by = "canonical_drug", all.x = TRUE)
rep_detail[, pair_label := paste(drug_name, TF, sep = " | ")]
rep_detail[, effect := c_index - 0.5]
write_table(rep_detail, file.path(tables, "Figure5_replicated_pair_dataset_details.tsv"))

panel_5f <- ggplot(rep_detail, aes(TF, interaction(drug_name, pathway_or_moa, sep = " — "),
                                  size = -log10(pmax(FDR, 1e-300)), fill = effect)) +
  geom_point(shape = 21, color = "#333333", stroke = 0.35) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-0.25, 0.25), oob = squish) +
  facet_wrap(~dataset, ncol = 1) +
  labs(title = "F  Pathway context of replicated associations",
       subtitle = "Negative concordance effect indicates higher TF activity associated with resistance",
       x = NULL, y = NULL, size = expression(-log[10]("BH FDR")), fill = "C-index - 0.5") +
  theme_pub(7.7) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot("Figure5F_pathway_TF_bubble", panel_5f, 7.2, 5.6)

# All contributing cell-line scatterplots for the five replicated pairs.
scatter_plots <- list()
scatter_receipts <- list()
for (i in seq_len(nrow(rep_detail))) {
  z <- rep_detail[i]
  ds <- cell_inputs$datasets[[z$dataset]]
  response <- copy(ds$paired)[eligible_luad == TRUE & canonical_drug == z$canonical_drug]
  response <- response[, .(aac = median(aac)), by = depmap_id]
  response <- response[depmap_id %in% colnames(cell_activity)]
  response[, activity := as.numeric(cell_activity[z$TF, depmap_id])]
  response <- response[is.finite(activity) & is.finite(aac)]
  rho <- if (nrow(response) >= 3L) cor(response$activity, response$aac, method = "spearman") else NA_real_
  scatter_receipts[[length(scatter_receipts) + 1L]] <- data.table(
    dataset = z$dataset, canonical_drug = z$canonical_drug, drug_name = z$drug_name,
    TF = z$TF, n = nrow(response), spearman_rho_recomputed = rho, reported_FDR = z$FDR
  )
  scatter_plots[[length(scatter_plots) + 1L]] <- ggplot(response, aes(activity, aac)) +
    geom_point(size = 1.5, alpha = 0.78, color = pal[["LUAD"]]) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, color = "#333333") +
    labs(title = sprintf("%s | %s", z$drug_name, z$TF),
         subtitle = sprintf("%s; n=%d; Spearman rho=%.2f; FDR=%.2g", z$dataset, nrow(response), rho, z$FDR),
         x = "TF activity (NES)", y = "Drug response (AAC)") +
    theme_pub(6.5) + theme(legend.position = "none")
}
scatter_receipts <- rbindlist(scatter_receipts, fill = TRUE)
write_table(scatter_receipts, file.path(tables, "SupplementaryFigure7B_scatter_receipts.tsv"))
if (length(scatter_plots)) {
  scatter_grid <- wrap_plots(scatter_plots, ncol = 3) +
    plot_annotation(title = "Supplementary Figure 7B | Replicated association examples",
                    subtitle = "Every contributing dataset is shown; response is AAC")
  save_plot("SupplementaryFigure7B_replicated_pair_scatters", scatter_grid, 12.6,
            max(5.8, ceiling(length(scatter_plots) / 3) * 3.0))
}

# Full significant matrix, not only the replicated subset.
all_sig <- cell_all[significant == TRUE]
all_sig <- merge(all_sig, drug_annotation,
                 by = "canonical_drug", all.x = TRUE)
drug_rank <- all_sig[, .(best_FDR = min(FDR), association_n = .N), by = .(canonical_drug, drug_name, pathway_or_moa)]
setorder(drug_rank, best_FDR, -association_n)
drug_rank[, display_label := paste0(drug_name, " — ", fifelse(is.na(pathway_or_moa), "Unclassified", pathway_or_moa))]
all_sig <- merge(all_sig, drug_rank[, .(canonical_drug, display_label)], by = "canonical_drug")
all_sig[, display_label := factor(display_label, levels = rev(unique(drug_rank$display_label)))]
all_sig[, effect := c_index - 0.5]
write_table(all_sig, file.path(tables, "SupplementaryFigure7C_all_significant_associations.tsv.gz"))
panel_s7c <- ggplot(all_sig, aes(TF, display_label, size = -log10(pmax(FDR, 1e-300)), fill = effect)) +
  geom_point(shape = 21, color = "#333333", stroke = 0.25, alpha = 0.85) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(title = "Supplementary Figure 7C | Complete significant drug-TF matrix",
       subtitle = "All within-dataset BH FDR<=0.05 associations; rows ordered by best FDR",
       x = NULL, y = NULL, size = expression(-log[10]("BH FDR")), fill = "C-index - 0.5") +
  theme_pub(6.1) + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                         strip.text.y = element_text(angle = 0))
full_height <- min(36, max(12, 0.11 * uniqueN(all_sig$display_label) + 4))
save_plot("SupplementaryFigure7C_complete_bubble_matrix", panel_s7c, 14.2, full_height, dpi = 200)

# PDX validation funnel and complete null screen replace blank placeholders.
pdx_scores <- fread(file.path(fig5, "results", "tables", "Figure5_PDXE_v2_lung_state_scores.tsv"))
pdx_coverage <- fread(file.path(fig5, "results", "tables", "SupplementaryFigure9A_PDXE_drug_coverage.tsv"))
pdx_tests <- fread(file.path(fig5, "results", "tables", "SupplementaryFigure9B_PDXE_all_concordance_tests.tsv"))
pdx_validated <- fread(file.path(fig5, "results", "tables", "Figure5GH_PDX_validated_pairs.tsv"))
funnel <- data.table(
  step = factor(c("NSCLC PDX classified", "LUAD-like PDX", "Drugs with >=10 PDX", "Drug-TF tests", "Validated replicated pairs"),
                levels = rev(c("NSCLC PDX classified", "LUAD-like PDX", "Drugs with >=10 PDX", "Drug-TF tests", "Validated replicated pairs"))),
  n = c(nrow(pdx_scores), sum(pdx_scores$prediction == "LUAD-like"),
        sum(pdx_coverage$PDX_n >= 10), nrow(pdx_tests), nrow(pdx_validated))
)
write_table(funnel, file.path(tables, "Figure5G_PDX_validation_funnel.tsv"))
panel_5g <- ggplot(funnel, aes(n, step, fill = step)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.15, fontface = "bold", size = 3) +
  scale_fill_manual(values = rev(c("#A6CEE3", "#6BAED6", "#4292C6", "#2171B5", "#B2182B"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "G  In-vivo validation funnel",
       subtitle = "Matched public PDX response data do not validate any replicated cell-line pair",
       x = "Models, drugs, or tests", y = NULL) +
  theme_pub(8.3) + theme(legend.position = "none")
save_plot("Figure5G_PDX_validation_funnel", panel_5g, 5.7, 3.7)

if (nrow(pdx_tests)) {
  pdx_tests[, status := fifelse(FDR <= 0.05 & cellline_replicated_same_direction == TRUE,
                                "Validated", fifelse(!is.na(dataset_n), "Cell-line replicated pair", "Other PDX screen"))]
  panel_5h <- ggplot(pdx_tests, aes(c_index - 0.5, -log10(pmax(FDR, 1e-300)), color = status)) +
    geom_vline(xintercept = 0, color = "#B8BDC2", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#777777", linewidth = 0.35) +
    geom_point(size = 1.5, alpha = 0.78) +
    geom_text_repel(data = pdx_tests[FDR <= 0.05 | !is.na(dataset_n)],
                    aes(label = paste(drug, TF, sep = " | ")), size = 2.0, max.overlaps = 8) +
    scale_color_manual(values = c(Validated = "#B2182B", `Cell-line replicated pair` = "#E69F00",
                                  `Other PDX screen` = "#B8BDC2")) +
    labs(title = "H  Complete PDX concordance screen",
         subtitle = sprintf("%d tests; %d validated replicated pairs", nrow(pdx_tests), nrow(pdx_validated)),
         x = "PDX concordance index - 0.5", y = expression(-log[10]("BH FDR")), color = NULL) +
    theme_pub(8.0)
} else {
  panel_5h <- ggplot() +
    annotate("rect", xmin = -1, xmax = 1, ymin = -1, ymax = 1, fill = "#F3F4F5", color = "#B8BDC2") +
    annotate("text", 0, 0.15, label = "No drug met the frozen >=10 LUAD-like PDX requirement", fontface = "bold") +
    annotate("text", 0, -0.18, label = "This is a data-coverage failure, not evidence of no biological association", size = 3.1) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + ggtitle("H  Complete PDX concordance screen")
}
save_plot("Figure5H_PDX_complete_screen", panel_5h, 6.4, 3.7)

drug_sets <- unique(rbindlist(lapply(names(cell_inputs$datasets), function(dataset) {
  p <- cell_inputs$datasets[[dataset]]$paired
  p <- p[eligible_luad == TRUE]
  data.table(dataset = dataset, canonical_drug = unique(p$canonical_drug))
})))
drug_membership <- dcast(drug_sets, canonical_drug ~ dataset, fun.aggregate = length, value.var = "dataset")
for (d in c("GDSC2", "CTRPv2", "PRISM")) {
  if (!d %in% names(drug_membership)) drug_membership[, (d) := 0L]
  drug_membership[, (d) := get(d) > 0]
}
write_table(drug_membership, file.path(tables, "SupplementaryFigure7A_drug_membership.tsv"))
panel_s7a <- upset(as.data.frame(drug_membership), intersect = c("GDSC2", "CTRPv2", "PRISM"),
                   name = "Dataset", min_size = 1, width_ratio = 0.23,
                   base_annotations = list(`Eligible drugs` = intersection_size(text = list(size = 2.8)))) +
  ggtitle("Supplementary Figure 7A | Eligible-drug overlap") +
  theme(plot.title = element_text(face = "bold", size = 10))
save_plot("SupplementaryFigure7A_drug_overlap", panel_s7a, 6.2, 4.2)

# ---------------------------------------------------------------------------
# Supplementary Data 1-9: full-resolution tables and provenance.
# ---------------------------------------------------------------------------

# Data 1: cohorts and annotations.
for (entry in list(
  list(file.path(fig1, "data", "processed", "tcga_manifest.tsv"), "TCGA discovery cohort manifest"),
  list(file.path(fig1, "data", "processed", "gse81089_manifest.tsv"), "GSE81089 validation cohort manifest"),
  list(file.path(fig1, "data", "processed", "gse41271_manifest.tsv"), "GSE41271 microarray validation manifest"),
  list(file.path(fig1, "results", "tables", "SupplementaryFigure1_inclusion_counts.tsv"), "Patient inclusion counts"),
  list(file.path(fig1, "data", "processed", "pdmr_manifest.tsv"), "PDMR PDX manifest"),
  list(file.path(fig1, "data", "processed", "depmap_manifest.tsv"), "DepMap cell-line manifest")
)) add_copy(entry[[1]], 1, entry[[2]])

# Data 2: networks and TF discovery.
for (entry in list(
  list(file.path(fig1, "results", "tcga", "TCGA_msviper.tsv"), "TCGA differential TF activity"),
  list(file.path(fig1, "results", "gse81089", "GSE81089_msviper.tsv"), "GSE81089 differential TF activity"),
  list(file.path(fig1, "results", "gse41271", "GSE41271_msviper.tsv"), "GSE41271 differential TF activity"),
  list(file.path(fig1, "results", "tables", "SupplementaryFigure2B_TF_overlap_counts.tsv"), "Independent-cohort TF overlap"),
  list(file.path(fig1, "results", "tables", "SupplementaryFigure2F_replicated_LUAD_TF_effect_zscores.tsv"), "Replicated LUAD TF effects")
)) add_copy(entry[[1]], 2, entry[[2]])

# Data 3: frozen TF selection and network.
for (entry in list(
  list(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"), "Promoter and activity gate for all discovery TFs"),
  list(file.path(fig2, "results", "promoter_gate", "system_tf_promoter_activity_summary.tsv"), "System-level promoter and activity evidence"),
  list(file.path(fig3, "results", "tables", "Figure3_HC_TF_signed_regulon_edges.tsv.gz"), "Signed 31-HC-TF regulon edges"),
  list(file.path(tables, "Figure3B_target_gene_network_nodes.tsv"), "Shared-target community nodes"),
  list(file.path(tables, "Figure3B_target_gene_network_edges.tsv"), "Shared-target community edges"),
  list(file.path(tables, "Figure3B_target_gene_community_GO.tsv"), "Shared-target community GO enrichment")
)) add_copy(entry[[1]], 3, entry[[2]])

# Data 4: ATAC QC.
for (entry in list(
  list(file.path(fig2, "audit", "manifests", "figure2_atomic", "figure2_atomic_sample_manifest.tsv"), "Figure 2 atomic sample manifest"),
  list(file.path(fig2, "audit", "raw_atac_qc", "figure2_raw_atac_qc.tsv"), "Raw ATAC QC receipt"),
  list(file.path(fig2, "results", "supplementary_figure3", "saturation_1000_permutations.tsv"), "Peak saturation permutations"),
  list(file.path(fig2, "results", "supplementary_figure3", "saturation_summary.tsv"), "Peak saturation summary"),
  list(file.path(fig2, "results", "supplementary_figure3", "sample_accessibility_pearson_correlations.tsv"), "Pairwise accessibility correlations"),
  list(file.path(fig2, "results", "supplementary_figure3", "peak_genomic_annotation.tsv"), "Peak genomic annotations")
)) add_copy(entry[[1]], 4, entry[[2]])

# Data 5: promoter and motif evidence.
for (entry in list(
  list(file.path(fig2, "results", "promoter_gate", "sample_tf_promoter_accessibility.tsv"), "Sample-level promoter accessibility"),
  list(file.path(fig2, "results", "motif_primary", "anchor_consensus_author", "anchor_consensus_sample_tf_best_motif.tsv"), "Sample-level primary motif results"),
  list(file.path(fig2, "results", "motif_primary", "anchor_consensus_author", "anchor_consensus_system_tf_motif_summary.tsv"), "System-level primary motif summary"),
  list(file.path(fig2, "results", "motif_primary", "anchor_consensus_author", "anchor_consensus_hc_tf_triple_system_summary.tsv"), "HC-TF three-system motif summary"),
  list(threshold_file, "Motif threshold sensitivity")
)) add_copy(entry[[1]], 5, entry[[2]])

# Data 6: clinical outcome.
for (entry in list(
  list(file.path(fig4, "results", "tables", "Figure4_all_Cox_models.tsv"), "All uni- and multivariable Cox models"),
  list(file.path(fig4, "results", "tables", "Figure4FH_KM_receipts.tsv"), "Kaplan-Meier receipts"),
  list(file.path(fig4, "results", "tables", "Figure4EG_endpoint_overlap.tsv"), "Endpoint-overlap table"),
  list(file.path(fig4, "results", "tables", "SupplementaryFigure7_permutation_summary.tsv"), "Permutation enrichment summary"),
  list(file.path(fig4, "results", "tables", "SupplementaryFigure7_permutation_null_counts.tsv.gz"), "Five-thousand-permutation counts")
)) add_copy(entry[[1]], 6, entry[[2]])

# Data 7: pharmacogenomics.
for (entry in list(
  list(file.path(fig5, "results", "tables", "Figure5_dataset_audit.tsv"), "Pharmacogenomic dataset audit"),
  list(file.path(fig5, "results", "tables", "Figure5_drug_identity_map.tsv"), "Drug identity mapping"),
  list(file.path(fig5, "results", "tables", "Figure5ABC_cellline_concordance_all.tsv.gz"), "All cell-line concordance tests"),
  list(file.path(fig5, "results", "tables", "Figure5D_replicated_drug_TF_pairs.tsv"), "Replicated drug-TF pairs"),
  list(file.path(tables, "Figure5_replicated_pair_dataset_details.tsv"), "Dataset-level details for replicated pairs"),
  list(file.path(tables, "SupplementaryFigure7C_all_significant_associations.tsv.gz"), "Complete significant association matrix")
)) add_copy(entry[[1]], 7, entry[[2]])

# Data 8: PDX response validation, including the complete null result.
for (entry in list(
  list(file.path(fig5, "results", "tables", "Figure5_PDXE_v2_lung_state_scores.tsv"), "PDX LUAD-state classifier scores"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure9A_PDXE_drug_coverage.tsv"), "PDX drug coverage"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure9B_PDXE_all_concordance_tests.tsv"), "All PDX concordance tests"),
  list(file.path(fig5, "results", "tables", "Figure5GH_PDX_validated_pairs.tsv"), "Validated PDX pairs; zero rows is the scientific result"),
  list(file.path(tables, "Figure5G_PDX_validation_funnel.tsv"), "PDX validation funnel")
)) add_copy(entry[[1]], 8, entry[[2]])

# Data 9: model-identity and classifier provenance.
for (entry in list(
  list(file.path(fig5, "results", "tables", "SupplementaryFigure8A_cell_mapping.tsv"), "Cell-line identifier mapping"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure8B_drug_eligibility.tsv"), "Drug eligibility audit"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure8C_lung_classifier_validation.tsv"), "Classifier validation"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure8E_PDMR_classifier_scores.tsv"), "PDX classifier scores"),
  list(file.path(fig5, "results", "tables", "SupplementaryFigure8F_DepMap_classifier_scores.tsv"), "Cell-line classifier scores")
)) add_copy(entry[[1]], 9, entry[[2]])

copy_index <- rbindlist(copy_rows, fill = TRUE)
write_table(copy_index, file.path(tables, "Supplementary_Data_1-9_file_index.tsv"))

# Machine-readable result and plotting receipt.
receipt <- as.data.table(data.frame(
  key = c(
    "frozen_discovery_TFs", "frozen_HC_TFs", "strict_three_system_motif_TFs",
    "sensitivity_q0.05_prev0.5_TFs", "target_network_nodes", "target_network_edges",
    "MKI67_cohorts", "ESTIMATE_rows", "replicated_cellline_pairs", "validated_PDX_pairs",
    "supplementary_data_files"
  ),
  value = c(
    158, length(hc), paste(strict_tf, collapse = ";"), paste(sensitivity_tf, collapse = ";"),
    nrow(target_nodes), nrow(target_pairs), uniqueN(mki67_rows$cohort), nrow(estimate_rows),
    nrow(rep_pairs), nrow(pdx_validated), nrow(copy_index)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
))
write_table(receipt, file.path(audit_root, "publication_rebuild_receipt.tsv"))
writeLines(capture.output(sessionInfo()), file.path(audit_root, "R_sessionInfo.txt"))

message("Publication rebuild atomic panels and Supplementary Data packages completed: ", out)

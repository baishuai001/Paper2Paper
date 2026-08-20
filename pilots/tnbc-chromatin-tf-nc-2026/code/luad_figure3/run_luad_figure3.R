#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(moments)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: run_luad_figure3.R FIGURE1_ROOT FIGURE2_ROOT OUTPUT_ROOT")
fig1 <- normalizePath(args[[1]], mustWork = TRUE)
fig2 <- normalizePath(args[[2]], mustWork = TRUE)
out <- normalizePath(args[[3]], mustWork = FALSE)
tables <- file.path(out, "results", "tables")
figures <- file.path(out, "results", "figures")
audit <- file.path(out, "audit")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
set.seed(20260820)

read_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  ids <- x[[1L]]
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

promoter <- fread(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
hc <- sort(promoter[
  analysis_definition == "PRIMARY_anchor_promoter_tcga_cpm1_both" &
    HC_TF_promoter_activity_definition == TRUE,
  unique(TF)
])
if (length(hc) != 31L) stop(sprintf("Expected 31 frozen HC-TFs, found %d", length(hc)))

motif_sensitivity <- fread(file.path(
  fig2, "results", "motif_threshold_sensitivity",
  "motif_threshold_sensitivity_triple_system_tfs.tsv"
))
triple_motif <- sort(unique(motif_sensitivity[scenario == "primary_anchor", TF]))
if (!identical(triple_motif, c("FOXA3", "NFATC4", "XBP1"))) {
  stop(sprintf("Unexpected primary three-system motif set: %s", paste(triple_motif, collapse = ",")))
}

regulon <- readRDS(file.path(fig1, "results", "tcga", "TCGA_regulon.rds"))
activity <- read_matrix(file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz"))
expression <- read_matrix(file.path(fig1, "data", "processed", "tcga_aracne_expression.tsv"))
manifest <- fread(file.path(fig1, "data", "processed", "tcga_manifest.tsv"))
manifest <- manifest[match(colnames(activity), sample_id)]
if (anyNA(manifest$sample_id) || !identical(manifest$sample_id, colnames(activity))) {
  stop("TCGA activity/manifest mismatch")
}
missing_regulon <- setdiff(hc, names(regulon))
missing_activity <- setdiff(hc, rownames(activity))
if (length(missing_regulon) || length(missing_activity)) {
  stop(sprintf("HC-TF missing from regulon/activity: %s / %s",
               paste(missing_regulon, collapse = ","), paste(missing_activity, collapse = ",")))
}

edge_list <- rbindlist(lapply(hc, function(tf) {
  object <- regulon[[tf]]
  targets <- names(object$tfmode)
  likelihood <- object$likelihood[targets]
  data.table(
    TF = tf,
    target = targets,
    tfmode = as.numeric(object$tfmode[targets]),
    likelihood = as.numeric(likelihood),
    direction = ifelse(object$tfmode[targets] >= 0, "Activated", "Repressed")
  )
}))
edge_list <- edge_list[is.finite(tfmode) & !is.na(target) & target != ""]
edge_list[, target_HC_TF_count := uniqueN(TF), by = target]
edge_list[, target_class := ifelse(target_HC_TF_count >= 2L, "Shared", "Private")]
edge_list[, motif_three_system := TF %in% triple_motif]
fwrite(edge_list, file.path(tables, "Figure3_HC_TF_signed_regulon_edges.tsv.gz"), sep = "\t")

composition <- edge_list[, .(targets = uniqueN(target)),
                         by = .(TF, direction, target_class, motif_three_system)]
complete_grid <- CJ(TF = hc, direction = c("Activated", "Repressed"),
                    target_class = c("Shared", "Private"), unique = TRUE)
composition <- merge(complete_grid, composition,
                     by = c("TF", "direction", "target_class"), all.x = TRUE)
composition[is.na(targets), targets := 0L]
composition[, motif_three_system := TF %in% triple_motif]
composition[, component := factor(paste(direction, target_class, sep = " / "),
                                  levels = c("Activated / Shared", "Activated / Private",
                                             "Repressed / Shared", "Repressed / Private"))]
totals <- composition[, .(regulon_targets = sum(targets)), by = TF]
setorder(totals, regulon_targets)
composition[, TF := factor(TF, levels = totals$TF)]
fwrite(composition, file.path(tables, "Figure3A_regulon_composition.tsv"), sep = "\t")

panel_a <- ggplot(composition, aes(TF, targets, fill = component)) +
  geom_col(width = 0.78) +
  coord_flip() +
  scale_fill_manual(values = c(
    `Activated / Shared` = "#A50F15", `Activated / Private` = "#FB6A4A",
    `Repressed / Shared` = "#08519C", `Repressed / Private` = "#6BAED6"
  )) +
  geom_point(data = unique(composition[motif_three_system == TRUE, .(TF)]),
             aes(TF, y = 0), inherit.aes = FALSE, shape = 8, size = 2.3, color = "#111111") +
  labs(title = "A  Signed regulon composition", subtitle = "Star: three-system motif support",
       x = NULL, y = "Unique ARACNe3 targets", fill = NULL) +
  theme_classic(base_size = 8.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(figures, "Figure3A_regulon_composition.pdf"), panel_a,
       width = 6.4, height = 8.2, useDingbats = FALSE)
ggsave(file.path(figures, "Figure3A_regulon_composition.png"), panel_a,
       width = 6.4, height = 8.2, dpi = 220)

activated <- split(unique(edge_list[direction == "Activated", .(TF, target)])$target,
                   unique(edge_list[direction == "Activated", .(TF, target)])$TF)
activated <- lapply(hc, function(tf) unique(edge_list[TF == tf & direction == "Activated", target])) |>
  setNames(hc)
pairs <- as.data.table(t(combn(hc, 2L)))
setnames(pairs, c("TF1", "TF2"))
pairs[, `:=`(
  shared_targets = lengths(Map(intersect, activated[TF1], activated[TF2])),
  union_targets = lengths(Map(union, activated[TF1], activated[TF2]))
)]
pairs[, jaccard := fifelse(union_targets > 0, shared_targets / union_targets, 0)]
fwrite(pairs, file.path(tables, "Figure3_TF_activated_target_pairwise_overlap.tsv"), sep = "\t")

positive_pairs <- pairs[shared_targets > 0]
graph <- graph_from_data_frame(
  positive_pairs[, .(from = TF1, to = TF2, weight = shared_targets, jaccard)],
  directed = FALSE, vertices = data.frame(name = hc)
)
if (ecount(graph) > 0L) {
  communities <- cluster_louvain(graph, weights = E(graph)$weight)
  membership_table <- data.table(TF = names(membership(communities)),
                                 module = as.integer(membership(communities)))
} else {
  membership_table <- data.table(TF = hc, module = seq_along(hc))
}
setorder(membership_table, module, TF)
fwrite(membership_table, file.path(tables, "Figure3B_HC_TF_modules.tsv"), sep = "\t")

universe <- unique(edge_list$target)
module_enrichment <- rbindlist(lapply(unique(membership_table$module), function(module_id) {
  module_tfs <- membership_table[module == module_id, TF]
  genes <- unique(edge_list[TF %in% module_tfs & direction == "Activated", target])
  enriched <- tryCatch(
    suppressMessages(enrichGO(
      gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
      readable = FALSE
    )), error = function(e) NULL
  )
  answer <- if (is.null(enriched)) data.table() else as.data.table(enriched@result)
  if (!nrow(answer)) {
    return(data.table(module = module_id, ID = NA_character_, Description = "No GO BP term",
                      GeneRatio = NA_character_, BgRatio = NA_character_, pvalue = NA_real_,
                      p.adjust = NA_real_, qvalue = NA_real_, geneID = NA_character_, Count = 0L,
                      module_TFs = length(module_tfs), module_activated_targets = length(genes)))
  }
  answer[, module := module_id]
  answer[, `:=`(module_TFs = length(module_tfs), module_activated_targets = length(genes))]
  answer
}), fill = TRUE)
setorder(module_enrichment, module, p.adjust, -Count)
fwrite(module_enrichment, file.path(tables, "Figure3B_module_GO_BP_enrichment.tsv"), sep = "\t")
top_go <- module_enrichment[, head(.SD, 3L), by = module]
top_go[, label := ifelse(nchar(Description) > 54L,
                         paste0(substr(Description, 1L, 51L), "..."), Description)]
top_go[, logFDR := -log10(pmax(p.adjust, 1e-300))]
panel_b <- ggplot(top_go[!is.na(ID)], aes(logFDR, reorder(label, logFDR),
                                          size = Count, color = factor(module))) +
  geom_point(alpha = 0.85) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  labs(title = "B  HC-TF communities and GO programs", x = expression(-log[10]("BH FDR")),
       y = NULL, color = "Module", size = "Genes") +
  theme_classic(base_size = 8.2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        strip.background = element_blank())
ggsave(file.path(figures, "Figure3B_module_GO_enrichment.pdf"), panel_b,
       width = 7.0, height = 6.2, useDingbats = FALSE)
ggsave(file.path(figures, "Figure3B_module_GO_enrichment.png"), panel_b,
       width = 7.0, height = 6.2, dpi = 220)

luad_samples <- manifest[group == "LUAD", sample_id]
lusc_samples <- manifest[group == "LUSC", sample_id]
luad_activity <- activity[hc, luad_samples, drop = FALSE]
cor_matrix <- cor(t(luad_activity), method = "pearson", use = "pairwise.complete.obs")
cor_long <- as.data.table(as.table(cor_matrix))
setnames(cor_long, c("TF1", "TF2", "pearson_r"))
fwrite(cor_long, file.path(tables, "Figure3C_HC_TF_activity_correlations.tsv.gz"), sep = "\t")
cor_heatmap <- Heatmap(
  cor_matrix, name = "Pearson r",
  col = colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE,
  show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6),
  column_title = sprintf("C  HC-TF activity correlation in TCGA-LUAD (n=%d)", length(luad_samples)),
  column_title_gp = grid::gpar(fontface = "bold", fontsize = 10)
)
cairo_pdf(file.path(figures, "Figure3C_HC_TF_activity_correlation.pdf"), width = 7.4, height = 7.2)
draw(cor_heatmap)
dev.off()
png(file.path(figures, "Figure3C_HC_TF_activity_correlation.png"), width = 1628, height = 1584, res = 220)
draw(cor_heatmap)
dev.off()

membership_table[, color := grDevices::hcl.colors(uniqueN(module), "Dark 3")[module]]
V(graph)$color <- membership_table$color[match(V(graph)$name, membership_table$TF)]
V(graph)$size <- 7 + 2.4 * log1p(totals$regulon_targets[match(V(graph)$name, totals$TF)])
if (ecount(graph) > 0L) {
  display_cut <- as.numeric(quantile(E(graph)$weight, 0.80, names = FALSE, type = 1))
  display_graph <- delete_edges(graph, E(graph)[weight < display_cut])
} else {
  display_cut <- NA_real_
  display_graph <- graph
}
network_layout <- layout_with_fr(display_graph, weights = if (ecount(display_graph)) E(display_graph)$weight else NULL)
plot_network <- function(path, png_device = FALSE) {
  if (png_device) png(path, width = 1600, height = 1400, res = 220) else cairo_pdf(path, width = 7.2, height = 6.4)
  par(mar = c(0.5, 0.5, 2.2, 0.5))
  plot(display_graph, layout = network_layout, vertex.label = V(display_graph)$name,
       vertex.label.cex = 0.55, vertex.label.color = "#202020",
       vertex.frame.color = "white", vertex.color = V(display_graph)$color,
       vertex.size = V(display_graph)$size,
       edge.color = adjustcolor("#526D82", alpha.f = 0.35),
       edge.width = if (ecount(display_graph)) 0.5 + 2.5 * E(display_graph)$weight / max(E(display_graph)$weight) else 1)
  title("D  Shared activated-target collaboration network", adj = 0, font.main = 2, cex.main = 1.0)
  # Keep the panel itself uncluttered; the display threshold and the fact that
  # communities use the full positive-overlap graph are reported in the legend.
  # The previous bottom margin note was clipped when the vector panel was
  # assembled into the composite Figure 3.
  dev.off()
}
plot_network(file.path(figures, "Figure3D_TF_collaboration_network.pdf"), FALSE)
plot_network(file.path(figures, "Figure3D_TF_collaboration_network.png"), TRUE)

partners <- rbind(
  positive_pairs[, .(TF = TF1, partner = TF2)],
  positive_pairs[, .(TF = TF2, partner = TF1)]
)[, .(partners = uniqueN(partner)), by = TF]
partners <- merge(data.table(TF = hc), partners, by = "TF", all.x = TRUE)
partners[is.na(partners), partners := 0L]
partners <- merge(partners, totals, by = "TF")
partners <- merge(partners, membership_table[, .(TF, module)], by = "TF")
fwrite(partners, file.path(tables, "Figure3E_regulon_size_and_partners.tsv"), sep = "\t")
panel_e <- ggplot(partners, aes(regulon_targets, partners, color = factor(module), label = TF)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#444444", linewidth = 0.45) +
  geom_text_repel(size = 2.2, max.overlaps = 12) +
  labs(title = "E  Regulon breadth and TF collaboration", x = "Regulon targets",
       y = "HC-TF partners", color = "Module") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(figures, "Figure3E_regulon_size_partners.pdf"), panel_e,
       width = 5.3, height = 4.5, useDingbats = FALSE)
ggsave(file.path(figures, "Figure3E_regulon_size_partners.png"), panel_e,
       width = 5.3, height = 4.5, dpi = 220)

skew_table <- rbindlist(lapply(hc, function(tf) {
  values <- as.numeric(luad_activity[tf, ])
  sk <- moments::skewness(values)
  z <- sk / sqrt(6 / length(values))
  data.table(TF = tf, n = length(values), mean_activity = mean(values),
             sd_activity = sd(values), skewness = sk, z = z,
             p_value = 2 * pnorm(-abs(z)))
}))
skew_table[, FDR := p.adjust(p_value, method = "BH")]
skew_table[, strong_skew := FDR <= 0.05 & abs(skewness) > 1]
skew_table[, motif_three_system := TF %in% triple_motif]
fwrite(skew_table, file.path(tables, "Figure3F_TCGA_LUAD_HC_TF_skewness.tsv"), sep = "\t")
panel_f <- ggplot(skew_table, aes(mean_activity, skewness, color = strong_skew, shape = motif_three_system)) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "#969696") +
  geom_point(size = 2.4) +
  geom_text_repel(data = skew_table[strong_skew == TRUE | motif_three_system == TRUE],
                  aes(label = TF), size = 2.4, max.overlaps = Inf) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#969696")) +
  labs(title = "F  Intertumour HC-TF activity heterogeneity", x = "Mean VIPER activity",
       y = "Skewness", color = "Strong skew", shape = "3-system motif") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(figures, "Figure3F_HC_TF_activity_skewness.pdf"), panel_f,
       width = 5.3, height = 4.5, useDingbats = FALSE)
ggsave(file.path(figures, "Figure3F_HC_TF_activity_skewness.png"), panel_f,
       width = 5.3, height = 4.5, dpi = 220)

representatives <- unique(c(
  skew_table[order(-skewness), head(TF, 2L)],
  skew_table[order(skewness), head(TF, 2L)]
))
density_table <- rbindlist(lapply(representatives, function(tf) {
  data.table(TF = tf, sample_id = manifest$sample_id, histology = manifest$group,
             activity = as.numeric(activity[tf, manifest$sample_id]))
}))
fwrite(density_table, file.path(tables, "Figure3G_representative_activity_distributions.tsv.gz"), sep = "\t")
panel_g <- ggplot(density_table, aes(activity, color = histology, fill = histology)) +
  geom_density(alpha = 0.16, linewidth = 0.65) +
  facet_wrap(~TF, scales = "free", ncol = 2) +
  scale_color_manual(values = c(LUAD = "#B2182B", LUSC = "#2166AC")) +
  scale_fill_manual(values = c(LUAD = "#B2182B", LUSC = "#2166AC")) +
  labs(title = "G  Representative skewed TF activity distributions", x = "VIPER activity",
       y = "Density", color = NULL, fill = NULL) +
  theme_classic(base_size = 8.5) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(figures, "Figure3G_representative_distributions.pdf"), panel_g,
       width = 5.8, height = 4.8, useDingbats = FALSE)
ggsave(file.path(figures, "Figure3G_representative_distributions.png"), panel_g,
       width = 5.8, height = 4.8, dpi = 220)

# Supplementary heterogeneity: proliferation and model-system skewness.
if (!"MKI67" %in% rownames(expression)) stop("MKI67 missing from TCGA expression matrix")
mki67 <- as.numeric(expression["MKI67", luad_samples])
mki67_cor <- rbindlist(lapply(hc, function(tf) {
  test <- cor.test(as.numeric(luad_activity[tf, ]), mki67, method = "pearson")
  data.table(TF = tf, pearson_r = unname(test$estimate), p_value = test$p.value)
}))
mki67_cor[, FDR := p.adjust(p_value, method = "BH")]
fwrite(mki67_cor, file.path(tables, "SupplementaryFigure6A_MKI67_correlations.tsv"), sep = "\t")
mki67_cor[, TF := factor(TF, levels = TF[order(pearson_r)])]
supp_a <- ggplot(mki67_cor, aes(TF, pearson_r, fill = FDR <= 0.05)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#BDBDBD")) +
  labs(title = "A  HC-TF activity association with MKI67", x = NULL,
       y = "Pearson r", fill = "BH FDR <= 0.05") +
  theme_classic(base_size = 8.5) + theme(legend.position = "bottom")
ggsave(file.path(figures, "SupplementaryFigure6A_MKI67_correlations.pdf"), supp_a,
       width = 5.6, height = 6.5, useDingbats = FALSE)

model_skew <- list(TCGA_LUAD = skew_table[, .(system = "TCGA_LUAD", TF, n, skewness, FDR)])
model_specs <- list(
  PDMR_PDX = list(activity = file.path(fig1, "results", "pdmr", "PDMR_viper_activity.tsv.gz"),
                  manifest = file.path(fig1, "data", "processed", "pdmr_manifest.tsv")),
  DepMap = list(activity = file.path(fig1, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
                manifest = file.path(fig1, "data", "processed", "depmap_manifest.tsv"))
)
for (system_name in names(model_specs)) {
  matrix <- read_matrix(model_specs[[system_name]]$activity)
  model_manifest <- fread(model_specs[[system_name]]$manifest)
  model_manifest <- model_manifest[match(colnames(matrix), sample_id)]
  ids <- model_manifest[group == "LUAD", sample_id]
  model_skew[[system_name]] <- rbindlist(lapply(hc, function(tf) {
    values <- as.numeric(matrix[tf, ids])
    sk <- moments::skewness(values)
    z <- sk / sqrt(6 / length(values))
    data.table(system = system_name, TF = tf, n = length(values), skewness = sk,
               p_value = 2 * pnorm(-abs(z)))
  }))
  model_skew[[system_name]][, FDR := p.adjust(p_value, method = "BH")]
}
model_skew_table <- rbindlist(model_skew, fill = TRUE)
fwrite(model_skew_table, file.path(tables, "SupplementaryFigure6B_cross_system_skewness.tsv"), sep = "\t")
model_skew_table[, TF := factor(TF, levels = skew_table[order(skewness), TF])]
supp_b <- ggplot(model_skew_table, aes(system, TF, fill = skewness)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "B  HC-TF activity skewness across LUAD systems", x = NULL, y = NULL) +
  theme_classic(base_size = 8.5) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(figures, "SupplementaryFigure6B_cross_system_skewness.pdf"), supp_b,
       width = 5.3, height = 7.0, useDingbats = FALSE)

receipt <- list(
  status = "passed",
  frozen_HC_TFs = length(hc),
  three_system_motif_TFs = triple_motif,
  signed_regulon_edges = nrow(edge_list),
  unique_regulon_targets = uniqueN(edge_list$target),
  activated_overlap_edges_all = nrow(positive_pairs),
  activated_overlap_display_threshold = display_cut,
  activated_overlap_edges_displayed = ecount(display_graph),
  modules = uniqueN(membership_table$module),
  TCGA_LUAD_samples = length(luad_samples),
  TCGA_LUSC_reference_samples = length(lusc_samples),
  strong_skew_TFs = skew_table[strong_skew == TRUE, TF],
  package_versions = list(
    R = as.character(getRversion()), data.table = as.character(packageVersion("data.table")),
    igraph = as.character(packageVersion("igraph")), viper = as.character(packageVersion("viper")),
    clusterProfiler = as.character(packageVersion("clusterProfiler"))
  )
)
write_json(receipt, file.path(audit, "figure3_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

#!/usr/bin/env Rscript

# Render the LUAD Figure 3 network panels with the graph construction and
# visual mappings used by the TNBC Code Ocean v1.0 Figure3/05 script.
#
# This wrapper consumes only the frozen LUAD Figure 3 tables. It does not rerun
# ARACNe3, VIPER, or any upstream statistical analysis.

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)
  library(dplyr)
  library(ggforce)
  library(ggplot2)
  library(ggrepel)
  library(ggraph)
  library(grid)
  library(igraph)
  library(patchwork)
  library(tidygraph)
  library(viridisLite)
})

argv_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", argv_full, value = TRUE)
if (!length(script_arg)) {
  stop("Cannot resolve the script directory from commandArgs().", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
source(file.path(script_dir, "tnbc_anchor_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: render_anchor_figure3.R FIGURE3_ROOT FIGURE1_ROOT FIGURE2_ROOT OUTPUT_ROOT",
      "Writes PDFs/PNGs to OUTPUT_ROOT/anchor_style/atomic."
    ),
    call. = FALSE
  )
}

fig3_root <- normalizePath(args[[1L]], mustWork = TRUE)
fig1_root <- normalizePath(args[[2L]], mustWork = TRUE)
fig2_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_root <- args[[4L]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)

atomic_dir <- file.path(output_root, "anchor_style", "atomic")
table_dir <- file.path(output_root, "anchor_style", "tables")
audit_dir <- file.path(output_root, "anchor_style", "audit")
for (path in c(atomic_dir, table_dir, audit_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

set.seed(2026)

save_anchor_atomic <- function(plot, stem, width, height) {
  pdf_path <- file.path(atomic_dir, paste0(stem, ".pdf"))
  png_path <- file.path(atomic_dir, paste0(stem, ".png"))
  anchor_save_plot(plot, pdf_path, width = width, height = height)
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
  data.table(
    artifact = basename(pdf_path),
    pdf_path = pdf_path,
    png_path = png_path,
    width_in = width,
    height_in = height
  )
}

save_anchor_heatmap <- function(heatmap, stem, width, height, tag = NULL) {
  pdf_path <- file.path(atomic_dir, paste0(stem, ".pdf"))
  png_path <- file.path(atomic_dir, paste0(stem, ".png"))
  draw_once <- function() {
    draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
    if (!is.null(tag)) {
      grid.text(
        tag, x = unit(0.01, "npc"), y = unit(0.995, "npc"),
        just = c("left", "top"),
        gp = gpar(family = anchor_font_family, fontsize = 16, fontface = "bold")
      )
    }
  }
  anchor_cairo_pdf(pdf_path, width = width, height = height)
  draw_once()
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 400)
  draw_once()
  dev.off()
  data.table(
    artifact = basename(pdf_path), pdf_path = pdf_path, png_path = png_path,
    width_in = width, height_in = height
  )
}

read_numeric_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  ids <- x[[1L]]
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

signed_edge_path <- file.path(
  fig3_root, "results", "tables", "Figure3_HC_TF_signed_regulon_edges.tsv.gz"
)
if (!file.exists(signed_edge_path)) {
  stop(sprintf("Frozen Figure 3 edge table not found: %s", signed_edge_path), call. = FALSE)
}

signed_edges <- fread(signed_edge_path)
anchor_assert_columns(
  signed_edges,
  c("TF", "target", "direction", "target_class"),
  "signed_edges"
)

hc_tfs <- sort(unique(signed_edges$TF))
if (length(hc_tfs) != 31L) {
  stop(sprintf("Expected 31 frozen LUAD HC-TFs; found %d.", length(hc_tfs)), call. = FALSE)
}

# The official script first restricts to Activated Shared Targets. A shared
# target is present in at least two HC-TF regulons; the pairwise matrix below
# then determines how many HC-TFs activate both members of each target pair.
activated_shared <- unique(
  signed_edges[direction == "Activated" & target_class == "Shared", .(TF, target)]
)
if (!nrow(activated_shared)) {
  stop("No activated shared targets were recovered.", call. = FALSE)
}

binary_incidence <- function(edges, row_ids, column_ids, row_col, column_col) {
  answer <- matrix(
    0L,
    nrow = length(row_ids),
    ncol = length(column_ids),
    dimnames = list(row_ids, column_ids)
  )
  row_index <- match(edges[[row_col]], row_ids)
  column_index <- match(edges[[column_col]], column_ids)
  answer[cbind(row_index, column_index)] <- 1L
  answer
}

# ---------------------------------------------------------------------------
# Figure 3B: official-equivalent target-gene network.
# ---------------------------------------------------------------------------

target_ids <- sort(unique(activated_shared$target))
tf_by_target <- binary_incidence(
  activated_shared,
  row_ids = hc_tfs,
  column_ids = target_ids,
  row_col = "TF",
  column_col = "target"
)
shared_tf_counts <- t(tf_by_target) %*% tf_by_target
target_edge_index <- which(
  shared_tf_counts >= 3L & row(shared_tf_counts) < col(shared_tf_counts),
  arr.ind = TRUE
)
target_edges <- data.table(
  from = rownames(shared_tf_counts)[target_edge_index[, 1L]],
  to = colnames(shared_tf_counts)[target_edge_index[, 2L]],
  weight = as.integer(shared_tf_counts[target_edge_index])
)

# Deliberately omit a vertices argument, matching the official script: only
# endpoints of qualifying edges belong to the plotted network. Targets with no
# qualifying edge remain auditable in the frozen table but are not scattered
# across the network canvas as artificial isolates.
target_graph <- graph_from_data_frame(target_edges, directed = FALSE)
if (vcount(target_graph) != 6L || ecount(target_graph) != 3L) {
  stop(
    sprintf(
      "Frozen anchor-equivalent Figure 3B should contain 6 nodes and 3 edges; found %d/%d.",
      vcount(target_graph), ecount(target_graph)
    ),
    call. = FALSE
  )
}

target_community <- cluster_louvain(target_graph, weights = E(target_graph)$weight)
target_degree <- degree(target_graph)
target_membership <- membership(target_community)
target_tbl <- as_tbl_graph(target_graph) |>
  activate(nodes) |>
  mutate(
    degree = unname(target_degree[name]),
    community = factor(unname(target_membership[name]))
  )

set.seed(2026)
panel_3b_audit <- ggraph(target_tbl, layout = "fr") +
  geom_mark_hull(
    aes(x = x, y = y, group = community, fill = community),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.15,
    radius = unit(2, "mm"),
    show.legend = FALSE,
    colour = NA
  ) +
  geom_edge_link(aes(alpha = weight), show.legend = FALSE) +
  geom_node_point(
    aes(colour = community, size = degree),
    show.legend = c(colour = TRUE, size = TRUE)
  ) +
  geom_node_text(
    aes(label = name),
    family = anchor_font_family,
    size = 3.5,
    repel = TRUE,
    max.overlaps = 20
  ) +
  scale_colour_viridis_d(option = "turbo", name = "Community") +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  labs(size = "Number of\nConnections", tag = "B") +
  theme_graph(base_family = anchor_font_family) +
  anchor_panel_tag_theme +
  theme(legend.position = "right")

artifact_rows <- list(save_anchor_atomic(
  panel_3b_audit,
  "Figure3B_anchor_equivalent_target_network",
  width = 10,
  height = 8
))

# The strict graph contains three disconnected dyads and therefore cannot
# support GO/community labels. This compact companion makes the negative
# topology interpretable without lowering the frozen >=3 shared-TF threshold.
target_regulators <- split(activated_shared$TF, activated_shared$target)
dyad_rows <- rbindlist(lapply(seq_len(nrow(target_edges)), function(index) {
  left <- target_edges$from[[index]]
  right <- target_edges$to[[index]]
  shared_tfs <- sort(intersect(target_regulators[[left]], target_regulators[[right]]))
  data.table(
    dyad = factor(index, levels = rev(seq_len(nrow(target_edges)))),
    left_target = left,
    right_target = right,
    shared_tfs = paste(shared_tfs, collapse = "  |  "),
    shared_tf_n = length(shared_tfs)
  )
}))

panel_3b_dyad <- ggplot(dyad_rows, aes(y = dyad, colour = dyad)) +
  geom_segment(
    aes(x = 0, xend = 4, yend = dyad),
    linewidth = 1.05,
    colour = "grey70",
    lineend = "round"
  ) +
  geom_point(aes(x = 0), size = 5.2, show.legend = FALSE) +
  geom_point(aes(x = 4), size = 5.2, show.legend = FALSE) +
  geom_text(
    aes(x = 0, label = left_target),
    family = anchor_font_family,
    colour = "black",
    fontface = "bold",
    size = 3.7,
    nudge_y = 0.22
  ) +
  geom_text(
    aes(x = 4, label = right_target),
    family = anchor_font_family,
    colour = "black",
    fontface = "bold",
    size = 3.7,
    nudge_y = 0.22
  ) +
  geom_label(
    aes(x = 2, label = paste0("Shared HC-TFs:  ", shared_tfs)),
    family = anchor_font_family,
    colour = "#303236",
    fill = "white",
    label.size = 0.25,
    label.padding = unit(0.16, "lines"),
    size = 3.15,
    show.legend = FALSE
  ) +
  scale_colour_viridis_d(option = "turbo", guide = "none") +
  scale_x_continuous(limits = c(-0.65, 4.65), expand = expansion(mult = 0)) +
  labs(
    title = "Anchor-threshold target-gene dyads",
    subtitle = "Three qualifying pairs; no connected component contains three or more genes",
    tag = "B",
    x = NULL,
    y = NULL
  ) +
  anchor_theme_void(base_size = 12) +
  anchor_panel_tag_theme +
  theme(
    plot.title = element_text(
      family = anchor_font_family, face = "bold", size = 13,
      margin = margin(b = 2, l = 28)
    ),
    plot.subtitle = element_text(
      family = anchor_font_family, size = 10, colour = "#55595C",
      margin = margin(b = 8, l = 28)
    )
  )

artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3b_dyad,
  "Figure3B_compact_dyad_explanation",
  width = 10,
  height = 4.6
)

# ---------------------------------------------------------------------------
# Figure 3D: official-equivalent HC-TF collaboration network.
# ---------------------------------------------------------------------------

target_ids_for_tf_graph <- sort(unique(activated_shared$target))
tf_by_target_for_tf_graph <- binary_incidence(
  activated_shared,
  row_ids = hc_tfs,
  column_ids = target_ids_for_tf_graph,
  row_col = "TF",
  column_col = "target"
)
shared_target_counts <- tf_by_target_for_tf_graph %*% t(tf_by_target_for_tf_graph)
tf_edge_index <- which(
  shared_target_counts >= 1L & row(shared_target_counts) < col(shared_target_counts),
  arr.ind = TRUE
)
tf_edges <- data.table(
  from = rownames(shared_target_counts)[tf_edge_index[, 1L]],
  to = colnames(shared_target_counts)[tf_edge_index[, 2L]],
  weight = as.integer(shared_target_counts[tf_edge_index])
)
tf_graph <- graph_from_data_frame(tf_edges, directed = FALSE)

if (vcount(tf_graph) != 31L || ecount(tf_graph) != 134L) {
  stop(
    sprintf(
      "Frozen anchor-equivalent Figure 3D should contain 31 nodes and 134 edges; found %d/%d.",
      vcount(tf_graph), ecount(tf_graph)
    ),
    call. = FALSE
  )
}
if (components(tf_graph)$no != 1L || any(degree(tf_graph) == 0L)) {
  stop("The frozen Figure 3D graph must be connected and contain no isolates.", call. = FALSE)
}

tf_degree <- degree(tf_graph)
tf_tbl <- as_tbl_graph(tf_graph) |>
  activate(nodes) |>
  mutate(degree = unname(tf_degree[name]))

set.seed(2026)
panel_3d <- ggraph(tf_tbl, layout = "fr") +
  geom_edge_link(
    aes(width = weight),
    colour = "grey70",
    alpha = 0.6
  ) +
  geom_node_point(
    aes(size = degree),
    colour = unname(anchor_colours$regulon[["activated_shared"]])
  ) +
  geom_node_text(
    aes(label = name),
    family = anchor_font_family,
    colour = "black",
    size = 3.5,
    repel = TRUE,
    max.overlaps = 20
  ) +
  scale_edge_width(
    range = c(0.2, 2),
    name = "Shared\nTargets"
  ) +
  scale_size(
    range = c(1, 6),
    name = "Number of\nConnections"
  ) +
  labs(tag = "D") +
  theme_graph(base_family = anchor_font_family) +
  anchor_panel_tag_theme +
  theme(legend.position = "right")

artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3d,
  "Figure3D_anchor_equivalent_TF_collaboration_network",
  width = 10,
  height = 8
)

# ---------------------------------------------------------------------------
# Figure 3A: signed regulon composition and net effect.
# ---------------------------------------------------------------------------

composition_edges <- unique(signed_edges[, .(TF, target, direction, target_class)])
composition_edges[, label := fcase(
  direction == "Activated" & target_class == "Shared", "Activated Shared Targets",
  direction == "Activated" & target_class != "Shared", "Activated Non-Shared Targets",
  direction == "Repressed" & target_class == "Shared", "Repressed Shared Targets",
  default = "Repressed Non-Shared Targets"
)]
composition_levels <- c(
  "Activated Shared Targets", "Activated Non-Shared Targets",
  "Repressed Shared Targets", "Repressed Non-Shared Targets"
)
composition_edges[, label := factor(label, levels = rev(composition_levels))]
composition_counts <- composition_edges[, .(targets = uniqueN(target)), by = .(TF, label)]
composition_total <- composition_counts[, .(regulon_size = sum(targets)), by = TF]
setorder(composition_total, regulon_size, TF)
tf_order_a <- composition_total$TF
composition_counts[, TF := factor(TF, levels = tf_order_a)]

motif_system_path <- file.path(
  fig2_root, "results", "motif_primary", "anchor_consensus_author",
  "anchor_consensus_system_tf_motif_summary.tsv"
)
if (!file.exists(motif_system_path)) {
  stop(sprintf("Figure 2 motif-system table not found: %s", motif_system_path), call. = FALSE)
}
motif_system <- fread(motif_system_path)
anchor_assert_columns(
  motif_system,
  c("system", "TF", "support_fixed_original_denominator"),
  "motif_system"
)
motif_system[, support_fixed_original_denominator :=
  tolower(as.character(support_fixed_original_denominator)) %in% c("true", "t", "1")]

# This exactly mirrors the anchor Figure3/04 logic: the motif tile marks HC-TFs
# whose motif passes the >=50%-of-samples rule in at least one of the three
# systems. Triple-system conservation remains a separate Figure 2 annotation.
motif_annotation <- motif_system[TF %in% hc_tfs, .(
  motif = any(support_fixed_original_denominator, na.rm = TRUE)
), by = TF]
motif_annotation <- merge(
  data.table(TF = hc_tfs), motif_annotation, by = "TF", all.x = TRUE
)
motif_annotation[is.na(motif), motif := FALSE]
motif_annotation[, TF := factor(TF, levels = tf_order_a)]
net_effect <- composition_edges[, .(
  net_score = uniqueN(target[direction == "Activated"]) / uniqueN(target) -
    uniqueN(target[direction == "Repressed"]) / uniqueN(target)
), by = TF]
net_effect[, TF := factor(TF, levels = tf_order_a)]

motif_tile <- ggplot(motif_annotation, aes("Motif", TF, fill = motif)) +
  geom_tile() +
  scale_fill_manual(values = c(`TRUE` = "#070926", `FALSE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12, base_family = anchor_font_family) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 4.5, colour = "black"),
    panel.grid = element_blank(), plot.margin = margin(5.5, 0, 5.5, 5.5)
  )

bar_limit <- max(10, ceiling(max(composition_total$regulon_size) / 25) * 25)
regulon_bar <- ggplot(composition_counts, aes(targets, TF, fill = label)) +
  geom_col(width = 0.84) +
  geom_text(
    aes(label = ifelse(targets > 0, targets, "")),
    position = position_stack(vjust = 0.5), family = anchor_font_family,
    colour = "white", size = 2
  ) +
  scale_fill_manual(values = c(
    `Activated Shared Targets` = unname(anchor_colours$regulon[["activated_shared"]]),
    `Activated Non-Shared Targets` = unname(anchor_colours$regulon[["activated_private"]]),
    `Repressed Shared Targets` = unname(anchor_colours$regulon[["repressed_shared"]]),
    `Repressed Non-Shared Targets` = unname(anchor_colours$regulon[["repressed_private"]])
  ), drop = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, bar_limit)) +
  labs(tag = "A", x = "Number of interactions in the LUAD regulatory network", y = NULL, fill = NULL) +
  anchor_theme_bw(10) +
  anchor_panel_tag_theme +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8, colour = "black"),
    legend.position = "right", plot.margin = margin(5.5, 0, 5.5, 0)
  )

net_tile <- ggplot(net_effect, aes("Net effect", TF, fill = net_score)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2C3E6B", mid = "white", high = "#8B1A1A",
    midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
    labels = c("Repressing", "Equal", "Activating"), name = "Net effect"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = anchor_font_family) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

panel_3a <- motif_tile + regulon_bar + net_tile +
  plot_layout(widths = c(0.07, 1, 0.08), guides = "collect") &
  theme(legend.position = "right")
artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3a, "Figure3A_anchor_signed_regulon_composition", 5.75, 5.25
)

# ---------------------------------------------------------------------------
# Figure 3C: HC-TF activity correlation and top/bottom TF pairs.
# ---------------------------------------------------------------------------

correlation_path <- file.path(
  fig3_root, "results", "tables", "Figure3C_HC_TF_activity_correlations.tsv.gz"
)
correlation_long <- fread(correlation_path)
anchor_assert_columns(correlation_long, c("TF1", "TF2", "pearson_r"), "correlation_long")
correlation_wide <- dcast(correlation_long, TF1 ~ TF2, value.var = "pearson_r")
correlation_ids <- correlation_wide$TF1
correlation_matrix <- as.matrix(correlation_wide[, -"TF1"])
storage.mode(correlation_matrix) <- "double"
rownames(correlation_matrix) <- correlation_ids
common_ids <- intersect(rownames(correlation_matrix), colnames(correlation_matrix))
correlation_matrix <- correlation_matrix[common_ids, common_ids, drop = FALSE]
if (!identical(sort(common_ids), sort(hc_tfs))) {
  stop("Figure 3C correlation matrix does not contain the frozen 31 HC-TFs.", call. = FALSE)
}

correlation_heatmap <- Heatmap(
  correlation_matrix,
  name = "NES Pearson\ncorrelation",
  col = viridisLite::magma(3, begin = 0, end = 1, direction = 1),
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontfamily = anchor_font_family, fontsize = 4.5),
  column_names_gp = gpar(fontfamily = anchor_font_family, fontsize = 4),
  heatmap_legend_param = list(title_gp = gpar(fontfamily = anchor_font_family))
)

pair_table <- copy(correlation_long[TF1 != TF2])
pair_table[, pair_key := ifelse(TF1 < TF2, paste(TF1, TF2, sep = "__"), paste(TF2, TF1, sep = "__"))]
pair_table <- pair_table[order(pair_key)][, .SD[1L], by = pair_key]
top_pairs <- head(pair_table[order(-pearson_r)], 10L)
bottom_pairs <- head(pair_table[order(pearson_r)], 10L)
top_bottom_pairs <- rbind(top_pairs, bottom_pairs)
top_bottom_pairs[, label := paste(TF1, TF2, sep = "–")]
top_bottom_pairs <- unique(top_bottom_pairs, by = "pair_key")
top_bottom_pairs[, label := factor(label, levels = label[order(pearson_r)])]

pair_bar <- ggplot(top_bottom_pairs, aes(pearson_r, label, fill = pearson_r)) +
  geom_col() +
  scale_fill_gradientn(
    colours = c("#000004FF", "#B63679FF", "#FCFDBFFF"),
    limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1),
    name = "NES Pearson\ncorrelation"
  ) +
  labs(x = "NES Pearson correlation", y = "TF pair (top 10, bottom 10)") +
  anchor_theme_bw(10) +
  theme(
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(size = 9, colour = "black"),
    legend.position = "none"
  )

heatmap_grob <- grid.grabExpr({
  draw(correlation_heatmap, heatmap_legend_side = "right")
  grid.text(
    "C", x = unit(0.01, "npc"), y = unit(0.995, "npc"),
    just = c("left", "top"),
    gp = gpar(family = anchor_font_family, fontsize = 16, fontface = "bold")
  )
})
panel_3c <- wrap_elements(full = heatmap_grob) / pair_bar + plot_layout(heights = c(3.6, 1.4))
artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3c, "Figure3C_anchor_activity_correlation", 5.7, 7.0
)

# ---------------------------------------------------------------------------
# Figure 3E: regulon size versus graph degree, using the same degree as 3D.
# ---------------------------------------------------------------------------

partners_path <- file.path(
  fig3_root, "results", "tables", "Figure3E_regulon_size_and_partners.tsv"
)
partners_table <- fread(partners_path)
anchor_assert_columns(partners_table, c("TF", "partners", "regulon_targets"), "partners_table")
degree_check <- data.table(TF = names(tf_degree), graph_degree = as.integer(tf_degree))
partners_table <- merge(partners_table, degree_check, by = "TF", all = TRUE)
if (anyNA(partners_table$partners) || any(partners_table$partners != partners_table$graph_degree)) {
  stop("Figure 3E partner counts are not identical to the complete Figure 3D graph degree.", call. = FALSE)
}
safe_mad_z <- function(x) {
  scale <- mad(x)
  if (!is.finite(scale) || scale == 0) return(rep(0, length(x)))
  (x - median(x)) / scale
}
partners_table[, `:=`(
  z_partners = safe_mad_z(partners),
  z_regulon = safe_mad_z(regulon_targets)
)]
partners_table[, is_outlier := z_partners > 2 | z_regulon > 2]
partners_table[, label := fifelse(is_outlier, TF, NA_character_)]
panel_3e <- ggplot(partners_table, aes(regulon_targets, partners)) +
  geom_point(aes(colour = is_outlier), size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = label), family = anchor_font_family, size = 3, max.overlaps = 40) +
  scale_colour_manual(values = c(
    `FALSE` = "grey70",
    `TRUE` = unname(anchor_colours$regulon[["activated_shared"]])
  ), guide = "none") +
  labs(tag = "E", x = "Regulon size (# target genes)", y = "Number of HC-TF partners") +
  anchor_theme_bw(10) + anchor_panel_tag_theme +
  theme(axis.text = element_text(colour = "black", size = 8))
artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3e, "Figure3E_anchor_regulon_size_vs_partners", 5, 2.25
)

# ---------------------------------------------------------------------------
# Figure 3F: intertumour skewness, matching the anchor axis and mappings.
# ---------------------------------------------------------------------------

skew_path <- file.path(
  fig3_root, "results", "tables", "Figure3F_TCGA_LUAD_HC_TF_skewness.tsv"
)
skew_table <- fread(skew_path)
anchor_assert_columns(skew_table, c("TF", "mean_activity", "skewness", "FDR"), "skew_table")
skew_table[, is_skewed := FDR < 0.05]
skew_table[, label := fifelse(is_skewed & abs(skewness) > 1, TF, NA_character_)]
skew_table[, display_fdr := pmax(FDR, .Machine$double.xmin)]
panel_3f <- ggplot(skew_table, aes(skewness, mean_activity)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_point(aes(colour = is_skewed, size = -log10(display_fdr)), alpha = 0.8) +
  geom_text_repel(aes(label = label), family = anchor_font_family, size = 3, max.overlaps = 40) +
  scale_colour_manual(values = c(
    `FALSE` = "grey70",
    `TRUE` = unname(anchor_colours$regulon[["activated_shared"]])
  ), name = "FDR < 0.05") +
  scale_size_continuous(name = "-log10(FDR)", range = c(1, 5)) +
  labs(tag = "F", x = "Skewness (LUAD samples)", y = "Mean HC-TF activity (NES)") +
  anchor_theme_bw(10) + anchor_panel_tag_theme +
  theme(axis.text = element_text(colour = "black", size = 8))
artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
  panel_3f, "Figure3F_anchor_TCGA_LUAD_skewness", 5, 2.75
)

# ---------------------------------------------------------------------------
# Figure 3G: representative distributions for the four highest-degree TFs.
# ---------------------------------------------------------------------------

if (!is.na(fig1_root)) {
  activity_path <- file.path(fig1_root, "results", "tcga", "TCGA_viper_activity.tsv.gz")
  manifest_path <- file.path(fig1_root, "data", "processed", "tcga_manifest.tsv")
  activity <- read_numeric_matrix(activity_path)
  manifest <- fread(manifest_path)
  activity_ids <- intersect(colnames(activity), manifest$sample_id)
  manifest <- manifest[match(activity_ids, sample_id)]
  top_degree_tfs <- names(sort(tf_degree, decreasing = TRUE))[seq_len(min(4L, length(tf_degree)))]
  density_table <- rbindlist(lapply(top_degree_tfs, function(tf) {
    data.table(
      TF = tf, sample_id = activity_ids, histology = manifest$group,
      activity = as.numeric(activity[tf, activity_ids])
    )
  }))
  density_table <- density_table[histology %in% c("LUAD", "LUSC") & is.finite(activity)]
  density_stats <- density_table[, .(
    p_value = tryCatch(wilcox.test(activity ~ histology)$p.value, error = function(e) NA_real_)
  ), by = TF]
  density_table <- merge(density_table, density_stats, by = "TF", all.x = TRUE)
  density_table[, p_label := paste0("p = ", format.pval(p_value, digits = 2))]
  density_table[, TF := factor(TF, levels = top_degree_tfs)]

  panel_3g <- ggplot(density_table, aes(activity, fill = histology)) +
    geom_density(alpha = 0.70, colour = "black", linewidth = 0.45) +
    geom_text(
      data = unique(density_table[, .(TF, p_label)]),
      aes(x = Inf, y = Inf, label = p_label), inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.4, family = anchor_font_family, size = 3.2
    ) +
    facet_wrap(~TF, scales = "free", ncol = 2) +
    scale_fill_manual(values = anchor_colours$sample_category[c("LUAD", "LUSC")]) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(tag = "G", x = "TF activity (NES)", y = "Density", fill = "Histology") +
    anchor_theme_bw(10) + anchor_panel_tag_theme +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold"), legend.position = "right"
    )
  artifact_rows[[length(artifact_rows) + 1L]] <- save_anchor_atomic(
    panel_3g, "Figure3G_anchor_high_degree_activity_distributions", 6, 4.2
  )
  fwrite(density_table, file.path(table_dir, "Figure3G_anchor_high_degree_distributions.tsv"), sep = "\t")
}

# Machine-readable graph objects make the visual mappings auditable.
target_node_table <- data.table(
  target = V(target_graph)$name,
  degree = as.integer(degree(target_graph)[V(target_graph)$name]),
  community = as.integer(target_membership[V(target_graph)$name])
)
tf_node_table <- data.table(
  TF = V(tf_graph)$name,
  degree = as.integer(degree(tf_graph)[V(tf_graph)$name])
)
setorder(target_node_table, community, target)
setorder(tf_node_table, -degree, TF)
setorder(target_edges, from, to)
setorder(tf_edges, -weight, from, to)

fwrite(target_node_table, file.path(table_dir, "Figure3B_anchor_equivalent_nodes.tsv"), sep = "\t")
fwrite(target_edges, file.path(table_dir, "Figure3B_anchor_equivalent_edges.tsv"), sep = "\t")
fwrite(dyad_rows, file.path(table_dir, "Figure3B_dyad_shared_TFs.tsv"), sep = "\t")
fwrite(tf_node_table, file.path(table_dir, "Figure3D_anchor_equivalent_nodes.tsv"), sep = "\t")
fwrite(tf_edges, file.path(table_dir, "Figure3D_anchor_equivalent_edges.tsv"), sep = "\t")

artifact_table <- rbindlist(artifact_rows)
artifact_table[, pdf_md5 := unname(tools::md5sum(pdf_path))]
artifact_table[, png_md5 := unname(tools::md5sum(png_path))]
fwrite(artifact_table, file.path(audit_dir, "Figure3_anchor_style_artifacts.tsv"), sep = "\t")

receipt <- data.table(
  key = c(
    "HC_TFs",
    "Figure3B_nodes",
    "Figure3B_edges",
    "Figure3B_components",
    "Figure3B_min_shared_TFs",
    "Figure3D_nodes",
    "Figure3D_edges",
    "Figure3D_components",
    "Figure3D_isolates",
    "Figure3D_min_shared_targets",
    "Figure3D_degree_range",
    "Figure3D_edge_weight_range"
  ),
  value = c(
    length(hc_tfs),
    vcount(target_graph),
    ecount(target_graph),
    components(target_graph)$no,
    3L,
    vcount(tf_graph),
    ecount(tf_graph),
    components(tf_graph)$no,
    sum(degree(tf_graph) == 0L),
    1L,
    paste(range(degree(tf_graph)), collapse = "-"),
    paste(range(E(tf_graph)$weight), collapse = "-")
  )
)
fwrite(receipt, file.path(audit_dir, "Figure3_anchor_style_receipt.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "Figure3_anchor_style_sessionInfo.txt"))

message("Anchor-style Figure 3 network renderer is ready: ", atomic_dir)

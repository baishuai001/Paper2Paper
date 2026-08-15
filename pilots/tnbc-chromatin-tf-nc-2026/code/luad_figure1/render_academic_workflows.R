#!/usr/bin/env Rscript

# Publication-style, render-only workflow figures for the LUAD Figure 1 transfer.
# This script reads frozen manifests/results and never changes numeric analyses.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: render_academic_workflows.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
processed <- file.path(root, "data", "processed")
results <- file.path(root, "results")
tables <- file.path(results, "tables")
figures <- file.path(results, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

open_device <- function(path, width, height, dpi = 300) {
  if (grepl("[.]pdf$", path, ignore.case = TRUE)) {
    cairo_pdf(path, width = width, height = height, family = "sans")
  } else {
    png(path, width = round(width * dpi), height = round(height * dpi),
        res = dpi, type = "cairo")
  }
}

COL <- list(
  ink = "#18232D",
  sub = "#63717D",
  line = "#AEB9C2",
  border = "#D7DEE4",
  paper = "#FAFBFC",
  blue = "#3C6485",
  blue_dark = "#294F70",
  blue_light = "#EAF1F6",
  teal = "#4F7F7A",
  teal_light = "#EAF3F1",
  rust = "#A95E52",
  rust_light = "#F5EAE7",
  lusc = "#B57A58",
  lusc_light = "#F4ECE7"
)

arrow_line <- function(x0, y0, x1, y1, dashed = FALSE, col = COL$line, lwd = 0.85) {
  grid.lines(
    unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
    arrow = grid::arrow(length = unit(0.052, "inches"), type = "closed"),
    gp = gpar(col = col, lwd = lwd, lty = if (dashed) 2 else 1)
  )
}

connector <- function(x, y, dashed = FALSE, col = COL$line, arrow = FALSE, lwd = 0.8) {
  grid.lines(
    unit(x, "npc"), unit(y, "npc"),
    arrow = if (arrow) grid::arrow(length = unit(0.052, "inches"), type = "closed") else NULL,
    gp = gpar(col = col, lwd = lwd, lty = if (dashed) 2 else 1)
  )
}

section_header <- function(x0, x1, y, label, accent) {
  grid.text(label, unit(x0, "npc"), unit(y, "npc"), just = "left",
            gp = gpar(fontsize = 7.4, fontface = "bold", col = accent,
                      fontfamily = "sans"))
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y - 0.025, y - 0.025), "npc"),
             gp = gpar(col = "#D9E0E5", lwd = 0.7))
}

round_box <- function(x, y, width, height, label, fill = COL$paper, col = COL$border,
                      text_col = COL$ink, fontsize = 7.8, fontface = "plain", lwd = 0.7) {
  grid.roundrect(
    unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
    r = unit(2.0, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd)
  )
  grid.text(label, unit(x, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = fontsize, fontface = fontface, col = text_col,
                      lineheight = 1.06, fontfamily = "sans"))
}

dataset_card <- function(x, y, width, height, title, subtitle, accent = COL$blue) {
  grid.roundrect(
    unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
    r = unit(2.2, "mm"), gp = gpar(fill = COL$paper, col = COL$border, lwd = 0.7)
  )
  grid.roundrect(
    unit(x - width / 2 + 0.008, "npc"), unit(y, "npc"),
    unit(0.010, "npc"), unit(height - 0.018, "npc"), r = unit(0.9, "mm"),
    gp = gpar(fill = accent, col = NA)
  )
  grid.text(title, unit(x - width / 2 + 0.026, "npc"),
            unit(y + height * 0.19, "npc"), just = "left",
            gp = gpar(fontsize = 8.7, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text(subtitle, unit(x - width / 2 + 0.026, "npc"),
            unit(y - height * 0.19, "npc"), just = "left",
            gp = gpar(fontsize = 6.9, col = COL$sub, lineheight = 1.03,
                      fontfamily = "sans"))
}

network_mark <- function(x, y, scale = 1, accent = COL$blue) {
  grid.circle(unit(x, "npc"), unit(y, "npc"), r = unit(0.052 * scale, "npc"),
              gp = gpar(fill = "#F7F9FA", col = "#E1E6EA", lwd = 0.6))
  theta <- seq(0, 2 * pi, length.out = 13)[-13]
  px <- x + 0.038 * scale * cos(theta)
  py <- y + 0.038 * scale * sin(theta)
  edges <- rbind(
    c(1, 4), c(1, 8), c(2, 5), c(2, 10), c(3, 7), c(3, 11),
    c(4, 8), c(4, 12), c(5, 9), c(6, 10), c(7, 11), c(8, 12),
    c(1, 10), c(2, 8), c(5, 12)
  )
  for (i in seq_len(nrow(edges))) {
    a <- edges[i, 1]
    b <- edges[i, 2]
    grid.lines(unit(c(px[a], px[b]), "npc"), unit(c(py[a], py[b]), "npc"),
               gp = gpar(col = "#C9D2D9", lwd = 0.5))
  }
  grid.points(unit(px, "npc"), unit(py, "npc"), pch = 16,
              size = unit(1.45, "mm"), gp = gpar(col = accent))
  grid.points(unit(x, "npc"), unit(y, "npc"), pch = 16,
              size = unit(2.0, "mm"), gp = gpar(col = COL$rust))
}

activity_mark <- function(x, y, width = 0.095, height = 0.027) {
  cols <- colorRampPalette(c("#356A91", "#E9EEF1", "#B75F53"))(36)
  for (i in seq_along(cols)) {
    grid.rect(
      unit(x - width / 2 + width * (i - 0.5) / length(cols), "npc"), unit(y, "npc"),
      unit(width / length(cols), "npc"), unit(height, "npc"),
      gp = gpar(fill = cols[i], col = NA)
    )
  }
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
                 r = unit(0.7, "mm"), gp = gpar(fill = NA, col = "#D2D9DF", lwd = 0.45))
}

number_callout <- function(x, y, number, label, accent = COL$blue) {
  grid.text(fmt_n(number), unit(x, "npc"), unit(y + 0.022, "npc"),
            gp = gpar(fontsize = 18, fontface = "bold", col = accent,
                      fontfamily = "sans"))
  grid.text(label, unit(x, "npc"), unit(y - 0.032, "npc"),
            gp = gpar(fontsize = 7.3, col = COL$sub, lineheight = 1.03,
                      fontfamily = "sans"))
}

sample_chip <- function(x, y, title, value, accent = COL$blue) {
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(0.103, "npc"), unit(0.064, "npc"),
                 r = unit(2.1, "mm"), gp = gpar(fill = "white", col = "#DCE3E8", lwd = 0.55))
  grid.circle(unit(x - 0.038, "npc"), unit(y, "npc"), r = unit(0.006, "npc"),
              gp = gpar(fill = accent, col = NA))
  grid.text(title, unit(x - 0.025, "npc"), unit(y + 0.012, "npc"), just = "left",
            gp = gpar(fontsize = 6.2, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text(value, unit(x - 0.025, "npc"), unit(y - 0.014, "npc"), just = "left",
            gp = gpar(fontsize = 6.2, col = COL$sub, fontfamily = "sans"))
}

cohort_node <- function(x, y, width, height, title, n, subtitle = NULL,
                        accent = COL$blue, final = FALSE) {
  fill <- if (final) COL$blue_light else COL$paper
  border <- if (final) accent else COL$border
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
                 r = unit(2.2, "mm"),
                 gp = gpar(fill = fill, col = border, lwd = if (final) 1.0 else 0.65))
  grid.roundrect(unit(x - width / 2 + 0.008, "npc"), unit(y, "npc"),
                 unit(0.009, "npc"), unit(height - 0.018, "npc"), r = unit(0.8, "mm"),
                 gp = gpar(fill = accent, col = NA))
  grid.text(title, unit(x - width / 2 + 0.027, "npc"), unit(y + height * 0.19, "npc"),
            just = "left", gp = gpar(fontsize = 8.1, fontface = "bold",
                                      col = COL$ink, fontfamily = "sans"))
  grid.text(sprintf("n = %s", fmt_n(n)), unit(x + width / 2 - 0.024, "npc"),
            unit(y + height * 0.19, "npc"), just = "right",
            gp = gpar(fontsize = 8.7, fontface = "bold", col = accent,
                      fontfamily = "sans"))
  if (!is.null(subtitle)) {
    grid.text(subtitle, unit(x - width / 2 + 0.027, "npc"),
              unit(y - height * 0.19, "npc"), just = "left",
              gp = gpar(fontsize = 6.8, col = COL$sub, lineheight = 1.02,
                        fontfamily = "sans"))
  }
}

filter_pill <- function(x, y, label, width = 0.23) {
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(0.042, "npc"),
                 r = unit(2.6, "mm"), gp = gpar(fill = "#EEF3F7", col = NA))
  grid.text(label, unit(x, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = 7.0, fontface = "bold", col = "#526575",
                      fontfamily = "sans"))
}

exclusion_note <- function(x, y, title, n) {
  grid.circle(unit(x, "npc"), unit(y, "npc"), r = unit(0.010, "npc"),
              gp = gpar(fill = COL$rust_light, col = NA))
  grid.text("-", unit(x, "npc"), unit(y + 0.001, "npc"),
            gp = gpar(fontsize = 8.2, fontface = "bold", col = COL$rust,
                      fontfamily = "sans"))
  grid.text(sprintf("%s\nexcluded, n = %s", title, fmt_n(n)),
            unit(x + 0.024, "npc"), unit(y, "npc"), just = "left",
            gp = gpar(fontsize = 6.9, col = "#6B5A56", lineheight = 1.04,
                      fontfamily = "sans"))
}

subtype_chip <- function(x, y, label, n, fill, text_col = "white") {
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(0.205, "npc"), unit(0.060, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = fill, col = NA))
  grid.text(sprintf("%s  n = %s", label, fmt_n(n)), unit(x, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = 7.6, fontface = "bold", col = text_col,
                      fontfamily = "sans"))
}

tcga_manifest <- fread(file.path(processed, "tcga_manifest.tsv"))
gse_manifest <- fread(file.path(processed, "gse81089_manifest.tsv"))
pdmr_manifest <- fread(file.path(processed, "pdmr_manifest.tsv"))
depmap_manifest <- fread(file.path(processed, "depmap_manifest.tsv"))
selected <- fread(file.path(tables, "tcga_specific_TFs.tsv"))
replicated <- fread(file.path(tables, "externally_replicated_TFs.tsv"))
flow_receipt <- fread(file.path(tables, "Figure1_anchor_style_render_receipt.tsv"))

flow_value <- function(metric_name) {
  keep <- flow_receipt$section == "SupplementaryFigure1" &
    flow_receipt$metric == metric_name
  hit <- flow_receipt$value[keep]
  if (length(hit) != 1) stop(sprintf("Missing or duplicated flow metric: %s", metric_name))
  as.integer(hit)
}

count <- list(
  tcga_luad = sum(tcga_manifest$group == "LUAD"),
  tcga_lusc = sum(tcga_manifest$group == "LUSC"),
  gse_luad = sum(gse_manifest$group == "LUAD"),
  gse_lusc = sum(gse_manifest$group == "LUSC"),
  pdx = nrow(pdmr_manifest),
  cell_line = nrow(depmap_manifest),
  luad_tf = sum(selected$discovery_category == "LUAD"),
  replicated_tf = sum(replicated$discovery_category == "LUAD"),
  tcga_all = flow_value("tcga_all"),
  tcga_primary = flow_value("tcga_primary"),
  tcga_non_primary = flow_value("tcga_non_primary"),
  tcga_unique = flow_value("tcga_unique"),
  tcga_duplicate = flow_value("tcga_duplicate"),
  gse_all = flow_value("gse_all"),
  gse_tumor = flow_value("gse_tumor"),
  gse_normal = flow_value("gse_normal"),
  gse_target_histology = flow_value("gse_target_histology"),
  gse_non_target_histology = flow_value("gse_non_target_histology"),
  gse_target_with_expression = flow_value("gse_target_with_expression"),
  gse_unmatched = flow_value("gse_unmatched")
)

stopifnot(
  count$tcga_luad == 516L, count$tcga_lusc == 501L,
  count$gse_luad == 106L, count$gse_lusc == 67L,
  count$luad_tf == 158L, count$replicated_tf == 97L,
  count$tcga_unique == nrow(tcga_manifest),
  count$gse_target_with_expression == nrow(gse_manifest)
)

render_figure1a <- function(path) {
  open_device(path, 10.8, 4.45)
  grid.newpage()
  grid.text("A", unit(0.016, "npc"), unit(0.978, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink, fontfamily = "sans"))

  section_header(0.055, 0.704, 0.925, "DISCOVERY AND SINGLE-SAMPLE PROJECTION", COL$blue)
  section_header(0.758, 0.958, 0.925, "INDEPENDENT REPLICATION", COL$teal)
  grid.lines(unit(c(0.731, 0.731), "npc"), unit(c(0.085, 0.885), "npc"),
             gp = gpar(col = "#D7DEE4", lwd = 0.65, lty = 3))

  dataset_card(
    0.140, 0.700, 0.190, 0.145,
    "TCGA lung RNA-seq",
    sprintf("%s tumors\n%s LUAD | %s LUSC",
            fmt_n(nrow(tcga_manifest)), fmt_n(count$tcga_luad), fmt_n(count$tcga_lusc)),
    COL$blue
  )
  arrow_line(0.240, 0.700, 0.276, 0.700)
  network_mark(0.330, 0.700, 0.95, COL$blue)
  grid.text("ARACNe3", unit(0.330, "npc"), unit(0.615, "npc"),
            gp = gpar(fontsize = 8.0, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("lung cancer interactome", unit(0.330, "npc"), unit(0.585, "npc"),
            gp = gpar(fontsize = 6.6, col = COL$sub, fontfamily = "sans"))

  arrow_line(0.386, 0.700, 0.424, 0.700)
  activity_mark(0.475, 0.708, 0.092, 0.026)
  grid.text("VIPER / msVIPER", unit(0.475, "npc"), unit(0.655, "npc"),
            gp = gpar(fontsize = 8.0, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("TF activity", unit(0.475, "npc"), unit(0.625, "npc"),
            gp = gpar(fontsize = 6.6, col = COL$sub, fontfamily = "sans"))

  connector(c(0.525, 0.548, 0.548), c(0.700, 0.700, 0.755))
  connector(c(0.548, 0.548), c(0.700, 0.520))
  arrow_line(0.548, 0.755, 0.574, 0.755)
  arrow_line(0.548, 0.520, 0.574, 0.520)
  round_box(0.628, 0.755, 0.108, 0.095, "Cohort-level\ncontrast",
            fill = COL$blue_light, col = NA, text_col = COL$blue_dark,
            fontsize = 7.5, fontface = "bold")
  round_box(0.628, 0.520, 0.108, 0.095, "Single-sample\nprojection",
            fill = "#F1F4F6", col = NA, text_col = "#465763",
            fontsize = 7.5, fontface = "bold")

  number_callout(0.628, 0.645, count$luad_tf, "LUAD-specific TFs", COL$blue)
  sample_chip(0.535, 0.365, "PATIENTS", sprintf("n = %s", fmt_n(nrow(tcga_manifest))), COL$blue)
  sample_chip(0.643, 0.365, "PDX", sprintf("n = %s", fmt_n(count$pdx)), COL$teal)
  sample_chip(0.589, 0.275, "CELL LINES", sprintf("n = %s", fmt_n(count$cell_line)), COL$rust)
  connector(c(0.628, 0.628, 0.589), c(0.470, 0.420, 0.397), col = "#C1CBD2")

  round_box(
    0.430, 0.125, 0.370, 0.095,
    sprintf("%s LUAD-specific TFs carried forward to chromatin filtering (Figure 2)",
            fmt_n(count$luad_tf)),
    fill = COL$blue_light, col = COL$blue, text_col = COL$blue_dark,
    fontsize = 8.0, fontface = "bold", lwd = 0.9
  )
  connector(c(0.670, 0.700, 0.700, 0.615), c(0.620, 0.620, 0.175, 0.175),
            col = COL$blue, arrow = TRUE)

  dataset_card(
    0.858, 0.755, 0.190, 0.135,
    "GSE81089 RNA-seq",
    sprintf("%s tumors\n%s LUAD | %s LUSC",
            fmt_n(nrow(gse_manifest)), fmt_n(count$gse_luad), fmt_n(count$gse_lusc)),
    COL$teal
  )
  arrow_line(0.858, 0.682, 0.858, 0.625)
  network_mark(0.858, 0.555, 0.80, COL$teal)
  grid.text("ARACNe3 interactome", unit(0.858, "npc"), unit(0.480, "npc"),
            gp = gpar(fontsize = 7.1, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  arrow_line(0.858, 0.452, 0.858, 0.410)
  activity_mark(0.858, 0.372, 0.084, 0.024)
  grid.text("VIPER / msVIPER", unit(0.858, "npc"), unit(0.327, "npc"),
            gp = gpar(fontsize = 7.2, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  arrow_line(0.858, 0.300, 0.858, 0.260)
  round_box(
    0.858, 0.190, 0.184, 0.098,
    sprintf("%s TFs replicated", fmt_n(count$replicated_tf)),
    fill = COL$teal_light, col = COL$teal, text_col = "#315F5B",
    fontsize = 8.2, fontface = "bold", lwd = 0.9
  )
  grid.text(sprintf("%s of %s discovery TFs", fmt_n(count$replicated_tf), fmt_n(count$luad_tf)),
            unit(0.858, "npc"), unit(0.112, "npc"),
            gp = gpar(fontsize = 6.6, col = COL$teal, fontfamily = "sans"))
  dev.off()
}

render_tcga_flow <- function(path) {
  open_device(path, 6.7, 6.6)
  grid.newpage()
  grid.text("A", unit(0.025, "npc"), unit(0.978, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("TCGA lung cancer discovery cohort", unit(0.095, "npc"), unit(0.945, "npc"),
            just = "left", gp = gpar(fontsize = 11.3, fontface = "bold",
                                      col = COL$ink, fontfamily = "sans"))
  grid.text("RNA-seq sample derivation", unit(0.095, "npc"), unit(0.905, "npc"),
            just = "left", gp = gpar(fontsize = 7.4, col = COL$sub, fontfamily = "sans"))
  grid.lines(unit(c(0.095, 0.905), "npc"), unit(c(0.875, 0.875), "npc"),
             gp = gpar(col = "#D9E0E5", lwd = 0.7))

  x <- 0.405
  w <- 0.485
  h <- 0.105
  cohort_node(x, 0.795, w, h, "Available expression columns", count$tcga_all,
              "Counts and TPM present", COL$blue)
  connector(c(x, x), c(0.742, 0.690), col = COL$line)
  filter_pill(x, 0.665, "PRIMARY TUMOR ONLY", 0.235)
  connector(c(x, x), c(0.640, 0.602), col = COL$line, arrow = TRUE)
  connector(c(0.523, 0.690), c(0.665, 0.665), col = "#C9B4AF")
  exclusion_note(0.716, 0.665, "Non-primary samples", count$tcga_non_primary)

  cohort_node(x, 0.545, w, h, "Primary-tumor columns", count$tcga_primary,
              "TCGA sample type 01", COL$blue)
  connector(c(x, x), c(0.492, 0.440), col = COL$line)
  filter_pill(x, 0.415, "ONE TUMOR PER PATIENT", 0.265)
  connector(c(x, x), c(0.390, 0.352), col = COL$line, arrow = TRUE)
  connector(c(0.538, 0.690), c(0.415, 0.415), col = "#C9B4AF")
  exclusion_note(0.716, 0.415, "Duplicate primary aliquots", count$tcga_duplicate)

  cohort_node(x, 0.295, w, h, "Final analytic cohort", count$tcga_unique,
              "Unique patients used for ARACNe3 and VIPER", COL$blue, final = TRUE)
  connector(c(x, x), c(0.242, 0.183), col = COL$blue)
  connector(c(0.405, 0.295, 0.295), c(0.183, 0.155, 0.132), col = COL$blue)
  connector(c(0.405, 0.515, 0.515), c(0.183, 0.155, 0.132), col = COL$blue)
  subtype_chip(0.295, 0.092, "LUAD", count$tcga_luad, COL$blue_dark)
  subtype_chip(0.515, 0.092, "LUSC", count$tcga_lusc, COL$lusc)
  dev.off()
}

render_gse_flow <- function(path) {
  open_device(path, 6.7, 7.1)
  grid.newpage()
  grid.text("B", unit(0.025, "npc"), unit(0.978, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("GSE81089 independent cohort", unit(0.095, "npc"), unit(0.945, "npc"),
            just = "left", gp = gpar(fontsize = 11.3, fontface = "bold",
                                      col = COL$ink, fontfamily = "sans"))
  grid.text("RNA-seq sample derivation", unit(0.095, "npc"), unit(0.907, "npc"),
            just = "left", gp = gpar(fontsize = 7.4, col = COL$sub, fontfamily = "sans"))
  grid.lines(unit(c(0.095, 0.905), "npc"), unit(c(0.878, 0.878), "npc"),
             gp = gpar(col = "#D9E0E5", lwd = 0.7))

  x <- 0.405
  w <- 0.485
  h <- 0.094
  cohort_node(x, 0.812, w, h, "Available RNA-seq samples", count$gse_all,
              "featureCounts matrix", COL$teal)
  connector(c(x, x), c(0.765, 0.730), col = COL$line)
  filter_pill(x, 0.707, "TUMOR TISSUE ONLY", 0.218)
  connector(c(x, x), c(0.684, 0.654), col = COL$line, arrow = TRUE)
  connector(c(0.515, 0.690), c(0.707, 0.707), col = "#C9B4AF")
  exclusion_note(0.716, 0.707, "Matched normal tissue", count$gse_normal)

  cohort_node(x, 0.607, w, h, "Tumor samples", count$gse_tumor,
              "All reported lung tumor histologies", COL$teal)
  connector(c(x, x), c(0.560, 0.525), col = COL$line)
  filter_pill(x, 0.502, "LUAD OR LUSC HISTOLOGY", 0.275)
  connector(c(x, x), c(0.479, 0.449), col = COL$line, arrow = TRUE)
  connector(c(0.543, 0.690), c(0.502, 0.502), col = "#C9B4AF")
  exclusion_note(0.716, 0.502, "Large-cell or NOS tumors", count$gse_non_target_histology)

  cohort_node(x, 0.402, w, h, "LUAD/LUSC by metadata", count$gse_target_histology,
              "Histology label retained", COL$teal)
  connector(c(x, x), c(0.355, 0.320), col = COL$line)
  filter_pill(x, 0.297, "EXPRESSION ID MATCH", 0.235)
  connector(c(x, x), c(0.274, 0.244), col = COL$line, arrow = TRUE)
  connector(c(0.523, 0.690), c(0.297, 0.297), col = "#C9B4AF")
  exclusion_note(0.716, 0.297, "Unmatched sample IDs", count$gse_unmatched)

  cohort_node(x, 0.197, w, h, "Final replication cohort", count$gse_target_with_expression,
              "Used for independent ARACNe3 and VIPER", COL$teal, final = TRUE)
  connector(c(x, x), c(0.150, 0.125), col = COL$teal)
  connector(c(0.405, 0.295, 0.295), c(0.125, 0.105, 0.082), col = COL$teal)
  connector(c(0.405, 0.515, 0.515), c(0.125, 0.105, 0.082), col = COL$teal)
  subtype_chip(0.295, 0.052, "LUAD", count$gse_luad, COL$blue_dark)
  subtype_chip(0.515, 0.052, "LUSC", count$gse_lusc, COL$lusc)
  dev.off()
}

render_figure1a(file.path(figures, "Figure1A_workflow.pdf"))
render_figure1a(file.path(figures, "Figure1A_workflow.png"))
render_tcga_flow(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.pdf"))
render_tcga_flow(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.png"))
render_gse_flow(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.pdf"))
render_gse_flow(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.png"))

fwrite(
  data.table(
    artifact = c("Figure1A", "SupplementaryFigure1A", "SupplementaryFigure1B"),
    style_version = "academic-v2",
    numeric_source = "frozen Figure1 manifests and Figure1_anchor_style_render_receipt.tsv"
  ),
  file.path(tables, "academic_workflow_render_receipt.tsv"),
  sep = "\t"
)

cat("Rendered academic-v2 Figure 1A and Supplementary Figure 1A/B without changing numeric inputs.\n")

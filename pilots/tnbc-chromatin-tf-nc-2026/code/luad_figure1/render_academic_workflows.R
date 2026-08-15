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
  # Anchor-faithful evidence logic: discovery defines the TF set; single-sample
  # analysis projects it; the independent cohort supports but does not shrink it.
  open_device(path, 10.2, 5.65)
  grid.newpage()
  grid.text("A", unit(0.016, "npc"), unit(0.982, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink, fontfamily = "sans"))

  panel_border <- function(x, y, width, height) {
    grid.rect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
              gp = gpar(fill = "white", col = "#7C8891", lwd = 0.75))
  }
  analysis_box <- function(x, y, label, width = 0.135) {
    grid.rect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(0.092, "npc"),
              gp = gpar(fill = "#FBFCFD", col = "#78858E", lwd = 0.7))
    grid.text(label, unit(x, "npc"), unit(y, "npc"),
              gp = gpar(fontsize = 7.4, col = COL$ink, lineheight = 1.03,
                        fontfamily = "sans"))
  }
  cohort_label <- function(x, y, total, first, second, accent) {
    grid.circle(unit(x - 0.012, "npc"), unit(y + 0.006, "npc"), r = unit(0.008, "npc"),
                gp = gpar(fill = accent, col = NA))
    grid.text(sprintf("%s tumors", fmt_n(total)), unit(x + 0.008, "npc"), unit(y + 0.014, "npc"),
              just = "left", gp = gpar(fontsize = 8.0, fontface = "bold",
                                        col = COL$ink, fontfamily = "sans"))
    grid.text(sprintf("%s LUAD  |  %s LUSC", fmt_n(first), fmt_n(second)),
              unit(x + 0.008, "npc"), unit(y - 0.020, "npc"), just = "left",
              gp = gpar(fontsize = 6.7, col = COL$sub, fontfamily = "sans"))
  }
  sample_line <- function(x, y, label, n, accent) {
    grid.circle(unit(x, "npc"), unit(y, "npc"), r = unit(0.006, "npc"),
                gp = gpar(fill = accent, col = NA))
    grid.text(sprintf("%s  %s", fmt_n(n), label), unit(x + 0.016, "npc"), unit(y, "npc"),
              just = "left", gp = gpar(fontsize = 6.9, col = COL$ink, fontfamily = "sans"))
  }

  # Two same-height cohort modules; width ratio follows the TNBC anchor (~2:1).
  panel_border(0.350, 0.565, 0.600, 0.650)
  panel_border(0.825, 0.565, 0.300, 0.650)
  grid.text("RNA-seq (TCGA-LUAD/LUSC)", unit(0.350, "npc"), unit(0.852, "npc"),
            gp = gpar(fontsize = 10.0, col = COL$ink, fontfamily = "sans"))
  grid.text("RNA-seq (GSE81089)", unit(0.825, "npc"), unit(0.852, "npc"),
            gp = gpar(fontsize = 10.0, col = COL$ink, fontfamily = "sans"))

  cohort_label(0.135, 0.748, nrow(tcga_manifest), count$tcga_luad, count$tcga_lusc, COL$blue)
  arrow_line(0.205, 0.705, 0.205, 0.655, col = "#5F6B73", lwd = 0.9)
  grid.text("ARACNe3", unit(0.205, "npc"), unit(0.625, "npc"),
            gp = gpar(fontsize = 8.5, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  network_mark(0.205, 0.492, 1.18, COL$blue)
  grid.text("Lung cancer interactome", unit(0.205, "npc"), unit(0.352, "npc"),
            gp = gpar(fontsize = 7.2, col = COL$ink, fontfamily = "sans"))
  arrow_line(0.270, 0.535, 0.335, 0.620, col = "#5F6B73", lwd = 0.9)
  grid.text("VIPER / msVIPER", unit(0.445, "npc"), unit(0.675, "npc"),
            gp = gpar(fontsize = 8.5, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  activity_mark(0.445, 0.632, 0.205, 0.027)
  grid.text("Differential TF activity signature", unit(0.445, "npc"), unit(0.596, "npc"),
            gp = gpar(fontsize = 6.9, col = COL$sub, fontfamily = "sans"))
  analysis_box(0.375, 0.505, "Multi-sample\nanalysis")
  analysis_box(0.535, 0.505, "Single-sample\nanalysis")
  grid.text(sprintf("%s LUSC-specific TFs\n%s LUAD-specific TFs",
                    fmt_n(sum(selected$discovery_category == "LUSC")), fmt_n(count$luad_tf)),
            unit(0.305, "npc"), unit(0.389, "npc"), just = "left",
            gp = gpar(fontsize = 7.0, col = COL$ink, lineheight = 1.08,
                      fontfamily = "sans"))
  sample_line(0.480, 0.408, "TCGA tumors", nrow(tcga_manifest), COL$blue)
  sample_line(0.480, 0.365, "PDXs", count$pdx, COL$teal)
  sample_line(0.480, 0.322, "cell lines", count$cell_line, COL$rust)

  cohort_label(0.760, 0.748, nrow(gse_manifest), count$gse_luad, count$gse_lusc, COL$teal)
  arrow_line(0.825, 0.705, 0.825, 0.660, col = "#5F6B73", lwd = 0.9)
  grid.text("ARACNe3", unit(0.825, "npc"), unit(0.625, "npc"),
            gp = gpar(fontsize = 8.5, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("Independent lung\ninteractome", unit(0.825, "npc"), unit(0.570, "npc"),
            gp = gpar(fontsize = 7.1, col = COL$sub, lineheight = 1.03,
                      fontfamily = "sans"))
  arrow_line(0.825, 0.520, 0.825, 0.480, col = "#5F6B73", lwd = 0.9)
  grid.text("VIPER / msVIPER", unit(0.825, "npc"), unit(0.445, "npc"),
            gp = gpar(fontsize = 8.5, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  analysis_box(0.750, 0.350, "Multi-sample\nanalysis", width = 0.125)
  analysis_box(0.900, 0.350, "Single-sample\nanalysis", width = 0.125)
  grid.text("LUSC-specific TFs\nLUAD-specific TFs", unit(0.690, "npc"), unit(0.267, "npc"),
            just = "left", gp = gpar(fontsize = 6.8, col = COL$ink, lineheight = 1.06,
                                      fontfamily = "sans"))

  # Only multi-sample outputs converge; single-sample projection is not a filter.
  connector(c(0.350, 0.350, 0.375), c(0.345, 0.155, 0.155),
            col = COL$blue, arrow = TRUE, lwd = 0.9)
  connector(c(0.750, 0.750, 0.655), c(0.242, 0.155, 0.155),
            dashed = TRUE, col = COL$teal, arrow = TRUE, lwd = 0.9)
  round_box(0.515, 0.155, 0.280, 0.092,
            sprintf("%s LUAD-specific TFs", fmt_n(count$luad_tf)),
            fill = COL$blue_light, col = COL$blue, text_col = COL$blue_dark,
            fontsize = 8.2, fontface = "plain", lwd = 0.9)
  dev.off()
}
render_figure1a_stacked <- function(path) {
  # Original stacked composition. Evidence semantics are identical to Figure1A:
  # TCGA defines the discovery set; single-sample analyses project it; GSE81089
  # independently recapitulates it without reducing the Figure 2 input set.
  open_device(path, 8.7, 7.2)
  grid.newpage()
  grid.text("A", unit(0.018, "npc"), unit(0.982, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink, fontfamily = "sans"))
  grid.text("Transcriptional regulator discovery and cross-system validation",
            unit(0.075, "npc"), unit(0.952, "npc"), just = "left",
            gp = gpar(fontsize = 11.0, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text("LUAD versus LUSC", unit(0.075, "npc"), unit(0.920, "npc"), just = "left",
            gp = gpar(fontsize = 7.1, col = COL$sub, fontfamily = "sans"))

  lane <- function(x, y, width, height, fill, accent) {
    grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(width, "npc"), unit(height, "npc"),
                   r = unit(3.0, "mm"), gp = gpar(fill = fill, col = COL$border, lwd = 0.75))
    grid.roundrect(unit(x - width / 2 + 0.012, "npc"), unit(y, "npc"),
                   unit(0.010, "npc"), unit(height - 0.030, "npc"), r = unit(1.1, "mm"),
                   gp = gpar(fill = accent, col = NA))
  }
  step_tag <- function(x, y, label, accent) {
    grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(0.170, "npc"), unit(0.038, "npc"),
                   r = unit(2.3, "mm"), gp = gpar(fill = accent, col = NA))
    grid.text(label, unit(x, "npc"), unit(y, "npc"),
              gp = gpar(fontsize = 6.7, fontface = "bold", col = "white",
                        fontfamily = "sans"))
  }
  analysis_card <- function(x, y, title, body, accent, fill, height = 0.100) {
    grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(0.215, "npc"), unit(height, "npc"),
                   r = unit(2.2, "mm"), gp = gpar(fill = fill, col = accent, lwd = 0.75))
    grid.text(title, unit(x, "npc"), unit(y + 0.020, "npc"),
              gp = gpar(fontsize = 6.6, fontface = "bold", col = accent,
                        fontfamily = "sans"))
    grid.text(body, unit(x, "npc"), unit(y - 0.020, "npc"),
              gp = gpar(fontsize = 6.4, col = COL$ink, lineheight = 1.03,
                        fontfamily = "sans"))
  }

  # Upper lane: TCGA discovery and cross-system single-sample projection.
  lane(0.515, 0.715, 0.900, 0.360, "#F6F9FB", COL$blue)
  step_tag(0.160, 0.862, "01  DISCOVERY COHORT", COL$blue)
  dataset_card(
    0.180, 0.720, 0.180, 0.112,
    "TCGA RNA-seq",
    sprintf("%s tumors\n%s LUAD | %s LUSC", fmt_n(nrow(tcga_manifest)),
            fmt_n(count$tcga_luad), fmt_n(count$tcga_lusc)),
    COL$blue
  )
  arrow_line(0.275, 0.720, 0.320, 0.720, col = COL$line, lwd = 0.8)
  network_mark(0.370, 0.720, 0.78, COL$blue)
  grid.text("ARACNe3", unit(0.370, "npc"), unit(0.657, "npc"),
            gp = gpar(fontsize = 7.2, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text("lung interactome", unit(0.370, "npc"), unit(0.633, "npc"),
            gp = gpar(fontsize = 6.1, col = COL$sub, fontfamily = "sans"))
  arrow_line(0.420, 0.720, 0.475, 0.720, col = COL$line, lwd = 0.8)
  activity_mark(0.545, 0.724, 0.125, 0.024)
  grid.text("VIPER / msVIPER", unit(0.545, "npc"), unit(0.772, "npc"),
            gp = gpar(fontsize = 7.3, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text("TF activity", unit(0.545, "npc"), unit(0.682, "npc"),
            gp = gpar(fontsize = 6.1, col = COL$sub, fontfamily = "sans"))
  connector(c(0.610, 0.650, 0.650), c(0.720, 0.720, 0.775), col = COL$line)
  connector(c(0.650, 0.650), c(0.720, 0.625), col = COL$line)
  arrow_line(0.650, 0.775, 0.677, 0.775, col = COL$line, lwd = 0.75)
  arrow_line(0.650, 0.625, 0.677, 0.625, col = COL$line, lwd = 0.75)
  analysis_card(0.795, 0.775, "SINGLE-SAMPLE PROJECTION",
                sprintf("%s tumors | %s PDXs | %s cell lines",
                        fmt_n(nrow(tcga_manifest)), fmt_n(count$pdx), fmt_n(count$cell_line)),
                "#657580", "white")
  analysis_card(0.795, 0.625, "MULTI-SAMPLE CONTRAST",
                sprintf("%s LUAD | %s LUSC TFs", fmt_n(count$luad_tf),
                        fmt_n(sum(selected$discovery_category == "LUSC"))),
                COL$blue_dark, COL$blue_light, height = 0.105)

  # Lower lane: the independent cohort repeats network inference and TF analysis.
  lane(0.515, 0.275, 0.900, 0.310, "#F5F9F8", COL$teal)
  step_tag(0.185, 0.407, "02  INDEPENDENT REPLICATION", COL$teal)
  dataset_card(
    0.180, 0.275, 0.180, 0.108,
    "GSE81089 RNA-seq",
    sprintf("%s tumors\n%s LUAD | %s LUSC", fmt_n(nrow(gse_manifest)),
            fmt_n(count$gse_luad), fmt_n(count$gse_lusc)),
    COL$teal
  )
  arrow_line(0.275, 0.275, 0.320, 0.275, col = COL$line, lwd = 0.8)
  network_mark(0.370, 0.275, 0.74, COL$teal)
  grid.text("ARACNe3", unit(0.370, "npc"), unit(0.217, "npc"),
            gp = gpar(fontsize = 7.2, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text("independent interactome", unit(0.370, "npc"), unit(0.193, "npc"),
            gp = gpar(fontsize = 6.0, col = COL$sub, fontfamily = "sans"))
  arrow_line(0.420, 0.275, 0.475, 0.275, col = COL$line, lwd = 0.8)
  activity_mark(0.545, 0.279, 0.125, 0.024)
  grid.text("VIPER / msVIPER", unit(0.545, "npc"), unit(0.327, "npc"),
            gp = gpar(fontsize = 7.3, fontface = "bold", col = COL$ink,
                      fontfamily = "sans"))
  grid.text("TF activity", unit(0.545, "npc"), unit(0.237, "npc"),
            gp = gpar(fontsize = 6.1, col = COL$sub, fontfamily = "sans"))
  connector(c(0.610, 0.650, 0.650), c(0.275, 0.275, 0.340), col = COL$line)
  connector(c(0.650, 0.650), c(0.275, 0.210), col = COL$line)
  arrow_line(0.650, 0.340, 0.677, 0.340, col = COL$line, lwd = 0.75)
  arrow_line(0.650, 0.210, 0.677, 0.210, col = COL$line, lwd = 0.75)
  analysis_card(0.795, 0.340, "MULTI-SAMPLE REPLICATION",
                sprintf("%s/%s LUAD TFs recapitulated", fmt_n(count$replicated_tf),
                        fmt_n(count$luad_tf)),
                "#315F5B", COL$teal_light, height = 0.100)
  analysis_card(0.795, 0.210, "SINGLE-SAMPLE ANALYSIS",
                sprintf("%s independent tumors", fmt_n(nrow(gse_manifest))),
                "#657580", "white", height = 0.090)

  # The central result is defined by discovery; solid and dashed arrows distinguish
  # discovery output from independent support. Projection branches are not filters.
  arrow_line(0.795, 0.570, 0.795, 0.523, col = COL$blue, lwd = 0.9)
  arrow_line(0.795, 0.392, 0.795, 0.447, dashed = TRUE, col = COL$teal, lwd = 0.9)
  round_box(0.795, 0.485, 0.255, 0.070,
            sprintf("%s LUAD-specific TFs", fmt_n(count$luad_tf)),
            fill = COL$blue_light, col = COL$blue, text_col = COL$blue_dark,
            fontsize = 8.5, fontface = "plain", lwd = 0.9)
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
render_figure1a_stacked(file.path(figures, "Figure1A_workflow_stacked.pdf"))
render_figure1a_stacked(file.path(figures, "Figure1A_workflow_stacked.png"))
render_tcga_flow(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.pdf"))
render_tcga_flow(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.png"))
render_gse_flow(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.pdf"))
render_gse_flow(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.png"))

fwrite(
  data.table(
    artifact = c("Figure1A", "Figure1A_stacked", "SupplementaryFigure1A", "SupplementaryFigure1B"),
    style_version = c("anchor-faithful-v3", "stacked-original-v1", "academic-v2", "academic-v2"),
    numeric_source = "frozen Figure1 manifests and Figure1_anchor_style_render_receipt.tsv"
  ),
  file.path(tables, "academic_workflow_render_receipt.tsv"),
  sep = "\t"
)

cat("Rendered anchor-faithful-v3 and stacked-original-v1 Figure 1A layouts with unchanged numeric inputs.\n")

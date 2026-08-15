#!/usr/bin/env Rscript

# Render-only pass for the LUAD Figure 1 transfer.  This script deliberately
# reads the frozen Figure 1 matrices and manifests and never rewrites them.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: render_anchor_style_figure1.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
raw <- file.path(root, "data", "raw")
processed <- file.path(root, "data", "processed")
results <- file.path(root, "results")
tables <- file.path(results, "tables")
figures <- file.path(results, "figures")
annotations_dir <- file.path(processed, "annotations")
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

fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

open_device <- function(path, width, height, dpi = 240) {
  if (grepl("[.]pdf$", path, ignore.case = TRUE)) {
    cairo_pdf(path, width = width, height = height)
  } else {
    png(path, width = round(width * dpi), height = round(height * dpi), res = dpi)
  }
}

draw_arrow <- function(x0, y0, x1, y1, dashed = FALSE, lwd = 0.9) {
  grid.lines(
    x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
    arrow = grid::arrow(length = unit(0.075, "inches"), type = "closed"),
    gp = gpar(col = "#222222", lwd = lwd, lty = if (dashed) 2 else 1)
  )
}

draw_poly_arrow <- function(x, y, dashed = FALSE, lwd = 0.9) {
  grid.lines(
    x = unit(x, "npc"), y = unit(y, "npc"),
    arrow = grid::arrow(length = unit(0.075, "inches"), type = "closed"),
    gp = gpar(col = "#222222", lwd = lwd, lty = if (dashed) 2 else 1)
  )
}

draw_box <- function(x, y, width, height, label, fill = "white", fontsize = 8.0,
                     col = "#222222", text_col = "#111111", lwd = 0.9) {
  grid.rect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(width, "npc"), height = unit(height, "npc"),
    gp = gpar(fill = fill, col = col, lwd = lwd)
  )
  grid.text(
    label, x = unit(x, "npc"), y = unit(y, "npc"),
    gp = gpar(fontsize = fontsize, col = text_col, lineheight = 1.05)
  )
}

draw_diamond <- function(x, y, width, height, label, fontsize = 7.7) {
  grid.polygon(
    x = unit(c(x, x + width / 2, x, x - width / 2), "npc"),
    y = unit(c(y + height / 2, y, y - height / 2, y), "npc"),
    gp = gpar(fill = "white", col = "#222222", lwd = 0.9)
  )
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
            gp = gpar(fontsize = fontsize, col = "#111111", lineheight = 1.0))
}

draw_person <- function(x, y, scale = 1, accent = "#C44E52") {
  grid.circle(unit(x, "npc"), unit(y + 0.030 * scale, "npc"),
              r = unit(0.010 * scale, "npc"), gp = gpar(fill = "#BDBDBD", col = NA))
  grid.lines(unit(c(x, x), "npc"), unit(c(y + 0.018 * scale, y - 0.014 * scale), "npc"),
             gp = gpar(col = "#BDBDBD", lwd = 3.0 * scale))
  grid.lines(unit(c(x - 0.018 * scale, x, x + 0.018 * scale), "npc"),
             unit(c(y + 0.010 * scale, y + 0.016 * scale, y + 0.010 * scale), "npc"),
             gp = gpar(col = "#BDBDBD", lwd = 2.0 * scale))
  grid.lines(unit(c(x - 0.010 * scale, x, x + 0.010 * scale), "npc"),
             unit(c(y - 0.044 * scale, y - 0.014 * scale, y - 0.044 * scale), "npc"),
             gp = gpar(col = "#BDBDBD", lwd = 2.0 * scale))
  grid.circle(unit(x + 0.005 * scale, "npc"), unit(y + 0.006 * scale, "npc"),
              r = unit(0.004 * scale, "npc"), gp = gpar(fill = accent, col = NA))
}

draw_mouse <- function(x, y, scale = 1) {
  grid.ellipse <- function(cx, cy, w, h, fill) {
    grid.circle(unit(cx, "npc"), unit(cy, "npc"), r = unit(w / 2, "npc"),
                gp = gpar(fill = fill, col = "#9E9E9E", lwd = 0.45))
  }
  draw_mouse_body <- "#D7D0C8"
  grid.ellipse(x, y, 0.030 * scale, 0.017 * scale, draw_mouse_body)
  grid.circle(unit(x + 0.016 * scale, "npc"), unit(y + 0.004 * scale, "npc"),
              r = unit(0.007 * scale, "npc"), gp = gpar(fill = draw_mouse_body, col = "#9E9E9E", lwd = 0.45))
  grid.circle(unit(x + 0.016 * scale, "npc"), unit(y + 0.012 * scale, "npc"),
              r = unit(0.003 * scale, "npc"), gp = gpar(fill = "#D7D0C8", col = "#9E9E9E", lwd = 0.4))
  grid.lines(unit(c(x - 0.016 * scale, x - 0.025 * scale), "npc"),
             unit(c(y, y + 0.008 * scale), "npc"), gp = gpar(col = "#B08C7A", lwd = 0.7))
}

draw_dish <- function(x, y, scale = 1) {
  grid.circle(unit(x, "npc"), unit(y, "npc"), r = unit(0.014 * scale, "npc"),
              gp = gpar(fill = "#D86666", col = "#8C8C8C", lwd = 0.7))
  grid.circle(unit(x, "npc"), unit(y + 0.003 * scale, "npc"), r = unit(0.011 * scale, "npc"),
              gp = gpar(fill = "#E58A8A", col = "white", lwd = 0.45))
}

draw_network <- function(x, y, width = 0.14, height = 0.18) {
  angles <- seq(0, 2 * pi, length.out = 17)[-17]
  px <- x + width * (0.45 * cos(angles) + 0.08 * cos(3 * angles))
  py <- y + height * (0.40 * sin(angles) + 0.06 * sin(2 * angles))
  for (i in seq_along(px)) {
    for (j in seq_along(px)) {
      if (j > i && ((i * 7 + j * 11) %% 13) < 3) {
        grid.lines(unit(c(px[i], px[j]), "npc"), unit(c(py[i], py[j]), "npc"),
                   gp = gpar(col = "#D9D9D9", lwd = 0.35))
      }
    }
  }
  grid.points(unit(px, "npc"), unit(py, "npc"), pch = 16,
              size = unit(1.2, "mm"), gp = gpar(col = "#B2182B"))
}

draw_gradient <- function(x, y, width, height) {
  cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D7191C"))(60)
  for (i in seq_along(cols)) {
    grid.rect(
      x = unit(x - width / 2 + width * (i - 0.5) / length(cols), "npc"),
      y = unit(y, "npc"), width = unit(width / length(cols), "npc"),
      height = unit(height, "npc"), gp = gpar(fill = cols[i], col = NA)
    )
  }
}

tcga_manifest <- fread(file.path(processed, "tcga_manifest.tsv"))
gse_manifest <- fread(file.path(processed, "gse81089_manifest.tsv"))
pdmr_manifest <- fread(file.path(processed, "pdmr_manifest.tsv"))
depmap_manifest <- fread(file.path(processed, "depmap_manifest.tsv"))
selected <- fread(file.path(tables, "tcga_specific_TFs.tsv"))
replicated <- fread(file.path(tables, "externally_replicated_TFs.tsv"))

cohort_counts <- list(
  tcga_luad = sum(tcga_manifest$group == "LUAD"),
  tcga_lusc = sum(tcga_manifest$group == "LUSC"),
  gse_luad = sum(gse_manifest$group == "LUAD"),
  gse_lusc = sum(gse_manifest$group == "LUSC"),
  pdx = nrow(pdmr_manifest),
  cell_line = nrow(depmap_manifest),
  luad_tf = sum(selected$discovery_category == "LUAD"),
  lusc_tf = sum(selected$discovery_category == "LUSC"),
  replicated_luad_tf = sum(replicated$discovery_category == "LUAD")
)
stopifnot(
  cohort_counts$tcga_luad == 516L, cohort_counts$tcga_lusc == 501L,
  cohort_counts$gse_luad == 106L, cohort_counts$gse_lusc == 67L,
  cohort_counts$luad_tf == 158L, cohort_counts$lusc_tf == 193L
)

render_figure1a <- function(path) {
  open_device(path, 10.0, 4.7)
  grid.newpage()
  grid.text("A", unit(0.015, "npc"), unit(0.975, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold"))

  grid.rect(unit(0.34, "npc"), unit(0.55, "npc"), unit(0.61, "npc"), unit(0.67, "npc"),
            gp = gpar(fill = "white", col = "#222222", lwd = 0.9))
  grid.rect(unit(0.815, "npc"), unit(0.55, "npc"), unit(0.285, "npc"), unit(0.67, "npc"),
            gp = gpar(fill = "white", col = "#222222", lwd = 0.9))
  grid.text("RNA-seq (TCGA lung cancer)", unit(0.34, "npc"), unit(0.855, "npc"),
            gp = gpar(fontsize = 10.5))
  grid.text("RNA-seq (GSE81089)", unit(0.815, "npc"), unit(0.855, "npc"),
            gp = gpar(fontsize = 10.5))

  draw_person(0.065, 0.772, 0.85)
  grid.text(sprintf("%s tumors\n%s LUAD / %s LUSC",
                    fmt_n(cohort_counts$tcga_luad + cohort_counts$tcga_lusc),
                    fmt_n(cohort_counts$tcga_luad), fmt_n(cohort_counts$tcga_lusc)),
            unit(0.09, "npc"), unit(0.77, "npc"), just = "left",
            gp = gpar(fontsize = 8.2, lineheight = 1.05))
  draw_arrow(0.145, 0.715, 0.145, 0.665)
  grid.text("ARACNe3", unit(0.145, "npc"), unit(0.635, "npc"), gp = gpar(fontsize = 9.2))
  draw_network(0.145, 0.47, 0.12, 0.15)
  grid.text("Lung cancer interactome", unit(0.145, "npc"), unit(0.345, "npc"), gp = gpar(fontsize = 8.1))

  draw_arrow(0.225, 0.51, 0.285, 0.56)
  grid.text("VIPER", unit(0.395, "npc"), unit(0.655, "npc"), gp = gpar(fontsize = 9.2))
  draw_gradient(0.395, 0.605, 0.19, 0.028)
  grid.text("Differential TF activity signature", unit(0.395, "npc"), unit(0.57, "npc"), gp = gpar(fontsize = 7.3))
  draw_box(0.335, 0.475, 0.12, 0.09, "Multi-sample\nanalysis", fontsize = 7.7)
  draw_box(0.48, 0.475, 0.12, 0.09, "Single-sample\nanalysis", fontsize = 7.7)
  grid.text("LUAD-specific TFs\nLUSC-specific TFs", unit(0.275, "npc"), unit(0.365, "npc"),
            just = "left", gp = gpar(fontsize = 7.4, lineheight = 1.0))
  draw_person(0.445, 0.37, 0.60)
  draw_mouse(0.505, 0.367, 0.85)
  draw_dish(0.565, 0.368, 0.85)
  grid.text(fmt_n(nrow(tcga_manifest)), unit(0.445, "npc"), unit(0.318, "npc"), gp = gpar(fontsize = 6.8))
  grid.text(sprintf("%s PDXs", fmt_n(cohort_counts$pdx)), unit(0.505, "npc"), unit(0.318, "npc"), gp = gpar(fontsize = 6.8))
  grid.text(sprintf("%s lines", fmt_n(cohort_counts$cell_line)), unit(0.565, "npc"), unit(0.318, "npc"), gp = gpar(fontsize = 6.8))

  draw_person(0.695, 0.772, 0.80)
  grid.text(sprintf("%s tumors\n%s LUAD / %s LUSC",
                    fmt_n(cohort_counts$gse_luad + cohort_counts$gse_lusc),
                    fmt_n(cohort_counts$gse_luad), fmt_n(cohort_counts$gse_lusc)),
            unit(0.72, "npc"), unit(0.77, "npc"), just = "left",
            gp = gpar(fontsize = 7.8, lineheight = 1.05))
  draw_arrow(0.815, 0.715, 0.815, 0.67)
  grid.text("ARACNe3", unit(0.815, "npc"), unit(0.64, "npc"), gp = gpar(fontsize = 9.0))
  grid.text("Independent lung\ninteractome", unit(0.815, "npc"), unit(0.57, "npc"),
            gp = gpar(fontsize = 7.8, lineheight = 1.0))
  draw_arrow(0.815, 0.52, 0.815, 0.485)
  grid.text("VIPER", unit(0.815, "npc"), unit(0.455, "npc"), gp = gpar(fontsize = 9.0))
  draw_box(0.75, 0.36, 0.105, 0.09, "Multi-sample\nanalysis", fontsize = 7.2)
  draw_box(0.88, 0.36, 0.105, 0.09, "Single-sample\nanalysis", fontsize = 7.2)

  draw_box(
    0.52, 0.105, 0.34, 0.095,
    sprintf("%s LUAD-specific TFs (Figure 2 input)\n%s independently replicated in GSE81089",
            fmt_n(cohort_counts$luad_tf), fmt_n(cohort_counts$replicated_luad_tf)),
    fill = "#EAF3F8", fontsize = 8.3, lwd = 1.0
  )
  draw_poly_arrow(c(0.34, 0.34, 0.35), c(0.215, 0.105, 0.105))
  draw_poly_arrow(c(0.815, 0.815, 0.69), c(0.315, 0.105, 0.105), dashed = TRUE)
  dev.off()
}

render_figure1a(file.path(figures, "Figure1A_workflow.pdf"))
render_figure1a(file.path(figures, "Figure1A_workflow.png"))

unavailable <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x %in% c("NA", "N/A", "Not Reported", "not reported", "Unknown", "-")] <- "Unavailable"
  x
}

make_colors <- function(values, preferred = NULL) {
  levels <- sort(unique(as.character(values)))
  palette <- setNames(grDevices::hcl.colors(max(3, length(levels)), "Dark 3")[seq_along(levels)], levels)
  if (!is.null(preferred)) {
    hit <- intersect(names(preferred), names(palette))
    palette[hit] <- preferred[hit]
  }
  if ("Unavailable" %in% names(palette)) palette[["Unavailable"]] <- "#D9D9D9"
  palette
}

preferred <- list(
  Histology = c(LUAD = "#B2182B", LUSC = "#2166AC"),
  Sex = c(Female = "#CC79A7", Male = "#0072B2"),
  Age = c(`<=60` = "#F1EEF6", `61-70` = "#BDC9E1", `>70` = "#74A9CF"),
  `Vital status` = c(Dead = "#8E0152", Alive = "#F7F7F7"),
  Stage = c(I = "#EDF8E9", II = "#BAE4B3", III = "#74C476", IV = "#238B45"),
  `Pack-years` = c(`0` = "#F7FCF5", `>0-20` = "#C7E9C0", `>20-40` = "#74C476", `>40` = "#238B45"),
  `Expression subtype` = c(
    `LUAD: TRU` = "#1B9E77", `LUAD: PI` = "#D95F02", `LUAD: PP` = "#7570B3",
    `LUSC: Basal` = "#66C2A5", `LUSC: Classical` = "#FC8D62",
    `LUSC: Primitive` = "#8DA0CB", `LUSC: Secretory` = "#E78AC3"
  )
)

render_figure1c <- function(path) {
  activity <- read_matrix(file.path(tables, "Figure1C_TCGA_activity_matrix.tsv.gz"))
  activity <- activity[apply(activity, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  display_matrix <- t(scale(t(activity)))
  display_matrix[!is.finite(display_matrix)] <- 0

  row_mean <- rowMeans(display_matrix)
  row_sd <- apply(display_matrix, 1, sd)
  if (max(abs(row_mean)) > 1e-8 || max(abs(row_sd - 1)) > 1e-8) {
    stop("Figure 1C display matrix is not row-standardized")
  }

  anno <- fread(file.path(annotations_dir, "tcga_figure1_annotations.tsv"))
  anno <- anno[match(colnames(display_matrix), sample_id)]
  if (anyNA(anno$sample_id)) stop("Figure 1C annotation/sample mismatch")
  anno_df <- data.frame(
    Histology = unavailable(anno$group),
    Sex = unavailable(anno$sex),
    Age = unavailable(anno$age_group),
    `Vital status` = unavailable(anno$vital_status),
    Stage = unavailable(anno$stage),
    `Pack-years` = unavailable(anno$pack_year_group),
    `Expression subtype` = unavailable(anno$published_expression_subtype),
    Oncogene = unavailable(anno$oncogene_mutation),
    Suppressor = unavailable(anno$suppressor_pathway_mutation),
    check.names = FALSE
  )
  rownames(anno_df) <- colnames(display_matrix)
  anno_colors <- setNames(lapply(names(anno_df), function(name) make_colors(anno_df[[name]], preferred[[name]])), names(anno_df))
  legend_params <- setNames(lapply(names(anno_df), function(name) list(
    title = name, title_gp = gpar(fontsize = 7.2, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.3), grid_width = unit(2.6, "mm"), grid_height = unit(2.6, "mm")
  )), names(anno_df))
  top <- HeatmapAnnotation(
    df = anno_df, col = anno_colors, na_col = "#D9D9D9",
    simple_anno_size = unit(2.1, "mm"),
    show_annotation_name = TRUE, annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 6.8), gap = unit(0.25, "mm"),
    annotation_legend_param = legend_params
  )

  row_group <- selected$discovery_category[match(rownames(display_matrix), selected$TF)]
  if (anyNA(row_group)) stop("Figure 1C TF category mismatch")
  category_colors <- c(LUAD = "#455F7D", LUSC = "#9EC3E6")
  left <- rowAnnotation(
    `VIPER TF category` = factor(row_group, levels = c("LUAD", "LUSC")),
    col = list(`VIPER TF category` = category_colors),
    simple_anno_size = unit(3.0, "mm"), show_annotation_name = FALSE,
    annotation_legend_param = list(`VIPER TF category` = list(
      title = "VIPER TF category", title_gp = gpar(fontsize = 7.2, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5), grid_width = unit(2.8, "mm"), grid_height = unit(2.8, "mm")
    ))
  )

  activity_colors <- colorRamp2(
    c(-4, -2, 0, 2, 4),
    c("#244B8A", "#74A9CF", "#F7F7F7", "#F4A582", "#B2182B")
  )
  column_tree <- hclust(dist(t(display_matrix)), method = "complete")
  heatmap <- Heatmap(
    display_matrix, name = "TF activity", col = activity_colors,
    cluster_rows = FALSE, cluster_columns = column_tree,
    show_row_names = FALSE, show_column_names = FALSE,
    use_raster = TRUE, raster_quality = 3,
    top_annotation = top, left_annotation = left,
    column_title = sprintf("n = %s", fmt_n(ncol(display_matrix))),
    column_title_side = "bottom", column_title_gp = gpar(fontsize = 8.2),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7.2, fontface = "bold"), labels_gp = gpar(fontsize = 6.5),
      at = c(-4, -2, 0, 2, 4), legend_height = unit(24, "mm")
    )
  )

  open_device(path, 12.4, 7.2)
  grid.newpage()
  draw(
    heatmap, merge_legends = FALSE,
    heatmap_legend_side = "right", annotation_legend_side = "right",
    padding = unit(c(7, 6, 7, 7), "mm")
  )
  grid.text("C", unit(1.5, "mm"), unit(1, "npc") - unit(1.5, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
  dev.off()

  data.table(
    metric = c("samples", "TFs", "LUAD_specific_TFs", "LUSC_specific_TFs",
               "max_abs_row_mean", "max_abs_row_sd_minus_1"),
    value = c(ncol(display_matrix), nrow(display_matrix), sum(row_group == "LUAD"), sum(row_group == "LUSC"),
              max(abs(row_mean)), max(abs(row_sd - 1)))
  )
}

heatmap_receipt <- render_figure1c(file.path(figures, "Figure1C_TCGA_TF_activity_annotated.pdf"))
invisible(render_figure1c(file.path(figures, "Figure1C_TCGA_TF_activity_annotated.png")))

matrix_samples <- function(path) names(fread(path, nrows = 0, check.names = FALSE))[-1]
tcga_counts_samples <- unlist(lapply(c("LUAD", "LUSC"), function(histology) {
  matrix_samples(file.path(raw, "tcga", sprintf("TCGA-%s.star_counts.tsv.gz", histology)))
}), use.names = FALSE)
tcga_tpm_samples <- unlist(lapply(c("LUAD", "LUSC"), function(histology) {
  matrix_samples(file.path(raw, "tcga", sprintf("TCGA-%s.star_tpm.tsv.gz", histology)))
}), use.names = FALSE)
tcga_common <- intersect(tcga_counts_samples, tcga_tpm_samples)
tcga_primary <- tcga_common[substr(tcga_common, 14, 15) == "01"]
tcga_unique_patients <- unique(substr(tcga_primary, 1, 12))

parse_gse_soft_audit <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1] - 1L, length(lines))
  rbindlist(lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    title_line <- block[grep("^!Sample_title = ", block)[1]]
    hist_line <- block[grep("^!Sample_characteristics_ch1 = histology: ", block)[1]]
    data.table(
      sample_id = sub("^!Sample_title = ", "", title_line),
      histology_code = if (length(hist_line)) sub("^!Sample_characteristics_ch1 = histology: ", "", hist_line) else NA_character_
    )
  }), fill = TRUE)
}

gse_soft <- parse_gse_soft_audit(file.path(raw, "gse81089", "GSE81089_family.soft.gz"))
gse_expression_samples <- matrix_samples(file.path(raw, "gse81089", "GSE81089_readcounts_featurecounts.tsv.gz"))
gse_soft[, tumor := grepl("T(_|$)", sample_id)]
gse_soft[, target_histology := tumor & histology_code %in% c("1", "2")]
gse_target <- gse_soft[target_histology == TRUE, sample_id]
gse_target_with_expression <- intersect(gse_target, gse_expression_samples)

flow_counts <- list(
  tcga_all = length(tcga_common),
  tcga_primary = length(tcga_primary),
  tcga_non_primary = length(tcga_common) - length(tcga_primary),
  tcga_unique = length(tcga_unique_patients),
  tcga_duplicate = length(tcga_primary) - length(tcga_unique_patients),
  gse_all = nrow(gse_soft),
  gse_tumor = sum(gse_soft$tumor),
  gse_normal = sum(!gse_soft$tumor),
  gse_target_histology = sum(gse_soft$target_histology),
  gse_non_target_histology = sum(gse_soft$tumor & !gse_soft$target_histology),
  gse_target_with_expression = length(gse_target_with_expression),
  gse_unmatched = length(gse_target) - length(gse_target_with_expression)
)
stopifnot(
  flow_counts$tcga_all == 1141L, flow_counts$tcga_primary == 1029L,
  flow_counts$tcga_non_primary == 112L, flow_counts$tcga_unique == 1017L,
  flow_counts$tcga_duplicate == 12L,
  flow_counts$gse_all == 218L, flow_counts$gse_tumor == 199L,
  flow_counts$gse_normal == 19L, flow_counts$gse_target_histology == 175L,
  flow_counts$gse_non_target_histology == 24L,
  flow_counts$gse_target_with_expression == 173L, flow_counts$gse_unmatched == 2L
)

render_supplementary_1a <- function(path) {
  open_device(path, 7.6, 7.0)
  grid.newpage()
  grid.text("A", unit(0.02, "npc"), unit(0.98, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold"))
  draw_box(0.50, 0.91, 0.50, 0.105,
           sprintf("Total TCGA lung RNA-seq columns\nwith count and TPM data\nn = %s", fmt_n(flow_counts$tcga_all)), fontsize = 8.8)
  draw_arrow(0.50, 0.855, 0.50, 0.805)
  draw_diamond(0.50, 0.73, 0.36, 0.15, "Primary-tumor\nsample-type filter", fontsize = 8.2)
  draw_poly_arrow(c(0.32, 0.20, 0.20), c(0.73, 0.66, 0.61))
  draw_box(0.20, 0.52, 0.28, 0.13,
           sprintf("Excluded non-primary\nsamples\nn = %s", fmt_n(flow_counts$tcga_non_primary)), fontsize = 8.1)
  draw_poly_arrow(c(0.68, 0.72, 0.72), c(0.73, 0.66, 0.61))
  draw_box(0.72, 0.52, 0.30, 0.13,
           sprintf("Primary-tumor columns\nn = %s", fmt_n(flow_counts$tcga_primary)), fontsize = 8.3)
  draw_arrow(0.72, 0.455, 0.72, 0.405)
  draw_diamond(0.72, 0.34, 0.34, 0.13, "One primary tumor\nper patient", fontsize = 8.0)
  draw_poly_arrow(c(0.55, 0.33, 0.33), c(0.34, 0.28, 0.23))
  draw_box(0.33, 0.16, 0.28, 0.12,
           sprintf("Excluded duplicate\nprimary aliquots\nn = %s", fmt_n(flow_counts$tcga_duplicate)), fontsize = 8.0)
  draw_arrow(0.72, 0.275, 0.72, 0.225)
  draw_box(0.72, 0.16, 0.30, 0.12,
           sprintf("Final analysis cohort\nn = %s", fmt_n(nrow(tcga_manifest))), fontsize = 8.4)
  draw_poly_arrow(c(0.67, 0.58, 0.58), c(0.10, 0.075, 0.055))
  draw_poly_arrow(c(0.77, 0.86, 0.86), c(0.10, 0.075, 0.055))
  draw_box(0.57, 0.035, 0.22, 0.07, sprintf("LUAD  n = %s", fmt_n(cohort_counts$tcga_luad)),
           fill = "#455F7D", text_col = "white", fontsize = 8.0)
  draw_box(0.86, 0.035, 0.22, 0.07, sprintf("LUSC  n = %s", fmt_n(cohort_counts$tcga_lusc)),
           fill = "#9EC3E6", fontsize = 8.0)
  dev.off()
}

render_supplementary_1b <- function(path) {
  open_device(path, 7.6, 7.0)
  grid.newpage()
  grid.text("B", unit(0.02, "npc"), unit(0.98, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 18, fontface = "bold"))
  draw_box(0.50, 0.92, 0.50, 0.105,
           sprintf("Total GSE81089 RNA-seq samples\nwith featureCounts data\nn = %s", fmt_n(flow_counts$gse_all)), fontsize = 8.8)
  draw_arrow(0.50, 0.865, 0.50, 0.815)
  draw_diamond(0.50, 0.75, 0.34, 0.13, "Tumor versus\nmatched normal", fontsize = 8.2)
  draw_poly_arrow(c(0.33, 0.20, 0.20), c(0.75, 0.68, 0.63))
  draw_box(0.20, 0.56, 0.26, 0.12,
           sprintf("Matched normal tissue\nn = %s", fmt_n(flow_counts$gse_normal)), fontsize = 8.1)
  draw_poly_arrow(c(0.67, 0.70, 0.70), c(0.75, 0.68, 0.63))
  draw_box(0.70, 0.56, 0.28, 0.12,
           sprintf("Tumor samples\nn = %s", fmt_n(flow_counts$gse_tumor)), fontsize = 8.2)
  draw_arrow(0.70, 0.50, 0.70, 0.455)
  draw_diamond(0.70, 0.39, 0.34, 0.13, "Histology filter\nLUAD or LUSC", fontsize = 8.1)
  draw_poly_arrow(c(0.53, 0.31, 0.31), c(0.39, 0.33, 0.285))
  draw_box(0.31, 0.22, 0.29, 0.115,
           sprintf("Excluded large-cell/NOS\ntumors\nn = %s", fmt_n(flow_counts$gse_non_target_histology)), fontsize = 7.9)
  draw_arrow(0.70, 0.325, 0.70, 0.285)
  draw_box(0.70, 0.22, 0.29, 0.115,
           sprintf("LUAD/LUSC by metadata\nn = %s", fmt_n(flow_counts$gse_target_histology)), fontsize = 8.1)
  draw_arrow(0.70, 0.16, 0.70, 0.13)
  draw_diamond(0.70, 0.09, 0.32, 0.08, "Expression-column match", fontsize = 7.5)
  draw_poly_arrow(c(0.54, 0.39, 0.39), c(0.09, 0.055, 0.035))
  grid.text(sprintf("Unmatched IDs\nn = %s", fmt_n(flow_counts$gse_unmatched)),
            unit(0.34, "npc"), unit(0.025, "npc"), gp = gpar(fontsize = 7.2, lineheight = 1.0))
  draw_poly_arrow(c(0.72, 0.72), c(0.05, 0.035))
  draw_box(0.67, 0.025, 0.20, 0.05, sprintf("LUAD  n = %s", fmt_n(cohort_counts$gse_luad)),
           fill = "#455F7D", text_col = "white", fontsize = 7.5)
  draw_poly_arrow(c(0.80, 0.88, 0.88), c(0.09, 0.055, 0.035))
  draw_box(0.88, 0.025, 0.20, 0.05, sprintf("LUSC  n = %s", fmt_n(cohort_counts$gse_lusc)),
           fill = "#9EC3E6", fontsize = 7.5)
  dev.off()
}

render_supplementary_1a(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.pdf"))
render_supplementary_1a(file.path(figures, "SupplementaryFigure1A_TCGA_inclusion.png"))
render_supplementary_1b(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.pdf"))
render_supplementary_1b(file.path(figures, "SupplementaryFigure1B_GSE81089_inclusion.png"))

flow_receipt <- data.table(
  metric = names(flow_counts),
  value = as.numeric(unlist(flow_counts, use.names = FALSE))
)
receipt <- rbind(
  data.table(section = "Figure1C", heatmap_receipt),
  data.table(section = "SupplementaryFigure1", flow_receipt),
  fill = TRUE
)
fwrite(receipt, file.path(tables, "Figure1_anchor_style_render_receipt.tsv"), sep = "\t")

cat(sprintf(
  "Rendered anchor-style Figure 1A, Figure 1C, Supplementary Figure 1A and 1B. Figure 1C: %d TFs x %d samples.\n",
  as.integer(heatmap_receipt[metric == "TFs", value]),
  as.integer(heatmap_receipt[metric == "samples", value])
))

#!/usr/bin/env Rscript

# Flat, journal-style schematics for panels for which the TNBC capsule did not
# contain an authoring script.  The information flow follows the frozen LUAD
# manifests, while the visual grammar follows the final TNBC figure: thin black
# rules, square boxes, restrained cohort colours, no dashboard cards, pills,
# page headers, subtitles, or decorative rails.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

argv_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", argv_full, value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script directory.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
source(file.path(script_dir, "tnbc_anchor_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: render_anchor_schematics.R FIGURE1_ROOT FIGURE2_ROOT OUTPUT_ROOT", call. = FALSE)
}

figure1_root <- normalizePath(args[[1L]], mustWork = TRUE)
figure2_root <- normalizePath(args[[2L]], mustWork = TRUE)
output_root <- normalizePath(args[[3L]], mustWork = FALSE)
atomic_dir <- file.path(output_root, "anchor_style", "atomic")
audit_dir <- file.path(output_root, "anchor_style", "audit")
dir.create(atomic_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop(sprintf("Required input is missing: %s", path), call. = FALSE)
  fread(path)
}

processed <- file.path(figure1_root, "data", "processed")
f1_tables <- file.path(figure1_root, "results", "tables")
tcga <- read_required(file.path(processed, "tcga_manifest.tsv"))
gse <- read_required(file.path(processed, "gse81089_manifest.tsv"))
micro <- read_required(file.path(processed, "gse41271_manifest.tsv"))
micro_all <- read_required(file.path(processed, "gse41271_all_samples.tsv"))
pdmr <- read_required(file.path(processed, "pdmr_manifest.tsv"))
depmap <- read_required(file.path(processed, "depmap_manifest.tsv"))
selected <- read_required(file.path(f1_tables, "tcga_specific_TFs.tsv"))
replication <- read_required(file.path(f1_tables, "figure1_cross_platform_replication_summary.tsv"))

counts <- list(
  tcga_luad = sum(tcga$group == "LUAD"),
  tcga_lusc = sum(tcga$group == "LUSC"),
  gse_luad = sum(gse$group == "LUAD"),
  gse_lusc = sum(gse$group == "LUSC"),
  micro_luad = sum(micro$group == "LUAD"),
  micro_lusc = sum(micro$group == "LUSC"),
  micro_other = sum(micro_all$group == "OTHER"),
  pdx = nrow(pdmr),
  cell_line = nrow(depmap),
  luad_tf = sum(selected$discovery_category == "LUAD"),
  gse_replicated = replication[specificity == "LUAD", TCGA_replicated_GSE81089],
  micro_replicated = replication[specificity == "LUAD", TCGA_replicated_GSE41271],
  both_replicated = replication[specificity == "LUAD", TCGA_replicated_both]
)

stopifnot(
  counts$tcga_luad == 516L,
  counts$tcga_lusc == 501L,
  counts$gse_luad == 108L,
  counts$gse_lusc == 67L,
  counts$micro_luad == 183L,
  counts$micro_lusc == 80L,
  counts$micro_other == 12L,
  counts$luad_tf == 187L
)

ink <- "#111111"
muted <- "#5C6165"
rule <- "#2C2C2C"
light_rule <- "#B9BDC0"
luad_col <- unname(anchor_colours$phenotype[["LUAD"]])
lusc_col <- unname(anchor_colours$phenotype[["LUSC"]])
patient_col <- unname(anchor_colours$atac_system[["Patient"]])
pdx_col <- unname(anchor_colours$atac_system[["PDX"]])
cell_col <- unname(anchor_colours$atac_system[["CellLine"]])

open_device <- function(path, width, height) {
  anchor_cairo_pdf(path, width = width, height = height)
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = anchor_font_family)))
}

close_device <- function() {
  popViewport()
  dev.off()
}

box <- function(x, y, w, h, label = NULL, fill = "white", border = rule,
                fontsize = 8.2, bold = FALSE, lineheight = 1.05) {
  grid.rect(unit(x, "npc"), unit(y, "npc"), unit(w, "npc"), unit(h, "npc"),
            gp = gpar(fill = fill, col = border, lwd = 0.75))
  if (!is.null(label)) {
    grid.text(label, unit(x, "npc"), unit(y, "npc"),
              gp = gpar(col = ink, fontsize = fontsize,
                        fontface = if (bold) "bold" else "plain",
                        lineheight = lineheight))
  }
}

txt <- function(label, x, y, fontsize = 8, bold = FALSE, col = ink,
                just = "centre", rot = 0, lineheight = 1.05) {
  grid.text(label, unit(x, "npc"), unit(y, "npc"), just = just, rot = rot,
            gp = gpar(col = col, fontsize = fontsize,
                      fontface = if (bold) "bold" else "plain",
                      lineheight = lineheight))
}

arrow_segment <- function(x0, y0, x1, y1, col = rule, dashed = FALSE) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             arrow = arrow(length = unit(0.065, "in"), type = "closed"),
             gp = gpar(col = col, lwd = 0.8, lty = if (dashed) 2 else 1))
}

plain_segment <- function(x0, y0, x1, y1, col = rule, dashed = FALSE, lwd = 0.7) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             gp = gpar(col = col, lwd = lwd, lty = if (dashed) 2 else 1))
}

panel_letter <- function(letter) {
  txt(letter, 0.015, 0.98, fontsize = 17, bold = TRUE, just = c("left", "top"))
}

save_png_copy <- function(pdf_path, png_path) {
  # Raster previews are generated later by the PDF QA pipeline.  Record the
  # intended path here rather than introducing a second rendering engine.
  invisible(c(pdf = pdf_path, png_preview = png_path))
}

render_figure1a <- function(path) {
  open_device(path, 9.2, 5.2)
  panel_letter("A")

  # Discovery arm.  All objects remain inside this frame; the earlier version
  # allowed the VIPER box to cross the validation-frame boundary.
  box(0.30, 0.56, 0.54, 0.70)
  txt("RNA-seq discovery", 0.30, 0.865, fontsize = 10.5, bold = TRUE)
  txt(sprintf("TCGA primary tumours  n = %d", nrow(tcga)), 0.30, 0.79,
      fontsize = 8.5)
  txt(sprintf("LUAD %d    LUSC %d", counts$tcga_luad, counts$tcga_lusc),
      0.30, 0.745, fontsize = 7.9, col = muted)
  arrow_segment(0.30, 0.705, 0.30, 0.645)
  box(0.19, 0.575, 0.18, 0.12, "ARACNe3\ninteractome", fontsize = 8.1)
  arrow_segment(0.285, 0.575, 0.35, 0.575)
  box(0.45, 0.575, 0.18, 0.12, "VIPER /\nmsVIPER", fontsize = 8.1)
  plain_segment(0.45, 0.515, 0.45, 0.475)
  plain_segment(0.19, 0.475, 0.45, 0.475)
  arrow_segment(0.19, 0.475, 0.19, 0.435)
  arrow_segment(0.45, 0.475, 0.45, 0.435)
  box(0.19, 0.38, 0.19, 0.105, "Multi-sample\nanalysis", fontsize = 7.6)
  box(0.45, 0.38, 0.19, 0.105, "Single-sample\nprojection", fontsize = 7.6)
  txt(sprintf("%d LUAD-specific TFs", counts$luad_tf), 0.19, 0.27,
      fontsize = 9.4, bold = TRUE, col = luad_col)
  txt(sprintf("%d tumours\n%d PDXs\n%d cell lines", nrow(tcga), counts$pdx, counts$cell_line),
      0.45, 0.27, fontsize = 7.7, lineheight = 1.12)

  # Independent validation arm: one RNA-seq and one microarray network.
  box(0.79, 0.56, 0.37, 0.70)
  txt("Independent patient validation", 0.79, 0.865, fontsize = 10.5, bold = TRUE)
  txt(sprintf("GSE81089 RNA-seq\n%d LUAD + %d LUSC", counts$gse_luad, counts$gse_lusc),
      0.69, 0.755, fontsize = 7.8)
  txt(sprintf("GSE41271 microarray\n%d LUAD + %d LUSC", counts$micro_luad, counts$micro_lusc),
      0.89, 0.755, fontsize = 7.8)
  arrow_segment(0.69, 0.69, 0.69, 0.625)
  arrow_segment(0.89, 0.69, 0.89, 0.625)
  box(0.69, 0.565, 0.16, 0.11, "Independent\nARACNe3 + VIPER", fontsize = 7.3)
  box(0.89, 0.565, 0.16, 0.11, "Independent\nARACNe3 + VIPER", fontsize = 7.3)
  txt(sprintf("%d/%d replicated", counts$gse_replicated, counts$luad_tf),
      0.69, 0.43, fontsize = 8.2, bold = TRUE, col = luad_col)
  txt(sprintf("%d/%d replicated", counts$micro_replicated, counts$luad_tf),
      0.89, 0.43, fontsize = 8.2, bold = TRUE, col = luad_col)
  plain_segment(0.69, 0.385, 0.69, 0.345)
  plain_segment(0.89, 0.385, 0.89, 0.345)
  plain_segment(0.69, 0.345, 0.89, 0.345)
  arrow_segment(0.79, 0.345, 0.79, 0.30)
  txt(sprintf("%d TFs supported by both cohorts", counts$both_replicated),
      0.79, 0.255, fontsize = 8.3, bold = TRUE)
  txt("Independent replication is reported separately; it does not filter Figure 2 input.",
      0.79, 0.19, fontsize = 7.1, col = muted)

  arrow_segment(0.19, 0.215, 0.19, 0.145)
  box(0.30, 0.09, 0.48, 0.09,
      sprintf("Frozen Figure 2 input: all %d TCGA-discovered LUAD TFs", counts$luad_tf),
      fill = "#F4F6F7", border = luad_col, fontsize = 8.5, bold = TRUE)
  close_device()
}

render_inclusion <- function(path, letter, title, start_n, steps, final_labels) {
  open_device(path, 5.2, 5.5)
  panel_letter(letter)
  txt(title, 0.51, 0.92, fontsize = 11, bold = TRUE)
  box(0.40, 0.80, 0.46, 0.095, sprintf("Available expression profiles\nn = %d", start_n),
      fontsize = 8.5)
  step_centres <- seq(0.64, 0.39, length.out = length(steps))
  previous_bottom <- 0.752
  for (index in seq_along(steps)) {
    step <- steps[[index]]
    centre <- step_centres[[index]]
    arrow_segment(0.40, previous_bottom, 0.40, centre + 0.047)
    box(0.40, centre, 0.46, 0.085, step$keep, fontsize = 7.8)
    if (!is.null(step$exclude)) {
      plain_segment(0.63, centre, 0.72, centre)
      txt(step$exclude, 0.75, centre, fontsize = 7.1, col = muted,
          just = "left", lineheight = 1.05)
    }
    previous_bottom <- centre - 0.047
  }
  final_centre <- 0.225
  arrow_segment(0.40, previous_bottom, 0.40, final_centre + 0.052)
  box(0.40, final_centre, 0.46, 0.095, "Frozen analytic cohort", fill = "#F4F6F7",
      border = luad_col, fontsize = 8.6, bold = TRUE)
  split_y <- 0.105
  plain_segment(0.40, final_centre - 0.052, 0.40, split_y + 0.035)
  plain_segment(0.27, split_y + 0.035, 0.53, split_y + 0.035)
  for (i in seq_along(final_labels)) {
    x <- seq(0.27, 0.53, length.out = length(final_labels))[[i]]
    plain_segment(x, split_y + 0.035, x, split_y)
    box(x, split_y - 0.035, 0.21, 0.070, final_labels[[i]],
        fill = if (i == 1L) luad_col else lusc_col, border = NA,
        fontsize = 7.7, bold = TRUE)
  }
  close_device()
}

render_figure2a <- function(path) {
  promoter_summary <- read_required(file.path(
    figure2_root, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"
  ))
  primary <- promoter_summary[analysis_definition == "PRIMARY_anchor_promoter_tcga_cpm1_both"]
  manifest <- read_required(file.path(
    figure2_root, "audit", "manifests", "figure2_atomic", "figure2_atomic_sample_manifest.tsv"
  ))
  system_n <- manifest[, .N, by = system]
  n_lookup <- setNames(system_n$N, system_n$system)

  # Figure 2A is a wide schematic in the published anchor.  A square canvas
  # forced the two evidence branches to cross labels in the previous build.
  open_device(path, 7.5, 3.5)
  panel_letter("A")
  txt("ATAC-seq (LUAD)", 0.68, 0.91, fontsize = 11, bold = TRUE)
  txt(sprintf("%d primary tumours", n_lookup[["patient"]]), 0.08, 0.78,
      fontsize = 8.4, col = patient_col, just = "left")
  txt(sprintf("%d LUAD PDX models", n_lookup[["PDX"]]), 0.08, 0.65,
      fontsize = 8.4, col = pdx_col, just = "left")
  txt(sprintf("%d LUAD cell lines", n_lookup[["cell_line"]]), 0.08, 0.52,
      fontsize = 8.4, col = cell_col, just = "left")
  arrow_segment(0.38, 0.65, 0.53, 0.65)
  box(0.69, 0.65, 0.25, 0.24, "Accessible chromatin\nprofiling", fontsize = 9)

  box(
    0.24, 0.27, 0.38, 0.14,
    sprintf("VIPER-identified LUAD-specific TFs  (n = %d)", nrow(primary)),
    fill = "#F4F6F7", border = rule, fontsize = 8.4, bold = TRUE
  )

  # Both candidate identity and accessible regions are required for the two
  # downstream readouts; the converging junction makes that logic explicit.
  plain_segment(0.43, 0.27, 0.52, 0.27)
  plain_segment(0.69, 0.53, 0.69, 0.43)
  arrow_segment(0.69, 0.43, 0.52, 0.34)
  plain_segment(0.52, 0.18, 0.52, 0.40)
  arrow_segment(0.52, 0.40, 0.65, 0.40)
  arrow_segment(0.52, 0.18, 0.65, 0.18)
  txt("Promoter accessibility", 0.79, 0.40, fontsize = 8.8)
  txt("HOMER genome-wide enrichment\nof TF binding-site motifs", 0.79, 0.18, fontsize = 8.8)
  close_device()
}

render_supp4a <- function(path) {
  promoter <- read_required(file.path(
    figure2_root, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"
  ))
  promoter <- promoter[analysis_definition == "PRIMARY_anchor_promoter_tcga_cpm1_both"]
  inventory <- read_required(file.path(
    figure2_root, "results", "motif_gate", "hc_tf_motif_inventory.tsv"
  ))
  motif <- read_required(file.path(
    figure2_root, "results", "motif_primary", "anchor_consensus_author",
    "anchor_consensus_hc_tf_triple_system_summary.tsv"
  ))
  n_open <- promoter[triple_system_promoter_accessible == "TRUE", .N]
  n_hc <- promoter[HC_TF_promoter_activity_definition == "TRUE", .N]
  n_testable <- inventory[motif_testable == "TRUE", uniqueN(TF)]
  motif_col <- intersect(c("triple_system_fixed_original_denominator", "triple_system_support"), names(motif))
  n_triple <- if (length(motif_col)) motif[get(motif_col[[1L]]) == "TRUE", uniqueN(TF)] else NA_integer_

  open_device(path, 6.2, 3.6)
  panel_letter("A")
  stages <- list(
    c(sprintf("Figure 1 LUAD-specific TFs\nn = %d", nrow(promoter)), 0.13),
    c(sprintf("Promoter open in at least half\nof samples in all three systems\nn = %d", n_open), 0.37),
    c(sprintf("Exclude only TFs with negative\nmean NES in all three systems\nn = %d HC-TFs", n_hc), 0.63),
    c(sprintf("Known motif available\nn = %d", n_testable), 0.86)
  )
  for (i in seq_along(stages)) {
    label <- stages[[i]][[1L]]
    x <- as.numeric(stages[[i]][[2L]])
    box(x, 0.58, 0.20, 0.26, label,
        fill = if (i == 3L) "#F3E8EC" else "white",
        border = if (i == 3L) unname(anchor_colours$regulon[["activated_shared"]]) else rule,
        fontsize = 7.7, bold = i == 3L)
    if (i < length(stages)) arrow_segment(x + 0.11, 0.58, as.numeric(stages[[i + 1L]][[2L]]) - 0.11, 0.58)
  }
  txt(sprintf("Three-system motif support: %s TFs", ifelse(is.na(n_triple), "not resolved", n_triple)),
      0.86, 0.30, fontsize = 8.2, bold = TRUE)
  txt("Motif support annotates the HC-TF set; it does not redefine downstream eligibility.",
      0.50, 0.13, fontsize = 7.5, col = muted)
  close_device()
}

render_figure1a(file.path(atomic_dir, "Figure1A_anchor_schematic.pdf"))
render_inclusion(
  file.path(atomic_dir, "SupplementaryFigure1A_TCGA_inclusion.pdf"),
  "A", "TCGA-LUAD/LUSC discovery cohort", nrow(tcga),
  list(
    list(keep = "Retain one primary tumour per patient", exclude = "non-primary / duplicate aliquots"),
    list(keep = "Require complete LUAD/LUSC histology", exclude = "ambiguous histology")
  ),
  c(sprintf("LUAD  n = %d", counts$tcga_luad), sprintf("LUSC  n = %d", counts$tcga_lusc))
)
render_inclusion(
  file.path(atomic_dir, "SupplementaryFigure1B_GSE81089_inclusion.pdf"),
  "B", "GSE81089 independent RNA-seq cohort", 218L,
  list(
    list(keep = "Retain tumour expression columns", exclude = "19 matched normal samples"),
    list(keep = "Canonicalize two GEO expression aliases", exclude = "no biological sample removed"),
    list(keep = "Require LUAD or LUSC histology", exclude = "24 large-cell / NOS tumours")
  ),
  c(sprintf("LUAD  n = %d", counts$gse_luad), sprintf("LUSC  n = %d", counts$gse_lusc))
)
render_inclusion(
  file.path(atomic_dir, "SupplementaryFigure1C_GSE41271_inclusion.pdf"),
  "C", "GSE41271 cross-platform microarray cohort", nrow(micro_all),
  list(list(keep = "Require exact LUAD or LUSC histology", exclude = sprintf("%d other histologies", counts$micro_other))),
  c(sprintf("LUAD  n = %d", counts$micro_luad), sprintf("LUSC  n = %d", counts$micro_lusc))
)
render_figure2a(file.path(atomic_dir, "Figure2A_anchor_schematic.pdf"))
render_supp4a(file.path(atomic_dir, "SupplementaryFigure4A_anchor_gate_schematic.pdf"))

receipt <- data.table(
  metric = names(counts),
  value = as.character(unlist(counts, use.names = FALSE))
)
fwrite(receipt, file.path(audit_dir, "anchor_schematic_counts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "anchor_schematic_sessionInfo.txt"))

message("Anchor-style schematics written to: ", atomic_dir)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: render_figure1_revision.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
processed <- file.path(root, "data", "processed")
tables <- file.path(root, "results", "tables")
figures <- file.path(root, "results", "figures")
audit <- file.path(root, "audit", "figure1_complete")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

tcga <- fread(file.path(processed, "tcga_manifest.tsv"))
gse <- fread(file.path(processed, "gse81089_manifest.tsv"))
micro <- fread(file.path(processed, "gse41271_manifest.tsv"))
pdmr <- fread(file.path(processed, "pdmr_manifest.tsv"))
depmap <- fread(file.path(processed, "depmap_manifest.tsv"))
selected <- fread(file.path(tables, "tcga_specific_TFs.tsv"))
replicated_gse <- fread(file.path(tables, "externally_replicated_TFs.tsv"))
replication <- fread(file.path(tables, "figure1_cross_platform_replication_summary.tsv"))
micro_all <- fread(file.path(processed, "gse41271_all_samples.tsv"))

count <- list(
  tcga_luad = sum(tcga$group == "LUAD"), tcga_lusc = sum(tcga$group == "LUSC"),
  gse_luad = sum(gse$group == "LUAD"), gse_lusc = sum(gse$group == "LUSC"),
  micro_luad = sum(micro$group == "LUAD"), micro_lusc = sum(micro$group == "LUSC"),
  micro_other = sum(micro_all$group == "OTHER"),
  pdx = nrow(pdmr), cell_line = nrow(depmap),
  luad_tf = sum(selected$discovery_category == "LUAD"),
  lusc_tf = sum(selected$discovery_category == "LUSC"),
  gse_luad_tf = sum(replicated_gse$discovery_category == "LUAD"),
  micro_luad_tf = replication[specificity == "LUAD", TCGA_replicated_GSE41271],
  both_luad_tf = replication[specificity == "LUAD", TCGA_replicated_both]
)
stopifnot(
  count$tcga_luad == 516L, count$tcga_lusc == 501L,
  count$gse_luad == 108L, count$gse_lusc == 67L,
  count$micro_luad == 183L, count$micro_lusc == 80L, count$micro_other == 12L,
  count$luad_tf == 158L
)

COL <- list(
  ink = "#17252D", sub = "#5D6A72", line = "#7C8991", border = "#CED7DC",
  blue = "#386A8C", blue_light = "#EAF2F7", teal = "#397D72",
  teal_light = "#E9F4F1", gold = "#A97429", gold_light = "#FBF3E6",
  rust = "#A04B3D", lusc = "#315A88"
)

open_device <- function(path, width, height) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) cairo_pdf(path, width = width, height = height)
  else png(path, width = width * 220, height = height * 220, res = 220)
}
arrow_line <- function(x0, y0, x1, y1, col = COL$line, dashed = FALSE) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             arrow = arrow(length = unit(0.075, "inches"), type = "closed"),
             gp = gpar(col = col, lwd = 0.85, lty = ifelse(dashed, 2, 1)))
}
lane <- function(y, accent, fill, step, cohort, modality, n_luad, n_lusc,
                 validation_text, projection = NULL) {
  grid.roundrect(unit(0.50, "npc"), unit(y, "npc"), unit(0.91, "npc"), unit(0.225, "npc"),
                 r = unit(3, "mm"), gp = gpar(fill = fill, col = COL$border, lwd = 0.8))
  grid.roundrect(unit(0.065, "npc"), unit(y, "npc"), unit(0.010, "npc"), unit(0.190, "npc"),
                 r = unit(1, "mm"), gp = gpar(fill = accent, col = NA))
  grid.roundrect(unit(0.160, "npc"), unit(y + 0.078, "npc"), unit(0.190, "npc"), unit(0.036, "npc"),
                 r = unit(2, "mm"), gp = gpar(fill = accent, col = NA))
  grid.text(step, unit(0.160, "npc"), unit(y + 0.078, "npc"),
            gp = gpar(fontsize = 6.5, fontface = "bold", col = "white"))
  grid.text(cohort, unit(0.155, "npc"), unit(y + 0.017, "npc"),
            gp = gpar(fontsize = 7.8, fontface = "bold", col = COL$ink))
  grid.text(sprintf("%s\n%d LUAD | %d LUSC", modality, n_luad, n_lusc),
            unit(0.155, "npc"), unit(y - 0.045, "npc"),
            gp = gpar(fontsize = 6.4, col = COL$sub, lineheight = 1.05))
  arrow_line(0.260, y, 0.320, y)
  grid.roundrect(unit(0.390, "npc"), unit(y, "npc"), unit(0.125, "npc"), unit(0.095, "npc"),
                 r = unit(2, "mm"), gp = gpar(fill = "white", col = accent, lwd = 0.7))
  grid.text("ARACNe3\nindependent regulon", unit(0.390, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = 6.4, col = COL$ink, lineheight = 1.03))
  arrow_line(0.455, y, 0.510, y)
  grid.roundrect(unit(0.580, "npc"), unit(y, "npc"), unit(0.125, "npc"), unit(0.095, "npc"),
                 r = unit(2, "mm"), gp = gpar(fill = "white", col = accent, lwd = 0.7))
  grid.text("VIPER / msVIPER\nexpression contrast", unit(0.580, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = 6.4, col = COL$ink, lineheight = 1.03))
  arrow_line(0.645, y, 0.700, y)
  grid.roundrect(unit(0.805, "npc"), unit(y, "npc"), unit(0.185, "npc"), unit(0.105, "npc"),
                 r = unit(2, "mm"), gp = gpar(fill = "white", col = accent, lwd = 0.8))
  grid.text(validation_text, unit(0.805, "npc"), unit(y, "npc"),
            gp = gpar(fontsize = 6.5, col = COL$ink, lineheight = 1.06))
  if (!is.null(projection)) {
    grid.text(projection, unit(0.590, "npc"), unit(y - 0.082, "npc"),
              gp = gpar(fontsize = 5.9, col = COL$sub))
  }
}

render_workflow <- function(path) {
  open_device(path, 8.7, 8.4)
  grid.newpage()
  grid.text("A", unit(0.018, "npc"), unit(0.985, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink))
  grid.text("LUAD transcriptional regulator discovery and independent validation",
            unit(0.070, "npc"), unit(0.957, "npc"), just = "left",
            gp = gpar(fontsize = 10.8, fontface = "bold", col = COL$ink))
  grid.text("Separate networks preserve cohort and platform independence; validation does not shrink the Figure 2 input",
            unit(0.070, "npc"), unit(0.928, "npc"), just = "left",
            gp = gpar(fontsize = 6.6, col = COL$sub))
  lane(0.765, COL$blue, "#F4F8FB", "01  DISCOVERY", "TCGA-LUAD/LUSC",
       sprintf("RNA-seq; n=%d", nrow(tcga)), count$tcga_luad, count$tcga_lusc,
       sprintf("%d LUAD-specific TFs\n%d LUSC-specific TFs", count$luad_tf, count$lusc_tf),
       sprintf("Single-sample projection: %d tumors | %d PDXs | %d cell lines",
               nrow(tcga), count$pdx, count$cell_line))
  lane(0.500, COL$teal, "#F3F9F7", "02  RNA-SEQ REPLICATION", "GSE81089",
       sprintf("RNA-seq; n=%d", nrow(gse)), count$gse_luad, count$gse_lusc,
       sprintf("%d/%d LUAD TFs\ndirection-concordant", count$gse_luad_tf, count$luad_tf))
  lane(0.235, COL$gold, "#FCF8EF", "03  CROSS-PLATFORM REPLICATION", "GSE41271",
       sprintf("Illumina microarray; n=%d", nrow(micro)), count$micro_luad, count$micro_lusc,
       sprintf("%d/%d LUAD TFs\n%d supported by both cohorts",
               count$micro_luad_tf, count$luad_tf, count$both_luad_tf))
  arrow_line(0.805, 0.700, 0.805, 0.640, COL$blue)
  arrow_line(0.805, 0.435, 0.805, 0.375, COL$teal, TRUE)
  arrow_line(0.805, 0.170, 0.805, 0.125, COL$gold, TRUE)
  grid.roundrect(unit(0.805, "npc"), unit(0.080, "npc"), unit(0.250, "npc"), unit(0.055, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = COL$blue_light, col = COL$blue, lwd = 0.9))
  grid.text(sprintf("Frozen Figure 2 input: all %d TCGA LUAD-specific TFs", count$luad_tf),
            unit(0.805, "npc"), unit(0.080, "npc"), gp = gpar(fontsize = 7.1, col = COL$ink))
  dev.off()
}

render_micro_flow <- function(path) {
  open_device(path, 6.7, 6.8)
  grid.newpage()
  grid.text("C", unit(0.025, "npc"), unit(0.978, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 17, fontface = "bold", col = COL$ink))
  grid.text("GSE41271 cross-platform cohort", unit(0.095, "npc"), unit(0.945, "npc"),
            just = "left", gp = gpar(fontsize = 11.2, fontface = "bold", col = COL$ink))
  grid.text("Illumina HumanWG-6 v3 microarray sample derivation", unit(0.095, "npc"), unit(0.905, "npc"),
            just = "left", gp = gpar(fontsize = 7.2, col = COL$sub))
  node <- function(y, title, n, note, accent, final = FALSE) {
    grid.roundrect(unit(0.405, "npc"), unit(y, "npc"), unit(0.50, "npc"), unit(0.112, "npc"),
                   r = unit(2.5, "mm"), gp = gpar(fill = ifelse(final, "#EEF7F4", "white"),
                                                   col = accent, lwd = ifelse(final, 1.2, 0.75)))
    grid.text(title, unit(0.275, "npc"), unit(y + 0.018, "npc"), just = "left",
              gp = gpar(fontsize = 7.2, fontface = "bold", col = COL$ink))
    grid.text(note, unit(0.275, "npc"), unit(y - 0.023, "npc"), just = "left",
              gp = gpar(fontsize = 6.2, col = COL$sub))
    grid.text(sprintf("n = %d", n), unit(0.570, "npc"), unit(y, "npc"),
              gp = gpar(fontsize = 8.2, fontface = "bold", col = accent))
  }
  node(0.790, "Available primary lung tumours", 275, "GEO Series Matrix and SOFT metadata", COL$gold)
  arrow_line(0.405, 0.730, 0.405, 0.660, COL$line)
  grid.roundrect(unit(0.405, "npc"), unit(0.625, "npc"), unit(0.260, "npc"), unit(0.046, "npc"),
                 r = unit(2, "mm"), gp = gpar(fill = COL$gold_light, col = COL$gold))
  grid.text("EXACT LUAD OR LUSC HISTOLOGY", unit(0.405, "npc"), unit(0.625, "npc"),
            gp = gpar(fontsize = 6.3, fontface = "bold", col = COL$gold))
  grid.text(sprintf("Excluded other histologies  n = %d", count$micro_other),
            unit(0.745, "npc"), unit(0.625, "npc"), just = "left",
            gp = gpar(fontsize = 6.2, col = COL$sub))
  arrow_line(0.405, 0.590, 0.405, 0.515, COL$line)
  node(0.455, "Frozen analytic cohort", nrow(micro), "One GEO accession per tumour", COL$gold, TRUE)
  grid.lines(unit(c(0.405, 0.295, 0.295), "npc"), unit(c(0.395, 0.325, 0.270), "npc"),
             gp = gpar(col = COL$gold, lwd = 0.9))
  grid.lines(unit(c(0.405, 0.515, 0.515), "npc"), unit(c(0.395, 0.325, 0.270), "npc"),
             gp = gpar(col = COL$gold, lwd = 0.9))
  grid.roundrect(unit(0.295, "npc"), unit(0.225, "npc"), unit(0.205, "npc"), unit(0.065, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = COL$blue, col = NA))
  grid.text(sprintf("LUAD  n = %d", count$micro_luad), unit(0.295, "npc"), unit(0.225, "npc"),
            gp = gpar(fontsize = 7.5, fontface = "bold", col = "white"))
  grid.roundrect(unit(0.515, "npc"), unit(0.225, "npc"), unit(0.205, "npc"), unit(0.065, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = COL$lusc, col = NA))
  grid.text(sprintf("LUSC  n = %d", count$micro_lusc), unit(0.515, "npc"), unit(0.225, "npc"),
            gp = gpar(fontsize = 7.5, fontface = "bold", col = "white"))
  grid.text("Probe-to-gene rule: highest-MAD single-symbol probe; network inferred independently",
            unit(0.405, "npc"), unit(0.115, "npc"),
            gp = gpar(fontsize = 6.2, col = COL$sub))
  dev.off()
}

render_workflow(file.path(figures, "Figure1A_workflow_final.pdf"))
render_workflow(file.path(figures, "Figure1A_workflow_final.png"))
render_micro_flow(file.path(figures, "SupplementaryFigure1C_GSE41271_inclusion.pdf"))
render_micro_flow(file.path(figures, "SupplementaryFigure1C_GSE41271_inclusion.png"))

receipt <- data.table(
  metric = names(count),
  value = as.integer(unlist(count)),
  source = "frozen manifests and cross-platform replication summary"
)
fwrite(receipt, file.path(audit, "figure1_complete_counts.tsv"), sep = "\t")
cat(sprintf("Rendered final Figure 1A with %d/%d GSE81089 and %d/%d GSE41271 LUAD TF replication.\n",
            count$gse_luad_tf, count$luad_tf, count$micro_luad_tf, count$luad_tf))


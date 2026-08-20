#!/usr/bin/env Rscript

# Render LUAD Figure 1B-E, Supplementary Figure 2A-F, and Figure 2B-E with
# the plotting grammar used by the TNBC Code Ocean v1 capsule.
#
# This is a rendering wrapper only.  It consumes frozen LUAD results and does
# not rerun ARACNe3, VIPER, ATAC-seq processing, promoter tests, or HOMER.
# Missing scientific inputs are never reconstructed from panel summaries.
# Every missing/partial input is written to anchor_style/audit/input_gaps.tsv.
#
# CLI
# ----
# Rscript render_anchor_figure12.R \
#   FIGURE1_ROOT FIGURE2_ROOT PUBLICATION_ROOT OUTPUT_ROOT [CASE_TFS]
#
# FIGURE1_ROOT
#   Frozen LUAD Figure 1 run, normally execution/luad-figure1-v2.
# FIGURE2_ROOT
#   Frozen LUAD Figure 2 run, normally execution/luad-figure2.
# PUBLICATION_ROOT
#   Publication rebuild run, normally execution/luad-publication-rebuild.
# OUTPUT_ROOT
#   Normally execution/luad-publication-rebuild/results.  Files are written
#   below OUTPUT_ROOT/anchor_style/{atomic,tables,audit}.
# CASE_TFS
#   Optional comma-separated, predeclared TFs for the Figure 2E case panel.
#   The default AUTO deterministically uses all triple-system primary hits,
#   then fills to at most six TFs by overall median LOR.  The exact choice and
#   rule are written to Figure2E_case_examples_used.tsv.
#
# Frozen input inventory (relative to the roots above)
# ----------------------------------------------------
# Figure 1:
#   results/tables/figure1b_tcga_limma_msviper.tsv
#   results/tables/gse81089_limma_msviper.tsv
#   results/tables/tcga_specific_TFs.tsv
#   results/tables/Figure1{C,D,E}_*_activity_matrix.tsv.gz
#   data/processed/annotations/{tcga,pdmr,depmap}_figure1_annotations.tsv
#   results/gse81089/GSE81089_viper_activity.tsv.gz  [needed for fourth S2F]
# Figure 2 (direct results first; publication supplementary-data copies are
# accepted as immutable fallbacks):
#   results/promoter_gate/{system_tf_promoter_activity_summary,
#     tf_promoter_gate_summary}.tsv
#   results/motif_primary/anchor_consensus_author/
#     anchor_consensus_{sample_tf_best_motif,
#     system_tf_motif_summary,hc_tf_triple_system_summary}.tsv
#   results/motif_gate/hc_tf_motif_inventory.tsv
#   publication results/supplementary_data/Data_05/
#     sample_tf_promoter_accessibility.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

argv_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", argv_full, value = TRUE)
if (!length(script_arg)) {
  stop("Cannot resolve the script directory from commandArgs().", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
source(file.path(script_dir, "tnbc_anchor_theme.R"))

usage <- paste(
  "Usage: render_anchor_figure12.R FIGURE1_ROOT FIGURE2_ROOT",
  "PUBLICATION_ROOT OUTPUT_ROOT [CASE_TFS]"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
if (length(args) < 4L || length(args) > 5L) stop(usage, call. = FALSE)

fig1_root <- normalizePath(args[[1L]], mustWork = TRUE)
fig2_root <- normalizePath(args[[2L]], mustWork = TRUE)
publication_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_root <- args[[4L]]
case_tf_arg <- if (length(args) == 5L) args[[5L]] else "AUTO"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
atomic_dir <- file.path(output_root, "anchor_style", "atomic")
table_dir <- file.path(output_root, "anchor_style", "tables")
audit_dir <- file.path(output_root, "anchor_style", "audit")
for (path in c(atomic_dir, table_dir, audit_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

set.seed(2026)

# -------------------------------------------------------------------------
# Input discovery and auditable failure handling.
# -------------------------------------------------------------------------

f1_tables <- file.path(fig1_root, "results", "tables")
f1_processed <- file.path(fig1_root, "data", "processed")
pub_supp <- file.path(publication_root, "results", "supplementary_data")

input_specs <- list(
  f1_tcga_stats = file.path(f1_tables, "figure1b_tcga_limma_msviper.tsv"),
  f1_gse_stats = file.path(f1_tables, "gse81089_limma_msviper.tsv"),
  f1_tf_membership = file.path(f1_tables, "tcga_specific_TFs.tsv"),
  f1_tcga_activity = file.path(f1_tables, "Figure1C_TCGA_activity_matrix.tsv.gz"),
  f1_pdx_activity = file.path(f1_tables, "Figure1D_PDMR_PDX_activity_matrix.tsv.gz"),
  f1_cell_activity = file.path(f1_tables, "Figure1E_DepMap_22Q2_activity_matrix.tsv.gz"),
  f1_gse_activity = c(
    file.path(fig1_root, "results", "gse81089", "GSE81089_viper_activity.tsv.gz"),
    file.path(f1_tables, "GSE81089_viper_activity.tsv.gz")
  ),
  f1_tcga_annotation = file.path(f1_processed, "annotations", "tcga_figure1_annotations.tsv"),
  f1_pdx_annotation = file.path(f1_processed, "annotations", "pdmr_figure1_annotations.tsv"),
  f1_cell_annotation = file.path(f1_processed, "annotations", "depmap_figure1_annotations.tsv"),
  f1_tcga_manifest = c(
    file.path(f1_processed, "tcga_manifest.tsv"),
    file.path(pub_supp, "Data_01", "tcga_manifest.tsv")
  ),
  f1_gse_manifest = c(
    file.path(f1_processed, "gse81089_manifest.tsv"),
    file.path(pub_supp, "Data_01", "gse81089_manifest.tsv")
  ),
  f1_pdx_manifest = c(
    file.path(f1_processed, "pdmr_manifest.tsv"),
    file.path(pub_supp, "Data_01", "pdmr_manifest.tsv")
  ),
  f1_cell_manifest = c(
    file.path(f1_processed, "depmap_manifest.tsv"),
    file.path(pub_supp, "Data_01", "depmap_manifest.tsv")
  ),
  f2_promoter_sample = c(
    file.path(fig2_root, "results", "promoter_gate", "sample_tf_promoter_accessibility.tsv"),
    file.path(pub_supp, "Data_05", "sample_tf_promoter_accessibility.tsv")
  ),
  f2_promoter_system = c(
    file.path(fig2_root, "results", "promoter_gate", "system_tf_promoter_activity_summary.tsv"),
    file.path(pub_supp, "Data_03", "system_tf_promoter_activity_summary.tsv")
  ),
  f2_promoter_gate = c(
    file.path(fig2_root, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"),
    file.path(pub_supp, "Data_03", "tf_promoter_gate_summary.tsv")
  ),
  f2_motif_sample = c(
    file.path(fig2_root, "results", "motif_primary", "anchor_consensus_author",
              "anchor_consensus_sample_tf_best_motif.tsv"),
    file.path(pub_supp, "Data_05", "anchor_consensus_sample_tf_best_motif.tsv")
  ),
  f2_motif_system = c(
    file.path(fig2_root, "results", "motif_primary", "anchor_consensus_author",
              "anchor_consensus_system_tf_motif_summary.tsv"),
    file.path(pub_supp, "Data_05", "anchor_consensus_system_tf_motif_summary.tsv")
  ),
  f2_motif_triple = c(
    file.path(fig2_root, "results", "motif_primary", "anchor_consensus_author",
              "anchor_consensus_hc_tf_triple_system_summary.tsv"),
    file.path(pub_supp, "Data_05", "anchor_consensus_hc_tf_triple_system_summary.tsv")
  ),
  f2_motif_inventory = c(
    file.path(fig2_root, "results", "motif_gate", "hc_tf_motif_inventory.tsv"),
    file.path(pub_supp, "Data_05", "hc_tf_motif_inventory.tsv")
  )
)

input_inventory <- rbindlist(lapply(names(input_specs), function(key) {
  data.table(
    input_key = key,
    candidate_rank = seq_along(input_specs[[key]]),
    path = input_specs[[key]],
    exists = file.exists(input_specs[[key]])
  )
}))
fwrite(input_inventory, file.path(audit_dir, "expected_input_inventory.tsv"), sep = "\t")

find_input <- function(key, required = TRUE) {
  candidates <- input_specs[[key]]
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(normalizePath(hit[[1L]], mustWork = TRUE))
  if (required) {
    stop(
      sprintf("Missing input '%s'; checked: %s", key, paste(candidates, collapse = " | ")),
      call. = FALSE
    )
  }
  NA_character_
}

gap_rows <- list()
status_rows <- list()

record_gap <- function(panel, requirement, reason, severity = "PARTIAL") {
  gap_rows[[length(gap_rows) + 1L]] <<- data.table(
    panel = panel,
    severity = severity,
    requirement = requirement,
    reason = reason
  )
  invisible(NULL)
}

run_panel <- function(panel, code) {
  started <- Sys.time()
  artifacts <- tryCatch(
    {
      value <- force(code)
      value <- as.character(value)
      status_rows[[length(status_rows) + 1L]] <<- data.table(
        panel = panel,
        status = "generated",
        artifacts = paste(basename(value), collapse = ";"),
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
        message = ""
      )
      value
    },
    error = function(e) {
      record_gap(panel, "panel inputs or rendering dependency", conditionMessage(e), "ERROR")
      status_rows[[length(status_rows) + 1L]] <<- data.table(
        panel = panel,
        status = "failed",
        artifacts = "",
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
        message = conditionMessage(e)
      )
      character()
    }
  )
  invisible(artifacts)
}

require_panel_package <- function(package, panel) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("%s requires R package '%s'.", panel, package), call. = FALSE)
  }
  invisible(TRUE)
}

save_anchor_pdf <- function(plot, stem, width, height) {
  path <- file.path(atomic_dir, paste0(stem, ".pdf"))
  anchor_save_plot(plot, path, width = width, height = height)
  path
}

save_anchor_heatmap <- function(heatmap, stem, width, height) {
  path <- file.path(atomic_dir, paste0(stem, ".pdf"))
  anchor_cairo_pdf(path, width = width, height = height)
  tryCatch(
    ComplexHeatmap::draw(
      heatmap,
      merge_legends = TRUE,
      padding = unit(c(2, 2, 2, 2), "mm")
    ),
    finally = dev.off()
  )
  path
}

read_activity_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  if (ncol(x) < 2L) stop(sprintf("Activity matrix has fewer than two columns: %s", path))
  ids <- as.character(x[[1L]])
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

row_zscore <- function(x) {
  answer <- t(scale(t(x)))
  answer[!is.finite(answer)] <- 0
  answer
}

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

clean_annotation <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | tolower(x) %in% c("na", "n/a", "unknown", "not reported", "-")] <-
    "Unavailable"
  x
}

pick_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) return(rep("Unavailable", nrow(x)))
  clean_annotation(x[[hit[[1L]]]])
}

age_group <- function(x) {
  raw <- suppressWarnings(as.numeric(as.character(x)))
  answer <- ifelse(
    is.na(raw), "Unavailable",
    ifelse(raw <= 60, "<=60", ifelse(raw <= 70, "61-70", ">70"))
  )
  answer
}

annotation_pool <- c(
  "#49a12e", "#2836a0", "#85a6d7", "#a02774", "#f2a039", "#53bab0",
  "#875591", "#b9b037", "#D99AB1", "#8C3048", "#798774", "#D4947D",
  "#96B9D9", "#C23E34", "#73a12f", "#afd7a5", "#745570", "#D9BBA9"
)

fixed_annotation_colours <- list(
  Histology = c(anchor_colours$phenotype, Unavailable = "#e1e1dd"),
  Published_subtype = c(
    `LUAD: TRU` = "#49a12e", `LUAD: PI` = "#f2a039", `LUAD: PP` = "#875591",
    `LUSC: Basal` = "#53bab0", `LUSC: Classical` = "#D4947D",
    `LUSC: Primitive` = "#85a6d7", `LUSC: Secretory` = "#D99AB1",
    Unavailable = "#e1e1dd"
  ),
  Sex = c(Female = "#D99AB1", Male = "#2836a0", Unavailable = "#e1e1dd"),
  Age = c(`<=60` = "#d7afd2", `61-70` = "#a34d9d", `>70` = "#560f55",
          Unavailable = "#e1e1dd"),
  Stage = c(I = "#afd7a5", II = "#73a12f", III = "#49a146", IV = "#2b6519",
            Unavailable = "#e1e1dd"),
  Vital_status = c(Alive = "#D99AB1", Dead = "#8C3048", Unavailable = "#e1e1dd"),
  Smoking = c(Current = "#8C3048", Former = "#D99AB1", Never = "#85a6d7",
              Yes = "#8C3048", No = "#85a6d7", Unavailable = "#e1e1dd"),
  Primary_or_metastasis = c(Primary = "#beb7ff", Metastasis = "#b20037",
                            Unavailable = "#e1e1dd"),
  Molecular_record = c(Reported = "#49a12e", Unavailable = "#e1e1dd")
)

palette_for_annotation <- function(values, name, panel) {
  levels_present <- sort(unique(clean_annotation(values)))
  predefined <- fixed_annotation_colours[[name]]
  answer <- if (is.null(predefined)) character() else predefined[names(predefined) %in% levels_present]
  missing_levels <- setdiff(levels_present, names(answer))
  if (length(missing_levels)) {
    if (length(missing_levels) > length(annotation_pool)) {
      record_gap(
        panel,
        paste0("annotation palette: ", name),
        sprintf("%d unrecognised levels exceed the fixed muted palette; overflow is grey.",
                length(missing_levels)),
        "PARTIAL"
      )
    }
    assigned <- rep("#b0b0b0", length(missing_levels))
    take <- seq_len(min(length(missing_levels), length(annotation_pool)))
    assigned[take] <- annotation_pool[take]
    names(assigned) <- missing_levels
    answer <- c(answer, assigned)
  }
  answer
}

build_annotation <- function(cohort, sample_ids, panel) {
  key <- switch(
    cohort,
    TCGA = "f1_tcga_annotation",
    PDX = "f1_pdx_annotation",
    CellLine = "f1_cell_annotation"
  )
  detail_path <- find_input(key, required = FALSE)
  if (!is.na(detail_path)) {
    x <- fread(detail_path)
    id_candidates <- c("sample_id", "Sample", "depmap_id")
    id_col <- id_candidates[id_candidates %in% names(x)][1L]
    if (is.na(id_col)) stop(sprintf("No sample identifier in %s", detail_path))
    matched <- match(sample_ids, x[[id_col]])
    if (anyNA(matched)) {
      record_gap(
        panel,
        paste0(cohort, " annotation/sample matching"),
        sprintf("%d/%d activity samples have no detailed annotation row.",
                sum(is.na(matched)), length(matched)),
        "PARTIAL"
      )
    }
    x <- x[matched]
  } else {
    manifest_key <- switch(
      cohort,
      TCGA = "f1_tcga_manifest",
      PDX = "f1_pdx_manifest",
      CellLine = "f1_cell_manifest"
    )
    manifest_path <- find_input(manifest_key)
    x <- fread(manifest_path)
    id_candidates <- c("sample_id", "Sample", "depmap_id")
    id_col <- id_candidates[id_candidates %in% names(x)][1L]
    if (is.na(id_col)) stop(sprintf("No sample identifier in %s", manifest_path))
    matched <- match(sample_ids, x[[id_col]])
    if (anyNA(matched)) {
      record_gap(
        panel,
        paste0(cohort, " manifest/sample matching"),
        sprintf("%d/%d activity samples have no manifest row.",
                sum(is.na(matched)), length(matched)),
        "PARTIAL"
      )
    }
    x <- x[matched]
    record_gap(
      panel,
      paste0(cohort, " detailed annotations"),
      sprintf("Detailed annotation file is absent; rendered only fields present in %s.", manifest_path),
      "PARTIAL"
    )
  }

  if (cohort == "TCGA") {
    answer <- data.frame(
      Histology = pick_column(x, c("group", "Histology")),
      Published_subtype = pick_column(x, c("published_expression_subtype", "Published_subtype")),
      Sex = pick_column(x, c("sex", "Sex")),
      Age = pick_column(x, c("age_group", "Age")),
      Vital_status = pick_column(x, c("vital_status", "Vital_status")),
      Stage = pick_column(x, c("stage", "Stage")),
      Pack_years = pick_column(x, c("pack_year_group", "Pack_years")),
      Oncogene_mutation = pick_column(x, c("oncogene_mutation", "Oncogene_mutation")),
      Suppressor_mutation = pick_column(x, c("suppressor_pathway_mutation", "Suppressor_mutation")),
      check.names = FALSE
    )
  } else if (cohort == "PDX") {
    answer <- data.frame(
      Histology = pick_column(x, c("group", "Histology")),
      Sex = pick_column(x, c("sex", "Sex")),
      Age = pick_column(x, c("age_group", "Age")),
      Metastatic_disease = pick_column(x, c("known_metastatic_disease", "Metastatic_disease")),
      Smoking = pick_column(x, c("has_smoked_100_cigarettes", "Smoking")),
      Biopsy_site = pick_column(x, c("biopsy_site", "sample_collection_site", "Biopsy_site")),
      Tissue_type = pick_column(x, c("tissue_type", "pdm_type", "Tissue_type")),
      Molecular_record = pick_column(x, c("molecular_record", "Molecular_record")),
      check.names = FALSE
    )
  } else {
    answer <- data.frame(
      Histology = pick_column(x, c("group", "Histology")),
      Sex = pick_column(x, c("sex", "Sex")),
      Age = if ("age_group" %in% names(x)) clean_annotation(x$age_group) else age_group(x$age),
      Primary_or_metastasis = pick_column(x, c("primary_or_metastasis", "Primary_or_metastasis")),
      Collection_site = pick_column(x, c("sample_collection_site", "Collection_site")),
      Growth = pick_column(x, c("growth_pattern", "default_growth_pattern", "Growth")),
      Oncogene_mutation = pick_column(x, c("oncogene_mutation", "Oncogene_mutation")),
      Suppressor_mutation = pick_column(x, c("suppressor_pathway_mutation", "Suppressor_mutation")),
      check.names = FALSE
    )
  }

  answer[] <- lapply(answer, clean_annotation)
  keep <- vapply(answer, function(z) any(z != "Unavailable"), logical(1))
  if (!any(keep)) {
    stop(sprintf("No usable %s annotations matched the activity matrix.", cohort))
  }
  answer <- answer[, keep, drop = FALSE]
  rownames(answer) <- sample_ids
  colours <- setNames(
    lapply(names(answer), function(name) palette_for_annotation(answer[[name]], name, panel)),
    names(answer)
  )
  list(df = answer, colours = colours)
}

tf_membership <- NULL
get_tf_membership <- function() {
  if (is.null(tf_membership)) {
    tf_membership <<- fread(find_input("f1_tf_membership"))
    anchor_assert_columns(tf_membership, c("TF", "discovery_category"), "tf_membership")
  }
  tf_membership
}

luad_tfs <- function() {
  answer <- unique(get_tf_membership()[discovery_category == "LUAD", as.character(TF)])
  if (length(answer) != 158L) {
    stop(sprintf("Expected 158 frozen LUAD-specific TFs; found %d.", length(answer)))
  }
  answer
}

# -------------------------------------------------------------------------
# Figure 1 and Supplementary Figure 2.
# -------------------------------------------------------------------------

discovery_colour_map <- c(
  LUAD = unname(anchor_colours$discovery["LUAD"]),
  LUSC = unname(anchor_colours$discovery["LUSC"]),
  `Non-Significant` = unname(anchor_colours$discovery["nonsignificant"]),
  `LogFC Non-Significant` = unname(anchor_colours$discovery["expression_only"]),
  `LogFC Significant` = unname(anchor_colours$discovery["expression_only"]),
  `VIPER Significant` = unname(anchor_colours$discovery["activity_only"])
)

discovery_scatter <- function(x, add_counts = FALSE) {
  anchor_assert_columns(x, c("TF", "NES", "logFC", "category"), "discovery statistics")
  x <- copy(x)
  x[, label := fifelse(category %in% c("LUAD", "LUSC"), TF, NA_character_)]
  pearson <- cor.test(x$NES, x$logFC, method = "pearson")
  x_pos <- unname(quantile(x$NES, 0.70, na.rm = TRUE))
  y_pos <- unname(quantile(x$logFC, 0.06, na.rm = TRUE))

  p <- ggplot(x, aes(NES, logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "lightgrey", linewidth = 0.4) +
    geom_point(aes(colour = category), size = 0.75) +
    ggrepel::geom_text_repel(
      aes(label = label), family = anchor_font_family, size = 3,
      min.segment.length = 0.1, max.overlaps = 40, seed = 2026
    ) +
    scale_colour_manual(values = discovery_colour_map, drop = FALSE) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
    annotate(
      "text", x = x_pos, y = y_pos, hjust = 0,
      label = sprintf("Pearson r = %.2f\np = %.2g", unname(pearson$estimate), pearson$p.value),
      family = anchor_font_family, size = 3.2
    ) +
    labs(x = "NES", y = expression(log[2] * "FC gene expression")) +
    anchor_theme_classic(12) +
    theme(legend.position = "none")

  if (add_counts) {
    counts <- x[category %in% c("LUAD", "LUSC"), .(
      NES = min(NES, na.rm = TRUE),
      logFC = max(logFC, na.rm = TRUE),
      n_TFs = .N
    ), by = category]
    p <- p + geom_text(
      data = counts,
      aes(NES, logFC, label = n_TFs, colour = category),
      inherit.aes = FALSE,
      family = anchor_font_family,
      fontface = "bold",
      size = 3.5
    )
  }
  p
}

run_panel("Figure1B", {
  require_panel_package("ggrepel", "Figure1B")
  p <- discovery_scatter(fread(find_input("f1_tcga_stats")), add_counts = FALSE)
  save_anchor_pdf(p, "Figure1B_TCGA_logFC_VIPER_anchor", 3.5, 4)
})

plot_activity_heatmap <- function(matrix_key, cohort, panel, cluster_rows, width) {
  require_panel_package("ComplexHeatmap", panel)
  matrix <- read_activity_matrix(find_input(matrix_key))
  wanted <- get_tf_membership()[discovery_category %in% c("LUAD", "LUSC"), TF]
  missing_tfs <- setdiff(wanted, rownames(matrix))
  if (length(missing_tfs)) {
    record_gap(
      panel, "activity rows",
      sprintf("%d discovery TFs are absent and were not plotted: %s",
              length(missing_tfs), paste(head(missing_tfs, 12L), collapse = ", ")),
      "PARTIAL"
    )
  }
  matrix <- matrix[intersect(wanted, rownames(matrix)), , drop = FALSE]
  matrix <- row_zscore(matrix)
  annotation <- build_annotation(cohort, colnames(matrix), panel)

  top <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation$df,
    col = annotation$colours,
    annotation_name_side = "left",
    which = "column",
    simple_anno_size = unit(1, "mm"),
    annotation_name_gp = gpar(fontfamily = anchor_font_family, fontsize = 7),
    annotation_legend_param = setNames(
      lapply(names(annotation$df), function(x) {
        list(
          title_gp = gpar(fontfamily = anchor_font_family, fontsize = 8, fontface = "bold"),
          labels_gp = gpar(fontfamily = anchor_font_family, fontsize = 7)
        )
      }),
      names(annotation$df)
    )
  )

  row_group <- get_tf_membership()$discovery_category[
    match(rownames(matrix), get_tf_membership()$TF)
  ]
  left <- ComplexHeatmap::rowAnnotation(
    Specificity = factor(row_group, levels = c("LUAD", "LUSC")),
    col = list(Specificity = anchor_colours$phenotype),
    simple_anno_size = unit(2, "mm"),
    annotation_name_gp = gpar(fontfamily = anchor_font_family, fontsize = 8),
    annotation_legend_param = list(
      Specificity = list(
        title_gp = gpar(fontfamily = anchor_font_family, fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontfamily = anchor_font_family, fontsize = 7)
      )
    )
  )

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "TF activity",
    cluster_rows = cluster_rows,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    use_raster = ncol(matrix) > 200L,
    raster_quality = 2,
    left_annotation = left,
    top_annotation = top,
    heatmap_legend_param = list(
      title = "TF activity",
      title_gp = gpar(fontfamily = anchor_font_family, fontsize = 9),
      labels_gp = gpar(fontfamily = anchor_font_family, fontsize = 8)
    ),
    column_title_gp = gpar(fontfamily = anchor_font_family, fontsize = 10),
    column_names_gp = gpar(fontfamily = anchor_font_family, fontsize = 8)
  )
  save_anchor_heatmap(heatmap, paste0(panel, "_", cohort, "_activity_anchor"), width, 6.5)
}

run_panel("Figure1C", plot_activity_heatmap(
  "f1_tcga_activity", "TCGA", "Figure1C", FALSE, 7.5
))
run_panel("Figure1D", plot_activity_heatmap(
  "f1_pdx_activity", "PDX", "Figure1D", TRUE, 7
))
run_panel("Figure1E", plot_activity_heatmap(
  "f1_cell_activity", "CellLine", "Figure1E", TRUE, 7
))

run_panel("SupplementaryFigure2A", {
  require_panel_package("ggrepel", "SupplementaryFigure2A")
  p <- discovery_scatter(fread(find_input("f1_gse_stats")), add_counts = TRUE)
  save_anchor_pdf(p, "SupplementaryFigure2A_GSE81089_logFC_VIPER_anchor", 3.5, 4)
})

run_panel("SupplementaryFigure2B", {
  require_panel_package("eulerr", "SupplementaryFigure2B")
  tcga <- fread(find_input("f1_tcga_stats"))
  gse <- fread(find_input("f1_gse_stats"))
  paths <- character()
  for (direction in c("LUAD", "LUSC")) {
    discovery <- unique(tcga[category == direction, TF])
    validation <- unique(gse[category == direction, TF])
    fit <- eulerr::euler(list(TCGA = discovery, GSE81089 = validation))
    path <- file.path(
      atomic_dir,
      sprintf("SupplementaryFigure2B_%s_TF_overlap_anchor.pdf", direction)
    )
    anchor_cairo_pdf(path, width = 3.5, height = 3)
    tryCatch(
      plot(
        fit,
        fills = list(fill = rep(unname(anchor_colours$phenotype[direction]), 2L), alpha = 0.7),
        labels = list(font = 2, family = anchor_font_family),
        quantities = TRUE,
        edges = list(col = "white")
      ),
      finally = dev.off()
    )
    paths <- c(paths, path)
  }
  paths
})

plot_correlation_heatmap <- function(matrix_key, cohort, panel, stem) {
  require_panel_package("ComplexHeatmap", panel)
  require_panel_package("circlize", panel)
  matrix <- read_activity_matrix(find_input(matrix_key))
  wanted <- get_tf_membership()[discovery_category %in% c("LUAD", "LUSC"), TF]
  matrix <- matrix[intersect(wanted, rownames(matrix)), , drop = FALSE]
  matrix <- row_zscore(matrix)
  corr <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
  annotation <- build_annotation(cohort, colnames(corr), panel)
  top <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation$df,
    col = annotation$colours,
    annotation_name_side = "left",
    which = "column",
    simple_anno_size = unit(1, "mm"),
    annotation_name_gp = gpar(fontfamily = anchor_font_family, fontsize = 7),
    annotation_legend_param = setNames(
      lapply(names(annotation$df), function(x) {
        list(
          title_gp = gpar(fontfamily = anchor_font_family, fontsize = 8, fontface = "bold"),
          labels_gp = gpar(fontfamily = anchor_font_family, fontsize = 7)
        )
      }),
      names(annotation$df)
    )
  )
  col_fun <- circlize::colorRamp2(
    c(0.5, 0.75, 1),
    hcl_palette = "Purple-Green",
    reverse = TRUE
  )
  heatmap <- ComplexHeatmap::Heatmap(
    corr,
    name = "Pearson correlation",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    use_raster = ncol(corr) > 200L,
    raster_quality = 2,
    top_annotation = top,
    heatmap_legend_param = list(
      title = "Pearson correlation",
      title_gp = gpar(fontfamily = anchor_font_family, fontsize = 9),
      labels_gp = gpar(fontfamily = anchor_font_family, fontsize = 8)
    )
  )
  save_anchor_heatmap(heatmap, stem, 7.5, 6)
}

run_panel("SupplementaryFigure2C", plot_correlation_heatmap(
  "f1_tcga_activity", "TCGA", "SupplementaryFigure2C",
  "SupplementaryFigure2C_TCGA_sample_correlations_anchor"
))
run_panel("SupplementaryFigure2D", plot_correlation_heatmap(
  "f1_pdx_activity", "PDX", "SupplementaryFigure2D",
  "SupplementaryFigure2D_PDX_sample_correlations_anchor"
))
run_panel("SupplementaryFigure2E", plot_correlation_heatmap(
  "f1_cell_activity", "CellLine", "SupplementaryFigure2E",
  "SupplementaryFigure2E_CellLine_sample_correlations_anchor"
))

activity_distribution_plot <- function(matrix, selected_tfs, selected_samples, tf_order, y_size) {
  present_tfs <- intersect(tf_order, intersect(selected_tfs, rownames(matrix)))
  present_samples <- intersect(selected_samples, colnames(matrix))
  if (!length(present_tfs) || !length(present_samples)) {
    stop("No TF-by-sample activity block remains after LUAD sample selection.")
  }
  sub <- matrix[present_tfs, present_samples, drop = FALSE]
  long <- as.data.table(as.table(sub))
  setnames(long, c("TF", "Sample", "NES"))
  long[, TF := factor(as.character(TF), levels = tf_order)]
  long[, sign := fifelse(NES >= 0, "Positive", "Negative")]
  ggplot(long[!is.na(TF)], aes(NES, TF, colour = sign)) +
    geom_point(size = 0.10) +
    scale_colour_manual(values = c(Positive = "darkred", Negative = "navy")) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
    labs(x = "NES", y = NULL) +
    anchor_theme_classic(12) +
    theme(
      axis.text.y = element_text(size = y_size, family = anchor_font_family),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}

run_panel("SupplementaryFigure2F", {
  tfs <- luad_tfs()
  tcga_matrix <- read_activity_matrix(find_input("f1_tcga_activity"))
  tcga_manifest <- fread(find_input("f1_tcga_manifest"))
  anchor_assert_columns(tcga_manifest, c("sample_id", "group"), "TCGA manifest")
  tcga_samples <- tcga_manifest[group == "LUAD", sample_id]
  tcga_sub <- tcga_matrix[
    intersect(tfs, rownames(tcga_matrix)),
    intersect(tcga_samples, colnames(tcga_matrix)),
    drop = FALSE
  ]
  medians <- apply(tcga_sub, 1L, median, na.rm = TRUE)
  tf_order <- names(sort(medians, decreasing = FALSE))

  specs <- list(
    TCGA = list(matrix_key = "f1_tcga_activity", manifest_key = "f1_tcga_manifest", y_size = 4.5),
    PDX = list(matrix_key = "f1_pdx_activity", manifest_key = "f1_pdx_manifest", y_size = 4.5),
    CellLine = list(matrix_key = "f1_cell_activity", manifest_key = "f1_cell_manifest", y_size = 4.5),
    GSE81089 = list(matrix_key = "f1_gse_activity", manifest_key = "f1_gse_manifest", y_size = 3.5)
  )
  paths <- character()
  for (cohort in names(specs)) {
    matrix_path <- find_input(specs[[cohort]]$matrix_key, required = FALSE)
    if (is.na(matrix_path)) {
      record_gap(
        "SupplementaryFigure2F",
        paste0(cohort, " sample-level VIPER activity"),
        sprintf("No matrix found for %s; this cohort panel was not synthesized from contrast NES.", cohort),
        "PARTIAL"
      )
      next
    }
    matrix <- read_activity_matrix(matrix_path)
    manifest <- fread(find_input(specs[[cohort]]$manifest_key))
    anchor_assert_columns(manifest, c("sample_id", "group"), paste0(cohort, " manifest"))
    samples <- manifest[group == "LUAD", sample_id]
    p <- activity_distribution_plot(matrix, tfs, samples, tf_order, specs[[cohort]]$y_size)
    paths <- c(
      paths,
      save_anchor_pdf(
        p,
        sprintf("SupplementaryFigure2F_%s_LUAD_TF_activity_anchor", cohort),
        3.75, 8
      )
    )
  }
  if (!length(paths)) stop("No Supplementary Figure 2F cohort could be rendered.")
  paths
})

# -------------------------------------------------------------------------
# Figure 2B-E: primary, original-denominator rule only.
# -------------------------------------------------------------------------

primary_definition <- "PRIMARY_anchor_promoter_tcga_cpm1_both"
system_levels_raw <- c("patient", "PDX", "cell_line")
system_labels <- c(patient = "Patient", PDX = "PDX", cell_line = "CellLine")
system_colours <- anchor_colours$atac_system[c("Patient", "PDX", "CellLine")]

read_primary_table <- function(key, required_columns) {
  x <- fread(find_input(key))
  anchor_assert_columns(x, required_columns, key)
  if ("analysis_definition" %in% names(x)) {
    x <- x[analysis_definition == primary_definition]
    if (!nrow(x)) stop(sprintf("%s has no rows for %s.", key, primary_definition))
  }
  x
}

run_panel("Figure2B", {
  require_panel_package("ComplexHeatmap", "Figure2B")
  require_panel_package("ComplexUpset", "Figure2B")
  gate <- read_primary_table(
    "f2_promoter_gate",
    c("TF", "patient_promoter_accessible", "PDX_promoter_accessible", "cell_line_promoter_accessible")
  )
  sets <- list(
    Patient = gate[as_flag(patient_promoter_accessible), TF],
    PDX = gate[as_flag(PDX_promoter_accessible), TF],
    CellLine = gate[as_flag(cell_line_promoter_accessible), TF]
  )
  membership <- ComplexHeatmap::list_to_matrix(sets)
  membership <- membership[, c("Patient", "PDX", "CellLine"), drop = FALSE]
  p <- ComplexUpset::upset(
    data = as.data.frame(membership),
    intersect = c("Patient", "PDX", "CellLine"),
    set_sizes = FALSE,
    queries = list(
      ComplexUpset::upset_query("Patient", color = system_colours[["Patient"]]),
      ComplexUpset::upset_query("PDX", color = system_colours[["PDX"]]),
      ComplexUpset::upset_query("CellLine", color = system_colours[["CellLine"]])
    ),
    name = "TFs with promoter accessibility\nin at least 50% of samples",
    themes = theme(text = element_text(
      size = 10, colour = "black", family = anchor_font_family
    )),
    base_annotations = list(
      "Number of TFs" = ComplexUpset::intersection_size() +
        labs(y = "Number of TFs") +
        anchor_theme_classic(10) +
        theme(
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    )
  ) + theme(text = element_text(colour = "black", family = anchor_font_family))
  save_anchor_pdf(p, "Figure2B_promoter_accessibility_UpSet_anchor", 7, 7)
})

run_panel("Figure2C", {
  require_panel_package("patchwork", "Figure2C")
  tfs <- luad_tfs()
  promoter <- read_primary_table(
    "f2_promoter_sample",
    c("system", "sample_id", "TF", "promoter_accessible")
  )[TF %in% tfs & system %in% system_levels_raw]
  activity <- read_primary_table(
    "f2_promoter_system",
    c("system", "TF", "mean_LUAD_NES")
  )[TF %in% tfs & system %in% system_levels_raw]
  gate <- read_primary_table(
    "f2_promoter_gate",
    c("TF", "HC_TF_promoter_activity_definition")
  )[TF %in% tfs]
  motif <- fread(find_input("f2_motif_sample"))
  anchor_assert_columns(
    motif,
    c("system", "sample_id", "TF", "q_value", "motif_enriched_q1e-5"),
    "motif sample table"
  )
  motif <- motif[TF %in% tfs & system %in% system_levels_raw]
  inventory <- fread(find_input("f2_motif_inventory"))
  anchor_assert_columns(inventory, c("TF", "motif_testable"), "motif inventory")

  expected_promoter_rows <- uniqueN(promoter[, .(system, sample_id)]) * length(tfs)
  if (nrow(unique(promoter[, .(system, sample_id, TF)])) != expected_promoter_rows) {
    stop(sprintf(
      "Figure2C requires a complete 158 TF x sample promoter grid; observed %d unique rows, expected %d.",
      nrow(unique(promoter[, .(system, sample_id, TF)])), expected_promoter_rows
    ))
  }
  if (uniqueN(gate$TF) != 158L || uniqueN(activity$TF) != 158L) {
    stop(sprintf(
      "Figure2C frozen TF coverage is incomplete: gate=%d, activity=%d, expected=158.",
      uniqueN(gate$TF), uniqueN(activity$TF)
    ))
  }
  if (nrow(unique(activity[, .(TF, system)])) != 158L * length(system_levels_raw)) {
    stop("Figure2C requires one mean NES value for every TF in every primary system.")
  }
  if (sum(as_flag(gate$HC_TF_promoter_activity_definition)) != 31L) {
    stop(sprintf(
      "Expected 31 frozen HC-TFs in Figure2C; found %d.",
      sum(as_flag(gate$HC_TF_promoter_activity_definition))
    ))
  }

  activity[, System := factor(system_labels[system], levels = unname(system_labels))]
  activity[, NES_bin := fifelse(
    mean_LUAD_NES < 0, "<0",
    fifelse(mean_LUAD_NES < 1, "0-1", fifelse(mean_LUAD_NES < 3, "1-3", ">3"))
  )]
  activity[, NES_bin := factor(NES_bin, levels = c(">3", "1-3", "0-1", "<0"))]

  testable_tfs <- inventory[as_flag(motif_testable), TF]
  average_nes <- activity[, .(overall_mean_NES = mean(mean_LUAD_NES, na.rm = TRUE)), by = TF]
  no_motif_order <- average_nes[!TF %in% testable_tfs][order(-overall_mean_NES), TF]
  with_motif_order <- average_nes[TF %in% testable_tfs][order(-overall_mean_NES), TF]
  tf_order <- c(rev(no_motif_order), rev(with_motif_order))
  if (length(tf_order) != 158L || uniqueN(tf_order) != 158L) {
    stop("Figure2C TF ordering did not retain exactly 158 unique TFs.")
  }

  gate[, TF := factor(TF, levels = tf_order)]
  gate[, HC_group := fifelse(
    as_flag(HC_TF_promoter_activity_definition), "HC-TF", "Not HC-TF"
  )]
  activity[, TF := factor(TF, levels = tf_order)]

  sample_order <- promoter[, unique(sample_id), by = system]
  sample_order[, system_rank := match(system, system_levels_raw)]
  setorder(sample_order, system_rank)
  sample_levels <- sample_order$V1
  promoter[, System := system_labels[system]]
  promoter[, Sample := factor(sample_id, levels = sample_levels)]
  promoter[, TF := factor(TF, levels = tf_order)]
  promoter[, accessibility_category := fifelse(
    as_flag(promoter_accessible), System, "NA"
  )]

  motif_flag_column <- if ("motif_enriched_q1e-5" %in% names(motif)) {
    "motif_enriched_q1e-5"
  } else if ("motif_enriched_q1e.5" %in% names(motif)) {
    "motif_enriched_q1e.5"
  } else {
    stop("The motif sample table lacks the frozen q<=1e-5 enrichment flag.")
  }
  motif_small <- motif[, .(
    motif_enriched = any(as_flag(get(motif_flag_column))),
    q_value = {
      z <- suppressWarnings(as.numeric(q_value))
      if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)
    }
  ), by = .(system, sample_id, TF)]

  motif_grid <- unique(promoter[, .(system, sample_id, TF = as.character(TF), System)])
  motif_grid <- merge(
    motif_grid,
    motif_small,
    by = c("system", "sample_id", "TF"),
    all.x = TRUE,
    sort = FALSE
  )
  motif_grid[, testable := TF %in% testable_tfs]
  missing_tested <- motif_grid[testable & is.na(motif_enriched), unique(TF)]
  if (length(missing_tested)) {
    stop(sprintf(
      "Figure2C lacks sample-level motif rows for %d testable TFs: %s",
      length(missing_tested), paste(missing_tested, collapse = ", ")
    ))
  }
  motif_grid[, motif_category := fifelse(
    !testable, "No Motif",
    fifelse(motif_enriched, System, "NA")
  )]
  motif_grid[, Sample := factor(sample_id, levels = sample_levels)]
  motif_grid[, TF := factor(TF, levels = tf_order)]

  cohort_fill <- c(system_colours, `NA` = "white", `No Motif` = "grey")
  p_hc <- ggplot(gate, aes("TF group", TF, fill = HC_group)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = c(`HC-TF` = "#0CABA8", `Not HC-TF` = "grey")) +
    labs(x = NULL, y = NULL, fill = "TF group") +
    theme_minimal(base_size = 12, base_family = anchor_font_family) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 8),
      axis.text.y = element_text(size = 3), panel.grid = element_blank()
    )
  p_nes <- ggplot(activity, aes(System, TF, fill = NES_bin)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = anchor_colours$nes_bin, drop = FALSE) +
    labs(x = NULL, y = NULL, fill = "NES") +
    theme_minimal(base_size = 12, base_family = anchor_font_family) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 8),
      axis.text.y = element_blank(), panel.grid = element_blank()
    )
  p_promoter <- ggplot(promoter, aes(Sample, TF, fill = accessibility_category)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = cohort_fill, guide = "none", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12, base_family = anchor_font_family) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 5),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank()
    )
  p_motif <- ggplot(motif_grid, aes(Sample, TF, fill = motif_category)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = cohort_fill, drop = FALSE) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 12, base_family = anchor_font_family) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 5),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank()
    )

  combined <- patchwork::wrap_plots(
    p_hc, p_nes, p_promoter, p_motif,
    nrow = 1,
    widths = c(0.05, 0.1, 1, 1),
    guides = "collect"
  ) & theme(legend.position = "right")
  save_anchor_pdf(
    combined,
    "Figure2C_158TF_sample_level_promoter_motif_anchor",
    12, 8
  )
})

run_panel("Figure2D", {
  require_panel_package("ComplexHeatmap", "Figure2D")
  require_panel_package("ComplexUpset", "Figure2D")
  system_summary <- fread(find_input("f2_motif_system"))
  anchor_assert_columns(
    system_summary,
    c("system", "TF", "support_fixed_original_denominator"),
    "motif system summary"
  )
  inventory <- fread(find_input("f2_motif_inventory"))
  anchor_assert_columns(
    inventory,
    c("TF", "motif_database_category", "motif_testable"),
    "motif inventory"
  )
  tfs <- luad_tfs()
  system_summary <- system_summary[TF %in% tfs & system %in% system_levels_raw]
  sets <- list(
    Patient = system_summary[system == "patient" & as_flag(support_fixed_original_denominator), TF],
    PDX = system_summary[system == "PDX" & as_flag(support_fixed_original_denominator), TF],
    CellLine = system_summary[system == "cell_line" & as_flag(support_fixed_original_denominator), TF]
  )
  membership <- ComplexHeatmap::list_to_matrix(sets)
  membership <- as.data.frame(membership[, c("Patient", "PDX", "CellLine"), drop = FALSE])
  membership$TF <- rownames(membership)
  db <- inventory$motif_database_category[match(membership$TF, inventory$TF)]
  membership$db_group <- factor(
    fifelse(
      db == "both", "Both",
      fifelse(db == "CIS-BP_only", "CISBP", fifelse(db == "JASPAR_only", "JASPAR", "Not HC-TF"))
    ),
    levels = c("Not HC-TF", "Both", "CISBP", "JASPAR")
  )
  motif_fill <- c(`Not HC-TF` = "grey", anchor_colours$motif_database)
  y_label <- sprintf("Number of TFs\n(n = %d/%d)", nrow(membership), sum(as_flag(inventory$motif_testable)))
  annotation_list <- setNames(
    list(
      ggplot(mapping = aes(x = intersection)) +
        geom_bar(aes(fill = db_group), stat = "count", position = "stack") +
        geom_text(
          stat = "count",
          aes(label = after_stat(count), group = db_group),
          position = position_stack(vjust = 0.5),
          colour = "black", size = 2, family = anchor_font_family
        ) +
        scale_fill_manual(values = motif_fill, drop = FALSE, name = "Motif database") +
        anchor_theme_classic(10) +
        theme(
          axis.title.y = element_text(size = 9), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 15, unit = "mm")
        )
    ),
    y_label
  )
  p <- ComplexUpset::upset(
    data = membership,
    intersect = c("Patient", "PDX", "CellLine"),
    set_sizes = FALSE,
    queries = list(
      ComplexUpset::upset_query("Patient", color = system_colours[["Patient"]]),
      ComplexUpset::upset_query("PDX", color = system_colours[["PDX"]]),
      ComplexUpset::upset_query("CellLine", color = system_colours[["CellLine"]])
    ),
    name = "TFs with genome-wide motif\nenrichment in at least 50% of samples",
    themes = theme(text = element_text(
      size = 10, colour = "black", family = anchor_font_family
    )),
    base_annotations = annotation_list
  )
  save_anchor_pdf(p, "Figure2D_motif_enrichment_UpSet_anchor", 7, 7)
})

prepare_lor_data <- function() {
  motif <- fread(find_input("f2_motif_sample"))
  anchor_assert_columns(
    motif,
    c("system", "sample_id", "TF", "q_value", "log2_odds_ratio"),
    "motif sample table"
  )
  tfs <- luad_tfs()
  motif <- motif[TF %in% tfs & system %in% system_levels_raw]
  motif[, LOR := suppressWarnings(as.numeric(log2_odds_ratio))]
  motif[, q_numeric := suppressWarnings(as.numeric(q_value))]
  motif[, System := factor(system_labels[system], levels = unname(system_labels))]
  motif[is.finite(LOR)]
}

run_panel("Figure2E_summary", {
  lor <- prepare_lor_data()
  support <- fread(find_input("f2_motif_system"))
  anchor_assert_columns(
    support,
    c("system", "TF", "support_fixed_original_denominator"),
    "motif system summary"
  )
  support[, System := system_labels[system]]
  support[, supported := as_flag(support_fixed_original_denominator)]
  summary <- lor[, .(
    Mean_LOR = mean(LOR, na.rm = TRUE),
    SD_LOR = sd(LOR, na.rm = TRUE),
    n_samples = uniqueN(sample_id)
  ), by = .(TF, System)]
  summary <- merge(
    summary,
    support[, .(TF, System, supported)],
    by = c("TF", "System"),
    all.x = TRUE
  )
  summary[is.na(supported), supported := FALSE]
  order_tfs <- summary[, .(median_LOR = median(Mean_LOR, na.rm = TRUE)), by = TF][
    order(median_LOR), TF
  ]
  summary[, TF_ordered := factor(TF, levels = order_tfs)]
  summary[, supported_factor := factor(supported, levels = c(FALSE, TRUE))]

  p <- ggplot(summary, aes(Mean_LOR, TF_ordered)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_point(
      aes(
        colour = System, size = SD_LOR, shape = supported_factor,
        alpha = supported_factor, stroke = supported_factor
      )
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 19),
      name = "Enriched in\n>=50% of samples",
      labels = c(`FALSE` = "No", `TRUE` = "Yes")
    ) +
    scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.9), guide = "none") +
    scale_discrete_manual(
      aesthetics = "stroke", values = c(`FALSE` = 1.2, `TRUE` = 0.3), guide = "none"
    ) +
    scale_colour_manual(values = system_colours, name = "Cohort") +
    scale_size_continuous(name = "SD of LOR", range = c(6, 1)) +
    labs(
      x = "Mean log2 odds ratio\n(target vs. background)",
      y = "Transcription factor"
    ) +
    anchor_theme_bw(12) +
    theme(
      axis.text = element_text(colour = "black", size = 8),
      axis.title = element_text(size = 10)
    )
  fwrite(summary, file.path(table_dir, "Figure2E_mean_LOR_summary.tsv"), sep = "\t")
  save_anchor_pdf(p, "Figure2E_motif_LOR_summary_anchor", 5, 8)
})

run_panel("Figure2E_cases", {
  require_panel_package("ggbeeswarm", "Figure2E_cases")
  lor <- prepare_lor_data()
  triple <- fread(find_input("f2_motif_triple"))
  anchor_assert_columns(
    triple,
    c("TF", "triple_system_fixed_original_denominator"),
    "triple-system motif summary"
  )
  ranked <- lor[, .(overall_median_LOR = median(LOR, na.rm = TRUE)), by = TF][
    order(-overall_median_LOR)
  ]
  if (toupper(case_tf_arg) == "AUTO") {
    primary <- triple[as_flag(triple_system_fixed_original_denominator), TF]
    candidates <- unique(c(primary, ranked$TF))
    chosen <- head(candidates[candidates %in% ranked$TF], 6L)
    selection_rule <- "AUTO: primary triple-system hits, then descending overall median LOR"
  } else {
    chosen <- trimws(strsplit(case_tf_arg, ",", fixed = TRUE)[[1L]])
    missing <- setdiff(chosen, unique(lor$TF))
    if (length(missing)) {
      stop(sprintf("Predeclared Figure2E case TFs lack sample LOR data: %s", paste(missing, collapse = ", ")))
    }
    selection_rule <- "CLI-predeclared"
  }
  if (!length(chosen)) stop("No TF qualified for the Figure2E case-example panel.")
  cases <- lor[TF %in% chosen]
  cases[, TF := factor(TF, levels = chosen)]
  cases[, significant := !is.na(q_numeric) & q_numeric < 0.05]
  selection <- data.table(
    rank = seq_along(chosen), TF = chosen, selection_rule = selection_rule
  )
  fwrite(selection, file.path(table_dir, "Figure2E_case_examples_used.tsv"), sep = "\t")

  p <- ggplot(cases, aes(System, LOR)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_boxplot(aes(colour = System), width = 0.65, outlier.shape = NA, linewidth = 0.6) +
    ggbeeswarm::geom_quasirandom(
      data = cases[significant], aes(colour = System),
      width = 0.12, varwidth = FALSE, size = 1.8, alpha = 0.9
    ) +
    ggbeeswarm::geom_quasirandom(
      data = cases[!significant], aes(colour = "Non significant"),
      width = 0.12, varwidth = FALSE, size = 1.8, alpha = 0.9
    ) +
    scale_colour_manual(
      values = c(system_colours, `Non significant` = "#F28E2B"),
      breaks = c(names(system_colours), "Non significant")
    ) +
    labs(x = NULL, y = "Log2 odds ratio (LOR)\n(target vs. background)", colour = NULL) +
    facet_grid(TF ~ ., scales = "free_y", space = "free_y") +
    coord_flip() +
    anchor_theme_bw(12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(colour = "black"),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(face = "bold", size = 12, angle = 0),
      strip.background = element_blank(),
      panel.spacing.y = unit(1.1, "lines")
    )
  save_anchor_pdf(p, "Figure2E_motif_LOR_case_examples_anchor", 5, 7.5)
})

# -------------------------------------------------------------------------
# Receipts.  A non-zero exit is intentional when a required panel failed;
# successfully rendered atomics remain available and the gap table explains
# exactly why the remaining panels were not generated.
# -------------------------------------------------------------------------

status <- if (length(status_rows)) rbindlist(status_rows, fill = TRUE) else data.table(
  panel = character(), status = character(), artifacts = character(),
  seconds = numeric(), message = character()
)
gaps <- if (length(gap_rows)) rbindlist(gap_rows, fill = TRUE) else data.table(
  panel = character(), severity = character(), requirement = character(), reason = character()
)
fwrite(status, file.path(audit_dir, "panel_render_status.tsv"), sep = "\t")
fwrite(gaps, file.path(audit_dir, "input_gaps.tsv"), sep = "\t")

receipt <- data.table(
  key = c(
    "figure1_root", "figure2_root", "publication_root", "output_root",
    "font_family", "primary_promoter_definition", "case_tf_argument",
    "generated_panels", "failed_panels", "partial_gaps"
  ),
  value = c(
    fig1_root, fig2_root, publication_root, output_root,
    anchor_font_family, primary_definition, case_tf_arg,
    sum(status$status == "generated"), sum(status$status == "failed"),
    sum(gaps$severity == "PARTIAL")
  )
)
fwrite(receipt, file.path(audit_dir, "anchor_figure12_render_receipt.tsv"), sep = "\t")

if (any(status$status == "failed")) {
  stop(
    sprintf(
      "%d anchor Figure 1/2 panel group(s) failed. See %s.",
      sum(status$status == "failed"),
      file.path(audit_dir, "input_gaps.tsv")
    ),
    call. = FALSE
  )
}

cat(sprintf(
  "Rendered %d anchor-style panel group(s); %d partial gap(s) recorded in %s.\n",
  sum(status$status == "generated"),
  sum(gaps$severity == "PARTIAL"),
  audit_dir
))

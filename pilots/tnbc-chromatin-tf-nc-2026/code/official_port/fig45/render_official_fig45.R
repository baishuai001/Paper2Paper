#!/usr/bin/env Rscript

# Strict rendering driver. Plot constructors are evaluated from materialized
# CodeOcean source blocks; this file only prepares objects and calls them.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(survival)
  library(survminer)
  library(eulerr)
  library(ComplexUpset)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args_cli <- commandArgs(trailingOnly = TRUE)
if (length(args_cli) != 3L) {
  stop("Usage: render_official_fig45.R GENERATED_DIR ADAPTED_INPUT_DIR OUTPUT_DIR", call. = FALSE)
}
generated <- normalizePath(args_cli[[1L]], mustWork = TRUE)
adapted <- normalizePath(args_cli[[2L]], mustWork = TRUE)
out <- normalizePath(args_cli[[3L]], mustWork = FALSE)
atomic <- file.path(out, "atomic")
audit <- file.path(out, "audit")
dir.create(atomic, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

receipt_path <- file.path(generated, "materialization_receipt.json")
if (!file.exists(receipt_path)) stop("Missing materialization receipt: ", receipt_path)
if (!file.exists(file.path(adapted, "adapter_receipt.json"))) {
  stop("Missing LUAD adapter receipt")
}

block_path <- function(id, patched = TRUE) {
  path <- file.path(generated, if (patched) "patched" else "exact", paste0(id, ".R"))
  if (!file.exists(path)) stop("Missing materialized block: ", path)
  path
}

source_block <- function(id, env = parent.frame(), patched = TRUE) {
  executed_blocks <<- rbind(
    executed_blocks,
    data.table(
      sequence = nrow(executed_blocks) + 1L,
      block_id = id,
      variant = ifelse(patched, "patched", "exact"),
      source_path = block_path(id, patched)
    )
  )
  sys.source(block_path(id, patched), envir = env, keep.source = TRUE)
  invisible(env)
}

read_tab <- function(path, required = character()) {
  if (!file.exists(path)) stop("Missing adapted input: ", path)
  x <- fread(path, na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(path, " missing: ", paste(missing, collapse = ","))
  x
}

status_rows <- list()
executed_blocks <- data.table(
  sequence = integer(), block_id = character(),
  variant = character(), source_path = character()
)
register <- function(panel, status, block_ids, output = "", note = "") {
  status_rows[[length(status_rows) + 1L]] <<- data.table(
    panel = panel,
    status = status,
    official_block_ids = paste(block_ids, collapse = ";"),
    output = output,
    note = note
  )
}

save_plot <- function(panel, plot, width, height, block_ids, note = "") {
  path <- file.path(atomic, paste0(panel, ".pdf"))
  cairo_pdf(path, width = width, height = height, onefile = FALSE)
  print(plot)
  dev.off()
  register(panel, "RENDERED_OFFICIAL_BLOCK", block_ids, path, note)
  invisible(path)
}

# Exact official palette definitions.
source_block("aesthetics", globalenv(), patched = FALSE)

# ---------------------------------------------------------------------------
# Figure 4A-D: official significant-only forest constructor.
# ---------------------------------------------------------------------------

forest_env <- new.env(parent = globalenv())
source_block("fig4_forest", forest_env, patched = FALSE)
cox <- read_tab(
  file.path(adapted, "figure4", "models_long.tsv"),
  c("cohort", "endpoint", "model", "TF", "HR", "lower95", "upper95", "p_value")
)
forest_defs <- data.table(
  panel = c("Figure4A", "Figure4B", "Figure4C", "Figure4D"),
  cohort = c("GSE41271_LUAD", "GSE41271_LUAD", "TCGA_LUAD", "TCGA_LUAD"),
  endpoint = c("OS", "RFS", "OS", "RFS"),
  title = c(
    "Multivariate GSE41271-LUAD OS", "Multivariate GSE41271-LUAD RFS",
    "Multivariate TCGA-LUAD OS", "Multivariate TCGA-LUAD DFS"
  )
)
for (i in seq_len(nrow(forest_defs))) {
  d <- forest_defs[i]
  z <- cox[cohort == d$cohort & endpoint == d$endpoint & model == "multivariable"]
  official <- z[, .(
    TR = TF, Hazard_Ratio = HR, HR_Lower_95 = lower95, HR_Upper_95 = upper95,
    Significant = fifelse(p_value <= 0.05, "Yes", "No")
  )][Significant == "Yes"]
  if (!nrow(official)) {
    register(d$panel, "SUPPORTED_ZERO", "fig4_forest", note = "No nominally significant multivariable TF")
    next
  }
  plot <- forest_env$create_forest_plot(
    official, plot_title = d$title, cohort_name = d$cohort,
    survival_type = ifelse(d$endpoint == "RFS" & d$cohort == "TCGA_LUAD", "DFS", d$endpoint)
  )
  save_plot(d$panel, plot, 2.5, 4, "fig4_forest", sprintf("%d significant TFs", nrow(official)))
}

# ---------------------------------------------------------------------------
# Figure 4E/G: official protective and risk Euler constructor.
# ---------------------------------------------------------------------------

euler_env <- new.env(parent = globalenv())
source_block("fig4_euler", euler_env, patched = FALSE)
for (cohort_name in c("GSE41271_LUAD", "TCGA_LUAD")) {
  tag <- ifelse(cohort_name == "GSE41271_LUAD", "E", "G")
  endpoint2 <- ifelse(cohort_name == "TCGA_LUAD", "DFS", "RFS")
  z <- cox[cohort == cohort_name & model == "multivariable" & p_value <= 0.05]
  for (direction_name in c("Risk", "Protective")) {
    want_risk <- direction_name == "Risk"
    selected <- z[(HR > 1) == want_risk]
    sets <- list(
      OS = selected[endpoint == "OS", TF],
      selected[endpoint == "RFS", TF]
    )
    names(sets)[[2L]] <- endpoint2
    panel <- sprintf("Figure4%s_%s", tag, tolower(direction_name))
    path <- file.path(atomic, paste0(panel, ".pdf"))
    cairo_pdf(path, width = 3.5, height = 2.75)
    print(euler_env$plot_venn(
      sets,
      title = ifelse(cohort_name == "TCGA_LUAD", "TCGA-LUAD Multivariate", "GSE41271-LUAD Multivariate"),
      group = ifelse(want_risk, "poor", "good")
    ))
    dev.off()
    register(panel, "RENDERED_OFFICIAL_BLOCK", "fig4_euler", path)
  }
}

# ---------------------------------------------------------------------------
# Figure 4F/H: official ggsurvplot constructor, rebuilt from frozen receipts.
# ---------------------------------------------------------------------------

km_env <- new.env(parent = globalenv())
source_block("fig4_km", km_env, patched = FALSE)
km_objects <- readRDS(file.path(adapted, "figure4", "km_objects.rds"))
for (key in names(km_objects)) {
  obj <- km_objects[[key]]
  rec <- obj$receipt
  dat <- obj$data
  fit <- survfit(Surv(time, event) ~ activity_group, data = dat)
  lr <- survdiff(Surv(time, event) ~ activity_group, data = dat)
  p <- pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  tolerance <- max(1e-10, abs(rec$logrank_p) * 1e-6)
  if (!is.finite(p) || abs(p - rec$logrank_p) > tolerance) {
    stop("KM log-rank receipt mismatch: ", key)
  }
  gg <- dat$activity_group
  model <- coxph(Surv(dat$time, dat$event) ~ gg)
  result <- list(
    km_fit = fit,
    data = dat,
    logrank_test = list(p_value = p),
    univariable = list(model = model)
  )
  result_list <- setNames(list(result), rec$TF)
  endpoint_print <- ifelse(rec$cohort == "TCGA_LUAD" & rec$endpoint == "RFS", "DFS", rec$endpoint)
  plot <- km_env$plot_km_single(
    rec$TF, result_list,
    cohort_name = ifelse(rec$cohort == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD"),
    survival_type = endpoint_print
  )
  tag <- ifelse(rec$cohort == "TCGA_LUAD", "H", "F")
  panel <- sprintf("Figure4%s_%s_%s", tag, rec$TF, endpoint_print)
  save_plot(panel, plot, 5, 5.25, "fig4_km", sprintf("Frozen cutpoint %.8g", rec$cutpoint))
}

# ---------------------------------------------------------------------------
# Supplementary Figure 6A-H: execute the eight official Cox volcano blocks.
# ---------------------------------------------------------------------------

supp6_defs <- data.table(
  panel = paste0("SupplementaryFigure6", LETTERS[1:8]),
  block = c(
    "supp6a_uni_gse_os", "supp6b_uni_gse_rfs", "supp6c_uni_tcga_os", "supp6d_uni_tcga_rfs",
    "supp6e_multi_gse_os", "supp6f_multi_gse_rfs", "supp6g_multi_tcga_os", "supp6h_multi_tcga_rfs"
  ),
  object = c(
    "GSE41271_LUAD_OS_Univariate_plt", "GSE41271_LUAD_RFS_Univariate_plt",
    "TCGA_LUAD_OS_Univariate_plt", "TCGA_LUAD_RFS_Univariate_plt",
    "GSE41271_LUAD_OS_Multivariate_plt", "GSE41271_LUAD_RFS_Multivariate_plt",
    "TCGA_LUAD_OS_Multivariate_plt", "TCGA_LUAD_RFS_Multivariate_plt"
  )
)
for (i in seq_len(nrow(supp6_defs))) {
  d <- supp6_defs[i]
  env <- new.env(parent = globalenv())
  env$args <- list(
    GSE41271_LUAD_dir = file.path(adapted, "figure4"),
    TCGA_LUAD_dir = file.path(adapted, "figure4")
  )
  source_block(d$block, env, patched = TRUE)
  if (!exists(d$object, envir = env, inherits = FALSE)) stop("Official block did not create ", d$object)
  save_plot(d$panel, get(d$object, envir = env), 5, 5, d$block)
}

# ---------------------------------------------------------------------------
# Figure 5A-C: official cell-line concordance-index screens.
# ---------------------------------------------------------------------------

fig5_datasets <- c(Figure5A = "PRISM", Figure5B = "GDSC2", Figure5C = "CTRPv2")
for (panel in names(fig5_datasets)) {
  dataset_name <- fig5_datasets[[panel]]
  env <- new.env(parent = globalenv())
  env$ConcordanceIndex <- as.data.frame(read_tab(file.path(
    adapted, "figure5", "by_dataset", paste0(dataset_name, "_drug_TR_concordance_results.tsv")
  )))
  env$args <- list(cohort = dataset_name)
  source_block("fig5_ci", env, patched = FALSE)
  save_plot(panel, env$CI_plot, 4, 3.75, "fig5_ci", sprintf("%s; n=%d tests", dataset_name, nrow(env$plot_df)))
}

# Figure 5D: official ComplexUpset block. The adapter uses canonical drug IDs
# in Drug_TR so overlap is identity-safe across resources.
upset_env <- new.env(parent = globalenv())
for (dataset_name in c("PRISM", "GDSC2", "CTRPv2")) {
  z <- read_tab(file.path(
    adapted, "figure5", "by_dataset", paste0(dataset_name, "_drug_TR_concordance_results.tsv")
  ))
  z <- z[FDR <= 0.05]
  # The published block calls the GDSC1000 object `GDSC_CI` (without the
  # resource-version suffix).  Preserve that exact interface while adapting
  # only the cohort represented by the object from GDSC1000 to GDSC2.
  official_prefix <- if (dataset_name == "GDSC2") "GDSC" else dataset_name
  assign(paste0(official_prefix, "_CI_sig"), as.data.frame(z), envir = upset_env)
  assign(paste0(official_prefix, "_CI"), as.data.frame(z), envir = upset_env)
}
source_block("fig5_upset", upset_env, patched = TRUE)
save_plot("Figure5D", upset_env$Drug_TR_upset, 4, 5, "fig5_upset")

# Figure 5E: LUAD counts, official 3x3 constructor unchanged.
grid_env <- new.env(parent = globalenv())
grid_env$conf_df <- as.data.frame(read_tab(file.path(adapted, "figure5", "figure5e_confusion.tsv")))
grid_env$conf_df$OS <- factor(
  grid_env$conf_df$OS,
  levels = c("Good prognosis", "Non-Prognostic", "Poor prognosis")
)
grid_env$conf_df$RFS <- factor(
  grid_env$conf_df$RFS,
  levels = c("Recurrence-Free", "Non-Prognostic", "Recurrence")
)
source_block("fig5e_grid", grid_env, patched = FALSE)
save_plot("Figure5E", grid_env$conf_df_plt, 5, 2.75, "fig5e_grid")

# Figure 5F and Supplementary Figure 7C: official direction-separated
# pathway/TF dot-matrix blocks. Empty direction/category strata are recorded,
# not replaced by another plot.
combined <- read_tab(file.path(adapted, "figure5", "combined_replicated_ci.tsv"))
combined[, drug := display_drug]
dot_env <- new.env(parent = globalenv())
dot_env$args <- list(visual_outdir = atomic)
dot_env$Combined_CI_low <- as.data.frame(combined[Category == "Negative_Correlated_CI_Significant"])
dot_env$Combined_CI_high <- as.data.frame(combined[Category == "Positive_Correlated_CI_Significant"])
main_labels <- c("PI3K/Akt/mTOR", "MAPK/ERK", "EGFR/ERBB2", "Kinase Inhibitor", "Antimitotic")
if (nrow(dot_env$Combined_CI_low) && any(dot_env$Combined_CI_low$Label %in% main_labels)) {
  source_block("fig5f_low", dot_env, patched = TRUE)
  register(
    "Figure5F_low", "RENDERED_OFFICIAL_BLOCK", "fig5f_low",
    file.path(atomic, "Figure5F_Concordance_Index_Low_Sig_PI3K_MAPK_EGFR_Kinase.pdf")
  )
} else {
  register("Figure5F_low", "SUPPORTED_ZERO", "fig5f_low", note = "No low-CI pair in an official main pathway class")
}
if (nrow(dot_env$Combined_CI_high) && any(dot_env$Combined_CI_high$Label %in% main_labels)) {
  source_block("fig5f_high", dot_env, patched = TRUE)
  register(
    "Figure5F_high", "RENDERED_OFFICIAL_BLOCK", "fig5f_high",
    file.path(atomic, "Figure5F_Concordance_Index_High_Sig_PI3K_MAPK_EGFR_Kinase.pdf")
  )
} else {
  register("Figure5F_high", "SUPPORTED_ZERO", "fig5f_high", note = "No high-CI pair in an official main pathway class")
}

# ---------------------------------------------------------------------------
# Supplementary Figure 5A-B.
# ---------------------------------------------------------------------------

skew_env <- new.env(parent = globalenv())
skew_env$args <- list(visuals_outdir = atomic)
skew_env$skew_results_all <- as.data.frame(read_tab(file.path(
  adapted, "supplementary5", "skewness_official.tsv"
)))
source_block("supp5a_skew", skew_env, patched = TRUE)
register(
  "SupplementaryFigure5A", "RENDERED_OFFICIAL_BLOCK", "supp5a_skew",
  file.path(atomic, "SupplementaryFigure5A_Skewness_PDX_CellLines.pdf")
)

estimate <- read_tab(file.path(adapted, "supplementary5", "estimate_official.tsv"))
volcano_env <- new.env(parent = globalenv())
source_block("supp5b_volcano", volcano_env, patched = FALSE)
cohort_order <- c("TCGA-LUAD", "GSE41271-LUAD")
score_order <- c("StromalScore", "ImmuneScore", "TumorPurity")
volcanoes <- list()
for (cohort_name in cohort_order) {
  for (score_name in score_order) {
    z <- as.data.frame(estimate[cohort == cohort_name & marker == score_name])
    volcanoes[[paste(cohort_name, score_name, sep = "__")]] <- volcano_env$make_volcano(
      z, sprintf("%s — %s", score_name, cohort_name), 0.4, 0.05, 40
    )
  }
}
save_plot(
  "SupplementaryFigure5B_volcano",
  wrap_plots(volcanoes, ncol = length(score_order), nrow = length(cohort_order)),
  15, 8, "supp5b_volcano"
)

heat_env <- new.env(parent = globalenv())
heat_env$cor_all <- as.data.frame(estimate)
heat_env$score_order <- score_order
heat_env$cohort_order <- cohort_order
heat_env$args <- list(
  min_cohorts = 1L,
  fdr_threshold = 0.05,
  r_threshold = 0.4,
  col_width_mm = 10,
  row_height_mm = 4,
  plot_title = "HC-TR activity correlation",
  out_name = "SupplementaryFigure5B_ESTIMATE_Combined_Heatmap_HC-TR",
  visuals_out = atomic
)
# Lines 87-247 of the official block perform its own significance filtering,
# constructor, dynamic device sizing, Helvetica rendering, FDR-star footer,
# and device closure.  No wrapper-side font or constructor change is allowed.
source_block("supp5b_heatmap", heat_env, patched = FALSE)
heat_path <- file.path(
  atomic, "SupplementaryFigure5B_ESTIMATE_Combined_Heatmap_HC-TR.pdf"
)
register("SupplementaryFigure5B_heatmap", "RENDERED_OFFICIAL_BLOCK", "supp5b_heatmap", heat_path)

# ---------------------------------------------------------------------------
# Supplementary Figure 7A-C.
# ---------------------------------------------------------------------------

membership <- read_tab(file.path(adapted, "figure5", "eligible_drug_membership.tsv"))
all_drugs <- sort(unique(membership$canonical_drug))
wide <- data.table(Standardized_Name = all_drugs)
for (dataset_name in c("PRISM", "CTRPv2", "GDSC2")) {
  present <- membership[dataset == dataset_name, canonical_drug]
  wide[, (dataset_name) := fifelse(Standardized_Name %in% present, Standardized_Name, NA_character_)]
}
drug_env <- new.env(parent = globalenv())
drug_env$drug_cohorts <- as.data.frame(wide)
drug_env$args <- list(visual_outdir = atomic)
drug_env$drug_dataset_colours <- c(
  PRISM = unname(drug_dataset_colours[["GRAY"]]),
  CTRPv2 = unname(drug_dataset_colours[["CTRPv2"]]),
  GDSC2 = unname(drug_dataset_colours[["GDSC1000"]])
)
source_block("supp7a_drugs", drug_env, patched = TRUE)
register(
  "SupplementaryFigure7A", "RENDERED_OFFICIAL_BLOCK", "supp7a_drugs",
  file.path(atomic, "SupplementaryFigure7A_CellLine_DrugSensitivity_Cohort_DrugIntersection.pdf")
)

examples_path <- file.path(adapted, "figure5", "supp7b_example_values.tsv")
examples <- read_tab(examples_path)
example_pairs <- read_tab(file.path(adapted, "figure5", "supp7b_example_pairs.tsv"))
if (!nrow(examples) || nrow(example_pairs) != 4L ||
    !setequal(example_pairs$case_block, paste0("supp7b_case", 1:4))) {
  stop("Strict Supplementary Figure 7B requires all four frozen official call slots")
}
aac_list <- list()
for (dataset_name in unique(examples$dataset)) {
  z <- examples[dataset == dataset_name]
  m <- dcast(z, display_drug ~ sample_id, value.var = "AAC", fun.aggregate = mean)
  rn <- m$display_drug
  m[, display_drug := NULL]
  mat <- as.matrix(m)
  rownames(mat) <- rn
  aac_list[[dataset_name]] <- mat
}
nes_wide <- dcast(examples, TF ~ sample_id, value.var = "NES", fun.aggregate = mean)
rn <- nes_wide$TF
nes_wide[, TF := NULL]
nes_mat <- as.matrix(nes_wide)
rownames(nes_mat) <- rn
ci_sub <- combined[
  canonical_drug %in% example_pairs$canonical_drug &
  TR %in% example_pairs$TR
]
ci_sub[, `:=`(drug = display_drug, Cohort = dataset)]

# Execute the four published call/save blocks.  Their default cairo_pdf()
# device is intentionally left at the official 7 x 7 inches; only frozen LUAD
# identifiers, cohort IDs, and the observed low-CI colour object are adapted
# in the materialized, audited diffs.
case_env <- new.env(parent = globalenv())
case_env$prep <- list(
  AAC_list = aac_list,
  NES_TR = nes_mat,
  combined_CI = as.data.frame(ci_sub)
)
case_env$args <- list(visuals_outdir = atomic)
case_env$AAC_concordance_colours_Low <- c(
  PRISM = unname(AAC_concordance_colours_Low[["GRAY"]]),
  GDSC2 = unname(AAC_concordance_colours_Low[["GDSC1000"]]),
  CTRPv2 = unname(AAC_concordance_colours_Low[["CTRPv2"]])
)
case_env$AAC_concordance_colours_High <- c(
  PRISM = unname(AAC_concordance_colours_High[["GRAY"]]),
  GDSC2 = unname(AAC_concordance_colours_High[["GDSC1000"]]),
  CTRPv2 = unname(AAC_concordance_colours_High[["CTRPv2"]])
)
source_block("aac_helper", case_env, patched = FALSE)
case_outputs <- data.table(
  block_id = paste0("supp7b_case", 1:4),
  panel = c(
    "SupplementaryFigure7B_case1_WWC2_Pha793887",
    "SupplementaryFigure7B_case2_ZNF254_5Fluorouracil",
    "SupplementaryFigure7B_case3_ZNF254_Dactolisib",
    "SupplementaryFigure7B_case4_ZNF254_Vorinostat"
  ),
  filename = c(
    "SupplementaryFigure7B_WWC2_Pha793887_AAC_NES.pdf",
    "SupplementaryFigure7B_ZNF254_5Fluorouracil_AAC_NES.pdf",
    "SupplementaryFigure7B_ZNF254_Dactolisib_AAC_NES.pdf",
    "SupplementaryFigure7B_ZNF254_Vorinostat_AAC_NES.pdf"
  )
)
for (i in seq_len(nrow(case_outputs))) {
  d <- case_outputs[i]
  source_block(d$block_id, case_env, patched = TRUE)
  case_path <- file.path(atomic, d$filename)
  if (!file.exists(case_path)) stop("Official Supplementary Figure 7B block did not write: ", case_path)
  register(
    d$panel, "RENDERED_OFFICIAL_BLOCK",
    c("aac_helper", d$block_id), case_path,
    note = "Official cairo_pdf default device: 7 x 7 inches"
  )
}

# Execute official Supplementary 7C blocks only for strata that have an
# official pathway class. All unmatched 'Other' rows remain in the table.
supp_labels <- c(
  "Small Molecule Inhibitor", "Anti-metabolite", "Bromodomain",
  "FGFR1/VEGFR inhibitor", "Telomerase inhibitor", "tyrosine kinase inhibitor",
  "DNA replication", "Apoptosis", "NAE inhibitor", "Survivin inhibitor",
  "inducer of DNA damage"
)
if (nrow(dot_env$Combined_CI_low) && any(dot_env$Combined_CI_low$Label %in% supp_labels)) {
  source_block("supp7c_low", dot_env, patched = TRUE)
  register(
    "SupplementaryFigure7C_low", "RENDERED_OFFICIAL_BLOCK", "supp7c_low",
    file.path(atomic, "SupplementaryFigure7C_Concordance_Index_Low_Sig_Supplementary.pdf")
  )
} else {
  register("SupplementaryFigure7C_low", "SUPPORTED_ZERO", "supp7c_low", note = "No low-CI pair in an official supplementary pathway class")
}
if (nrow(dot_env$Combined_CI_high) && any(!(dot_env$Combined_CI_high$Label %in% main_labels))) {
  source_block("supp7c_high", dot_env, patched = TRUE)
  register(
    "SupplementaryFigure7C_high", "RENDERED_OFFICIAL_BLOCK", "supp7c_high",
    file.path(atomic, "SupplementaryFigure7C_Concordance_Index_High_Sig_Supplementary.pdf")
  )
} else {
  register("SupplementaryFigure7C_high", "SUPPORTED_ZERO", "supp7c_high", note = "No high-CI pair in an official supplementary pathway class")
}

# Figure 5G/H and Supplementary 7D: no graphics are generated.
unsupported <- read_tab(file.path(adapted, "unsupported_panels.tsv"))
for (i in seq_len(nrow(unsupported))) {
  register(
    unsupported$panel[[i]], "UNSUPPORTED",
    switch(
      unsupported$panel[[i]],
      Figure5G = "fig5g_pdx",
      Figure5H = "fig5h_waterfall",
      SupplementaryFigure7D = "supp7d_pdx"
    ),
    note = unsupported$reason[[i]]
  )
}
fwrite(unsupported, file.path(out, "unsupported_panels.tsv"), sep = "\t")

manifest <- rbindlist(status_rows, fill = TRUE)
fwrite(manifest, file.path(out, "render_manifest.tsv"), sep = "\t")
fwrite(executed_blocks, file.path(audit, "executed_official_blocks.tsv"), sep = "\t")
cat(sprintf(
  "Rendered %d official-block PDF entries; %d panels are explicitly unsupported.\n",
  manifest[status == "RENDERED_OFFICIAL_BLOCK", .N],
  manifest[status == "UNSUPPORTED", .N]
))

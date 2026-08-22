#!/usr/bin/env Rscript

# Plot construction is sourced from the checksum-locked TNBC Code Ocean
# blocks. This wrapper supplies only LUAD HC45 results, labels and devices.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(survival)
  library(survminer)
  library(eulerr)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: render_figure4_official_logrank.R GENERATED_BLOCKS ANALYSIS_ROOT OUTPUT_ROOT")
}
generated <- normalizePath(args[[1L]], mustWork = TRUE)
analysis <- normalizePath(args[[2L]], mustWork = TRUE)
out <- normalizePath(args[[3L]], mustWork = FALSE)
atomic <- file.path(out, "atomic")
audit <- file.path(out, "audit")
dir.create(atomic, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

block_path <- function(id, patched = TRUE) {
  path <- file.path(generated, if (patched) "patched" else "exact", paste0(id, ".R"))
  if (!file.exists(path)) stop("Missing official block: ", path)
  path
}
executed <- data.table(sequence = integer(), block_id = character(), variant = character())
source_block <- function(id, env = parent.frame(), patched = TRUE) {
  executed <<- rbind(executed, data.table(
    sequence = nrow(executed) + 1L,
    block_id = id,
    variant = ifelse(patched, "patched", "exact")
  ))
  sys.source(block_path(id, patched), envir = env, keep.source = TRUE)
}
save_plot <- function(name, plot, width, height) {
  path <- file.path(atomic, paste0(name, ".pdf"))
  cairo_pdf(path, width = width, height = height, onefile = FALSE)
  print(plot)
  dev.off()
  path
}

source_block("aesthetics", globalenv(), patched = FALSE)
results <- fread(file.path(analysis, "results", "tables", "Figure4_logrank_all_HC45.tsv"))

# Figure 4A-D: official forest grammar; rows are selected by the user-frozen
# unadjusted log-rank rule and carry the matching binary-group HR/CI.
forest_env <- new.env(parent = globalenv())
source_block("fig4_forest", forest_env, patched = FALSE)
defs <- data.table(
  panel = paste0("Figure4", LETTERS[1:4]),
  cohort = c("GSE41271_LUAD", "GSE41271_LUAD", "TCGA_LUAD", "TCGA_LUAD"),
  endpoint = c("OS", "RFS", "OS", "RFS"),
  title = c(
    "GSE41271-LUAD OS",
    "GSE41271-LUAD RFS",
    "TCGA-LUAD OS",
    "TCGA-LUAD DFS"
  )
)
render_rows <- list()
for (i in seq_len(nrow(defs))) {
  d <- defs[i]
  z <- results[
    cohort == d$cohort & endpoint == d$endpoint & nominal_unadjusted_logrank_significant
  ]
  official <- z[, .(
    TR = TF, Hazard_Ratio = HR, HR_Lower_95 = lower95, HR_Upper_95 = upper95,
    Significant = "Yes"
  )]
  plot <- forest_env$create_forest_plot(
    official,
    plot_title = d$title,
    cohort_name = d$cohort,
    survival_type = ifelse(d$endpoint == "RFS" & d$cohort == "TCGA_LUAD", "DFS", d$endpoint)
  )
  height <- max(4.0, min(10.0, 1.7 + 0.19 * nrow(official)))
  path <- save_plot(d$panel, plot, 3.0, height)
  render_rows[[length(render_rows) + 1L]] <- data.table(
    panel = d$panel, official_block = "fig4_forest", rows = nrow(official), output = path,
    analysis = "unadjusted selected-cutpoint log-rank P<=0.05"
  )
}

# Figure 4E/G: exact official Euler constructor.
euler_env <- new.env(parent = globalenv())
source_block("fig4_euler", euler_env, patched = FALSE)
for (cohort_name in c("GSE41271_LUAD", "TCGA_LUAD")) {
  tag <- ifelse(cohort_name == "GSE41271_LUAD", "E", "G")
  endpoint2 <- ifelse(cohort_name == "TCGA_LUAD", "DFS", "RFS")
  z <- results[cohort == cohort_name & nominal_unadjusted_logrank_significant]
  for (direction_name in c("Risk", "Protective")) {
    want_risk <- direction_name == "Risk"
    selected <- z[(HR > 1) == want_risk]
    sets <- list(OS = selected[endpoint == "OS", TF], selected[endpoint == "RFS", TF])
    names(sets)[[2L]] <- endpoint2
    panel <- sprintf("Figure4%s_%s", tag, tolower(direction_name))
    plot <- euler_env$plot_venn(
      sets,
      # Keep the atom title inside the checksum-locked official device.  The
      # full primary rule is carried by the figure legend and audit receipt.
      title = ifelse(cohort_name == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD"),
      group = ifelse(want_risk, "poor", "good")
    )
    path <- save_plot(panel, plot, 3.5, 2.75)
    render_rows[[length(render_rows) + 1L]] <- data.table(
      panel = panel, official_block = "fig4_euler", rows = length(unique(unlist(sets))),
      output = path, analysis = "unadjusted selected-cutpoint log-rank P<=0.05"
    )
  }
}

# Figure 4F/H: exact official Kaplan-Meier constructor.
km_env <- new.env(parent = globalenv())
source_block("fig4_km", km_env, patched = FALSE)
km_objects <- readRDS(file.path(analysis, "results", "Figure4_km_objects.rds"))
for (key in names(km_objects)) {
  obj <- km_objects[[key]]
  rec <- obj$receipt
  dat <- obj$data
  fit <- survfit(Surv(time, event) ~ activity_group, data = dat)
  lr <- survdiff(Surv(time, event) ~ activity_group, data = dat)
  p <- pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  if (abs(p - rec$logrank_p) > max(1e-10, rec$logrank_p * 1e-6)) {
    stop("KM log-rank receipt mismatch: ", key)
  }
  # The checksum-locked TNBC KM constructor indexes the Cox confidence matrix
  # by the original author's model term name, `ggHigh`.  Preserve that input
  # contract without changing the LUAD grouping or survival estimates.
  gg <- dat$activity_group
  model <- coxph(Surv(dat$time, dat$event) ~ gg)
  result <- list(
    km_fit = fit, data = dat, logrank_test = list(p_value = p),
    univariable = list(model = model)
  )
  endpoint_print <- ifelse(rec$cohort == "TCGA_LUAD" & rec$endpoint == "RFS", "DFS", rec$endpoint)
  plot <- km_env$plot_km_single(
    rec$TF, setNames(list(result), rec$TF),
    cohort_name = ifelse(rec$cohort == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD"),
    survival_type = endpoint_print
  )
  tag <- ifelse(rec$cohort == "TCGA_LUAD", "H", "F")
  panel <- sprintf("Figure4%s_%s_%s", tag, rec$TF, endpoint_print)
  path <- save_plot(panel, plot, 5.0, 5.25)
  render_rows[[length(render_rows) + 1L]] <- data.table(
    panel = panel, official_block = "fig4_km", rows = nrow(dat), output = path,
    analysis = sprintf("cutpoint %.8g; unadjusted log-rank P %.4g", rec$cutpoint, p)
  )
}

# Supplementary Figure 6A-D: exact official univariable volcano blocks for the
# log-rank primary screen. E-H retain the 45-TF multivariable median-split
# sensitivity and are never used to decide the main Figure 4 display.
supp_defs <- data.table(
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
input_dir <- file.path(analysis, "results", "official_inputs")
for (i in seq_len(nrow(supp_defs))) {
  d <- supp_defs[i]
  env <- new.env(parent = globalenv())
  env$args <- list(GSE41271_LUAD_dir = input_dir, TCGA_LUAD_dir = input_dir)
  source_block(d$block, env, patched = TRUE)
  if (!exists(d$object, envir = env, inherits = FALSE)) stop("Official block did not create ", d$object)
  path <- save_plot(d$panel, get(d$object, envir = env), 5, 5)
  render_rows[[length(render_rows) + 1L]] <- data.table(
    panel = d$panel, official_block = d$block, rows = 45L, output = path,
    analysis = ifelse(i <= 4L, "primary unadjusted log-rank", "multivariable sensitivity")
  )
}

manifest <- rbindlist(render_rows, fill = TRUE)
fwrite(manifest, file.path(out, "render_manifest.tsv"), sep = "\t")
fwrite(executed, file.path(audit, "executed_official_blocks.tsv"), sep = "\t")
cat(sprintf("Rendered %d official-code Figure 4 atoms.\n", nrow(manifest)))

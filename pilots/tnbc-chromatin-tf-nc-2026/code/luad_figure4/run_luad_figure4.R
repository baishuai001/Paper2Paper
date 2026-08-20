#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: run_luad_figure4.R FIGURE1_ROOT FIGURE2_ROOT OUTPUT_ROOT")
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

hc_table <- fread(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
hc <- sort(hc_table[
  analysis_definition == "PRIMARY_anchor_promoter_tcga_cpm1_both" &
    HC_TF_promoter_activity_definition == TRUE,
  unique(TF)
])
if (length(hc) != 31L) stop(sprintf("Expected 31 HC-TFs, found %d", length(hc)))

read_json_pages <- function(directory, study, kind) {
  paths <- sort(Sys.glob(file.path(directory, sprintf("%s_%s_page*.json", study, kind))))
  if (!length(paths)) stop(sprintf("No cBioPortal pages for %s %s", study, kind))
  rbindlist(lapply(paths, function(path) as.data.table(fromJSON(path, flatten = TRUE))), fill = TRUE)
}
pivot_clinical <- function(tab) {
  tab <- tab[!is.na(patientId) & !is.na(clinicalAttributeId)]
  dcast(tab, patientId ~ clinicalAttributeId, value.var = "value",
        fun.aggregate = function(x) if (length(x)) x[[1L]] else NA_character_,
        fill = NA_character_)
}
parse_event <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fifelse(grepl("^1|DECEASED|RECUR|PROGRESS", x), 1L,
          fifelse(grepl("^0|LIVING|DISEASE.?FREE", x), 0L, NA_integer_))
}
stage_group <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE[[:space:]]*", "", x)
  fifelse(grepl("^IV", x), "IV", fifelse(grepl("^III", x), "III",
    fifelse(grepl("^II", x), "II", fifelse(grepl("^I", x), "I", NA_character_))))
}

tcga_activity <- read_matrix(file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz"))
tcga_manifest <- fread(file.path(fig1, "data", "processed", "tcga_manifest.tsv"))
tcga_annotation <- fread(file.path(fig1, "data", "processed", "annotations", "tcga_figure1_annotations.tsv"))
tcga_ids <- tcga_manifest[group == "LUAD", sample_id]
if (length(tcga_ids) != 516L || length(setdiff(hc, rownames(tcga_activity)))) {
  stop("TCGA survival input does not contain 516 LUADs and all 31 HC-TFs")
}
cbio <- pivot_clinical(read_json_pages(
  file.path(fig1, "data", "raw", "annotations", "cbioportal"),
  "luad_tcga_pan_can_atlas_2018", "clinical_patient"
))
tcga_clinical <- tcga_annotation[group == "LUAD"]
tcga_clinical <- tcga_clinical[match(tcga_ids, sample_id)]
tcga_clinical <- merge(tcga_clinical, cbio, by.x = "patient_id", by.y = "patientId",
                       all.x = TRUE, sort = FALSE, suffixes = c("", ".cbio"))
gdc_smoking <- fread(file.path(out, "data", "tcga_luad_smoking.tsv"))
tcga_clinical <- merge(tcga_clinical, gdc_smoking, by = "patient_id", all.x = TRUE, sort = FALSE)
tcga_clinical <- tcga_clinical[match(tcga_ids, sample_id)]
tcga_clinical[, `:=`(
  cohort = "TCGA_LUAD",
  OS_time = suppressWarnings(as.numeric(OS_MONTHS)),
  OS_event = parse_event(OS_STATUS),
  RFS_time = suppressWarnings(as.numeric(DFS_MONTHS)),
  RFS_event = parse_event(DFS_STATUS),
  age = suppressWarnings(as.numeric(age_at_diagnosis)),
  sex_model = factor(sex),
  stage_model = factor(stage_group(stage)),
  smoking_model = factor(fifelse(smoking_model %in% c("Ever", "Never"), smoking_model, NA_character_))
)]

gse_activity_independent <- read_matrix(file.path(fig1, "results", "gse41271", "GSE41271_viper_activity.tsv.gz"))
gse_activity <- read_matrix(file.path(
  out, "results", "gse41271_tcga_projection",
  "GSE41271_TCGA_PROJECTED_viper_activity.tsv.gz"
))
gse_manifest <- fread(file.path(fig1, "data", "processed", "gse41271_manifest.tsv"))
gse_ids <- gse_manifest[group == "LUAD", sample_id]
if (length(gse_ids) != 183L || length(setdiff(hc, rownames(gse_activity)))) {
  stop(sprintf("GSE41271 TCGA-regulon projection missing LUAD samples or HC-TFs: samples=%d, missing=%s",
               length(gse_ids), paste(setdiff(hc, rownames(gse_activity)), collapse = ",")))
}
gse_independent_available <- intersect(hc, rownames(gse_activity_independent))
gse_independent_missing <- setdiff(hc, gse_independent_available)
if (length(gse_independent_available) != 26L ||
    !identical(sort(gse_independent_missing), sort(c("CREBRF", "CSRNP1", "ZBTB18", "ZNF33B", "ZNF75D")))) {
  stop(sprintf("Unexpected GSE41271 independent-network HC-TF coverage: %d available; missing %s",
               length(gse_independent_available), paste(gse_independent_missing, collapse = ",")))
}
gse_clinical <- gse_manifest[group == "LUAD"]
gse_clinical <- gse_clinical[match(gse_ids, sample_id)]
gse_clinical[, `:=`(
  cohort = "GSE41271_LUAD",
  OS_time = os_days / 30.4375,
  OS_event = as.integer(os_event),
  RFS_time = rfs_days / 30.4375,
  RFS_event = as.integer(rfs_event),
  age = as.numeric(age_at_surgery),
  sex_model = factor(ifelse(toupper(sex) %in% c("M", "MALE"), "Male",
                            ifelse(toupper(sex) %in% c("F", "FEMALE"), "Female", NA))),
  stage_model = factor(stage_group(stage)),
  smoking_model = factor(ifelse(toupper(smoking) == "Y", "Ever",
                                ifelse(toupper(smoking) == "N", "Never", NA)))
)]

cohorts <- list(
  GSE41271_LUAD = list(clinical = gse_clinical, activity = gse_activity[, gse_ids, drop = FALSE]),
  TCGA_LUAD = list(clinical = tcga_clinical, activity = tcga_activity[, tcga_ids, drop = FALSE])
)

endpoint_data <- function(cohort_name, endpoint, tf) {
  object <- cohorts[[cohort_name]]
  clinical <- copy(object$clinical)
  clinical[, activity := as.numeric(object$activity[tf, sample_id])]
  time_col <- paste0(endpoint, "_time")
  event_col <- paste0(endpoint, "_event")
  clinical[, `:=`(time_raw = suppressWarnings(as.numeric(get(time_col))),
                  event_raw = suppressWarnings(as.integer(get(event_col))))]
  clinical <- clinical[is.finite(time_raw) & time_raw > 0 & event_raw %in% c(0L, 1L) &
                         is.finite(activity) & is.finite(age) & !is.na(sex_model) &
                         !is.na(stage_model) & !is.na(smoking_model)]
  clinical[, `:=`(time = pmin(time_raw, 120),
                  event = as.integer(event_raw == 1L & time_raw <= 120))]
  clinical[, activity_high := as.integer(activity > median(activity))]
  as.data.table(droplevels(as.data.frame(clinical)))
}

extract_cox <- function(fit, exposure = "activity_high") {
  summary_fit <- summary(fit)
  row <- summary_fit$coefficients[exposure, ]
  ci <- summary_fit$conf.int[exposure, ]
  data.table(
    beta = unname(row[["coef"]]), HR = unname(ci[["exp(coef)"]]),
    lower95 = unname(ci[["lower .95"]]), upper95 = unname(ci[["upper .95"]]),
    p_value = unname(row[["Pr(>|z|)"]]),
    model_likelihood_p = unname(summary_fit$logtest[["pvalue"]])
  )
}

model_rows <- list()
survival_inputs <- list()
for (cohort_name in names(cohorts)) {
  for (endpoint in c("OS", "RFS")) {
    reference <- endpoint_data(cohort_name, endpoint, hc[[1L]])
    wide <- reference[, .(sample_id, time, event, age, sex_model, stage_model, smoking_model)]
    activity_matrix <- cohorts[[cohort_name]]$activity[hc, wide$sample_id, drop = FALSE]
    activity_high <- t(apply(activity_matrix, 1L, function(x) as.integer(x > median(x))))
    rownames(activity_high) <- hc
    colnames(activity_high) <- wide$sample_id
    survival_inputs[[paste(cohort_name, endpoint, sep = "__")]] <- list(
      cohort = cohort_name, endpoint = endpoint, clinical = wide,
      activity_high = activity_high, HC_TFs = hc
    )
    for (tf in hc) {
      dat <- copy(wide)
      dat[, activity_high := as.integer(activity_high[tf, sample_id])]
      univ <- coxph(Surv(time, event) ~ activity_high, data = dat, ties = "efron", x = FALSE)
      multi <- coxph(Surv(time, event) ~ activity_high + age + sex_model + stage_model + smoking_model,
                     data = dat, ties = "efron", x = FALSE)
      model_rows[[length(model_rows) + 1L]] <- cbind(
        data.table(cohort = cohort_name, endpoint = endpoint, model = "univariable",
                   TF = tf, n = nrow(dat), events = sum(dat$event)), extract_cox(univ)
      )
      model_rows[[length(model_rows) + 1L]] <- cbind(
        data.table(cohort = cohort_name, endpoint = endpoint, model = "multivariable",
                   TF = tf, n = nrow(dat), events = sum(dat$event)), extract_cox(multi)
      )
    }
  }
}
models <- rbindlist(model_rows)
models[, FDR := p.adjust(p_value, method = "BH"), by = .(cohort, endpoint, model)]
models[, `:=`(log2_HR = log2(HR), log2_lower95 = log2(lower95), log2_upper95 = log2(upper95),
              direction = fifelse(HR < 1, "Favorable", "Adverse"),
              nominal_significant = p_value <= 0.05,
              FDR_significant = FDR <= 0.05)]
fwrite(models, file.path(tables, "Figure4_all_Cox_models.tsv"), sep = "\t")
saveRDS(survival_inputs, file.path(out, "results", "figure4_survival_inputs.rds"), compress = "xz")

# Sensitivity analysis: repeat the GSE41271 models for the 26 HC-TFs present in
# its independently inferred one-subnet regulon.  These estimates never replace
# the frozen 31-TF TCGA-regulon projection used by the main Figure 4 panels.
gse_sensitivity_rows <- list()
for (endpoint in c("OS", "RFS")) {
  for (tf in gse_independent_available) {
    dat <- copy(gse_clinical)
    dat[, activity := as.numeric(gse_activity_independent[tf, sample_id])]
    time_col <- paste0(endpoint, "_time")
    event_col <- paste0(endpoint, "_event")
    dat[, `:=`(time_raw = suppressWarnings(as.numeric(get(time_col))),
               event_raw = suppressWarnings(as.integer(get(event_col))))]
    dat <- dat[is.finite(time_raw) & time_raw > 0 & event_raw %in% c(0L, 1L) &
                 is.finite(activity) & is.finite(age) & !is.na(sex_model) &
                 !is.na(stage_model) & !is.na(smoking_model)]
    dat[, `:=`(time = pmin(time_raw, 120),
               event = as.integer(event_raw == 1L & time_raw <= 120),
               activity_high = as.integer(activity > median(activity)))]
    dat <- as.data.table(droplevels(as.data.frame(dat)))
    for (model_name in c("univariable", "multivariable")) {
      fit <- if (model_name == "univariable") {
        coxph(Surv(time, event) ~ activity_high, data = dat, ties = "efron")
      } else {
        coxph(Surv(time, event) ~ activity_high + age + sex_model + stage_model + smoking_model,
              data = dat, ties = "efron")
      }
      gse_sensitivity_rows[[length(gse_sensitivity_rows) + 1L]] <- cbind(
        data.table(cohort = "GSE41271_LUAD", endpoint = endpoint, model = model_name,
                   TF = tf, n = nrow(dat), events = sum(dat$event)), extract_cox(fit)
      )
    }
  }
}
gse_sensitivity <- rbindlist(gse_sensitivity_rows)
gse_sensitivity[, FDR := p.adjust(p_value, method = "BH"), by = .(endpoint, model)]
gse_sensitivity[, `:=`(log2_HR = log2(HR), direction = fifelse(HR < 1, "Favorable", "Adverse"),
                       nominal_significant = p_value <= 0.05, FDR_significant = FDR <= 0.05)]
fwrite(gse_sensitivity, file.path(tables, "SupplementaryFigure7_GSE41271_independent_network_Cox_sensitivity.tsv"), sep = "\t")
gse_definition_comparison <- merge(
  models[cohort == "GSE41271_LUAD" & TF %in% gse_independent_available,
         .(endpoint, model, TF, projected_log2_HR = log2_HR, projected_p = p_value)],
  gse_sensitivity[, .(endpoint, model, TF, independent_log2_HR = log2_HR,
                      independent_p = p_value)],
  by = c("endpoint", "model", "TF")
)
fwrite(gse_definition_comparison,
       file.path(tables, "SupplementaryFigure7_GSE41271_activity_definition_comparison.tsv"), sep = "\t")

panel_map <- data.table(
  panel = c("A", "B", "C", "D"),
  cohort = c("GSE41271_LUAD", "GSE41271_LUAD", "TCGA_LUAD", "TCGA_LUAD"),
  endpoint = c("OS", "RFS", "OS", "RFS"),
  title = c("GSE41271 OS", "GSE41271 RFS", "TCGA OS", "TCGA DFS")
)
for (i in seq_len(nrow(panel_map))) {
  info <- panel_map[i]
  tab <- models[cohort == info$cohort & endpoint == info$endpoint & model == "multivariable"]
  setorder(tab, log2_HR)
  tab[, TF_factor := factor(TF, levels = TF)]
  tab[, plot_class := fifelse(!nominal_significant, "Not significant", direction)]
  panel <- ggplot(tab, aes(log2_HR, TF_factor, color = plot_class)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#555555") +
    geom_errorbarh(aes(xmin = log2_lower95, xmax = log2_upper95), height = 0.15, linewidth = 0.45) +
    geom_point(shape = 18, size = 2.4) +
    scale_color_manual(values = c(Favorable = "#2166AC", Adverse = "#B2182B", `Not significant` = "#BDBDBD")) +
    labs(title = sprintf("%s  Multivariable %s", info$panel, info$title),
         subtitle = sprintf("n=%d; events=%d; median activity cutpoint", unique(tab$n), unique(tab$events)),
         x = expression(log[2]("hazard ratio")~"[95% CI]"), y = NULL, color = NULL) +
    theme_classic(base_size = 8.2) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  filename <- sprintf("Figure4%s_%s_%s_multivariable_forest", info$panel, info$cohort, info$endpoint)
  ggsave(file.path(figures, paste0(filename, ".pdf")), panel, width = 5.2, height = 7.2, useDingbats = FALSE)
  ggsave(file.path(figures, paste0(filename, ".png")), panel, width = 5.2, height = 7.2, dpi = 220)
}

significant <- models[model == "multivariable" & nominal_significant == TRUE]
fwrite(significant, file.path(tables, "Figure4_nominal_multivariable_prognostic_TFs.tsv"), sep = "\t")
overlap_rows <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  rbindlist(lapply(c("Favorable", "Adverse"), function(direction_name) {
    os <- significant[cohort == cohort_name & endpoint == "OS" & direction == direction_name, TF]
    rfs <- significant[cohort == cohort_name & endpoint == "RFS" & direction == direction_name, TF]
    data.table(
      cohort = cohort_name, direction = direction_name,
      OS_only = length(setdiff(os, rfs)), overlap = length(intersect(os, rfs)),
      RFS_only = length(setdiff(rfs, os)), overlap_TFs = paste(sort(intersect(os, rfs)), collapse = ";")
    )
  }))
}))
fwrite(overlap_rows, file.path(tables, "Figure4EG_endpoint_overlap.tsv"), sep = "\t")

render_overlap <- function(cohort_name, panel, path, png_device = FALSE) {
  if (png_device) png(path, width = 1500, height = 900, res = 220) else cairo_pdf(path, width = 7.0, height = 4.2)
  grid.newpage()
  grid.text(sprintf("%s  %s multivariable endpoint overlap", panel,
                    ifelse(cohort_name == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD")),
            unit(0.03, "npc"), unit(0.95, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 11, fontface = "bold"))
  for (j in seq_along(c("Favorable", "Adverse"))) {
    direction_name <- c("Favorable", "Adverse")[[j]]
    row <- overlap_rows[cohort == cohort_name & direction == direction_name]
    center <- c(0.27, 0.73)[[j]]
    color <- c("#2166AC", "#B2182B")[[j]]
    grid.circle(unit(center - 0.065, "npc"), unit(0.56, "npc"), r = unit(0.115, "npc"),
                gp = gpar(fill = adjustcolor(color, alpha.f = 0.14), col = color, lwd = 1.2))
    grid.circle(unit(center + 0.065, "npc"), unit(0.56, "npc"), r = unit(0.115, "npc"),
                gp = gpar(fill = adjustcolor(color, alpha.f = 0.14), col = color, lwd = 1.2))
    grid.text("OS", unit(center - 0.12, "npc"), unit(0.72, "npc"), gp = gpar(fontsize = 8))
    grid.text(ifelse(cohort_name == "TCGA_LUAD", "DFS", "RFS"),
              unit(center + 0.12, "npc"), unit(0.72, "npc"), gp = gpar(fontsize = 8))
    grid.text(row$OS_only, unit(center - 0.085, "npc"), unit(0.56, "npc"), gp = gpar(fontsize = 11))
    grid.text(row$overlap, unit(center, "npc"), unit(0.56, "npc"), gp = gpar(fontsize = 11, fontface = "bold"))
    grid.text(row$RFS_only, unit(center + 0.085, "npc"), unit(0.56, "npc"), gp = gpar(fontsize = 11))
    grid.text(direction_name, unit(center, "npc"), unit(0.30, "npc"), gp = gpar(fontsize = 8.5, fontface = "bold", col = color))
    label <- ifelse(row$overlap_TFs == "", "No shared TF", gsub(";", ", ", row$overlap_TFs))
    grid.text(label, unit(center, "npc"), unit(0.21, "npc"), gp = gpar(fontsize = 6.1, col = "#555555"))
  }
  dev.off()
}
render_overlap("GSE41271_LUAD", "E", file.path(figures, "Figure4E_GSE41271_endpoint_overlap.pdf"), FALSE)
render_overlap("GSE41271_LUAD", "E", file.path(figures, "Figure4E_GSE41271_endpoint_overlap.png"), TRUE)
render_overlap("TCGA_LUAD", "G", file.path(figures, "Figure4G_TCGA_endpoint_overlap.pdf"), FALSE)
render_overlap("TCGA_LUAD", "G", file.path(figures, "Figure4G_TCGA_endpoint_overlap.png"), TRUE)

choose_representatives <- function(cohort_name) {
  subset <- models[cohort == cohort_name & model == "multivariable"]
  x <- dcast(subset, TF ~ endpoint, value.var = "p_value")
  directions <- dcast(subset, TF ~ endpoint, value.var = "direction")
  setnames(directions, c("OS", "RFS"), c("OS_direction", "RFS_direction"))
  x <- merge(x, directions, by = "TF")
  x[, `:=`(worst_p = pmax(OS, RFS),
            consistent_direction = OS_direction == RFS_direction)]
  answer <- character()
  for (direction_name in c("Favorable", "Adverse")) {
    candidate <- x[consistent_direction == TRUE & OS_direction == direction_name][order(worst_p, TF), TF]
    if (length(candidate)) answer <- c(answer, candidate[[1L]])
  }
  if (length(answer) < 2L) answer <- unique(c(answer, x[order(worst_p, TF), TF]))[1:2]
  answer
}

km_receipts <- list()
for (cohort_name in names(cohorts)) {
  representatives <- choose_representatives(cohort_name)
  panel <- ifelse(cohort_name == "GSE41271_LUAD", "F", "H")
  for (tf in representatives) {
    for (endpoint in c("OS", "RFS")) {
      dat <- endpoint_data(cohort_name, endpoint, tf)
      cut <- surv_cutpoint(dat, time = "time", event = "event", variables = "activity", minprop = 0.1)
      cut_value <- as.numeric(cut$cutpoint["activity", "cutpoint"])
      dat[, activity_group := factor(ifelse(activity <= cut_value, "Low", "High"), levels = c("Low", "High"))]
      fit <- survfit(Surv(time, event) ~ activity_group, data = dat)
      logrank <- survdiff(Surv(time, event) ~ activity_group, data = dat)
      logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)
      title_endpoint <- ifelse(cohort_name == "TCGA_LUAD" & endpoint == "RFS", "DFS", endpoint)
      plot <- ggsurvplot(
        fit, data = dat, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
        xlim = c(0, 120), break.time.by = 30, xscale = "m_y",
        palette = c("#1B9E77", "#E7298A"),
        legend.title = sprintf("%s activity", tf), legend.labs = c("Low", "High"),
        title = sprintf("%s  %s %s", panel, tf, title_endpoint),
        xlab = "Time (years)", ylab = sprintf("%s probability", title_endpoint),
        risk.table.height = 0.23, ggtheme = theme_classic(base_size = 8.2)
      )
      filename <- sprintf("Figure4%s_%s_%s_%s_KM", panel, cohort_name, tf, endpoint)
      pdf(file.path(figures, paste0(filename, ".pdf")), width = 5.2, height = 4.8, useDingbats = FALSE)
      print(plot)
      dev.off()
      png(file.path(figures, paste0(filename, ".png")), width = 1144, height = 1056, res = 220)
      print(plot)
      dev.off()
      km_receipts[[length(km_receipts) + 1L]] <- data.table(
        cohort = cohort_name, panel = panel, TF = tf, endpoint = endpoint,
        cutpoint = cut_value, n = nrow(dat), events = sum(dat$event),
        low_n = sum(dat$activity_group == "Low"), high_n = sum(dat$activity_group == "High"),
        logrank_p = logrank_p, filename = filename
      )
    }
  }
}
km_table <- rbindlist(km_receipts)
fwrite(km_table, file.path(tables, "Figure4FH_KM_receipts.tsv"), sep = "\t")

endpoint_audit <- rbindlist(lapply(names(survival_inputs), function(name) {
  object <- survival_inputs[[name]]
  data.table(cohort = object$cohort, endpoint = object$endpoint,
             n_complete = nrow(object$clinical), events_10year = sum(object$clinical$event),
             age_available = sum(is.finite(object$clinical$age), na.rm = TRUE),
             sex_levels = nlevels(object$clinical$sex_model),
             stage_levels = nlevels(object$clinical$stage_model),
             smoking_levels = nlevels(object$clinical$smoking_model))
}))
fwrite(endpoint_audit, file.path(tables, "SupplementaryFigure7_endpoint_complete_case_audit.tsv"), sep = "\t")

receipt <- list(
  status = "passed", HC_TFs = length(hc),
  endpoint_audit = endpoint_audit,
  nominal_multivariable_counts = models[model == "multivariable", .(
    significant = sum(nominal_significant), FDR_significant = sum(FDR_significant)
  ), by = .(cohort, endpoint)],
  overlap = overlap_rows,
  KM_representatives = unique(km_table[, .(cohort, TF)]),
  cox_exposure = "activity above versus at-or-below cohort/endpoint median",
  censoring_months = 120,
  multivariable_covariates = c("age", "sex", "pathologic stage", "smoking history"),
  GSE41271_primary_activity_definition = "TCGA regulon projected to GSE41271 microarray expression; all 31 HC-TFs",
  GSE41271_independent_network_sensitivity_TFs = length(gse_independent_available),
  GSE41271_independent_network_missing_TFs = gse_independent_missing,
  package_versions = list(R = as.character(getRversion()), survival = as.character(packageVersion("survival")),
                          survminer = as.character(packageVersion("survminer")))
)
write_json(receipt, file.path(audit, "figure4_main_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

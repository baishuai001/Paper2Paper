#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(maxstat)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: run_figure4_nominal_logrank_hc45.R FIGURE1_ROOT FIGURE2_ROOT",
    "FIGURE4_SOURCE_ROOT HC45_PARAMETER_AUDIT_TSV OUTPUT_ROOT"
  ))
}
fig1 <- normalizePath(args[[1L]], mustWork = TRUE)
fig2 <- normalizePath(args[[2L]], mustWork = TRUE)
fig4_source <- normalizePath(args[[3L]], mustWork = TRUE)
multiv_audit_path <- normalizePath(args[[4L]], mustWork = TRUE)
out <- normalizePath(args[[5L]], mustWork = FALSE)
tables <- file.path(out, "results", "tables")
official_inputs <- file.path(out, "results", "official_inputs")
audit <- file.path(out, "audit")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(official_inputs, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

set.seed(20260822)
HORIZON_MONTHS <- 120
MINPROP <- 0.20

read_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  ids <- x[[1L]]
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

read_json_pages <- function(directory, study, kind) {
  paths <- sort(Sys.glob(file.path(directory, sprintf("%s_%s_page*.json", study, kind))))
  if (!length(paths)) stop(sprintf("No cBioPortal pages for %s %s", study, kind))
  rbindlist(lapply(paths, function(path) as.data.table(fromJSON(path, flatten = TRUE))), fill = TRUE)
}

pivot_clinical <- function(tab) {
  tab <- tab[!is.na(patientId) & !is.na(clinicalAttributeId)]
  dcast(
    tab, patientId ~ clinicalAttributeId, value.var = "value",
    fun.aggregate = function(x) if (length(x)) x[[1L]] else NA_character_,
    fill = NA_character_
  )
}

parse_event <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fifelse(
    grepl("^1|DECEASED|RECUR|PROGRESS", x), 1L,
    fifelse(grepl("^0|LIVING|DISEASE.?FREE", x), 0L, NA_integer_)
  )
}

gate <- fread(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
hc <- sort(unique(gate[as.logical(HC_TF_promoter_activity_definition), TF]))
if (length(hc) != 45L) stop(sprintf("Expected 45 HC-TFs, found %d", length(hc)))

tcga_activity <- read_matrix(file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz"))
tcga_manifest <- fread(file.path(fig1, "data", "processed", "tcga_manifest.tsv"))
tcga_annotation <- fread(file.path(
  fig1, "data", "processed", "annotations", "tcga_figure1_annotations.tsv"
))
tcga_ids <- tcga_manifest[group == "LUAD", sample_id]
if (length(tcga_ids) != 516L || length(setdiff(hc, rownames(tcga_activity)))) {
  stop("TCGA activity does not contain the complete LUAD HC45 input")
}
cbio <- pivot_clinical(read_json_pages(
  file.path(fig1, "data", "raw", "annotations", "cbioportal"),
  "luad_tcga_pan_can_atlas_2018", "clinical_patient"
))
tcga_clinical <- tcga_annotation[group == "LUAD"]
tcga_clinical <- tcga_clinical[match(tcga_ids, sample_id)]
tcga_clinical <- merge(
  tcga_clinical, cbio, by.x = "patient_id", by.y = "patientId",
  all.x = TRUE, sort = FALSE, suffixes = c("", ".cbio")
)
tcga_clinical <- tcga_clinical[match(tcga_ids, sample_id)]
tcga_clinical[, `:=`(
  cohort = "TCGA_LUAD",
  OS_time = suppressWarnings(as.numeric(OS_MONTHS)),
  OS_event = parse_event(OS_STATUS),
  RFS_time = suppressWarnings(as.numeric(DFS_MONTHS)),
  RFS_event = parse_event(DFS_STATUS)
)]

gse_activity <- read_matrix(file.path(
  fig4_source, "results", "gse41271_tcga_projection",
  "GSE41271_TCGA_PROJECTED_viper_activity.tsv.gz"
))
gse_manifest <- fread(file.path(fig1, "data", "processed", "gse41271_manifest.tsv"))
gse_ids <- gse_manifest[group == "LUAD", sample_id]
if (length(gse_ids) != 183L || length(setdiff(hc, rownames(gse_activity)))) {
  stop("GSE41271 projection does not contain the complete LUAD HC45 input")
}
gse_clinical <- gse_manifest[group == "LUAD"]
gse_clinical <- gse_clinical[match(gse_ids, sample_id)]
gse_clinical[, `:=`(
  cohort = "GSE41271_LUAD",
  OS_time = suppressWarnings(as.numeric(os_days) / 30.4375),
  OS_event = as.integer(os_event),
  RFS_time = suppressWarnings(as.numeric(rfs_days) / 30.4375),
  RFS_event = as.integer(rfs_event)
)]

cohorts <- list(
  GSE41271_LUAD = list(
    clinical = gse_clinical,
    activity = gse_activity[, gse_ids, drop = FALSE]
  ),
  TCGA_LUAD = list(
    clinical = tcga_clinical,
    activity = tcga_activity[, tcga_ids, drop = FALSE]
  )
)

endpoint_data <- function(cohort_name, endpoint, tf) {
  object <- cohorts[[cohort_name]]
  dat <- copy(object$clinical)
  dat[, activity := as.numeric(object$activity[tf, sample_id])]
  dat[, `:=`(
    time_raw = suppressWarnings(as.numeric(get(paste0(endpoint, "_time")))),
    event_raw = suppressWarnings(as.integer(get(paste0(endpoint, "_event"))))
  )]
  dat <- dat[
    is.finite(time_raw) & time_raw > 0 & event_raw %in% c(0L, 1L) & is.finite(activity)
  ]
  dat[, `:=`(
    time = pmin(time_raw, HORIZON_MONTHS),
    event = as.integer(event_raw == 1L & time_raw <= HORIZON_MONTHS)
  )]
  dat
}

fit_one <- function(cohort_name, endpoint, tf) {
  dat <- endpoint_data(cohort_name, endpoint, tf)
  fit_max <- maxstat.test(
    Surv(time, event) ~ activity,
    data = as.data.frame(dat),
    smethod = "LogRank", pmethod = "Lau92",
    minprop = MINPROP, maxprop = 1 - MINPROP
  )
  cutpoint <- as.numeric(fit_max$estimate)
  dat[, activity_group := factor(
    fifelse(activity <= cutpoint, "Low", "High"), levels = c("Low", "High")
  )]
  lr <- survdiff(Surv(time, event) ~ activity_group, data = dat)
  logrank_p <- pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  cox <- coxph(Surv(time, event) ~ activity_group, data = dat, ties = "efron", x = TRUE)
  s <- summary(cox)
  ci <- s$conf.int["activity_groupHigh", ]
  ph <- tryCatch(cox.zph(cox), error = function(e) NULL)
  ph_p <- if (is.null(ph)) NA_real_ else as.numeric(ph$table["activity_group", "p"])
  data.table(
    cohort = cohort_name,
    endpoint = endpoint,
    TF = tf,
    n = nrow(dat),
    events = sum(dat$event),
    cutpoint = cutpoint,
    cutpoint_percentile = mean(dat$activity <= cutpoint),
    low_n = sum(dat$activity_group == "Low"),
    high_n = sum(dat$activity_group == "High"),
    HR = unname(ci[["exp(coef)"]]),
    lower95 = unname(ci[["lower .95"]]),
    upper95 = unname(ci[["upper .95"]]),
    logrank_p = logrank_p,
    maxstat_Lau92_adjusted_p = as.numeric(fit_max$p.value),
    ph_group_p = ph_p
  )
}

rows <- list()
for (cohort_name in names(cohorts)) {
  for (endpoint in c("OS", "RFS")) {
    for (tf in hc) {
      rows[[length(rows) + 1L]] <- fit_one(cohort_name, endpoint, tf)
    }
  }
}
results <- rbindlist(rows)
results[, `:=`(
  logrank_BH_FDR = p.adjust(logrank_p, method = "BH"),
  maxstat_Lau92_BH_FDR = p.adjust(maxstat_Lau92_adjusted_p, method = "BH")
), by = .(cohort, endpoint)]
results[, `:=`(
  direction = fifelse(HR < 1, "Favorable", "Adverse"),
  nominal_unadjusted_logrank_significant = logrank_p <= 0.05,
  logrank_FDR_significant = logrank_BH_FDR <= 0.05,
  maxstat_adjusted_significant = maxstat_Lau92_adjusted_p <= 0.05,
  maxstat_FDR_significant = maxstat_Lau92_BH_FDR <= 0.05,
  PH_group_violation = is.finite(ph_group_p) & ph_group_p <= 0.05,
  horizon_months = HORIZON_MONTHS,
  minimum_group_proportion = MINPROP,
  primary_display_rule = "unadjusted log-rank P<=0.05 after max-selected cutpoint"
)]
fwrite(results, file.path(tables, "Figure4_logrank_all_HC45.tsv"), sep = "\t")
fwrite(
  results[(nominal_unadjusted_logrank_significant)],
  file.path(tables, "Figure4_nominal_unadjusted_logrank_P05.tsv"), sep = "\t"
)

overlap <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  rbindlist(lapply(c("Favorable", "Adverse"), function(direction_name) {
    selected <- results[
      cohort == cohort_name & direction == direction_name &
        (nominal_unadjusted_logrank_significant)
    ]
    os <- selected[endpoint == "OS", TF]
    rfs <- selected[endpoint == "RFS", TF]
    data.table(
      cohort = cohort_name,
      direction = direction_name,
      OS_only = length(setdiff(os, rfs)),
      overlap = length(intersect(os, rfs)),
      RFS_only = length(setdiff(rfs, os)),
      overlap_TFs = paste(sort(intersect(os, rfs)), collapse = ";")
    )
  }))
}))
fwrite(overlap, file.path(tables, "Figure4EG_endpoint_overlap_logrank.tsv"), sep = "\t")

choose_representatives <- function(cohort_name) {
  z <- dcast(
    results[cohort == cohort_name],
    TF + direction ~ endpoint,
    value.var = "logrank_p"
  )
  # Figure 4F/H illustrate TFs in the OS-RFS/DFS overlap.  A TF that is
  # nominally significant in only one endpoint must not be promoted to the
  # overlap example merely because it has the smallest worst P value.
  z <- z[is.finite(OS) & is.finite(RFS) & OS <= 0.05 & RFS <= 0.05]
  z[, worst_p := pmax(OS, RFS)]
  answer <- character()
  for (direction_name in c("Favorable", "Adverse")) {
    candidate <- z[direction == direction_name][order(worst_p, TF), TF]
    if (length(candidate)) answer <- c(answer, candidate[[1L]])
  }
  unique(answer)
}

km_objects <- list()
km_receipts <- list()
for (cohort_name in names(cohorts)) {
  representatives <- choose_representatives(cohort_name)
  for (tf in representatives) {
    for (endpoint_name in c("OS", "RFS")) {
      rec <- results[cohort == cohort_name & endpoint == endpoint_name & TF == tf]
      dat <- endpoint_data(cohort_name, endpoint_name, tf)
      dat[, activity_group := factor(
        fifelse(activity <= rec$cutpoint, "Low", "High"), levels = c("Low", "High")
      )]
      key <- paste(cohort_name, tf, endpoint_name, sep = "__")
      km_objects[[key]] <- list(data = dat, receipt = rec)
      km_receipts[[length(km_receipts) + 1L]] <- rec[, .(
        cohort, endpoint, TF, cutpoint, n, events, low_n, high_n,
        logrank_p, maxstat_Lau92_adjusted_p, logrank_BH_FDR,
        primary_display_rule
      )]
    }
  }
}
km_receipts <- rbindlist(km_receipts)
saveRDS(km_objects, file.path(out, "results", "Figure4_km_objects.rds"), compress = "xz")
fwrite(km_receipts, file.path(tables, "Figure4FH_KM_receipts.tsv"), sep = "\t")

write_official <- function(x, path, p_column, p_official_name) {
  z <- x[, .(
    TR = TF,
    Hazard_Ratio = HR,
    HR_Lower_95 = lower95,
    HR_Upper_95 = upper95,
    p_value = get(p_column),
    Significant = fifelse(get(p_column) <= 0.05, "Yes", "No")
  )]
  z[, LogRank_PVal := p_value]
  z[, CoxMultiPVal := p_value]
  z[, p_value := NULL]
  setcolorder(z, c(
    "TR", "Hazard_Ratio", "HR_Lower_95", "HR_Upper_95",
    "LogRank_PVal", "CoxMultiPVal", "Significant"
  ))
  fwrite(z, path, sep = "\t")
}

for (cohort_name in names(cohorts)) {
  for (endpoint_name in c("OS", "RFS")) {
    z <- results[cohort == cohort_name & endpoint == endpoint_name]
    write_official(
      z,
      file.path(official_inputs, sprintf("%s_%s_CoxUniVariate.tsv", cohort_name, endpoint_name)),
      "logrank_p", "LogRank_PVal"
    )
  }
}

multiv <- fread(multiv_audit_path)
multiv <- multiv[
  TF %chin% hc & horizon == "H120" & exposure == "median_binary" & model == "multivariable"
]
if (nrow(multiv) != 4L * length(hc)) {
  stop(sprintf("Expected %d HC45 multivariable sensitivity rows; found %d", 4L * length(hc), nrow(multiv)))
}
fwrite(multiv, file.path(tables, "SupplementaryFigure6_multivariable_HC45_sensitivity.tsv"), sep = "\t")
for (cohort_name in names(cohorts)) {
  for (endpoint_name in c("OS", "RFS")) {
    z <- multiv[cohort == cohort_name & endpoint == endpoint_name]
    write_official(
      z,
      file.path(official_inputs, sprintf("%s_%s_CoxMultiVariate.tsv", cohort_name, endpoint_name)),
      "p_value", "CoxMultiPVal"
    )
  }
}

counts <- results[, .(
  HC_TFs = .N,
  nominal_unadjusted_logrank_P05 = sum(nominal_unadjusted_logrank_significant),
  logrank_BH_FDR05 = sum(logrank_FDR_significant),
  maxstat_Lau92_P05 = sum(maxstat_adjusted_significant),
  maxstat_Lau92_BH_FDR05 = sum(maxstat_FDR_significant),
  PH_group_violations = sum(PH_group_violation)
), by = .(cohort, endpoint)]
fwrite(counts, file.path(tables, "Figure4_logrank_significance_counts.tsv"), sep = "\t")

receipt <- list(
  status = "PASS",
  HC_TFs = length(hc),
  horizon_months = HORIZON_MONTHS,
  minimum_group_proportion = MINPROP,
  cutpoint_estimator = "maxstat LogRank; cutpoint estimate from Lau92 call",
  primary_decision = "naive unadjusted selected-cutpoint log-rank P <= 0.05",
  warning = paste(
    "The primary p-value does not correct for cutpoint search or 45 TF tests.",
    "It is exploratory and must not be described as a validated prognostic classifier."
  ),
  counts = counts,
  representatives = unique(km_receipts[, .(cohort, TF)]),
  adjusted_columns_retained = c(
    "logrank_BH_FDR", "maxstat_Lau92_adjusted_p", "maxstat_Lau92_BH_FDR"
  ),
  multivariable_role = "Supplementary Figure 6 sensitivity only",
  packages = list(
    R = as.character(getRversion()),
    survival = as.character(packageVersion("survival")),
    maxstat = as.character(packageVersion("maxstat"))
  )
)
write_json(receipt, file.path(audit, "figure4_nominal_logrank_hc45_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

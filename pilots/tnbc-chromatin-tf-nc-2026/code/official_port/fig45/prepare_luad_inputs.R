#!/usr/bin/env Rscript

# Data-only adapter for the strict TNBC CodeOcean Figure 4/5 plotting port.
# It does not contain plotting code and does not refit any statistical model.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: prepare_luad_inputs.R FIGURE1_ROOT FIGURE4_ROOT FIGURE5_ROOT",
    "PUBLICATION_REBUILD_ROOT OUTPUT_DIR"
  ), call. = FALSE)
}

fig1 <- normalizePath(args[[1L]], mustWork = TRUE)
fig4 <- normalizePath(args[[2L]], mustWork = TRUE)
fig5 <- normalizePath(args[[3L]], mustWork = TRUE)
pub <- normalizePath(args[[4L]], mustWork = TRUE)
out <- normalizePath(args[[5L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figure4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figure5", "by_dataset"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "supplementary5"), recursive = TRUE, showWarnings = FALSE)

must_read <- function(path, required = character()) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  x <- fread(path, na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(path, " missing: ", paste(missing, collapse = ","))
  x
}

read_matrix <- function(path) {
  tab <- must_read(path)
  if (ncol(tab) < 2L) stop("Matrix has fewer than two columns: ", path)
  ids <- tab[[1L]]
  tab[[1L]] <- NULL
  mat <- as.matrix(tab)
  suppressWarnings(storage.mode(mat) <- "double")
  rownames(mat) <- ids
  mat
}

write_tab <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

input_rows <- list()
record_input <- function(role, path, rows = NA_integer_) {
  input_rows[[length(input_rows) + 1L]] <<- data.table(
    role = role,
    path = normalizePath(path, mustWork = TRUE),
    bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path)),
    rows = rows
  )
}

# ---------------------------------------------------------------------------
# Figure 4 and Supplementary Figure 6: adapt frozen Cox rows to the exact
# column names consumed by the upstream plotting scripts.
# ---------------------------------------------------------------------------

cox_path <- file.path(fig4, "results", "tables", "Figure4_all_Cox_models.tsv")
cox <- must_read(
  cox_path,
  c("cohort", "endpoint", "model", "TF", "HR", "lower95", "upper95",
    "p_value", "nominal_significant")
)
record_input("Figure4 frozen Cox models", cox_path, nrow(cox))

expected_cohorts <- c("GSE41271_LUAD", "TCGA_LUAD")
expected_endpoints <- c("OS", "RFS")
expected_models <- c("univariable", "multivariable")
stopifnot(setequal(unique(cox$cohort), expected_cohorts))
stopifnot(setequal(unique(cox$endpoint), expected_endpoints))
stopifnot(setequal(unique(cox$model), expected_models))

for (cohort_name in expected_cohorts) {
  for (endpoint_name in expected_endpoints) {
    for (model_name in expected_models) {
      z <- copy(cox[
        cohort == cohort_name & endpoint == endpoint_name & model == model_name
      ])
      if (!nrow(z)) stop("No Cox rows for ", cohort_name, "/", endpoint_name, "/", model_name)
      official <- z[, .(
        TR = TF,
        Hazard_Ratio = HR,
        HR_Lower_95 = lower95,
        HR_Upper_95 = upper95,
        LogRank_PVal = p_value,
        CoxMultiPVal = p_value,
        Significant = fifelse(nominal_significant, "Yes", "No")
      )]
      suffix <- if (model_name == "univariable") "CoxUniVariate" else "CoxMultiVariate"
      write_tab(
        official,
        file.path(out, "figure4", sprintf(
          "%s_%s_%s.tsv", cohort_name, endpoint_name, suffix
        ))
      )
    }
  }
}
write_tab(cox, file.path(out, "figure4", "models_long.tsv"))

overlap_path <- file.path(fig4, "results", "tables", "Figure4EG_endpoint_overlap.tsv")
overlap <- must_read(overlap_path, c("cohort", "direction", "OS_only", "overlap", "RFS_only"))
record_input("Figure4 endpoint overlap", overlap_path, nrow(overlap))
write_tab(overlap, file.path(out, "figure4", "endpoint_overlap.tsv"))

km_receipt_path <- file.path(fig4, "results", "tables", "Figure4FH_KM_receipts.tsv")
km_receipts <- must_read(
  km_receipt_path,
  c("cohort", "TF", "endpoint", "cutpoint", "n", "events", "low_n", "high_n", "logrank_p")
)
record_input("Figure4 KM receipts", km_receipt_path, nrow(km_receipts))

survival_path <- file.path(fig4, "results", "figure4_survival_inputs.rds")
if (!file.exists(survival_path)) stop("Missing frozen survival RDS: ", survival_path)
survival_inputs <- readRDS(survival_path)
record_input("Figure4 survival inputs", survival_path)

activity_paths <- c(
  GSE41271_LUAD = file.path(
    fig4, "results", "gse41271_tcga_projection",
    "GSE41271_TCGA_PROJECTED_viper_activity.tsv.gz"
  ),
  TCGA_LUAD = file.path(fig1, "results", "tcga", "TCGA_viper_activity.tsv.gz")
)
activity_mats <- lapply(activity_paths, read_matrix)
for (nm in names(activity_paths)) record_input(paste("KM activity", nm), activity_paths[[nm]])

km_objects <- list()
for (i in seq_len(nrow(km_receipts))) {
  rec <- km_receipts[i]
  key <- paste(rec$cohort, rec$endpoint, sep = "__")
  base <- survival_inputs[[key]]
  mat <- activity_mats[[rec$cohort]]
  if (is.null(base) || !rec$TF %in% rownames(mat)) stop("Cannot reconstruct KM receipt ", i)
  dat <- copy(base$clinical)
  dat <- dat[sample_id %in% colnames(mat)]
  dat[, activity := as.numeric(mat[rec$TF, sample_id])]
  dat <- dat[is.finite(time) & time > 0 & event %in% c(0L, 1L) & is.finite(activity)]
  dat[, activity_group := factor(
    fifelse(activity <= rec$cutpoint, "Low", "High"), levels = c("Low", "High")
  )]
  observed <- c(
    n = nrow(dat), events = sum(dat$event),
    low_n = sum(dat$activity_group == "Low"), high_n = sum(dat$activity_group == "High")
  )
  expected <- unlist(rec[, .(n, events, low_n, high_n)])
  if (!identical(as.numeric(observed), as.numeric(expected))) {
    stop("KM receipt count mismatch: ", rec$cohort, "/", rec$TF, "/", rec$endpoint)
  }
  km_objects[[paste(rec$cohort, rec$TF, rec$endpoint, sep = "__")]] <- list(
    receipt = rec,
    data = dat
  )
}
saveRDS(km_objects, file.path(out, "figure4", "km_objects.rds"), compress = "xz")
write_tab(km_receipts, file.path(out, "figure4", "km_receipts.tsv"))

# Cross-cohort three-state prognosis classification. This is the frozen LUAD
# rule already used by the analysis; it replaces only TNBC's hard-coded counts.
hc_tfs <- sort(unique(cox$TF))
support <- cox[
  model == "multivariable" & p_value <= 0.05,
  .(
    cohort_n = uniqueN(cohort),
    direction_n = uniqueN(ifelse(HR < 1, "Favorable", "Adverse")),
    direction = if (uniqueN(ifelse(HR < 1, "Favorable", "Adverse")) == 1L) {
      first(ifelse(HR < 1, "Favorable", "Adverse"))
    } else NA_character_
  ),
  by = .(TF, endpoint)
]
states <- CJ(TF = hc_tfs, endpoint = expected_endpoints, unique = TRUE)
states <- merge(states, support, by = c("TF", "endpoint"), all.x = TRUE)
states[, prognosis := fifelse(
  cohort_n == 2L & direction_n == 1L & direction == "Favorable", "Good prognosis",
  fifelse(cohort_n == 2L & direction_n == 1L & direction == "Adverse", "Poor prognosis", "Non-Prognostic")
)]
prognosis <- dcast(states, TF ~ endpoint, value.var = "prognosis")
setnames(prognosis, c("TF", "OS", "RFS"), c("HC_TR", "Overall_OS", "Overall_RFS"))
prognosis[, Overall_RFS := fifelse(
  Overall_RFS == "Good prognosis", "Recurrence-Free",
  fifelse(Overall_RFS == "Poor prognosis", "Recurrence", "Non-Prognostic")
)]
write_tab(prognosis, file.path(out, "figure5", "prognostic_metadata.tsv"))

conf_levels_os <- c("Good prognosis", "Non-Prognostic", "Poor prognosis")
conf_levels_rfs <- c("Recurrence-Free", "Non-Prognostic", "Recurrence")
conf <- CJ(OS = conf_levels_os, RFS = conf_levels_rfs, unique = TRUE)
counts <- prognosis[, .N, by = .(OS = Overall_OS, RFS = Overall_RFS)]
conf <- merge(conf, counts, by = c("OS", "RFS"), all.x = TRUE)
conf[is.na(N), N := 0L]
setnames(conf, "N", "n")
write_tab(conf, file.path(out, "figure5", "figure5e_confusion.tsv"))

# ---------------------------------------------------------------------------
# Figure 5A-F and Supplementary Figure 7A-C.
# ---------------------------------------------------------------------------

cell_path <- file.path(fig5, "results", "tables", "Figure5ABC_cellline_concordance_all.tsv.gz")
cell <- must_read(
  cell_path,
  c("dataset", "canonical_drug", "drug_name", "TF", "n", "c_index", "p_value", "FDR", "direction")
)
record_input("Figure5 all cell-line concordance tests", cell_path, nrow(cell))
stopifnot(setequal(unique(cell$dataset), c("GDSC2", "CTRPv2", "PRISM")))
cell[, Category := fifelse(
  c_index > 0.5, "Positive_Correlated_CI_Significant",
  "Negative_Correlated_CI_Significant"
)]
cell[, drug := canonical_drug]
cell[, TR := TF]
cell[, Drug_TR := paste(drug, TR, Category, sep = "_")]
cell[, Cohort := dataset]
for (dataset_name in c("PRISM", "GDSC2", "CTRPv2")) {
  write_tab(
    cell[dataset == dataset_name],
    file.path(out, "figure5", "by_dataset", paste0(dataset_name, "_drug_TR_concordance_results.tsv"))
  )
}

significant <- cell[FDR <= 0.05]
membership <- unique(significant[, .(Drug_TR, dataset)])
pair_n <- membership[, .(dataset_n = uniqueN(dataset)), by = Drug_TR]
replicated_ids <- pair_n[dataset_n >= 2L, Drug_TR]
replicated <- significant[Drug_TR %in% replicated_ids]
replicated[, display_drug := fifelse(is.na(drug_name) | drug_name == "", canonical_drug, drug_name)]

drug_map_path <- file.path(fig5, "results", "tables", "Figure5_drug_identity_map.tsv")
drug_map <- must_read(drug_map_path, c("dataset", "canonical_drug", "drug_name", "pathway_or_moa"))
record_input("Figure5 drug identity map", drug_map_path, nrow(drug_map))
drug_anno <- unique(drug_map[, .(
  dataset, canonical_drug, pathway_or_moa,
  mapped_drug_name = drug_name
)])
replicated <- merge(replicated, drug_anno, by = c("dataset", "canonical_drug"), all.x = TRUE)

pathway_label <- function(x) {
  y <- tolower(fifelse(is.na(x), "", x))
  fcase(
    grepl("pi3k|mtor|akt", y), "PI3K/Akt/mTOR",
    grepl("mapk|mek|erk", y), "MAPK/ERK",
    grepl("egfr|erbb", y), "EGFR/ERBB2",
    grepl("antimit|microtub|tubulin|taxane", y), "Antimitotic",
    grepl("kinase|cdk|fgfr|vegfr|aurora", y), "Kinase Inhibitor",
    grepl("anti[- ]?metabol|pyrimidine|thymidylate", y), "Anti-metabolite",
    grepl("bromodomain|brd[0-9]|bet inhibitor", y), "Bromodomain",
    grepl("telomerase", y), "Telomerase inhibitor",
    grepl("dna replication|dna synthesis", y), "DNA replication",
    grepl("apoptosis|bcl|caspase", y), "Apoptosis",
    grepl("nedd[- ]?8|nae inhibitor", y), "NAE inhibitor",
    grepl("survivin", y), "Survivin inhibitor",
    grepl("dna damage", y), "inducer of DNA damage",
    grepl("hdac|histone acetyl", y), "Small Molecule Inhibitor",
    default = "Other"
  )
}
replicated[, Label := pathway_label(pathway_or_moa)]
# A canonical drug must occupy one facet across resources.  When one source
# has a specific mechanism and another says "Other", retain the specific
# official Figure 5 pathway class rather than splitting the same drug.
label_priority <- c(
  "PI3K/Akt/mTOR", "MAPK/ERK", "EGFR/ERBB2", "Antimitotic",
  "Kinase Inhibitor", "Anti-metabolite", "Bromodomain",
  "Telomerase inhibitor", "DNA replication", "Apoptosis", "NAE inhibitor",
  "Survivin inhibitor", "inducer of DNA damage", "Small Molecule Inhibitor",
  "Other"
)
canonical_label <- replicated[, .(
  canonical_Label = Label[which.min(match(Label, label_priority))]
), by = canonical_drug]
replicated <- merge(replicated, canonical_label, by = "canonical_drug", all.x = TRUE)
replicated[, Label := canonical_Label]
replicated[, canonical_Label := NULL]
replicated[, label := paste(display_drug, TR, Cohort, sep = "_")]
replicated <- merge(
  replicated,
  prognosis[, .(TR = HC_TR, OS = Overall_OS, RFS = Overall_RFS)],
  by = "TR", all.x = TRUE
)
write_tab(replicated, file.path(out, "figure5", "combined_replicated_ci.tsv"))
write_tab(pair_n, file.path(out, "figure5", "drug_tf_membership_counts.tsv"))

drug_membership <- unique(cell[, .(dataset, canonical_drug)])
write_tab(drug_membership, file.path(out, "figure5", "eligible_drug_membership.tsv"))

# Freeze four LUAD pairs into the four published Supplementary Figure 7B
# call/save slots.  The renderer executes those official blocks; this table
# only verifies that each substituted identifier is supported in the frozen
# replicated-pair result and in exactly the declared cohorts.
example_pairs <- data.table(
  case_block = paste0("supp7b_case", 1:4),
  canonical_drug = c("CID:46191454", "CID:3385", "CID:11977753", "CID:5311"),
  TR = c("WWC2", "ZNF254", "ZNF254", "ZNF254"),
  display_drug = c("Pha-793887", "5-Fluorouracil", "Dactolisib", "Vorinostat"),
  expected_datasets = c("CTRPv2;PRISM", "CTRPv2;GDSC2", "GDSC2;PRISM", "GDSC2;PRISM")
)
for (i in seq_len(nrow(example_pairs))) {
  p <- example_pairs[i]
  hit <- replicated[canonical_drug == p$canonical_drug & TR == p$TR]
  expected_datasets <- strsplit(p$expected_datasets, ";", fixed = TRUE)[[1L]]
  if (!nrow(hit) || !setequal(unique(hit$dataset), expected_datasets)) {
    stop("Frozen Supplementary Figure 7B case no longer matches replicated evidence: ", p$case_block)
  }
  if (!identical(unique(hit$Category), "Negative_Correlated_CI_Significant")) {
    stop("Frozen Supplementary Figure 7B direction changed: ", p$case_block)
  }
}
write_tab(example_pairs, file.path(out, "figure5", "supp7b_example_pairs.tsv"))

cell_inputs_path <- file.path(fig5, "data", "processed", "figure5_cellline_inputs.rds")
if (!file.exists(cell_inputs_path)) stop("Missing cell-line input RDS: ", cell_inputs_path)
cell_inputs <- readRDS(cell_inputs_path)
record_input("Figure5 cell-line paired values", cell_inputs_path)
example_rows <- list()
for (i in seq_len(nrow(example_pairs))) {
  p <- example_pairs[i]
  for (dataset_name in replicated[canonical_drug == p$canonical_drug & TR == p$TR, unique(dataset)]) {
    obj <- cell_inputs$datasets[[dataset_name]]
    if (is.null(obj)) next
    z <- copy(obj$paired[canonical_drug == p$canonical_drug])
    z <- z[depmap_id %in% colnames(cell_inputs$activity)]
    if (!nrow(z) || !p$TR %in% rownames(cell_inputs$activity)) next
    z[, NES := as.numeric(cell_inputs$activity[p$TR, depmap_id])]
    z[, `:=`(
      TF = p$TR, dataset = dataset_name,
      display_drug = p$display_drug
    )]
    example_rows[[length(example_rows) + 1L]] <- z[, .(
      canonical_drug, display_drug, TF, dataset,
      sample_id = depmap_id, AAC = aac, NES
    )]
  }
}
if (length(example_rows)) {
  write_tab(rbindlist(example_rows), file.path(out, "figure5", "supp7b_example_values.tsv"))
} else {
  # Keep a parseable, typed zero-row table.  The renderer records a
  # SUPPORTED_ZERO stratum instead of attempting an unofficial substitute.
  write_tab(data.table(
    canonical_drug = character(), display_drug = character(), TF = character(),
    dataset = character(), sample_id = character(), AAC = numeric(), NES = numeric()
  ), file.path(out, "figure5", "supp7b_example_values.tsv"))
}

# Strict PDX boundary. Both Figure 5G and Figure 5H require at least one
# cell-line-replicated pair with matched PDX response. The frozen screen has 0.
pdx_path <- file.path(fig5, "results", "tables", "SupplementaryFigure9B_PDXE_all_concordance_tests.tsv")
pdx <- must_read(pdx_path, c("TF", "drug", "FDR", "cellline_replicated_same_direction"))
record_input("Figure5 complete PDX audit screen", pdx_path, nrow(pdx))
if (nrow(pdx) != 423L) stop("Frozen PDX audit screen changed from 423 rows; re-audit Figure 5G/H")
pdx_entry_n <- pdx[cellline_replicated_same_direction == TRUE, .N]
validated_path <- file.path(fig5, "results", "tables", "Figure5GH_PDX_validated_pairs.tsv")
validated <- must_read(validated_path)
record_input("Figure5 frozen PDX validated pairs", validated_path, nrow(validated))
if (pdx_entry_n != 0L || nrow(validated) != 0L) {
  stop("Frozen PDX boundary changed; re-audit Figure 5G/H before rendering")
}
unsupported <- data.table(
  panel = c("Figure5G", "Figure5H", "SupplementaryFigure7D"),
  status = "UNSUPPORTED",
  reason = c(
    "No cell-line-replicated drug-TF pair has matched PDX response (0/423 PDX tests enter validation)",
    "No validated LUAD PDX drug-TF pair; the official response waterfall has no eligible input",
    "No validated LUAD PDX drug-TF pair; the official PDX AAC example has no eligible input"
  ),
  replacement_plot_permitted = FALSE
)
write_tab(unsupported, file.path(out, "unsupported_panels.tsv"))

# ---------------------------------------------------------------------------
# Supplementary Figure 5A-B.
# ---------------------------------------------------------------------------

skew_path <- file.path(pub, "results", "tables", "SupplementaryFigure5A_cross_system_skewness.tsv")
skew <- must_read(skew_path, c("cohort", "cohort_label", "TF", "n", "skewness", "FDR"))
record_input("Supplementary Figure 5A skewness", skew_path, nrow(skew))

hetero_specs <- list(
  PDMR_LUAD = list(
    activity = file.path(fig1, "results", "pdmr", "PDMR_viper_activity.tsv.gz"),
    manifest = file.path(fig1, "data", "processed", "pdmr_manifest.tsv")
  ),
  DEPMAP_LUAD = list(
    activity = file.path(fig1, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
    manifest = file.path(fig1, "data", "processed", "depmap_manifest.tsv")
  )
)
mean_rows <- list()
for (cohort_name in names(hetero_specs)) {
  spec <- hetero_specs[[cohort_name]]
  mat <- read_matrix(spec$activity)
  man <- must_read(spec$manifest, c("sample_id", "group"))
  ids <- intersect(man[group == "LUAD", sample_id], colnames(mat))
  if (!length(ids)) stop("No LUAD samples for Supplementary Figure 5A: ", cohort_name)
  mean_rows[[length(mean_rows) + 1L]] <- data.table(
    cohort = cohort_name,
    TR_Source = rownames(mat),
    mean_NES = rowMeans(mat[, ids, drop = FALSE], na.rm = TRUE)
  )
  record_input(paste("Supplementary Figure 5A activity", cohort_name), spec$activity)
  record_input(paste("Supplementary Figure 5A manifest", cohort_name), spec$manifest, nrow(man))
}
supp5a <- merge(
  skew[cohort %in% names(hetero_specs), .(
    cohort, TR_Source = TF, n, skewness, skew_fdr = FDR,
    is_skewed = FDR < 0.05
  )],
  rbindlist(mean_rows), by = c("cohort", "TR_Source"), all = FALSE
)
supp5a[, Cohort := factor(
  fifelse(cohort == "PDMR_LUAD", "PDX", "CellLines"),
  levels = c("CellLines", "PDX")
)]
write_tab(supp5a, file.path(out, "supplementary5", "skewness_official.tsv"))

estimate_path <- file.path(pub, "results", "tables", "SupplementaryFigure5B_ESTIMATE_correlations.tsv")
estimate <- must_read(estimate_path, c("cohort", "cohort_label", "score", "TF", "pearson_r", "FDR"))
record_input("Supplementary Figure 5B ESTIMATE", estimate_path, nrow(estimate))
estimate_official <- estimate[, .(
  cohort = fifelse(cohort == "TCGA_LUAD", "TCGA-LUAD", "GSE41271-LUAD"),
  marker = score,
  TR = TF,
  r = pearson_r,
  fdr = FDR,
  sig = FDR < 0.05 & abs(pearson_r) > 0.4
)]
write_tab(estimate_official, file.path(out, "supplementary5", "estimate_official.tsv"))

input_manifest <- rbindlist(input_rows, fill = TRUE)
write_tab(input_manifest, file.path(out, "input_manifest.tsv"))
receipt <- list(
  status = "PASS",
  frozen_hc_tfs = length(hc_tfs),
  cox_rows = nrow(cox),
  cellline_tests = nrow(cell),
  replicated_drug_tf_ids = length(replicated_ids),
  pdx_tests = nrow(pdx),
  pdx_tests_entering_validation = pdx_entry_n,
  validated_pdx_pairs = nrow(validated),
  strict_main_supported = c(
    "Figure4A-H", "Figure5A-F", "SupplementaryFigure5A-B",
    "SupplementaryFigure6A-H", "SupplementaryFigure7A-C"
  ),
  unsupported = unsupported
)
write_json(receipt, file.path(out, "adapter_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, pretty = TRUE, auto_unbox = TRUE), "\n")

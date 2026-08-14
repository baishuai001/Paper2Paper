#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(TCGAbiolinks)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: build_figure1_annotations.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
raw <- file.path(root, "data", "raw")
processed <- file.path(root, "data", "processed")
out <- file.path(processed, "annotations")
audit <- file.path(root, "audit", "figure1_annotations")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

missing_string <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "N/A", "Not Reported", "not reported", "Unknown", "Unavailable", "--", "-")] <- NA_character_
  y
}

first_existing <- function(tab, candidates, default = NA_character_) {
  hit <- candidates[candidates %in% names(tab)]
  if (!length(hit)) return(rep(default, nrow(tab)))
  tab[[hit[[1]]]]
}

stage_group <- function(x) {
  x <- toupper(missing_string(x))
  fifelse(grepl("STAGE IV", x), "IV",
    fifelse(grepl("STAGE III", x), "III",
      fifelse(grepl("STAGE II", x), "II",
        fifelse(grepl("STAGE I", x), "I", NA_character_))))
}

age_group <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  fifelse(is.na(x), NA_character_,
    fifelse(x <= 60, "<=60", fifelse(x <= 70, "61-70", ">70")))
}

pack_year_group <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  fifelse(is.na(x), NA_character_,
    fifelse(x == 0, "0",
      fifelse(x <= 20, ">0-20", fifelse(x <= 40, ">20-40", ">40"))))
}

title_case_sex <- function(x) {
  x <- tolower(missing_string(x))
  fifelse(x == "female", "Female", fifelse(x == "male", "Male", NA_character_))
}

read_json_pages <- function(directory, study, kind) {
  paths <- sort(Sys.glob(file.path(directory, sprintf("%s_%s_page*.json", study, kind))))
  if (!length(paths)) stop(sprintf("No cBioPortal pages found for %s %s", study, kind))
  rbindlist(lapply(paths, function(path) as.data.table(fromJSON(path, flatten = TRUE))), fill = TRUE)
}

pivot_clinical <- function(tab, id_col) {
  tab <- tab[!is.na(get(id_col)) & !is.na(clinicalAttributeId)]
  dcast(
    tab,
    as.formula(sprintf("%s ~ clinicalAttributeId", id_col)),
    value.var = "value",
    fun.aggregate = function(x) if (length(x)) x[[1]] else NA_character_,
    fill = NA_character_
  )
}

mutation_summary <- function(tab, patients) {
  if (!nrow(tab)) return(data.table(patient_id = patients))
  gene_col <- if ("gene.hugoGeneSymbol" %in% names(tab)) "gene.hugoGeneSymbol" else "gene"
  tab[, gene_symbol := as.character(get(gene_col))]
  silent <- c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "lincRNA")
  tab <- tab[patientId %in% patients & gene_symbol != "" & !(mutationType %in% silent)]
  selected <- c("EGFR", "KRAS", "BRAF", "ERBB2", "MET", "PIK3CA", "ALK", "ROS1", "RET", "FGFR3", "DDR2",
                "TP53", "STK11", "KEAP1", "NFE2L2", "CDKN2A", "NF1", "PTEN")
  tab <- tab[gene_symbol %in% selected]
  collapse_unique <- function(x) paste(sort(unique(x)), collapse = ";")
  genes <- tab[, .(selected_mutations = collapse_unique(gene_symbol)), by = .(patient_id = patientId)]
  calls <- split(tab$gene_symbol, tab$patientId)
  oncogene_priority <- c("EGFR", "KRAS", "BRAF", "ERBB2", "MET", "PIK3CA", "ALK", "ROS1", "RET", "FGFR3", "DDR2")
  suppressor_set <- c("TP53", "STK11", "KEAP1", "NFE2L2", "CDKN2A", "NF1", "PTEN")
  choose_oncogene <- function(x) {
    hit <- oncogene_priority[oncogene_priority %in% x]
    if (length(hit)) hit[[1]] else "None detected"
  }
  choose_suppressor <- function(x) {
    hit <- unique(suppressor_set[suppressor_set %in% x])
    if (!length(hit)) "None detected" else if (length(hit) > 1) "Multiple" else hit[[1]]
  }
  summary <- data.table(
    patient_id = names(calls),
    oncogene_mutation = vapply(calls, choose_oncogene, character(1)),
    suppressor_pathway_mutation = vapply(calls, choose_suppressor, character(1))
  )
  answer <- merge(data.table(patient_id = patients), genes, by = "patient_id", all.x = TRUE)
  answer <- merge(answer, summary, by = "patient_id", all.x = TRUE)
  answer[is.na(selected_mutations), selected_mutations := "None detected"]
  answer[is.na(oncogene_mutation), oncogene_mutation := "None detected"]
  answer[is.na(suppressor_pathway_mutation), suppressor_pathway_mutation := "None detected"]
  answer
}

normalize_expression_subtype <- function(histology, subtype) {
  subtype <- tolower(missing_string(subtype))
  fifelse(histology == "LUAD" & subtype == "tru", "LUAD: TRU",
    fifelse(histology == "LUAD" & grepl("inflam", subtype), "LUAD: PI",
      fifelse(histology == "LUAD" & grepl("prolif", subtype), "LUAD: PP",
        fifelse(histology == "LUSC" & subtype %in% c("basal", "classical", "primitive", "secretory"),
          paste0("LUSC: ", tools::toTitleCase(subtype)), NA_character_))))
}

# Published TCGA expression-subtype calls are package data, not a new clustering.
tcga_subtype_dir <- file.path(raw, "annotations", "tcgabiolinks")
dir.create(tcga_subtype_dir, recursive = TRUE, showWarnings = FALSE)
luad_subtype <- as.data.table(TCGAquery_subtype("luad"))
lusc_subtype <- as.data.table(TCGAquery_subtype("lusc"))
fwrite(luad_subtype, file.path(tcga_subtype_dir, "TCGAquery_subtype_LUAD.tsv"), sep = "\t")
fwrite(lusc_subtype, file.path(tcga_subtype_dir, "TCGAquery_subtype_LUSC.tsv"), sep = "\t")
subtypes <- rbindlist(list(
  luad_subtype[, .(patient_id = patient, group = "LUAD", raw_expression_subtype = expression_subtype)],
  lusc_subtype[, .(patient_id = patient, group = "LUSC", raw_expression_subtype = Expression.Subtype)]
), fill = TRUE)
subtypes[, published_expression_subtype := normalize_expression_subtype(group, raw_expression_subtype)]

tcga_manifest <- fread(file.path(processed, "tcga_manifest.tsv"))
clinical <- rbindlist(lapply(c("LUAD", "LUSC"), function(project) {
  path <- file.path(raw, "tcga", sprintf("TCGA-%s.clinical.tsv.gz", project))
  x <- fread(path, check.names = FALSE)
  x[, group := project]
  x
}), fill = TRUE)
clinical[, sample_id := sample]
clinical[, patient_id := submitter_id]
clinical <- clinical[match(tcga_manifest$sample_id, sample_id)]

cbio_dir <- file.path(raw, "annotations", "cbioportal")
study_map <- c(LUAD = "luad_tcga_pan_can_atlas_2018", LUSC = "lusc_tcga_pan_can_atlas_2018")
cbio_patient <- rbindlist(lapply(study_map, function(study) {
  pivot_clinical(read_json_pages(cbio_dir, study, "clinical_patient"), "patientId")
}), fill = TRUE)
setnames(cbio_patient, "patientId", "patient_id")
cbio_sample <- rbindlist(lapply(study_map, function(study) {
  pivot_clinical(read_json_pages(cbio_dir, study, "clinical_sample"), "patientId")
}), fill = TRUE)
setnames(cbio_sample, "patientId", "patient_id")

mutation_tables <- rbindlist(lapply(study_map, function(study) {
  as.data.table(fromJSON(file.path(cbio_dir, sprintf("%s_driver_mutations_full.json", study)), flatten = TRUE))
}), fill = TRUE)
mutations <- mutation_summary(mutation_tables, unique(tcga_manifest$patient_id))

tcga <- data.table(
  sample_id = tcga_manifest$sample_id,
  patient_id = tcga_manifest$patient_id,
  group = tcga_manifest$group,
  sex = title_case_sex(first_existing(clinical, c("gender.demographic"))),
  age_at_diagnosis = suppressWarnings(as.numeric(first_existing(clinical, c("age_at_index.demographic", "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived")))),
  stage = stage_group(first_existing(clinical, c("ajcc_pathologic_stage.diagnoses"))),
  pathologic_t = missing_string(first_existing(clinical, c("ajcc_pathologic_t.diagnoses"))),
  pathologic_n = missing_string(first_existing(clinical, c("ajcc_pathologic_n.diagnoses"))),
  pathologic_m = missing_string(first_existing(clinical, c("ajcc_pathologic_m.diagnoses"))),
  pack_years = suppressWarnings(as.numeric(first_existing(clinical, c("pack_years_smoked.exposures")))),
  vital_status = tools::toTitleCase(tolower(missing_string(first_existing(clinical, c("vital_status.demographic"))))),
  biopsy_site = missing_string(first_existing(clinical, c("site_of_resection_or_biopsy.diagnoses", "tissue_or_organ_of_origin.diagnoses")))
)
tcga[, age_group := age_group(age_at_diagnosis)]
tcga[, pack_year_group := pack_year_group(pack_years)]
tcga <- merge(tcga, subtypes[, .(patient_id, published_expression_subtype)], by = "patient_id", all.x = TRUE, sort = FALSE)
tcga <- merge(tcga, mutations, by = "patient_id", all.x = TRUE, sort = FALSE)
tcga <- merge(tcga, cbio_patient[, intersect(names(cbio_patient), c("patient_id", "OS_STATUS", "OS_MONTHS", "SMOKING_HISTORY")), with = FALSE], by = "patient_id", all.x = TRUE, sort = FALSE)
tmb_col <- intersect(names(cbio_sample), c("TMB_NONSYNONYMOUS", "TMB_SCORE"))
if (length(tmb_col)) {
  tmb <- cbio_sample[, .(tmb_nonsynonymous = suppressWarnings(as.numeric(get(tmb_col[[1]])))), by = patient_id]
  tcga <- merge(tcga, tmb, by = "patient_id", all.x = TRUE, sort = FALSE)
} else {
  tcga[, tmb_nonsynonymous := NA_real_]
}
tcga <- tcga[match(tcga_manifest$sample_id, sample_id)]
tcga[, annotation_sources := "GDC/Xena clinical; TCGAbiolinks published subtype tables; cBioPortal PanCancer Atlas mutations"]
fwrite(tcga, file.path(out, "tcga_figure1_annotations.tsv"), sep = "\t")

gse_manifest <- fread(file.path(processed, "gse81089_manifest.tsv"))
stage_labels <- c(`1` = "IA", `2` = "IB", `3` = "IIA", `4` = "IIB", `5` = "IIIA", `6` = "IIIB", `7` = "IV")
smoking_labels <- c(`1` = "Current", `2` = "Former >1 year", `3` = "Never")
gse <- copy(gse_manifest)
gse[, sex := tools::toTitleCase(tolower(gender))]
gse[, age_at_diagnosis := suppressWarnings(as.numeric(age))]
gse[, age_group := age_group(age_at_diagnosis)]
gse[, stage := unname(stage_labels[as.character(stage_code)])]
gse[, smoking_status := unname(smoking_labels[as.character(smoking_code)])]
gse[, vital_status := fifelse(dead == 1, "Dead", fifelse(dead == 0, "Alive", NA_character_))]
gse[, annotation_sources := "GSE81089 GEO sample metadata and study codebook"]
fwrite(gse, file.path(out, "gse81089_figure1_annotations.tsv"), sep = "\t")

pdmr_path <- file.path(out, "pdmr_figure1_annotations.tsv")
if (!file.exists(pdmr_path)) stop("PDMR annotations must be fetched before building the complete annotation layer")
pdmr <- fread(pdmr_path)
pdmr[, age_group := age_group(fcoalesce(suppressWarnings(as.numeric(age_at_sampling)), suppressWarnings(as.numeric(age_at_diagnosis))))]
pdmr[, pack_year_group := pack_year_group(total_pack_years)]
pdmr[, molecular_record := fifelse(is.na(missing_string(known_molecular_ihc_data)), "Unavailable", "Reported")]
fwrite(pdmr, pdmr_path, sep = "\t")

depmap_manifest <- fread(file.path(processed, "depmap_manifest.tsv"))
mutation_path <- file.path(raw, "depmap", "EH7557_22Q2.RData.xz")
mutation_env <- new.env(parent = emptyenv())
loaded <- load(xzfile(mutation_path), envir = mutation_env)
if (length(loaded) != 1) stop("Expected exactly one object in DepMap mutation resource")
depmap_mut <- as.data.table(get(loaded[[1]], envir = mutation_env))
depmap_mut <- depmap_mut[depmap_id %in% depmap_manifest$depmap_id]
setnames(depmap_mut, c("depmap_id", "gene_name"), c("patient_id", "gene_symbol"))
depmap_summary <- mutation_summary(
  depmap_mut[, .(patientId = patient_id, gene.hugoGeneSymbol = gene_symbol, mutationType = var_class)],
  depmap_manifest$depmap_id
)
setnames(depmap_summary, "patient_id", "depmap_id")
depmap <- depmap_manifest[, .(
  sample_id, depmap_id, group,
  sex = title_case_sex(sex),
  age_at_diagnosis = suppressWarnings(as.numeric(age)),
  primary_or_metastasis = missing_string(primary_or_metastasis),
  sample_collection_site = missing_string(sample_collection_site),
  lineage_molecular_subtype = missing_string(lineage_molecular_subtype),
  growth_pattern = missing_string(default_growth_pattern)
)]
depmap[, age_group := age_group(age_at_diagnosis)]
depmap <- merge(depmap, depmap_summary, by = "depmap_id", all.x = TRUE, sort = FALSE)
depmap <- depmap[match(depmap_manifest$sample_id, sample_id)]
depmap[, annotation_sources := "DepMap 22Q2 model metadata and mutation calls"]
fwrite(depmap, file.path(out, "depmap_figure1_annotations.tsv"), sep = "\t")

coverage_fields <- list(
  TCGA = c("sex", "age_group", "stage", "pack_year_group", "published_expression_subtype", "oncogene_mutation", "suppressor_pathway_mutation"),
  GSE81089 = c("sex", "age_group", "stage", "smoking_status", "vital_status"),
  PDMR_PDX = c("sex", "age_group", "known_metastatic_disease", "has_smoked_100_cigarettes", "biopsy_site", "tissue_type", "molecular_record"),
  DepMap_22Q2 = c("sex", "age_group", "primary_or_metastasis", "sample_collection_site", "growth_pattern", "oncogene_mutation", "suppressor_pathway_mutation")
)
coverage_tables <- list(TCGA = tcga, GSE81089 = gse, PDMR_PDX = pdmr, DepMap_22Q2 = depmap)
coverage <- rbindlist(lapply(names(coverage_fields), function(cohort) {
  tab <- coverage_tables[[cohort]]
  rbindlist(lapply(coverage_fields[[cohort]], function(field) {
    value <- missing_string(tab[[field]])
    data.table(cohort = cohort, field = field, total_n = nrow(tab), available_n = sum(!is.na(value)), coverage = mean(!is.na(value)))
  }))
}))
fwrite(coverage, file.path(out, "figure1_annotation_coverage.tsv"), sep = "\t")

receipt <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  TCGAbiolinks_version = as.character(packageVersion("TCGAbiolinks")),
  published_subtype_counts = as.list(table(tcga$published_expression_subtype, useNA = "ifany")),
  sample_counts = list(TCGA = nrow(tcga), GSE81089 = nrow(gse), PDMR_PDX = nrow(pdmr), DepMap_22Q2 = nrow(depmap)),
  note = "Published TCGA subtype calls were imported; no new expression clustering was relabeled as a molecular subtype."
)
write_json(receipt, file.path(audit, "figure1_annotation_build_receipt.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")

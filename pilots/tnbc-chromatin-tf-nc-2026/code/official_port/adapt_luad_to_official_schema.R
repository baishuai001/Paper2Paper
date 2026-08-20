#!/usr/bin/env Rscript

# LUAD -> TNBC CodeOcean v1 plotting-schema adapter.
#
# This file contains no plotting code.  It only materialises the files and
# column names expected by the unmodified plotting constructors from commit
# edf5314ce5b9ee0e2f88b2310e7c2df5619ad888.

suppressPackageStartupMessages(library(data.table))

usage <- paste(
  "Usage: adapt_luad_to_official_schema.R",
  "FIGURE1_ROOT FIGURE2_ROOT PUBLICATION_ROOT ADAPTER_ROOT"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(usage, call. = FALSE)

figure1_root <- normalizePath(args[[1L]], mustWork = TRUE)
figure2_root <- normalizePath(args[[2L]], mustWork = TRUE)
publication_root <- normalizePath(args[[3L]], mustWork = TRUE)
adapter_root <- args[[4L]]
dir.create(adapter_root, recursive = TRUE, showWarnings = FALSE)
adapter_root <- normalizePath(adapter_root, mustWork = TRUE)

data_root <- file.path(adapter_root, "data")
results_root <- file.path(adapter_root, "results")
audit_root <- file.path(adapter_root, "audit")
for (path in c(data_root, results_root, audit_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

made <- list()
record <- function(adapter_object, source_path, destination_path, rows, columns,
                   semantic_translation = "") {
  made[[length(made) + 1L]] <<- data.table(
    adapter_object = adapter_object,
    source_path = normalizePath(source_path, mustWork = FALSE),
    destination_path = normalizePath(destination_path, mustWork = FALSE),
    rows = as.integer(rows),
    columns = as.integer(columns),
    semantic_translation = semantic_translation
  )
}

must_file <- function(...) {
  candidates <- c(...)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Required frozen input not found: ", paste(candidates, collapse = " | "),
         call. = FALSE)
  }
  normalizePath(hit[[1L]], mustWork = TRUE)
}

write_tsv <- function(x, destination, adapter_object, source_path,
                      semantic_translation = "", row_names = FALSE) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (row_names) {
    write.table(
      as.data.frame(x), destination, sep = "\t", quote = FALSE,
      # Official readers rely on read.table's "one fewer header field"
      # convention to promote the first data field to row names.
      row.names = TRUE, col.names = TRUE
    )
  } else {
    fwrite(as.data.table(x), destination, sep = "\t", quote = FALSE, na = "NA")
  }
  record(adapter_object, source_path, destination, nrow(x), ncol(x),
         semantic_translation)
  invisible(destination)
}

read_activity <- function(path) {
  x <- fread(path, check.names = FALSE)
  if (ncol(x) < 2L) stop("Activity matrix has fewer than two columns: ", path)
  ids <- as.character(x[[1L]])
  x[[1L]] <- NULL
  answer <- as.matrix(x)
  storage.mode(answer) <- "double"
  rownames(answer) <- ids
  answer
}

normalise_group <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fifelse(x %chin% c("LUAD", "ADENOCARCINOMA"), "LUAD",
          fifelse(x %chin% c("LUSC", "SQUAMOUS"), "LUSC", "Unavailable"))
}

normalise_sex <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  fifelse(raw %chin% c("f", "female"), "Female",
          fifelse(raw %chin% c("m", "male"), "Male", "Unavailable"))
}

normalise_vital <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  fifelse(grepl("dead|deceased|^1$|^d$", raw), "Dead",
          fifelse(grepl("alive|^0$|^a$", raw), "Alive", "NA"))
}

normalise_stage <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  fifelse(grepl("^IV", raw), "IV",
          fifelse(grepl("^III", raw), "III",
                  fifelse(grepl("^II", raw), "II",
                          fifelse(grepl("^I", raw), "I", "Unavailable"))))
}

normalise_smoking <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  fifelse(grepl("never|^no$|^3$", raw), "Never",
          fifelse(grepl("current|^yes$|^1$", raw), "Current",
                  fifelse(grepl("former|year|^2$", raw), "Former", "Unavailable")))
}

normalise_subtype <- function(x) {
  raw <- trimws(as.character(x))
  raw[is.na(raw) | raw == ""] <- "Unavailable"
  allowed <- c("LUAD: TRU", "LUAD: PI", "LUAD: PP")
  raw[!raw %chin% allowed] <- "Unavailable"
  raw
}

official_age <- function(x) {
  age <- suppressWarnings(as.numeric(as.character(x)))
  fifelse(is.na(age), "NA",
          fifelse(age <= 40, "<=40", fifelse(age < 70, "(40-70)", ">=70")))
}

first_present <- function(x, candidates, default = "Unavailable") {
  hit <- candidates[candidates %chin% names(x)]
  if (!length(hit)) return(rep(default, nrow(x)))
  value <- as.character(x[[hit[[1L]]]])
  value[is.na(value) | trimws(value) == ""] <- default
  value
}

matrix_as_official <- function(mat, destination, adapter_object, source_path,
                               semantic_translation = "") {
  write_tsv(mat, destination, adapter_object, source_path,
            semantic_translation = semantic_translation, row_names = TRUE)
}

copy_tree <- function(source, destination) {
  source <- normalizePath(source, mustWork = TRUE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(
    source,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = TRUE
  )
  if (!length(entries)) return(invisible(destination))
  relative <- substring(entries, nchar(source) + 2L)
  info <- file.info(entries)
  directories <- which(info$isdir)
  for (index in directories) {
    dir.create(
      file.path(destination, relative[[index]]),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  files <- which(!info$isdir)
  if (length(files)) {
    destinations <- file.path(destination, relative[files])
    for (directory in unique(dirname(destinations))) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    copied <- file.copy(entries[files], destinations, overwrite = TRUE)
    if (!all(copied)) {
      stop("Failed to materialise the isolated secondary-validation adapter tree.",
           call. = FALSE)
    }
  }
  invisible(destination)
}

figure1_tables <- file.path(figure1_root, "results", "tables")
figure1_annotations <- file.path(figure1_root, "data", "processed", "annotations")
publication_data01 <- file.path(
  publication_root, "results", "supplementary_data", "Data_01"
)

tcga_stats_path <- must_file(file.path(figure1_tables, "figure1b_tcga_limma_msviper.tsv"))
validation_stats_path <- must_file(file.path(figure1_tables, "gse41271_limma_msviper.tsv"))
secondary_validation_stats_path <- must_file(
  file.path(figure1_tables, "gse81089_limma_msviper.tsv")
)
tf_membership_path <- must_file(file.path(figure1_tables, "tcga_specific_TFs.tsv"))
tcga_activity_path <- must_file(file.path(figure1_tables, "Figure1C_TCGA_activity_matrix.tsv.gz"))
pdx_activity_path <- must_file(file.path(figure1_tables, "Figure1D_PDMR_PDX_activity_matrix.tsv.gz"))
cell_activity_path <- must_file(file.path(figure1_tables, "Figure1E_DepMap_22Q2_activity_matrix.tsv.gz"))
validation_activity_path <- must_file(
  file.path(figure1_root, "results", "gse41271", "GSE41271_viper_activity.tsv.gz"),
  file.path(figure1_tables, "GSE41271_viper_activity.tsv.gz")
)
secondary_validation_activity_path <- must_file(
  file.path(figure1_root, "results", "gse81089", "GSE81089_viper_activity.tsv.gz"),
  file.path(figure1_tables, "GSE81089_viper_activity.tsv.gz")
)

tcga_annotation_path <- must_file(
  file.path(figure1_annotations, "tcga_figure1_annotations.tsv"),
  file.path(publication_data01, "tcga_manifest.tsv")
)
pdx_annotation_path <- must_file(
  file.path(figure1_annotations, "pdmr_figure1_annotations.tsv"),
  file.path(publication_data01, "pdmr_manifest.tsv")
)
cell_annotation_path <- must_file(
  file.path(figure1_annotations, "depmap_figure1_annotations.tsv"),
  file.path(publication_data01, "depmap_manifest.tsv")
)
validation_annotation_path <- must_file(
  file.path(figure1_root, "data", "processed", "gse41271_manifest.tsv"),
  file.path(publication_data01, "gse41271_manifest.tsv")
)
secondary_validation_annotation_path <- must_file(
  file.path(figure1_annotations, "gse81089_figure1_annotations.tsv"),
  file.path(publication_data01, "gse81089_manifest.tsv")
)

tcga_stats <- fread(tcga_stats_path)
validation_stats <- fread(validation_stats_path)
secondary_validation_stats <- fread(secondary_validation_stats_path)
if ("TF" %chin% names(tcga_stats)) setnames(tcga_stats, "TF", "tf")
if ("NES" %chin% names(tcga_stats)) setnames(tcga_stats, "NES", "TCGA_NES")
if ("TF" %chin% names(validation_stats)) setnames(validation_stats, "TF", "tf")
if ("NES" %chin% names(validation_stats)) setnames(validation_stats, "NES", "METABRIC_NES")
if ("TF" %chin% names(secondary_validation_stats)) {
  setnames(secondary_validation_stats, "TF", "tf")
}
if ("NES" %chin% names(secondary_validation_stats)) {
  setnames(secondary_validation_stats, "NES", "METABRIC_NES")
}
required_tcga_stats <- c("tf", "TCGA_NES", "logFC", "category")
required_validation_stats <- c("tf", "METABRIC_NES", "logFC", "category")
if (!all(required_tcga_stats %chin% names(tcga_stats))) {
  stop("TCGA statistics do not satisfy official Figure1 schema.", call. = FALSE)
}
if (!all(required_validation_stats %chin% names(validation_stats))) {
  stop("GSE41271 statistics do not satisfy official METABRIC-analogue schema.", call. = FALSE)
}
if (!all(required_validation_stats %chin% names(secondary_validation_stats))) {
  stop("GSE81089 statistics do not satisfy secondary validation schema.", call. = FALSE)
}

f1_plot_dir <- file.path(data_root, "Figure1", "Tables_ForPlotting")
write_tsv(
  tcga_stats, file.path(f1_plot_dir, "TCGA_TR_activities_NES_viper.tsv"),
  "Figure1_TCGA_NES", tcga_stats_path,
  "TF->tf; NES->TCGA_NES; LUAD/LUSC categories retained"
)
write_tsv(
  validation_stats,
  file.path(f1_plot_dir, "METABRIC_TR_activities_NES_viper.tsv"),
  "Figure1_validation_NES", validation_stats_path,
  "GSE41271 microarray is the primary METABRIC-analogue; TF->tf; NES->METABRIC_NES"
)

tcga_activity_all <- read_activity(tcga_activity_path)
pdx_activity_all <- read_activity(pdx_activity_path)
cell_activity_all <- read_activity(cell_activity_path)
validation_activity <- read_activity(validation_activity_path)
secondary_validation_activity <- read_activity(secondary_validation_activity_path)

# The official cross-system heatmaps require identical TF rows because the PDX
# and cell-line blocks reuse the TCGA row annotation object.  Never fabricate a
# missing activity score.  Restrict only these Figure 1 heatmap adapter objects
# to TFs genuinely observed in all three systems; retain the full discovery
# statistics for Figure 1B and the external-validation overlap.
heatmap_tf_universe <- unique(toupper(rownames(tcga_activity_all)))
heatmap_tf_availability <- data.table(
  TF = heatmap_tf_universe,
  available_TCGA = heatmap_tf_universe %chin% toupper(rownames(tcga_activity_all)),
  available_PDX = heatmap_tf_universe %chin% toupper(rownames(pdx_activity_all)),
  available_CellLines = heatmap_tf_universe %chin% toupper(rownames(cell_activity_all))
)
heatmap_tf_availability[, display_in_cross_system_heatmaps :=
  available_TCGA & available_PDX & available_CellLines]
heatmap_tf_availability[, boundary_reason := fifelse(
  display_in_cross_system_heatmaps,
  "AVAILABLE_ALL_THREE_SYSTEMS",
  "UNAVAILABLE_ACTIVITY_NOT_ZERO_FILLED_OR_IMPUTED"
)]
fwrite(
  heatmap_tf_availability,
  file.path(audit_root, "figure1_cross_system_heatmap_tf_availability.tsv"),
  sep = "\t"
)
missing_heatmap_tfs <- heatmap_tf_availability[
  display_in_cross_system_heatmaps == FALSE, TF
]
if (length(missing_heatmap_tfs)) {
  projection_regulon_path <- must_file(
    file.path(figure1_root, "results", "tcga", "TCGA_regulon.rds")
  )
  depmap_expression_path <- must_file(
    file.path(
      figure1_root, "data", "processed", "depmap_viper_expression.tsv.gz"
    )
  )
  projection_regulon <- readRDS(projection_regulon_path)
  depmap_expression_genes <- toupper(
    fread(depmap_expression_path, select = 1L)[[1L]]
  )
  missing_heatmap_audit <- rbindlist(lapply(
    missing_heatmap_tfs,
    function(tf) {
      regulon_name <- names(projection_regulon)[
        match(tf, toupper(names(projection_regulon)))
      ]
      targets <- if (!is.na(regulon_name)) {
        names(projection_regulon[[regulon_name]]$tfmode)
      } else {
        character()
      }
      data.table(
        TF = tf,
        tcga_regulon_target_count = length(targets),
        depmap_expression_target_overlap = sum(
          toupper(targets) %chin% depmap_expression_genes
        ),
        activity_status = "NOT_ESTIMABLE_IN_DEPMAP",
        figure_boundary = "FIGURE1E_AVAILABLE_TF_SET_ONLY",
        resolution = "EXCLUDE_FROM_CROSS_SYSTEM_HEATMAP_NO_ZERO_FILL_NO_IMPUTATION",
        regulon_path = normalizePath(projection_regulon_path, mustWork = TRUE),
        depmap_expression_path = normalizePath(
          depmap_expression_path,
          mustWork = TRUE
        )
      )
    }
  ))
  fwrite(
    missing_heatmap_audit,
    file.path(audit_root, "figure1_model_heatmap_missing_tf_audit.tsv"),
    sep = "\t"
  )
}
heatmap_tfs <- heatmap_tf_availability[
  display_in_cross_system_heatmaps == TRUE, TF
]
tcga_activity <- tcga_activity_all[
  match(heatmap_tfs, toupper(rownames(tcga_activity_all))), , drop = FALSE
]
pdx_activity <- pdx_activity_all[
  match(heatmap_tfs, toupper(rownames(pdx_activity_all))), , drop = FALSE
]
cell_activity <- cell_activity_all[
  match(heatmap_tfs, toupper(rownames(cell_activity_all))), , drop = FALSE
]

model_heatmap_stats <- copy(tcga_stats)
model_heatmap_stats[
  category %chin% c("LUAD", "LUSC") & !toupper(tf) %chin% heatmap_tfs,
  category := "Non-Significant"
]
write_tsv(
  model_heatmap_stats,
  file.path(
    f1_plot_dir,
    "TCGA_TR_activities_NES_viper_ModelHeatmapsAvailable.tsv"
  ),
  "Figure1_model_heatmap_available_TFs", tcga_stats_path,
  "Only official PDX/cell-line heatmap row selection uses the 350 TFs observed in all systems; Figure1B remains full"
)

membership <- fread(tf_membership_path)
membership[, TF := toupper(as.character(TF))]
category_lookup <- setNames(as.character(membership$discovery_category), membership$TF)
ordered_tfs <- rownames(tcga_activity)
row_categories <- category_lookup[toupper(ordered_tfs)]
row_categories[is.na(row_categories)] <- "Non-Significant"

matrix_as_official(
  tcga_activity,
  file.path(f1_plot_dir, "TCGA_TR_Activity_TNBC150_NonTNBC155.tsv"),
  "Figure1C_activity", tcga_activity_path,
  "Official legacy filename retained only for schema compatibility; matrix is LUAD/LUSC"
)
row_annotation <- data.frame(VIPER_TR_Category = unname(row_categories),
                             row.names = ordered_tfs, check.names = FALSE)
matrix_as_official(
  as.matrix(row_annotation),
  file.path(f1_plot_dir, "TCGA_TR_Activity_TNBC150_NonTNBC155_TRAnnotation.tsv"),
  "Figure1_row_annotation", tf_membership_path,
  "discovery_category->VIPER_TR_Category; values LUAD/LUSC"
)

make_annotation <- function(path, cohort, sample_ids) {
  x <- fread(path, check.names = FALSE)
  sample_col <- c("sample_id", "Sample", "depmap_id", "geo_accession")
  sample_col <- sample_col[sample_col %chin% names(x)][1L]
  if (is.na(sample_col)) stop("No sample identifier in annotation: ", path)
  x[, Sample := as.character(get(sample_col))]
  x <- x[match(sample_ids, Sample)]
  if (anyNA(x$Sample)) {
    missing <- sample_ids[is.na(x$Sample)]
    stop(cohort, " annotation misses ", length(missing), " samples: ",
         paste(head(missing, 8L), collapse = ", "), call. = FALSE)
  }
  group <- normalise_group(first_present(x, c("group", "histology_normalized", "diagnosis")))
  sex <- normalise_sex(first_present(x, c("sex", "gender")))
  age_raw <- first_present(x, c("age_at_diagnosis", "age", "age_at_surgery", "age_at_sampling"), NA_character_)
  vital <- normalise_vital(first_present(x, c("vital_status", "dead", "OS_STATUS")))
  stage <- normalise_stage(first_present(x, c("stage", "stage_code")))
  smoking <- normalise_smoking(first_present(x, c("smoking_status", "tobacco_use_history", "smoking_code", "has_smoked_100_cigarettes")))
  subtype <- normalise_subtype(first_present(x, c("published_expression_subtype", "lineage_molecular_subtype", "diagnosis_subtype")))
  molecular_raw <- first_present(x, c("molecular_record", "selected_mutations", "oncogene_mutation"))
  molecular <- ifelse(molecular_raw %chin% c("Unavailable", "None detected"),
                      "Unavailable", "Reported")
  cancer_type <- first_present(x, c("primary_or_metastasis"), "Primary")
  cancer_type[!cancer_type %chin% c("Primary", "Metastasis", "Metastatic")] <- "Primary"
  cancer_type[cancer_type == "Metastasis"] <- "Metastatic"

  data.frame(
    Sample = sample_ids,
    Cancer_Type = cancer_type,
    HistologicalType = group,
    Sex = sex,
    Age = official_age(age_raw),
    Stage = stage,
    VitalSt = vital,
    Smoking = smoking,
    PublishedSubtype = subtype,
    MolecularRecord = molecular,
    # Compatibility-only fields used by non-plot side code in the official file.
    ER = "NA", PR = "NA", HER2 = "NA",
    PAM50 = ifelse(group == "LUAD", "Basal", "Normal"),
    SCMOD2 = ifelse(group == "LUAD", "Basal", "LumA"),
    IntClust = "iC1",
    group = group,
    check.names = FALSE
  )
}

tcga_anno <- make_annotation(tcga_annotation_path, "TCGA", colnames(tcga_activity))
pdx_anno <- make_annotation(pdx_annotation_path, "PDX", colnames(pdx_activity))
cell_anno <- make_annotation(cell_annotation_path, "CellLine", colnames(cell_activity))
validation_anno <- make_annotation(validation_annotation_path, "GSE41271", colnames(validation_activity))
secondary_validation_anno <- make_annotation(
  secondary_validation_annotation_path,
  "GSE81089",
  colnames(secondary_validation_activity)
)

write_tsv(
  tcga_anno,
  file.path(f1_plot_dir, "TCGA_TR_Activity_TNBC150_NonTNBC155_SampleAnnotation.tsv"),
  "Figure1C_annotation", tcga_annotation_path,
  "LUAD clinical/molecular annotations mapped to the unchanged HeatmapAnnotation constructor"
)

matrix_as_official(
  cor(tcga_activity, use = "pairwise.complete.obs"),
  file.path(f1_plot_dir, "TCGA_TR_Activity_TNBC150_NonTNBC155_PearsonCor.tsv"),
  "SupplementaryFigure2C_correlation", tcga_activity_path,
  "Sample-sample Pearson correlation from frozen LUAD activity matrix"
)
matrix_as_official(
  pdx_activity,
  file.path(data_root, "Figure1", "PDX", "PDX73_TR_activities_viper_usingTCGAregulons.tsv"),
  "Figure1D_activity", pdx_activity_path,
  "Official legacy filename retained only for schema compatibility"
)
write_tsv(
  pdx_anno,
  file.path(data_root, "Reference_Data", "PDX", "PDX_molecular_subtype.tsv"),
  "Figure1D_annotation", pdx_annotation_path,
  "LUAD PDMR annotations mapped to official annotation object"
)
matrix_as_official(
  cor(pdx_activity, use = "pairwise.complete.obs"),
  file.path(f1_plot_dir, "PDX_TR_Activity_TNBC150_NonTNBC154_PearsonCor.tsv"),
  "SupplementaryFigure2D_correlation", pdx_activity_path
)
matrix_as_official(
  cell_activity,
  file.path(data_root, "Figure1", "CellLines", "CellLines82_TR_activities_viper_usingTCGAregulons.tsv"),
  "Figure1E_activity", cell_activity_path,
  "Official legacy filename retained only for schema compatibility"
)
write_tsv(
  cell_anno,
  file.path(data_root, "Reference_Data", "CellLines_Marcotte_et_al_2016", "CellLines_molecular_subtype.tsv"),
  "Figure1E_annotation", cell_annotation_path,
  "DepMap LUAD/LUSC annotations mapped to official annotation object"
)
lehmann <- data.table(Sample = cell_anno$Sample, Annotation = cell_anno$HistologicalType)
write_tsv(
  lehmann,
  file.path(data_root, "Reference_Data", "CellLines_Marcotte_et_al_2016", "Lehmann_Annotation.tsv"),
  "Figure1E_secondary_annotation", cell_annotation_path,
  "LUAD/LUSC histology occupies the official secondary-annotation slot"
)
matrix_as_official(
  cor(cell_activity, use = "pairwise.complete.obs"),
  file.path(f1_plot_dir, "CellLines_TR_Activity_TNBC150_NonTNBC155_PearsonCor.tsv"),
  "SupplementaryFigure2E_correlation", cell_activity_path
)

matrix_as_official(
  tcga_activity_all,
  file.path(data_root, "Figure1", "TCGA", "TCGA_TR_activities_viper.tsv"),
  "Figure2C_TCGA_activity", tcga_activity_path
)
matrix_as_official(
  validation_activity,
  file.path(data_root, "Figure1", "METABRIC", "METABRIC_TR_activities_viper.tsv"),
  "SupplementaryFigure2F_validation_activity", validation_activity_path,
  "GSE41271 microarray occupies the primary official METABRIC-validation slot"
)
validation_meta <- validation_anno[, c("Sample", "PAM50", "SCMOD2", "group"), drop = FALSE]
write_tsv(
  validation_meta,
  file.path(data_root, "Reference_Data", "METABRIC", "METABRIC_molecular_subtype.tsv"),
  "SupplementaryFigure2F_validation_annotation", validation_annotation_path,
  "GSE41271 group mapped to primary official validation metadata slot"
)
write_tsv(
  data.table(SAMPLE_ID = validation_meta$Sample),
  file.path(data_root, "Reference_Data", "METABRIC", "METABRIC_TNBC233_clinical_info.tsv"),
  "SupplementaryFigure2F_validation_clinical", validation_annotation_path,
  "GSE41271 sample identifiers only"
)

# The official overlap constructor consumes cohort-specific categories.  Keep
# LUAD/LUSC values; the runtime source patch changes only the cancer labels.
# Build the overlap from the actual validation statistics in the active slot;
# the source membership flag predates the GSE41271 microarray validation.
comparison_for_validation <- function(validation_table) {
  tcga_luad <- unique(toupper(tcga_stats[category == "LUAD", tf]))
  tcga_lusc <- unique(toupper(tcga_stats[category == "LUSC", tf]))
  validation_luad <- unique(toupper(validation_table[category == "LUAD", tf]))
  validation_lusc <- unique(toupper(validation_table[category == "LUSC", tf]))
  rbindlist(list(
    data.table(Group = "TNBC_TCGA_Only", TR = setdiff(tcga_luad, validation_luad)),
    data.table(Group = "TNBC_METABRIC_Only", TR = setdiff(validation_luad, tcga_luad)),
    data.table(Group = "TNBC_METABRIC_TCGA_Shared", TR = intersect(tcga_luad, validation_luad)),
    data.table(Group = "NonTNBC_TCGA_Only", TR = setdiff(tcga_lusc, validation_lusc)),
    data.table(Group = "NonTNBC_METABRIC_Only", TR = setdiff(validation_lusc, tcga_lusc)),
    data.table(Group = "NonTNBC_METABRIC_TCGA_Shared", TR = intersect(tcga_lusc, validation_lusc))
  ), use.names = TRUE)
}
comparison <- comparison_for_validation(validation_stats)
write_tsv(
  comparison,
  file.path(f1_plot_dir, "TNBC_TCGA_METABRIC_comparision.tsv"),
  "Figure1_TF_comparison", validation_stats_path,
  "Legacy Group tokens are compatibility keys; overlap is TCGA versus primary GSE41271 LUAD/LUSC"
)

# Figure 2 adapter objects -------------------------------------------------
promoter_sample_path <- must_file(
  file.path(figure2_root, "results", "promoter_gate", "sample_tf_promoter_accessibility.tsv"),
  file.path(publication_root, "results", "supplementary_data", "Data_05", "sample_tf_promoter_accessibility.tsv")
)
promoter_gate_path <- must_file(
  file.path(figure2_root, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"),
  file.path(publication_root, "results", "supplementary_data", "Data_03", "tf_promoter_gate_summary.tsv")
)
motif_best_path <- must_file(
  file.path(figure2_root, "results", "motif_primary", "anchor_consensus_author",
            "anchor_consensus_sample_tf_best_motif.tsv"),
  file.path(publication_root, "results", "supplementary_data", "Data_05",
            "harmonized_sample_tf_best_motif.tsv")
)
motif_primary_all_path <- must_file(
  file.path(figure2_root, "results", "motif_primary", "anchor_consensus_author",
            "anchor_consensus_all_selected_motif_results.tsv"),
  file.path(publication_root, "results", "supplementary_data", "Data_05",
            "harmonized_all_selected_motif_results.tsv")
)
motif_primary_triple_summary_path <- must_file(
  file.path(
    figure2_root, "results", "motif_primary", "anchor_consensus_author",
    "anchor_consensus_hc_tf_triple_system_summary.tsv"
  )
)
motif_candidate_best_path <- must_file(
  file.path(
    figure2_root, "results", "motif_equivalence", "candidate_anchor_equivalent",
    "harmonized_sample_tf_best_motif.tsv"
  )
)
motif_candidate_all_path <- must_file(
  file.path(
    figure2_root, "results", "motif_equivalence", "candidate_anchor_equivalent",
    "harmonized_all_selected_motif_results.tsv"
  )
)
motif_candidate_triple_summary_path <- must_file(
  file.path(
    figure2_root, "results", "motif_equivalence", "candidate_anchor_equivalent",
    "harmonized_hc_tf_triple_system_summary.tsv"
  )
)
motif_inventory_path <- must_file(
  file.path(figure2_root, "results", "motif_gate", "hc_tf_motif_inventory.tsv"),
  file.path(publication_root, "results", "supplementary_data", "Data_05", "hc_tf_motif_inventory.tsv")
)
motif_all_path <- must_file(
  file.path(figure2_root, "results", "motif_equivalence",
            "harmonized_all_selected_motif_results.tsv"),
  file.path(figure2_root, "motif-equivalence", "results", "full1944", "harmonized_all_selected_motif_results.tsv")
)

system_label <- function(x) {
  fifelse(x == "patient", "TCGA", fifelse(x == "PDX", "PDX", "CellLines"))
}

promoter_sample <- fread(promoter_sample_path)
promoter_official <- promoter_sample[
  as.logical(promoter_accessible),
  .(geneID = toupper(TF), Sample = sample_id, Sample_Group = system_label(system))
]
promoter_out <- file.path(
  data_root, "Figure2", "Promoter_Accessibility",
  "Combined_Promoter_Accessibility_TNBC_TR_Genes.tsv"
)
write_tsv(
  promoter_official, promoter_out, "Figure2_promoter_accessibility",
  promoter_sample_path,
  "Only promoter_accessible=TRUE rows are materialised, matching the official intersect-file semantics"
)

# Official 06C derives TCGA RNA column selectors from patient-level ATAC
# barcodes (for example, TCGA-44-3918 -> TCGA.44.3918 after read.table).  Keep
# the Figure 1 aliquot-level activity object untouched and materialise a
# Figure2-only LUAD RNA adapter with one patient barcode per primary-tumour
# aliquot.  No values are aggregated, filled, or imputed.
tcga_activity_luad_columns <- tcga_anno$group == "LUAD"
tcga_figure2_activity <- tcga_activity_all[
  , tcga_activity_luad_columns, drop = FALSE
]
tcga_figure2_aliquot_ids <- colnames(tcga_figure2_activity)
tcga_figure2_patient_ids <- substr(tcga_figure2_aliquot_ids, 1L, 12L)
if (length(tcga_figure2_aliquot_ids) != 516L ||
    anyDuplicated(tcga_figure2_aliquot_ids) ||
    anyDuplicated(tcga_figure2_patient_ids)) {
  stop(
    "Figure2 TCGA RNA adapter must map 516 LUAD aliquots one-to-one to 516 patient barcodes.",
    call. = FALSE
  )
}
if (!setequal(
  tcga_figure2_aliquot_ids,
  tcga_anno$Sample[tcga_anno$group == "LUAD"]
)) {
  stop(
    "Figure2 TCGA RNA aliquots do not exactly match the LUAD annotation set.",
    call. = FALSE
  )
}
colnames(tcga_figure2_activity) <- tcga_figure2_patient_ids
tcga_figure2_activity_path <- file.path(
  data_root, "Figure2", "TCGA",
  "TCGA_TR_activities_viper_patient_barcode.tsv"
)
matrix_as_official(
  tcga_figure2_activity,
  tcga_figure2_activity_path,
  "Figure2C_TCGA_activity_patient_barcode",
  tcga_activity_path,
  paste(
    "Figure2-only 516-to-516 aliquot-to-patient schema mapping for official 06C;",
    "Figure1 activity remains aliquot-level; no aggregation, fill, or imputation"
  )
)

promoter_tcga_patient_ids <- unique(
  as.character(promoter_official[Sample_Group == "TCGA", Sample])
)
matched_activity_patient_ids <- tcga_figure2_patient_ids[
  tcga_figure2_patient_ids %chin% promoter_tcga_patient_ids
]
missing_promoter_rna_patient_ids <- setdiff(
  promoter_tcga_patient_ids,
  tcga_figure2_patient_ids
)
if (!identical(missing_promoter_rna_patient_ids, "TCGA-44-A47F")) {
  stop(
    "Frozen Figure2 paired boundary changed; expected only TCGA-44-A47F to lack RNA/VIPER activity.",
    call. = FALSE
  )
}
paired_patient_ids <- promoter_tcga_patient_ids[
  promoter_tcga_patient_ids %chin% tcga_figure2_patient_ids
]
tcga_figure2_paired_activity <- tcga_figure2_activity[
  , match(paired_patient_ids, colnames(tcga_figure2_activity)), drop = FALSE
]
tcga_figure2_paired_activity_path <- file.path(
  data_root, "Figure2", "TCGA",
  "TCGA_TR_activities_viper_Figure2C_paired21.tsv"
)
matrix_as_official(
  tcga_figure2_paired_activity,
  tcga_figure2_paired_activity_path,
  "Figure2C_TCGA_activity_paired_multiomic",
  tcga_activity_path,
  paste(
    "Figure2C-only paired multiomic intersection: 21 of 22 ATAC patients have",
    "RNA/VIPER activity; TCGA-44-A47F is not filled or imputed"
  )
)

promoter_figure2c_paired <- promoter_official[
  Sample_Group != "TCGA" | Sample %chin% paired_patient_ids
]
promoter_figure2c_paired_path <- file.path(
  data_root, "Figure2", "Promoter_Accessibility",
  "Combined_Promoter_Accessibility_TNBC_TR_Genes_Figure2C_PairedRNA.tsv"
)
write_tsv(
  promoter_figure2c_paired,
  promoter_figure2c_paired_path,
  "Figure2C_promoter_accessibility_paired_multiomic",
  promoter_sample_path,
  paste(
    "Figure2C-only patient/RNA intersection; Figure2B and the promoter gate retain",
    "all 22 ATAC patients; TCGA-44-A47F lacks RNA/VIPER and is not imputed"
  )
)

promoter_activity_exact_set_equal <- setequal(
  promoter_tcga_patient_ids,
  matched_activity_patient_ids
)
official_selector_exact_set_equal <- setequal(
  gsub("-", ".", paired_patient_ids),
  make.names(colnames(tcga_figure2_paired_activity))
)
paired_promoter_activity_exact_set_equal <- setequal(
  paired_patient_ids,
  colnames(tcga_figure2_paired_activity)
)
if (!paired_promoter_activity_exact_set_equal ||
    !official_selector_exact_set_equal) {
  stop(
    "Figure2C paired TCGA promoter barcodes do not map exactly to official 06C RNA selectors.",
    call. = FALSE
  )
}
tcga_figure2_id_audit <- data.table(
  aliquot_sample_id = tcga_figure2_aliquot_ids,
  patient_barcode = tcga_figure2_patient_ids,
  group = "LUAD",
  promoter_atac_patient = tcga_figure2_patient_ids %chin% promoter_tcga_patient_ids,
  mapping_status = "ONE_TO_ONE_NO_AGGREGATION"
)
fwrite(
  tcga_figure2_id_audit,
  file.path(audit_root, "figure2_tcga_rna_sample_id_mapping.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    n_input_luad_aliquots = length(tcga_figure2_aliquot_ids),
    n_unique_patient_barcodes = uniqueN(tcga_figure2_patient_ids),
    n_duplicate_patient_barcodes = sum(duplicated(tcga_figure2_patient_ids)),
    luad_activity_annotation_exact_set_equal = TRUE,
    n_promoter_atac_patients = length(promoter_tcga_patient_ids),
    n_promoter_patients_matched_in_activity = length(matched_activity_patient_ids),
    n_promoter_patients_missing_rna = length(missing_promoter_rna_patient_ids),
    missing_rna_patient_ids = paste(missing_promoter_rna_patient_ids, collapse = ";"),
    full_22_promoter_activity_exact_set_equal = promoter_activity_exact_set_equal,
    paired_21_promoter_activity_exact_set_equal =
      paired_promoter_activity_exact_set_equal,
    official_selector_exact_set_equal = official_selector_exact_set_equal,
    resolution = paste(
      "FIGURE2C_ONLY_PAIRED_N21;FIGURE2B_AND_PROMOTER_GATE_RETAIN_N22;",
      "NO_ZERO_FILL_NO_IMPUTATION",
      sep = ""
    ),
    legend_required_note = paste(
      "TCGA patient panel in Figure 2C uses the paired ATAC-RNA/VIPER",
      "intersection (n=21/22); TCGA-44-A47F lacked RNA/VIPER activity."
    )
  ),
  file.path(audit_root, "figure2_tcga_rna_sample_id_mapping_summary.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    patient_barcode = missing_promoter_rna_patient_ids,
    present_in_luad_atac = TRUE,
    present_in_luad_rna_viper_activity = FALSE,
    action_figure2b_and_promoter_gate = "RETAIN",
    action_figure2c_paired_panel = "EXCLUDE_NO_IMPUTATION",
    evidence_activity_path = normalizePath(tcga_activity_path, mustWork = TRUE),
    evidence_promoter_path = normalizePath(promoter_sample_path, mustWork = TRUE)
  ),
  file.path(audit_root, "figure2c_unpaired_tcga_patient_audit.tsv"),
  sep = "\t"
)

# Mean activity table consumed by the official Supplementary Figure 4B source
# (stored under the capsule's Figure3/02 script). The anchor uses the target-
# subtype RNA cohorts from Figure 1; its RNA and ATAC model collections are not
# sample-matched (for example, 73 RNA-profiled PDX versus 38 ATAC-profiled PDX).
# Preserve that design and average all LUAD-labelled RNA models in each system.
mean_activity_for_system <- function(mat, annotation, system_name) {
  selected <- which(annotation$group == "LUAD")
  if (!length(selected)) stop("No LUAD RNA activity columns in ", system_name, ".")
  data.table(
    geneID = toupper(rownames(mat)),
    Avg_NES_RNA = rowMeans(mat[, selected, drop = FALSE], na.rm = TRUE),
    Sample_Group = system_name
  )
}

mean_activity <- rbindlist(list(
  mean_activity_for_system(
    tcga_activity_all,
    tcga_anno,
    "TCGA"
  ),
  mean_activity_for_system(
    pdx_activity_all,
    pdx_anno,
    "PDX"
  ),
  mean_activity_for_system(
    cell_activity_all,
    cell_anno,
    "CellLines"
  ),
  data.table(
    geneID = toupper(rownames(validation_activity)),
    Avg_NES_RNA = rowMeans(
      validation_activity[, validation_anno$group == "LUAD", drop = FALSE],
      na.rm = TRUE
    ),
    Sample_Group = "METABRIC"
  )
), use.names = TRUE)
mean_activity <- mean_activity[geneID %chin% membership[discovery_category == "LUAD", toupper(TF)]]
mean_activity[, NES_Category := fcase(
  Avg_NES_RNA < 0, "<0",
  Avg_NES_RNA < 1, "0-1",
  Avg_NES_RNA < 3, "1-3",
  default = ">3"
)]
write_tsv(
  mean_activity,
  file.path(data_root, "Figure2", "Combined_RNA_NES_Heatmap.tsv"),
  "SupplementaryFigure4B_mean_activity", tcga_activity_path,
  "Anchor-equivalent mean TF activity across LUAD-labelled Figure 1 RNA cohorts; ATAC cohorts remain independent"
)

gate <- fread(promoter_gate_path)
hc_tfs <- unique(toupper(gate[as.logical(HC_TF_promoter_activity_definition), TF]))
hc_path <- file.path(data_root, "Figure2", "High_RNA_NES.tsv")
dir.create(dirname(hc_path), recursive = TRUE, showWarnings = FALSE)
write.table(hc_tfs, hc_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
record("Figure2_HC_TFs", promoter_gate_path, hc_path, length(hc_tfs), 1L,
       "31 LUAD HC-TFs from the frozen promoter/activity gate")

motif_best <- fread(motif_best_path)[grepl("HC_TF", TF_role)]
motif_primary_all <- fread(motif_primary_all_path)
motif_all <- fread(motif_all_path)
inventory <- fread(motif_inventory_path)

# The formal selected-motif table and all-motif count table are the paired
# outputs of the frozen author-priority run. A best table must never be joined
# to a different run/version because its testability states, q values and HOMER
# counts can differ. Use the complete five-field identity, including the
# filename-safe sample slug, and do not permit a cross-database fallback.
join_keys <- c("system", "sample_id", "sample_slug", "TF", "motif_token")
required_best <- c(
  join_keys, "TF_role", "motif_test_status", "q_value", "log2_odds_ratio"
)
required_all <- c(
  join_keys, "p_value", "target_with_motif", "target_total",
  "background_with_motif", "background_total", "q_value", "log2_odds_ratio"
)
missing_best <- setdiff(required_best, names(motif_best))
missing_all <- setdiff(required_all, names(motif_primary_all))
missing_complete_all <- setdiff(required_all, names(motif_all))
if (length(missing_best) || length(missing_all) || length(missing_complete_all)) {
  stop(
    "Motif schema mismatch. Missing best columns: ",
    paste(missing_best, collapse = ","),
    "; missing primary all-motif columns: ", paste(missing_all, collapse = ","),
    "; missing complete all-motif columns: ",
    paste(missing_complete_all, collapse = ","),
    call. = FALSE
  )
}

motif_best[, TF := toupper(TF)]
motif_primary_all[, TF := toupper(TF)]
motif_all[, TF := toupper(TF)]
all_key_lookup <- motif_primary_all[, .(
  all_match_count = .N,
  all_q_value = q_value[[1L]],
  all_log2_odds_ratio = log2_odds_ratio[[1L]]
), by = join_keys]
ambiguous_all_keys <- all_key_lookup[all_match_count != 1L]
if (nrow(ambiguous_all_keys)) {
  fwrite(
    ambiguous_all_keys,
    file.path(audit_root, "motif_all_ambiguous_join_keys.tsv"),
    sep = "\t"
  )
  stop("The harmonized all-motif table has non-unique motif keys.", call. = FALSE)
}

priority <- inventory[, .(
  TF = toupper(TF),
  motif_database_category,
  expected_primary_source = primary_source,
  motif_testable = as.logical(motif_testable)
)]
motif_best[, selected_database := fcase(
  grepl("^JASPAR", motif_token), "JASPAR2024",
  grepl("^CIS-BP", motif_token), "CIS-BP2.00",
  default = NA_character_
)]
motif_best <- merge(motif_best, priority, by = "TF", all.x = TRUE, sort = FALSE)
motif_best[, database_priority_ok := fifelse(
  motif_test_status == "TESTED",
  !is.na(selected_database) & selected_database == expected_primary_source,
  is.na(motif_token)
)]

motif_join_audit <- merge(
  motif_best,
  all_key_lookup,
  by = join_keys,
  all.x = TRUE,
  sort = FALSE
)
motif_join_audit[is.na(all_match_count), all_match_count := 0L]
motif_join_audit[, paired_q_value_match :=
  motif_test_status != "TESTED" |
  (!is.na(all_q_value) & abs(q_value - all_q_value) <= 1e-12)]
motif_join_audit[, paired_log2_odds_match :=
  motif_test_status != "TESTED" |
  (!is.na(all_log2_odds_ratio) &
     abs(log2_odds_ratio - all_log2_odds_ratio) <= 1e-12)]
motif_join_audit[, join_status := fcase(
  motif_test_status != "TESTED" & is.na(motif_token),
  "NOT_TESTABLE_NO_COUNT_REQUIRED",
  motif_test_status == "TESTED" & all_match_count == 1L &
    database_priority_ok & paired_q_value_match & paired_log2_odds_match,
  "EXACT_PAIRED_MATCH",
  motif_test_status == "TESTED" & all_match_count == 0L,
  "ERROR_TESTED_KEY_MISSING",
  motif_test_status == "TESTED" & !database_priority_ok,
  "ERROR_DATABASE_PRIORITY",
  motif_test_status == "TESTED" &
    (!paired_q_value_match | !paired_log2_odds_match),
  "ERROR_UNPAIRED_STATISTICS",
  default = "ERROR_UNEXPECTED_STATE"
)]
motif_join_audit_path <- file.path(audit_root, "motif_primary_count_join_audit.tsv")
fwrite(motif_join_audit, motif_join_audit_path, sep = "\t", na = "NA")
motif_join_summary <- motif_join_audit[, .N, by = .(
  system, motif_test_status, selected_database,
  expected_primary_source, database_priority_ok,
  paired_q_value_match, paired_log2_odds_match, join_status
)]
fwrite(
  motif_join_summary,
  file.path(audit_root, "motif_primary_count_join_summary.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  data.table(
    input_role = c(
      "formal_author_priority_best",
      "formal_author_priority_all_counts",
      "complete_two_database_inventory"
    ),
    source_path = c(
      normalizePath(motif_best_path, mustWork = TRUE),
      normalizePath(motif_primary_all_path, mustWork = TRUE),
      normalizePath(motif_all_path, mustWork = TRUE)
    ),
    rows = c(nrow(motif_best), nrow(motif_primary_all), nrow(motif_all)),
    deterministic_join_key = paste(join_keys, collapse = "+"),
    database_priority_rule =
      "separate database runs; JASPAR2024 when available; CIS-BP2.00 only for CIS-BP-only TFs; no cross-database fallback",
    permitted_use = c(
      "formal Figure2C-D motif identity and explicit non-testable states",
      "formal Figure2C-D HOMER statistics under JASPAR-first priority",
      "complete database availability for SupplementaryFigure4C only"
    )
  ),
  file.path(audit_root, "motif_paired_input_manifest.tsv"),
  sep = "\t"
)

motif_join_errors <- motif_join_audit[grepl("^ERROR_", join_status)]
if (nrow(motif_join_errors)) {
  stop(
    "Paired motif join/database-priority audit failed for ",
    nrow(motif_join_errors), " rows; see ", motif_join_audit_path,
    call. = FALSE
  )
}

# Rows explicitly marked non-testable remain in the system denominators but do
# not have (and must not be assigned) HOMER counts.
motif_best_tested <- motif_best[motif_test_status == "TESTED"]
motif_primary <- merge(
  motif_best_tested[, c(join_keys, "q_value", "log2_odds_ratio"), with = FALSE],
  motif_primary_all[, c(join_keys, "p_value", "target_with_motif", "target_total",
                "background_with_motif", "background_total"), with = FALSE],
  by = join_keys,
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(motif_primary$target_total)) {
  stop("Could not recover HOMER counts for tested primary motifs.", call. = FALSE)
}
decorate_motif <- function(x) {
  x <- copy(x)
  x[, Sample_Group := system_label(system)]
  x[, Sample := sample_id]
  x[, Motif.Name := toupper(TF)]
  x[, qvalue_Benjamini := q_value]
  x[, Number_TargetSequencesMotif := target_with_motif]
  x[, Number_BackgroundSequencesMotif := background_with_motif]
  x[, Percent_TargetSequencesMotif := paste0(100 * target_with_motif / target_total, "%")]
  x[, Percent_TargetSequencesMotif.1 := paste0(100 * background_with_motif / background_total, "%")]
  x
}
motif_primary <- decorate_motif(motif_primary)
# The formal Figure 2C/D input is exclusively the frozen author-priority run:
# JASPAR rows whenever a JASPAR motif is available, with CIS-BP used only for
# CIS-BP-only TFs. Complete two-database rows are written to separate files for
# Supplementary Figure 4C and never mixed into the formal Figure 2 input.
motif_formal_source <- motif_primary_all[toupper(TF) %chin% hc_tfs]
motif_formal_source[, TF := toupper(TF)]
motif_formal_source[, selected_database := fcase(
  grepl("^JASPAR", motif_token), "JASPAR2024",
  grepl("^CIS-BP", motif_token), "CIS-BP2.00",
  default = NA_character_
)]
motif_formal_source <- merge(
  motif_formal_source,
  priority[, .(TF, expected_primary_source)],
  by = "TF",
  all.x = TRUE,
  sort = FALSE
)
formal_database_priority_violations <- motif_formal_source[
  is.na(expected_primary_source) |
    is.na(selected_database) |
    selected_database != expected_primary_source
]
if (nrow(formal_database_priority_violations)) {
  fwrite(
    formal_database_priority_violations,
    file.path(audit_root, "motif_formal_database_priority_violations.tsv"),
    sep = "\t"
  )
  stop(
    "Formal author-priority motif input contains database-priority violations.",
    call. = FALSE
  )
}
motif_compiled <- decorate_motif(motif_formal_source)

motif_availability_source <- motif_all[toupper(TF) %chin% hc_tfs]
motif_availability_compiled <- decorate_motif(motif_availability_source)
official_motif_columns <- c(
  "Motif.Name", "qvalue_Benjamini", "p_value", "Sample", "Sample_Group",
  "Number_TargetSequencesMotif", "Number_BackgroundSequencesMotif",
  "Percent_TargetSequencesMotif", "Percent_TargetSequencesMotif.1"
)
motif_jaspar <- motif_compiled[grepl("^JASPAR", motif_token), ..official_motif_columns]
motif_cisbp <- motif_compiled[grepl("^CIS-BP", motif_token), ..official_motif_columns]
motif_availability_jaspar <- motif_availability_compiled[
  grepl("^JASPAR", motif_token), ..official_motif_columns
]
motif_availability_cisbp <- motif_availability_compiled[
  grepl("^CIS-BP", motif_token), ..official_motif_columns
]
motif_dir <- file.path(data_root, "Figure2", "Motif", "Combined_MotifEnrichment")
jaspar_path <- file.path(motif_dir, "JASPAR_Combined_MotifEnrichment_TNBC.tsv")
cisbp_path <- file.path(motif_dir, "CISBP_Combined_MotifEnrichment_TNBC.tsv")
jaspar_availability_path <- file.path(
  motif_dir, "JASPAR_Combined_MotifEnrichment_DatabaseAvailability.tsv"
)
cisbp_availability_path <- file.path(
  motif_dir, "CISBP_Combined_MotifEnrichment_DatabaseAvailability.tsv"
)
write_tsv(
  motif_jaspar, jaspar_path, "Figure2_JASPAR_motif_author_priority",
  motif_primary_all_path,
  "Formal author-priority JASPAR rows for Figure2C-D"
)
write_tsv(
  motif_cisbp, cisbp_path, "Figure2_CISBP_motif_author_fallback",
  motif_primary_all_path,
  "Formal CIS-BP rows only for TFs without JASPAR motifs, per author priority"
)
write_tsv(
  motif_availability_jaspar,
  jaspar_availability_path,
  "SupplementaryFigure4C_JASPAR_database_availability",
  motif_all_path,
  "Complete JASPAR availability; SupplementaryFigure4C only, never formal Figure2C-D"
)
write_tsv(
  motif_availability_cisbp,
  cisbp_availability_path,
  "SupplementaryFigure4C_CISBP_database_availability",
  motif_all_path,
  "Complete CIS-BP availability; SupplementaryFigure4C only, never formal Figure2C-D"
)

inventory <- inventory[as.logical(motif_testable)]
if (nrow(inventory) != 20L || uniqueN(toupper(inventory$TF)) != 20L) {
  stop("Frozen formal motif-testable HC-TF set must contain exactly 20 TFs.",
       call. = FALSE)
}
inventory[, Database := fcase(
  motif_database_category == "both", "Both",
  motif_database_category == "JASPAR_only", "JASPAR_only",
  motif_database_category == "CISBP_only", "CISBP_only",
  default = "Not HC-TR"
)]
compare <- inventory[, .(Database, TR = toupper(TF))]
compare_path <- file.path(data_root, "Figure2", "Motif", "JASPAR_CISBP_comparison_long.csv")
dir.create(dirname(compare_path), recursive = TRUE, showWarnings = FALSE)
fwrite(compare, compare_path, sep = ",", quote = FALSE)
record("Figure2_motif_database_map", motif_inventory_path, compare_path,
       nrow(compare), ncol(compare), "Frozen motif inventory category mapping")

# Preserve the frozen 22-patient/13-PDX/19-cell-line denominators.  Samples
# explicitly marked non-testable therefore contribute no enriched call rather
# than silently lowering the author's half-of-system threshold.
system_denominators <- motif_best[, .(
  n_system_samples = uniqueN(sample_id)
), by = .(Sample_Group = system_label(system))]
support_long <- motif_primary[, .(
  n_enriched = uniqueN(Sample[qvalue_Benjamini <= 1e-5]),
  n_observed = uniqueN(Sample)
), by = .(Motif.Name, Sample_Group)]
support_long <- merge(support_long, system_denominators, by = "Sample_Group", all.x = TRUE)
support_long[, supported := n_enriched >= ceiling(n_system_samples / 2)]
support <- dcast(
  support_long,
  Motif.Name ~ Sample_Group, value.var = "supported", fill = FALSE
)
for (column in c("TCGA", "PDX", "CellLines")) {
  if (!column %chin% names(support)) support[, (column) := FALSE]
}
support <- support[, .(TR = Motif.Name, TCGA, PDX, CellLines)]
support[, hc_group := compare$Database[match(TR, compare$TR)]]
support[is.na(hc_group), hc_group := "Not HC-TR"]

formal_triple_summary <- fread(motif_primary_triple_summary_path)
candidate_triple_summary <- fread(motif_candidate_triple_summary_path)
triple_column <- "triple_system_fixed_original_denominator"
if (!triple_column %chin% names(formal_triple_summary) ||
    !triple_column %chin% names(candidate_triple_summary)) {
  stop("Triple-system motif summary schema changed.", call. = FALSE)
}
formal_triple_tfs <- sort(unique(toupper(
  formal_triple_summary[as.logical(get(triple_column)), TF]
)))
candidate_triple_tfs <- sort(unique(toupper(
  candidate_triple_summary[as.logical(get(triple_column)), TF]
)))
support_triple_tfs <- sort(unique(toupper(
  support[TCGA & PDX & CellLines, TR]
)))
expected_formal_triple_tfs <- sort(c("FOXA3", "NFATC4", "XBP1"))
if (!identical(formal_triple_tfs, expected_formal_triple_tfs) ||
    !identical(support_triple_tfs, expected_formal_triple_tfs)) {
  stop(
    "Formal Figure2D triple-system TFs must be exactly FOXA3, NFATC4 and XBP1.",
    call. = FALSE
  )
}
if (length(candidate_triple_tfs) != 0L) {
  stop(
    "Candidate-anchor sensitivity boundary changed; expected zero triple-system TFs.",
    call. = FALSE
  )
}
fwrite(
  rbindlist(list(
    data.table(
      analysis_slot = "FORMAL_PRIMARY_AUTHOR_PRIORITY",
      best_motif_path = normalizePath(motif_best_path, mustWork = TRUE),
      all_counts_path = normalizePath(motif_primary_all_path, mustWork = TRUE),
      triple_summary_path = normalizePath(
        motif_primary_triple_summary_path, mustWork = TRUE
      ),
      n_motif_testable_hc_tfs = nrow(inventory),
      n_triple_system_tfs = length(formal_triple_tfs),
      triple_system_tfs = paste(formal_triple_tfs, collapse = ";"),
      database_priority_violations = 0L,
      permitted_use = "FORMAL_FIGURE2C_D_E"
    ),
    data.table(
      analysis_slot = "CANDIDATE_ANCHOR_VERSION_BOUNDARY",
      best_motif_path = normalizePath(motif_candidate_best_path, mustWork = TRUE),
      all_counts_path = normalizePath(motif_candidate_all_path, mustWork = TRUE),
      triple_summary_path = normalizePath(
        motif_candidate_triple_summary_path, mustWork = TRUE
      ),
      n_motif_testable_hc_tfs = 20L,
      n_triple_system_tfs = length(candidate_triple_tfs),
      triple_system_tfs = paste(candidate_triple_tfs, collapse = ";"),
      database_priority_violations = NA_integer_,
      permitted_use = "SENSITIVITY_AUDIT_ONLY_NOT_FORMAL_FIGURE"
    )
  ), use.names = TRUE),
  file.path(audit_root, "motif_analysis_version_boundary.tsv"),
  sep = "\t",
  na = "NA"
)
support_path <- file.path(data_root, "Figure2", "Motif", "JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv")
support_write <- as.data.frame(support[, .(TCGA, PDX, CellLines, hc_group)])
rownames(support_write) <- support$TR
write_tsv(
  support_write, support_path, "Figure2_motif_support_matrix", motif_best_path,
  "TF row names; system support and database category", row_names = TRUE
)

# Isolated secondary external-validation slot -----------------------------
# Reuse the same official blocks without changing their constructors.  Copy
# the complete schema tree, then replace only the METABRIC compatibility slot
# with GSE81089 objects.  Its outputs are written to a separate results root.
secondary_adapter_root <- file.path(
  adapter_root, "secondary_validation", "GSE81089"
)
secondary_data_root <- file.path(secondary_adapter_root, "data")
secondary_results_root <- file.path(secondary_adapter_root, "results")
copy_tree(data_root, secondary_data_root)
dir.create(secondary_results_root, recursive = TRUE, showWarnings = FALSE)

secondary_f1_plot_dir <- file.path(
  secondary_data_root, "Figure1", "Tables_ForPlotting"
)
write_tsv(
  secondary_validation_stats,
  file.path(secondary_f1_plot_dir, "METABRIC_TR_activities_NES_viper.tsv"),
  "Secondary_GSE81089_validation_NES", secondary_validation_stats_path,
  "GSE81089 RNA-seq occupies the isolated secondary official validation slot"
)
matrix_as_official(
  secondary_validation_activity,
  file.path(
    secondary_data_root, "Figure1", "METABRIC",
    "METABRIC_TR_activities_viper.tsv"
  ),
  "Secondary_GSE81089_validation_activity", secondary_validation_activity_path,
  "GSE81089 RNA-seq activity in the isolated secondary official validation slot"
)
secondary_validation_meta <- secondary_validation_anno[
  , c("Sample", "PAM50", "SCMOD2", "group"), drop = FALSE
]
write_tsv(
  secondary_validation_meta,
  file.path(
    secondary_data_root, "Reference_Data", "METABRIC",
    "METABRIC_molecular_subtype.tsv"
  ),
  "Secondary_GSE81089_validation_annotation",
  secondary_validation_annotation_path,
  "GSE81089 group mapped to the isolated secondary official metadata slot"
)
write_tsv(
  data.table(SAMPLE_ID = secondary_validation_meta$Sample),
  file.path(
    secondary_data_root, "Reference_Data", "METABRIC",
    "METABRIC_TNBC233_clinical_info.tsv"
  ),
  "Secondary_GSE81089_validation_clinical",
  secondary_validation_annotation_path,
  "GSE81089 sample identifiers only"
)
secondary_comparison <- comparison_for_validation(secondary_validation_stats)
write_tsv(
  secondary_comparison,
  file.path(
    secondary_f1_plot_dir, "TNBC_TCGA_METABRIC_comparision.tsv"
  ),
  "Secondary_GSE81089_TF_comparison", secondary_validation_stats_path,
  "Legacy Group tokens; overlap is TCGA versus secondary GSE81089 LUAD/LUSC"
)

validation_slot_manifest <- rbindlist(list(
  data.table(
    validation_role = "PRIMARY_METABRIC_ANALOGUE",
    dataset = "GSE41271",
    assay = "MICROARRAY",
    n_LUAD = sum(validation_anno$group == "LUAD"),
    n_LUSC = sum(validation_anno$group == "LUSC"),
    statistics_path = normalizePath(validation_stats_path, mustWork = TRUE),
    activity_path = normalizePath(validation_activity_path, mustWork = TRUE),
    annotation_path = normalizePath(validation_annotation_path, mustWork = TRUE),
    adapter_slot_root = normalizePath(adapter_root, mustWork = TRUE)
  ),
  data.table(
    validation_role = "SECONDARY_INDEPENDENT_VALIDATION",
    dataset = "GSE81089",
    assay = "RNA_SEQ",
    n_LUAD = sum(secondary_validation_anno$group == "LUAD"),
    n_LUSC = sum(secondary_validation_anno$group == "LUSC"),
    statistics_path = normalizePath(secondary_validation_stats_path, mustWork = TRUE),
    activity_path = normalizePath(secondary_validation_activity_path, mustWork = TRUE),
    annotation_path = normalizePath(
      secondary_validation_annotation_path,
      mustWork = TRUE
    ),
    adapter_slot_root = normalizePath(secondary_adapter_root, mustWork = TRUE)
  )
), use.names = TRUE)
fwrite(
  validation_slot_manifest,
  file.path(audit_root, "external_validation_slot_manifest.tsv"),
  sep = "\t"
)

fwrite(rbindlist(made, fill = TRUE), file.path(audit_root, "adapter_object_manifest.tsv"), sep = "\t")
cat("Adapter root:", adapter_root, "\n")
cat("Objects materialised:", length(made), "\n")
cat("HC-TFs:", length(hc_tfs), "\n")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: audit_luad_pharmacogenomics.R FIG1_ROOT FIG2_ROOT OUTPUT_ROOT")
}
fig1 <- normalizePath(args[[1]], mustWork = TRUE)
fig2 <- normalizePath(args[[2]], mustWork = TRUE)
out <- normalizePath(args[[3]], mustWork = TRUE)
raw_dir <- file.path(out, "data", "raw", "psets")
processed <- file.path(out, "data", "processed")
tables <- file.path(out, "results", "tables")
audit <- file.path(out, "audit")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)

norm_key <- function(x) toupper(gsub("[^A-Z0-9]", "", fifelse(is.na(x), "", as.character(x))))
first_present <- function(tab, candidates) {
  hit <- candidates[tolower(candidates) %in% tolower(names(tab))]
  if (!length(hit)) return(NULL)
  names(tab)[match(tolower(hit[[1]]), tolower(names(tab)))]
}
first_nonempty <- function(...) {
  values <- list(...)
  if (!length(values)) return(character())
  result <- rep(NA_character_, max(vapply(values, length, integer(1))))
  for (value in values) {
    value <- rep(as.character(value), length.out = length(result))
    take <- (is.na(result) | result == "") & !is.na(value) & value != ""
    result[take] <- value[take]
  }
  result
}

manifest <- fread(file.path(fig1, "data", "processed", "depmap_manifest.tsv"))
manifest <- manifest[group == "LUAD"]
activity_tab <- fread(file.path(fig1, "results", "depmap", "DEPMAP_viper_activity.tsv.gz"),
                      check.names = FALSE)
tf_col <- names(activity_tab)[1]
activity <- as.matrix(activity_tab[, -1, with = FALSE])
storage.mode(activity) <- "double"
rownames(activity) <- activity_tab[[tf_col]]
hc <- fread(file.path(fig2, "results", "promoter_gate", "tf_promoter_gate_summary.tsv"))
hc_col <- first_present(hc, c("TF", "tf", "regulator"))
flag_col <- first_present(hc, c("HC_TF_promoter_activity_definition", "HC_TF"))
if (is.null(hc_col) || is.null(flag_col)) stop("Cannot identify HC-TF columns")
hc_tfs <- unique(hc[get(flag_col) == TRUE, get(hc_col)])
if (length(hc_tfs) != 31L) stop("Figure 5 requires the frozen set of 31 HC-TFs")
if (!all(hc_tfs %in% rownames(activity))) stop("Some HC-TFs are absent from DepMap activity")
activity <- activity[hc_tfs, intersect(colnames(activity), manifest$depmap_id), drop = FALSE]
manifest <- manifest[match(colnames(activity), depmap_id)]

key_rows <- rbindlist(list(
  data.table(depmap_id = manifest$depmap_id, key_type = "depmap", lookup_key = norm_key(manifest$depmap_id)),
  data.table(depmap_id = manifest$depmap_id, key_type = "cellosaurus", lookup_key = norm_key(manifest$RRID)),
  data.table(depmap_id = manifest$depmap_id, key_type = "cosmic", lookup_key = norm_key(manifest$cosmic_id)),
  data.table(depmap_id = manifest$depmap_id, key_type = "name", lookup_key = norm_key(manifest$cell_line_name)),
  data.table(depmap_id = manifest$depmap_id, key_type = "name", lookup_key = norm_key(manifest$stripped_cell_line_name)),
  data.table(depmap_id = manifest$depmap_id, key_type = "name", lookup_key = norm_key(sub("_LUNG$", "", manifest$cell_line)))
), use.names = TRUE)
key_rows <- key_rows[lookup_key != "" & !is.na(lookup_key)]
key_rows[, key_n := uniqueN(depmap_id), by = .(key_type, lookup_key)]
key_rows <- unique(key_rows[key_n == 1L, .(depmap_id, key_type, lookup_key)])

map_cells <- function(cell, dataset) {
  cell <- as.data.table(cell, keep.rownames = "legacy_rowname")
  direct_col <- first_present(cell, c("depmap_id", "DepMap_ID", "depmapid"))
  rrid_col <- first_present(cell, c("Cellosaurus.Accession.id", "cellosaurus_accession", "rrid"))
  cosmic_col <- first_present(cell, c("cosmic_id", "cosmicid", "COSMIC_ID"))
  name_cols <- unique(na.omit(c(
    first_present(cell, c("cellid")), first_present(cell, c("ccl_name")),
    first_present(cell, c("CellLine", "cell_line_name")), "legacy_rowname"
  )))
  cell[, source_cell_id := if (!is.null(first_present(cell, c("cellid")))) {
    as.character(get(first_present(cell, c("cellid"))))
  } else as.character(legacy_rowname)]
  cell[, depmap_id := NA_character_]
  cell[, match_method := NA_character_]
  attempt <- function(values, type, method) {
    keys <- norm_key(values)
    lookup <- key_rows[key_type == type]
    ids <- lookup$depmap_id[match(keys, lookup$lookup_key)]
    take <- is.na(cell$depmap_id) & !is.na(ids)
    cell[take, `:=`(depmap_id = ids[take], match_method = method)]
  }
  if (!is.null(direct_col)) attempt(cell[[direct_col]], "depmap", "DepMap_ID")
  if (!is.null(rrid_col)) attempt(cell[[rrid_col]], "cellosaurus", "Cellosaurus")
  if (!is.null(cosmic_col)) attempt(cell[[cosmic_col]], "cosmic", "COSMIC")
  for (column in name_cols) attempt(cell[[column]], "name", paste0("name:", column))
  cell[, dataset := dataset]
  cell[, .(dataset, source_cell_id, depmap_id, match_method,
           eligible_luad = !is.na(depmap_id) & depmap_id %in% colnames(activity))]
}

drug_metadata <- function(drug, dataset) {
  drug <- as.data.table(drug, keep.rownames = "legacy_rowname")
  id_col <- first_present(drug, c("drugid", "drug.id", "PRISM.drugid", "master_cpd_id"))
  if (is.null(id_col)) id_col <- "legacy_rowname"
  name_cols <- unique(na.omit(c(
    first_present(drug, c("drugid", "drug.id")),
    first_present(drug, c("cpd_name", "standard.name", "drug_name", "compound")),
    "legacy_rowname"
  )))
  display <- do.call(first_nonempty, lapply(name_cols, function(z) drug[[z]]))
  cid_col <- first_present(drug, c("cid", "pubchem..CID", "pubchem_cid"))
  inchi_col <- first_present(drug, c("inchikey", "pubchem.InChIKey", "inchi_key"))
  cid <- if (!is.null(cid_col)) as.character(drug[[cid_col]]) else rep(NA_character_, nrow(drug))
  inchi <- if (!is.null(inchi_col)) as.character(drug[[inchi_col]]) else rep(NA_character_, nrow(drug))
  target_cols <- unique(na.omit(c(
    first_present(drug, c("TARGET", "target")),
    first_present(drug, c("gene_symbol_of_protein_target"))
  )))
  pathway_cols <- unique(na.omit(c(
    first_present(drug, c("TARGET_PATHWAY", "target_pathway")),
    first_present(drug, c("moa")),
    first_present(drug, c("target_or_activity_of_compound"))
  )))
  target <- if (length(target_cols)) do.call(first_nonempty, lapply(target_cols, function(z) drug[[z]])) else rep(NA_character_, nrow(drug))
  pathway <- if (length(pathway_cols)) do.call(first_nonempty, lapply(pathway_cols, function(z) drug[[z]])) else rep(NA_character_, nrow(drug))
  cid <- sub("\\.0$", "", cid)
  cid[!grepl("^[0-9]+$", cid)] <- NA_character_
  inchi[!grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", toupper(inchi))] <- NA_character_
  canonical <- ifelse(!is.na(cid), paste0("CID:", cid),
                      ifelse(!is.na(inchi), paste0("INCHI:", toupper(inchi)),
                             paste0("NAME:", norm_key(display))))
  data.table(dataset = dataset, source_drug_id = as.character(drug[[id_col]]),
             drug_name = display, canonical_drug = canonical, pubchem_cid = cid,
             inchikey = inchi, target = target, pathway_or_moa = pathway)
}

extract_dataset <- function(path, dataset) {
  object <- readRDS(path)
  attrs <- attributes(object)
  sens <- attrs[["sensitivity"]]
  if (is.null(sens)) stop(dataset, ": legacy sensitivity payload missing")
  info <- as.data.table(sens[["info"]], keep.rownames = "experiment_key")
  profiles <- as.data.table(sens[["profiles"]], keep.rownames = "experiment_key")
  if (!identical(info$experiment_key, profiles$experiment_key)) {
    profiles <- profiles[match(info$experiment_key, experiment_key)]
  }
  metric <- first_present(profiles, c("aac_recomputed", "AAC_recomputed", "aac_published", "AAC"))
  auc_metric <- first_present(profiles, c("auc_recomputed", "AUC_recomputed", "auc_published", "AUC"))
  if (is.null(metric) && is.null(auc_metric)) stop(dataset, ": no AAC/AUC column")
  aac <- if (!is.null(metric)) as.numeric(profiles[[metric]]) else 1 - as.numeric(profiles[[auc_metric]])
  if (!is.null(metric)) {
    fallback <- setdiff(c("aac_published", "AAC"), metric)
    fallback <- first_present(profiles, fallback)
    if (!is.null(fallback)) {
      missing <- !is.finite(aac)
      aac[missing] <- as.numeric(profiles[[fallback]])[missing]
    }
  }
  cell_col <- first_present(info, c("cellid", "cell.id", "sampleid", "sample_id"))
  drug_col <- first_present(info, c("drugid", "drug.id", "treatmentid", "treatment_id"))
  if (is.null(cell_col) || is.null(drug_col)) stop(dataset, ": response identifiers missing")
  mapping <- map_cells(attrs[["cell"]], dataset)
  drugs <- drug_metadata(attrs[["drug"]], dataset)
  response <- data.table(source_cell_id = as.character(info[[cell_col]]),
                         source_drug_id = as.character(info[[drug_col]]), aac = aac)
  response <- response[is.finite(aac) & aac >= 0 & aac <= 1]
  response <- response[, .(aac = median(aac), experiments = .N),
                       by = .(source_cell_id, source_drug_id)]
  response <- merge(response, mapping, by = "source_cell_id", all.x = TRUE)
  response <- merge(response, drugs[, .(source_drug_id, drug_name, canonical_drug,
                                        pubchem_cid, inchikey, target, pathway_or_moa)],
                    by = "source_drug_id", all.x = TRUE)
  response[, dataset := dataset]
  paired <- response[eligible_luad == TRUE & depmap_id %in% colnames(activity)]
  eligible <- paired[, .(paired_lines = uniqueN(depmap_id), max_aac = max(aac),
                         median_aac = median(aac), drug_name = first(na.omit(drug_name))),
                     by = .(dataset, canonical_drug)]
  eligible[, eligible_drug := paired_lines >= 10L & max_aac > 0.2]
  list(mapping = mapping, drugs = drugs, response = response, paired = paired,
       eligible = eligible, metric = if (!is.null(metric)) metric else paste0("1-", auc_metric),
       raw_cells = nrow(attrs[["cell"]]), raw_drugs = nrow(attrs[["drug"]]))
}

files <- c(GDSC2 = "GDSC2.rds", CTRPv2 = "PSet_CTRPv2.rds", PRISM = "PSet_PRISM.rds")
missing <- files[!file.exists(file.path(raw_dir, files))]
if (length(missing)) stop("Missing PSet files: ", paste(missing, collapse = ", "))
objects <- lapply(names(files), function(dataset) extract_dataset(file.path(raw_dir, files[[dataset]]), dataset))
names(objects) <- names(files)

mapping <- rbindlist(lapply(objects, `[[`, "mapping"), fill = TRUE)
eligible <- rbindlist(lapply(objects, `[[`, "eligible"), fill = TRUE)
drug_meta <- rbindlist(lapply(objects, `[[`, "drugs"), fill = TRUE)
summary <- rbindlist(lapply(names(objects), function(dataset) {
  object <- objects[[dataset]]
  data.table(dataset = dataset, raw_cells = object$raw_cells, raw_drugs = object$raw_drugs,
             mapped_LUAD_lines = uniqueN(object$paired$depmap_id),
             paired_cell_drug_rows = nrow(object$paired),
             eligible_drugs_n10_AACgt0.2 = object$eligible[eligible_drug == TRUE, .N],
             sensitivity_metric = object$metric)
}))

fwrite(mapping, file.path(tables, "SupplementaryFigure8A_cell_mapping.tsv"), sep = "\t")
fwrite(eligible, file.path(tables, "SupplementaryFigure8B_drug_eligibility.tsv"), sep = "\t")
fwrite(unique(drug_meta), file.path(tables, "Figure5_drug_identity_map.tsv"), sep = "\t")
fwrite(summary, file.path(tables, "Figure5_dataset_audit.tsv"), sep = "\t")
saveRDS(list(activity = activity, manifest = manifest, hc_tfs = hc_tfs, datasets = objects),
        file.path(processed, "figure5_cellline_inputs.rds"))

hashes <- tools::md5sum(file.path(raw_dir, files))
receipt <- list(
  status = "passed", frozen_HC_TFs = length(hc_tfs), eligible_LUAD_activity_lines = ncol(activity),
  dataset_summary = summary, raw_md5 = as.list(hashes),
  mapping_priority = c("DepMap_ID", "Cellosaurus", "COSMIC", "unambiguous normalized name"),
  drug_filter = "paired LUAD lines >=10 and at least one AAC >0.2"
)
write_json(receipt, file.path(audit, "figure5_cellline_input_receipt.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

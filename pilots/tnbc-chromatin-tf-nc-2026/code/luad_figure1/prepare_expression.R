#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: prepare_expression.R CLOUD_RUN_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
raw <- file.path(root, "data", "raw")
processed <- file.path(root, "data", "processed")
audit <- file.path(root, "audit")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)

collapse_symbols <- function(matrix, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "."
  tab <- as.data.table(matrix[keep, , drop = FALSE])
  tab[, gene := symbols[keep]]
  result <- tab[, lapply(.SD, mean), by = gene]
  genes <- result$gene
  result[, gene := NULL]
  matrix <- as.matrix(result)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  colnames(matrix) <- setdiff(names(result), "gene")
  matrix
}

collapse_counts <- function(matrix, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "."
  tab <- as.data.table(matrix[keep, , drop = FALSE])
  tab[, gene := symbols[keep]]
  # The anchor's deal_with_duplicates() uses column means even for raw counts.
  # Retain that behavior so limma input matches the published code path.
  result <- tab[, lapply(.SD, mean), by = gene]
  genes <- result$gene
  result[, gene := NULL]
  matrix <- as.matrix(result)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  colnames(matrix) <- setdiff(names(result), "gene")
  matrix
}

write_matrix <- function(matrix, path) {
  genes <- rownames(matrix)
  output <- as.data.table(matrix)
  output[, gene := genes]
  setcolorder(output, "gene")
  fwrite(output, path, sep = "\t", quote = FALSE)
}

select_primary <- function(samples) {
  primary <- samples[substr(samples, 14, 15) == "01"]
  patient <- substr(primary, 1, 12)
  primary[!duplicated(patient)]
}

probemap <- fread(file.path(raw, "tcga", "gencode.v36.annotation.gtf.gene.probemap"))
map_symbol <- setNames(probemap$gene, probemap$id)

read_tcga <- function(project) {
  counts_path <- file.path(raw, "tcga", sprintf("TCGA-%s.star_counts.tsv.gz", project))
  tpm_path <- file.path(raw, "tcga", sprintf("TCGA-%s.star_tpm.tsv.gz", project))
  counts_tab <- fread(counts_path, check.names = FALSE)
  tpm_tab <- fread(tpm_path, check.names = FALSE)
  ids <- counts_tab[[1]]
  if (!identical(ids, tpm_tab[[1]])) stop(sprintf("TCGA-%s count/TPM gene rows differ", project))
  counts_tab[[1]] <- NULL
  tpm_tab[[1]] <- NULL
  samples <- select_primary(intersect(names(counts_tab), names(tpm_tab)))
  log2_counts <- as.matrix(counts_tab[, ..samples])
  storage.mode(log2_counts) <- "double"
  counts <- round(pmax(2^log2_counts - 1, 0))
  log2_tpm <- as.matrix(tpm_tab[, ..samples])
  storage.mode(log2_tpm) <- "double"
  tpm <- pmax(2^log2_tpm - 0.001, 0)
  symbols <- unname(map_symbol[ids])
  counts <- collapse_counts(counts, symbols)
  # Xena star_tpm is stored as log2(TPM + 0.001). Collapse duplicate symbols
  # on the TPM scale, then log-transform; ARACNe applies the same copula/rank
  # transform as in the TNBC capsule.
  expression <- collapse_symbols(tpm, symbols)
  manifest <- data.frame(
    sample_id = samples,
    patient_id = substr(samples, 1, 12),
    group = project,
    system = "patient_discovery",
    source = sprintf("TCGA-%s", project),
    stringsAsFactors = FALSE
  )
  list(counts = counts, expression = expression, manifest = manifest)
}

luad <- read_tcga("LUAD")
lusc <- read_tcga("LUSC")
tcga_genes <- intersect(rownames(luad$counts), rownames(lusc$counts))
tcga_counts <- cbind(luad$counts[tcga_genes, ], lusc$counts[tcga_genes, ])
tcga_expression <- cbind(luad$expression[tcga_genes, ], lusc$expression[tcga_genes, ])
tcga_manifest <- rbind(luad$manifest, lusc$manifest)
stopifnot(identical(colnames(tcga_counts), tcga_manifest$sample_id))
write_matrix(tcga_counts, file.path(processed, "tcga_counts_symbols.tsv.gz"))
write_matrix(tcga_expression, file.path(processed, "tcga_aracne_expression.tsv"))
fwrite(tcga_manifest, file.path(processed, "tcga_manifest.tsv"), sep = "\t")

parse_soft <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  sample_ends <- c(sample_starts[-1] - 1L, length(lines))
  records <- lapply(seq_along(sample_starts), function(index) {
    block <- lines[sample_starts[index]:sample_ends[index]]
    value <- function(prefix) sub(prefix, "", block[grep(prefix, block)[1]])
    chars <- sub("^!Sample_characteristics_ch1 = ", "", block[grep("^!Sample_characteristics_ch1 = ", block)])
    pairs <- strsplit(chars, ": ", fixed = TRUE)
    fields <- setNames(vapply(pairs, function(x) paste(x[-1], collapse = ": "), character(1)), vapply(pairs, `[[`, character(1), 1))
    data.frame(
      geo_accession = value("^!Sample_geo_accession = "),
      sample_id = value("^!Sample_title = "),
      histology_code = unname(fields["histology"]),
      stage_code = unname(fields["stage tnm"]),
      age = unname(fields["age"]),
      gender = unname(fields["gender"]),
      dead = unname(fields["dead"]),
      smoking_code = unname(fields["smoking"]),
      stringsAsFactors = FALSE
    )
  })
  rbindlist(records, fill = TRUE)
}

gse_manifest <- parse_soft(file.path(raw, "gse81089", "GSE81089_family.soft.gz"))
# GSE81089 authors encode histology as: 1 squamous cell cancer,
# 2 adenocarcinoma unspecified, 3 large-cell/NOS. Only tumor samples ending in
# T and codes 1/2 are retained (official GEO series description).
gse_manifest <- gse_manifest[grepl("T(_|$)", sample_id) & histology_code %in% c("1", "2")]
gse_manifest[, group := ifelse(histology_code == "2", "LUAD", "LUSC")]
gse_manifest[, `:=`(system = "patient_validation", source = "GSE81089")]
gse_group_counts <- table(gse_manifest$group)
stopifnot(
  unname(gse_group_counts["LUAD"]) == 108L,
  unname(gse_group_counts["LUSC"]) == 67L
)
# GEO explicitly documents these processed-column aliases. Canonicalize them
# before matching expression columns to the sample-level SOFT manifest.
gse_processed_column_aliases <- c(L608T_2122 = "L608T", L771T_1 = "L771T")

tcga_map_base <- setNames(probemap$gene, sub("\\..*$", "", probemap$id))

read_gse <- function(filename) {
  tab <- fread(file.path(raw, "gse81089", filename), check.names = FALSE)
  ids <- sub("\\..*$", "", tab[[1]])
  tab[[1]] <- NULL
  processed_names <- names(tab)
  alias_index <- match(processed_names, names(gse_processed_column_aliases))
  has_alias <- !is.na(alias_index)
  processed_names[has_alias] <- unname(gse_processed_column_aliases[alias_index[has_alias]])
  if (anyDuplicated(processed_names)) {
    stop(sprintf("Duplicate GSE81089 sample names after alias canonicalization: %s", filename))
  }
  setnames(tab, processed_names)
  missing_samples <- setdiff(gse_manifest$sample_id, names(tab))
  if (length(missing_samples)) {
    stop(sprintf("GSE81089 target samples missing from %s: %s",
                 filename, paste(missing_samples, collapse = ", ")))
  }
  samples <- intersect(gse_manifest$sample_id, names(tab))
  matrix <- as.matrix(tab[, ..samples])
  storage.mode(matrix) <- "double"
  symbols <- unname(tcga_map_base[ids])
  if (grepl("readcounts", filename, fixed = TRUE)) {
    collapse_counts(matrix, symbols)
  } else {
    collapse_symbols(matrix, symbols)
  }
}

gse_counts <- read_gse("GSE81089_readcounts_featurecounts.tsv.gz")
gse_fpkm <- read_gse("GSE81089_FPKM_cufflinks.tsv.gz")
gse_genes <- intersect(rownames(gse_counts), rownames(gse_fpkm))
gse_samples <- intersect(colnames(gse_counts), colnames(gse_fpkm))
gse_counts <- gse_counts[gse_genes, gse_samples, drop = FALSE]
gse_expression <- log2(gse_fpkm[gse_genes, gse_samples, drop = FALSE] + 1)
gse_manifest <- as.data.frame(gse_manifest[match(gse_samples, sample_id)])
stopifnot(identical(gse_samples, gse_manifest$sample_id))
write_matrix(gse_counts, file.path(processed, "gse81089_counts_symbols.tsv.gz"))
write_matrix(gse_expression, file.path(processed, "gse81089_aracne_expression.tsv"))
fwrite(gse_manifest, file.path(processed, "gse81089_manifest.tsv"), sep = "\t")

pdmr_tab <- fread(file.path(processed, "pdmr_tpm_symbols.tsv"), check.names = FALSE)
pdmr_symbols <- pdmr_tab[[1]]
pdmr_tab[[1]] <- NULL
pdmr_matrix <- as.matrix(pdmr_tab)
storage.mode(pdmr_matrix) <- "double"
pdmr_expression <- collapse_symbols(pdmr_matrix, pdmr_symbols)
write_matrix(pdmr_expression, file.path(processed, "pdmr_viper_expression.tsv.gz"))

depmap_meta <- fread(file.path(raw, "depmap", "depmap_metadata_22Q2.tsv.gz"), check.names = FALSE)
strict_luad <- !is.na(depmap_meta$lineage) & !is.na(depmap_meta$subtype_disease) &
  depmap_meta$lineage == "lung" & depmap_meta$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma"
strict_lusc <- !is.na(depmap_meta$lineage) & !is.na(depmap_meta$subtype_disease) &
  depmap_meta$lineage == "lung" & depmap_meta$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma"
depmap_manifest <- depmap_meta[strict_luad | strict_lusc]
depmap_manifest[, group := ifelse(
  subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma",
  "LUAD", "LUSC"
)]
depmap_manifest[, `:=`(sample_id = depmap_id, system = "cell_line", source = "DepMap_22Q2")]

depmap_long <- fread(file.path(raw, "depmap", "depmap_TPM_22Q2.tsv.gz"), check.names = FALSE)
required_depmap <- c("depmap_id", "gene_name", "rna_expression")
if (!all(required_depmap %in% names(depmap_long))) {
  stop(sprintf("Unexpected DepMap TPM columns: %s", paste(names(depmap_long), collapse = ", ")))
}
depmap_long <- depmap_long[depmap_id %in% depmap_manifest$sample_id]
depmap_ids_with_expression <- unique(depmap_long$depmap_id)
missing_depmap_expression <- setdiff(depmap_manifest$sample_id, depmap_ids_with_expression)
depmap_manifest <- depmap_manifest[sample_id %in% depmap_ids_with_expression]
# EH7556 stores log2(TPM+1); invert it because the TNBC projection scripts
# supply TPM-scale PDX and cell-line matrices directly to VIPER.
depmap_long[, rna_expression := pmax(2^rna_expression - 1, 0)]
depmap_wide <- dcast(depmap_long, gene_name ~ depmap_id, value.var = "rna_expression")
depmap_genes <- depmap_wide$gene_name
depmap_wide[, gene_name := NULL]
depmap_expression <- as.matrix(depmap_wide)
storage.mode(depmap_expression) <- "double"
rownames(depmap_expression) <- depmap_genes
depmap_manifest <- as.data.frame(depmap_manifest[match(colnames(depmap_expression), sample_id)])
write_matrix(depmap_expression, file.path(processed, "depmap_viper_expression.tsv.gz"))
fwrite(depmap_manifest, file.path(processed, "depmap_manifest.tsv"), sep = "\t")

pango <- scan(file.path(root, "reference", "pango_regulators.txt"), what = "character", quiet = TRUE)
regulators <- sort(intersect(pango, rownames(tcga_expression)))
writeLines(regulators, file.path(processed, "aracne_regulators.txt"))

expression_scales <- data.table(
  cohort = c("TCGA_LUAD_LUSC", "GSE81089", "PDMR_PDX", "DepMap_22Q2"),
  file = c(
    "tcga_aracne_expression.tsv", "gse81089_aracne_expression.tsv",
    "pdmr_viper_expression.tsv.gz", "depmap_viper_expression.tsv.gz"
  ),
  genes = c(nrow(tcga_expression), nrow(gse_expression), nrow(pdmr_expression), nrow(depmap_expression)),
  samples = c(ncol(tcga_expression), ncol(gse_expression), ncol(pdmr_expression), ncol(depmap_expression)),
  scale = c(
    "raw_TPM", "log2_FPKM_plus_1", "RSEM_TPM",
    "TPM_inverted_from_EH7556_log2_TPM_plus_1"
  )
)
fwrite(expression_scales, file.path(audit, "expression_scale_receipt.tsv"), sep = "\t")

counts_by <- function(manifest) as.list(table(manifest$group))
receipt <- list(
  status = "passed",
  tcga = list(samples = ncol(tcga_expression), genes = nrow(tcga_expression), groups = counts_by(tcga_manifest)),
  gse81089 = list(samples = ncol(gse_expression), genes = nrow(gse_expression), groups = counts_by(gse_manifest)),
  pdmr = list(samples = ncol(pdmr_expression), genes = nrow(pdmr_expression), expression_scale = "RSEM_TPM"),
  depmap = list(
    metadata_eligible_samples = length(depmap_ids_with_expression) + length(missing_depmap_expression),
    samples_with_TPM = ncol(depmap_expression), genes = nrow(depmap_expression), expression_scale = "TPM_inverted_from_EH7556_log2_TPM_plus_1",
    groups = counts_by(depmap_manifest), missing_TPM_sample_ids = missing_depmap_expression
  ),
  pango_symbols = length(pango),
  aracne_regulators = length(regulators),
  package_versions = list(R = as.character(getRversion()), data.table = as.character(packageVersion("data.table")))
)
write_json(receipt, file.path(audit, "expression_preparation_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

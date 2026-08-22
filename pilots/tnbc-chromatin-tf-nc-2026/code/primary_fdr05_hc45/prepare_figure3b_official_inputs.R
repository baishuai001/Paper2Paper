#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: prepare_figure3b_official_inputs.R EDGES_TSV FIGURE2_ROOT OUTPUT_DIR")
}
edge_path <- normalizePath(args[[1L]], mustWork = TRUE)
fig2 <- normalizePath(args[[2L]], mustWork = TRUE)
out <- normalizePath(args[[3L]], mustWork = FALSE)
inputs <- file.path(out, "inputs")
dir.create(inputs, recursive = TRUE, showWarnings = FALSE)

edges <- fread(edge_path)
required <- c("TF", "target", "tfmode", "direction", "target_class")
if (length(setdiff(required, names(edges)))) {
  stop("Expanded regulon edge table lacks required columns")
}
hc <- sort(unique(as.character(edges$TF)))
if (length(hc) != 45L) stop(sprintf("Expected 45 HC-TFs; found %d", length(hc)))
writeLines(hc, file.path(inputs, "High_RNA_NES.tsv"), useBytes = TRUE)

official <- unique(edges[, .(
  TR_Source = as.character(TF),
  Target_Gene = as.character(target),
  VIPER_Weight = as.numeric(tfmode)
)])
if (anyNA(official$VIPER_Weight)) stop("Non-numeric VIPER weights")
fwrite(
  official,
  file.path(inputs, "TCGA_Regulon_TF_Target_Weights_Viper.tsv"),
  sep = "\t"
)

motif_path <- file.path(
  fig2, "results", "motif_primary", "anchor_consensus_author",
  "anchor_consensus_system_tf_motif_summary.tsv"
)
motif <- fread(motif_path)
motif_supported <- unique(motif[
  TF %chin% hc & as.logical(support_fixed_original_denominator),
  .(TF, hc_group = "Motif")
], by = "TF")
motif_table <- data.frame(
  hc_group = motif_supported$hc_group,
  row.names = motif_supported$TF,
  check.names = FALSE
)
write.table(
  motif_table,
  file.path(inputs, "JASPAR_CISBP_HC_TF_motif.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

receipt <- list(
  status = "PASS",
  HC_TFs = length(hc),
  regulon_edges = nrow(official),
  unique_targets = uniqueN(official$Target_Gene),
  any_system_motif_supported_HC_TFs = nrow(motif_supported),
  figure3B_minimum_shared_HC_TFs = 2L,
  input_edge_md5 = unname(tools::md5sum(edge_path))
)
write_json(receipt, file.path(out, "figure3b_input_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(receipt, auto_unbox = TRUE), "\n")

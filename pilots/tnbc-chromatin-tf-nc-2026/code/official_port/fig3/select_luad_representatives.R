#!/usr/bin/env Rscript

# Data-only semantic adapter for the official Figure 3G constructor. The TNBC
# script hard-codes four TNBC examples; the LUAD port selects four TFs from the
# complete Figure 3D edge table by degree, with alphabetical tie-breaking.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: select_luad_representatives.R FIGURE3D_EDGE_TSV OUTPUT_TXT RECEIPT_TSV",
    call. = FALSE
  )
}

edges <- fread(args[[1L]])
required <- c("from", "to", "weight")
if (length(setdiff(required, names(edges)))) {
  stop(sprintf("Official Figure 3D edge table lacks: %s",
               paste(setdiff(required, names(edges)), collapse = ", ")), call. = FALSE)
}
if (!nrow(edges) || any(edges$weight < 1L)) {
  stop("Official Figure 3D edge table is empty or contains non-positive weights.", call. = FALSE)
}

degree <- rbind(
  edges[, .(TF = from, partner = to)],
  edges[, .(TF = to, partner = from)]
)[, .(degree = uniqueN(partner)), by = TF]
setorder(degree, -degree, TF)
representatives <- degree[seq_len(min(4L, .N))]
if (nrow(representatives) != 4L) {
  stop(sprintf("Expected four connected LUAD TF representatives; found %d.", nrow(representatives)), call. = FALSE)
}

writeLines(representatives$TF, args[[2L]], useBytes = TRUE)
fwrite(representatives[, rank := seq_len(.N)], args[[3L]], sep = "\t")
message("LUAD Figure 3G representatives: ", paste(representatives$TF, collapse = ", "))

# Compare motif enrichment across TCGA, PDX and CellLines
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

## Inputs
parser$add_argument("-category", "--category", required=TRUE,
    help="Category (e.g., JASPAR, CISBP)")   
parser$add_argument("-compiled_motif_enrichment_JASPAR", "--compiled_motif_enrichment_JASPAR",
    default = "/data/Figure2/Motif/Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
    help="Compiled_motif enrichment")
parser$add_argument("-compiled_motif_enrichment_CISBP", "--compiled_motif_enrichment_CISBP",
    default = "/data/Figure2/Motif/Combined_MotifEnrichment/CISBP_Combined_MotifEnrichment_TNBC.tsv",
    help="Compiled_motif enrichment")
parser$add_argument("-JASPAR_CISBP_compare", "--JASPAR_CISBP_compare",
    default="/data/Figure2/Motif/JASPAR_CISBP_comparison_long.csv",
    help="Path to the file comparing JASPAR and CISBP TRs with motifs")
parser$add_argument("-TCGA_METABRIC_TR", "--TCGA_METABRIC_TR",
  default="/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
  help="Path to the TCGA and METABRIC TF file")
parser$add_argument("-HC_TR", "--HC_TR",
  default="/data/Figure2/High_RNA_NES.tsv",
  help="Path to the HC-TF file")

## Outputs
parser$add_argument("-data_outdir", "--data_outdir",
    default="/data/Figure2/Motif/",
    help="Path to output plot summarizing motif enrichment")
parser$add_argument("-visuals_outdir", "--visuals_outdir",
    default="/results/visuals/Figure2/",
    help="Path to output plot summarizing motif enrichment")

args <- parser$parse_args()
print(args)
##################################################################################################
#### For Testing
#args <- list(category = "JASPAR",
#    compiled_motif_enrichment_JASPAR = "data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/#JASPAR_Combined_MotifEnrichment_TNBC.tsv",
#    compiled_motif_enrichment_CISBP = "data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/#CISBP_Combined_MotifEnrichment_TNBC.tsv",
#    JASPAR_CISBP_compare = "data/Figure2_TR_Chromatin/Motif/JASPAR_CISBP_comparison_long.csv",
#    TCGA_METABRIC_TR = "data/Fig1_TRNetwork/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
#    HC_TR = "code/Brain_Dump/High_RNA_NES.tsv",
#    data_outdir = "data/Figure2_TR_Chromatin/Motif/",
#    visuals_outdir = "visuals/Figure2_TR_Chromatin/")

#print(args)
##################################################################################################
dir.create(args$data_outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(args$visuals_outdir, showWarnings=FALSE, recursive=TRUE)

# Read in the TCGA METABRIC to extract the TNBC-TR
TCGA_METABRIC_TR <- read.table(args$TCGA_METABRIC_TR, header=TRUE, sep='\t')
TCGA_METABRIC_TR <- TCGA_METABRIC_TR[TCGA_METABRIC_TR$Group %in% c("TNBC_TCGA_Only", "TNBC_METABRIC_TCGA_Shared"),]

# Read in the JASPAR vs. CISBP motifs called
JASPAR_CISBP_compare <- read.csv(args$JASPAR_CISBP_compare, header=TRUE, sep=',')
dim(JASPAR_CISBP_compare)

# Read in the HC-TRs
HC_TR <- read.table(args$HC_TR, header=FALSE, sep='\t')
print(head(HC_TR))
dim(HC_TR)

# Cutoffs
TCGA_Cutoff <- ceiling(7/2)
PDX_Cutoff <- ceiling(38/2)
CellLine_Cutoff <- ceiling(26/2)

# Read in the compiled motif enrichment
compiled_motif_enrich_JASPAR <- read.table(args$compiled_motif_enrichment_JASPAR, header=TRUE, sep='\t')
compiled_motif_enrich_CISBP <- read.table(args$compiled_motif_enrichment_CISBP, header=TRUE, sep='\t')
compiled_motif_enrich_CISBP$Motif.Name <- toupper(compiled_motif_enrich_CISBP$Motif.Name)
compiled_motif_enrich_CISBP_sub <- compiled_motif_enrich_CISBP[compiled_motif_enrich_CISBP$Motif.Name %in% 
  JASPAR_CISBP_compare[JASPAR_CISBP_compare$Database == "CISBP_only", "TR"],]

compiled_motif_enrich <- rbind(compiled_motif_enrich_JASPAR, compiled_motif_enrich_CISBP_sub)
compiled_motif_enrich$Motif.Name <- toupper(compiled_motif_enrich$Motif.Name)
compiled_motif_enrich <- compiled_motif_enrich[compiled_motif_enrich$Motif.Name %in% TCGA_METABRIC_TR$TR,]
compiled_motif_enrich <- compiled_motif_enrich[compiled_motif_enrich$Motif.Name %in% JASPAR_CISBP_compare$TR,]

setdiff(compiled_motif_enrich$Motif.Name, JASPAR_CISBP_compare$TR)

#### TCGA ####
TCGA_compiled_motif_enrich <- compiled_motif_enrich[compiled_motif_enrich$Sample_Group == "TCGA",]
TCGA_compiled_motif_enrich <- TCGA_compiled_motif_enrich[TCGA_compiled_motif_enrich$qvalue_Benjamini <= 0.00001,]
TCGA_freq_table <- table(TCGA_compiled_motif_enrich$Motif.Name, TCGA_compiled_motif_enrich$Sample)
TCGA_freq_table[TCGA_freq_table > 1] <- 1
TCGA_freq_table <- as.data.frame(TCGA_freq_table)

TCGA_freq_table <- TCGA_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq
  )

TCGA_freq_table$Var1 <- toupper(TCGA_freq_table$Var1)
TCGA_freq_table <- TCGA_freq_table[!duplicated(TCGA_freq_table$Var1),]

TCGA_freq_table_samp <- TCGA_freq_table

TCGA_freq_table$row_sum <- rowSums(as.matrix(TCGA_freq_table_samp[ , -1, drop = FALSE]))
TCGA_freq_table <- TCGA_freq_table[TCGA_freq_table$row_sum >= TCGA_Cutoff,]

#### PDX ####
PDX_compiled_motif_enrich <- compiled_motif_enrich[compiled_motif_enrich$Sample_Group == "PDX",]
PDX_compiled_motif_enrich <- PDX_compiled_motif_enrich[PDX_compiled_motif_enrich$qvalue_Benjamini <= 0.00001,]
PDX_freq_table <- table(PDX_compiled_motif_enrich$Motif.Name, PDX_compiled_motif_enrich$Sample)
PDX_freq_table[PDX_freq_table > 1] <- 1
PDX_freq_table <- as.data.frame(PDX_freq_table)

PDX_freq_table <- PDX_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq
  )

PDX_freq_table$Var1 <- toupper(PDX_freq_table$Var1)
PDX_freq_table <- PDX_freq_table[!duplicated(PDX_freq_table$Var1),]

PDX_freq_table_samp <- PDX_freq_table
colnames(PDX_freq_table_samp)[!grepl("Var1", colnames(PDX_freq_table_samp))] <- paste0("PDX_", colnames(PDX_freq_table_samp)[!grepl("Var1", colnames(PDX_freq_table_samp))])

PDX_freq_table$row_sum <- rowSums(as.matrix(PDX_freq_table[ , -1, drop = FALSE]))
PDX_freq_table <- PDX_freq_table[PDX_freq_table$row_sum >= PDX_Cutoff,]

#### CellLines ####
CellLines_compiled_motif_enrich <- compiled_motif_enrich[compiled_motif_enrich$Sample_Group == "CellLines",]
CellLines_compiled_motif_enrich <- CellLines_compiled_motif_enrich[CellLines_compiled_motif_enrich$qvalue_Benjamini <= 0.00001,]
CellLines_freq_table <- table(CellLines_compiled_motif_enrich$Motif.Name, CellLines_compiled_motif_enrich$Sample)
CellLines_freq_table[CellLines_freq_table > 1] <- 1
CellLines_freq_table <- as.data.frame(CellLines_freq_table)

CellLines_freq_table <- CellLines_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq
  )

CellLines_freq_table$Var1 <- toupper(CellLines_freq_table$Var1)
CellLines_freq_table <- CellLines_freq_table[!duplicated(CellLines_freq_table$Var1),]

CellLines_freq_table_samp <- CellLines_freq_table
colnames(CellLines_freq_table_samp)[!grepl("Var1", colnames(CellLines_freq_table_samp))] <- paste0("CellLines_", colnames(CellLines_freq_table_samp)[!grepl("Var1", colnames(CellLines_freq_table_samp))])

CellLines_freq_table$row_sum <- rowSums(as.matrix(CellLines_freq_table[ , -1, drop = FALSE]))
CellLines_freq_table <- CellLines_freq_table[CellLines_freq_table$row_sum >= CellLine_Cutoff,]

# Compile the TCGA, PDX and CellLine freq tables
freq_table_combined <- merge(TCGA_freq_table_samp, PDX_freq_table_samp, by = "Var1", all=TRUE)
freq_table_combined <- merge(freq_table_combined, CellLines_freq_table_samp, by = "Var1", all=TRUE)

colnames(freq_table_combined)[grepl("Var1", colnames(freq_table_combined))] <- "TNBC_TR"
freq_table_combined[is.na(freq_table_combined)] <- 0

not_acc <- setdiff(TCGA_METABRIC_TR$TR, freq_table_combined$TNBC_TR)
intersect(TCGA_METABRIC_TR$TR, freq_table_combined$TNBC_TR)

#write.table(freq_table_combined,
#  paste0("data/Figure2_TR_Chromatin/SupplementaryTables/MotifEnrichment_Across_Samples.tsv"),
#  col.names=TRUE,
#  row.names=FALSE,
#  sep='\t',
#  quote=FALSE)

# Now compare across TCGA, PDX and CellLines
# Create a grouping column: TRUE if intersection contains HC_TR
TR_list <- list(TCGA = TCGA_freq_table$Var1,
  PDX = PDX_freq_table$Var1,
  CellLines = CellLines_freq_table$Var1) 
TR_list_list_mat <- ComplexHeatmap::list_to_matrix(TR_list)

TR_list_list_mat <- as.data.frame(TR_list_list_mat)

TR_list_list_mat$hc_group <- ifelse(
  rownames(TR_list_list_mat) %in% JASPAR_CISBP_compare[JASPAR_CISBP_compare$Database == "CISBP_only", "TR"], 
  "CISBP only",
  ifelse(
    rownames(TR_list_list_mat) %in% JASPAR_CISBP_compare[JASPAR_CISBP_compare$Database == "JASPAR_only", "TR"],
    "JASPAR only",
    ifelse(
      rownames(TR_list_list_mat) %in% JASPAR_CISBP_compare[JASPAR_CISBP_compare$Database == "Both", "TR"], 
      "Both",
      "Not HC-TR"
    )
  )
)

TR_list_list_mat$hc_group <- factor(TR_list_list_mat$hc_group,
  levels=c("Not HC-TR", "Both", "CISBP only", "JASPAR only"))

yaxis_label <- paste0("Number of TRs", "\n", "(n = ", nrow(TR_list_list_mat), "/", nrow(JASPAR_CISBP_compare), ")")
cat(yaxis_label)

# Build annotation as a named list using setNames() so the variable value is used as the name
annotation_list <- setNames(
  list(
    ggplot(mapping = aes(x = intersection)) +
      geom_bar(aes(fill = hc_group), stat = "count", position = "stack") +
      geom_text(
        stat = "count",
        aes(label = after_stat(count), group = hc_group),
        position = position_stack(vjust = 0.5),
        color = "black",
        size = 2
      ) +
      scale_fill_manual(
        values = HC_TR_Motif,
        name   = "HC-TR present",
        labels = c("CISBP & JASPAR", "CISBP only", "JASPAR only")
      ) +
      theme_classic() +
      theme(
        axis.text.y  = element_text(colour = "black"),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 0, r = 0, b = 0, l = 15, unit = "mm")
      )
  ),
  yaxis_label  # variable value becomes the list name = y-axis label in ComplexUpset
)

TR_upset <- ComplexUpset::upset(
  data      = TR_list_list_mat,
  intersect = c("TCGA", "PDX", "CellLines"),
  set_sizes = FALSE,
  queries   = list(
    ComplexUpset::upset_query("TCGA",      color = cohort_colours["TCGA"]),
    ComplexUpset::upset_query("PDX",       color = cohort_colours["PDX"]),
    ComplexUpset::upset_query("CellLines", color = cohort_colours["CellLines"])
  ),
  name   = "TRs with genome-wide motif\nenrichment in at least 50% of samples",
  themes = theme(text = element_text(size = 10, colour = "black")),
  base_annotations = annotation_list
)

#write.table(TR_list_list_mat, paste0(args$data_outdir, "/", args$category, "_TR_MotifEnrichment_Intersection_Sample_Groups.tsv"),
#    col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

cairo_pdf(paste0(args$visuals_outdir, "/Figure2D_", args$category, "_TR_MotifEnrichment_Intersection_Sample_Groups.pdf"))
print(TR_upset)
dev.off()

sessionInfo()
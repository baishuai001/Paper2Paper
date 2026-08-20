# Compare promoter accessibility across cell-lines, PDX and primary tumours
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

## Inputs
parser$add_argument("-compiled_promoter_accessibility", "--compiled_promoter_accessibility",
    default = "/data/Figure2/Promoter_Accessibility/Combined_Promoter_Accessibility_TNBC_TR_Genes.tsv",
    help="Compiled promoter accessibility")              
parser$add_argument("-TCGA_METABRIC_TR", "-TCGA_METABRIC_TR",
    default = "/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
    help="Path to the TCGA and METABRIC TF file")

## Outputs
parser$add_argument("-data_outdir", "--data_outdir",
    default="/data/Figure2/Promoter_Accessibility/",
    help="Path to output plot summarizing genomic annotation")
parser$add_argument("-visuals_outdir", "--visuals_outdir",
    default="/results/visuals/Figure2/",
    help="Path to output plot summarizing genomic annotation")

args <- parser$parse_args()
print(args)
##################################################################################################
dir.create(args$data_outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(args$visuals_outdir, showWarnings=FALSE, recursive=TRUE)

# TCGA_METABRIC_TR
TCGA_METABRIC_TR <- read.table(args$TCGA_METABRIC_TR, header=TRUE, sep='\t')
TCGA_METABRIC_TR <- TCGA_METABRIC_TR[TCGA_METABRIC_TR$Group %in% c("TNBC_TCGA_Only", "TNBC_METABRIC_TCGA_Shared"),]

# Cutoffs
TCGA_Cutoff <- ceiling(7/2)
PDX_Cutoff <- ceiling(38/2)
CellLine_Cutoff <- ceiling(26/2)

# Read in the compiled motif enrichment
compiled_promoter_accessibility <- read.table(args$compiled_promoter_accessibility, header=TRUE, sep='\t')
print(head(compiled_promoter_accessibility))

#### TCGA ####
TCGA_compiled_promoter_accessibility <- compiled_promoter_accessibility[compiled_promoter_accessibility$Sample_Group == "TCGA",]
TCGA_freq_table <- table(TCGA_compiled_promoter_accessibility$geneID, TCGA_compiled_promoter_accessibility$Sample)
TCGA_freq_table[TCGA_freq_table > 1] <- 1
TCGA_freq_table <- as.data.frame(TCGA_freq_table)

TCGA_freq_table <- TCGA_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq)

TCGA_freq_table$Var1 <- toupper(TCGA_freq_table$Var1)
TCGA_freq_table_samp <- TCGA_freq_table
colnames(TCGA_freq_table_samp)[!grepl("Var1", colnames(TCGA_freq_table_samp))] <- paste0("TCGA_", colnames(TCGA_freq_table_samp)[!grepl("Var1", colnames(TCGA_freq_table_samp))])

TCGA_freq_table$row_sum <- rowSums(as.matrix(TCGA_freq_table[ , -1, drop = FALSE]))
TCGA_freq_table <- TCGA_freq_table[TCGA_freq_table$row_sum >= TCGA_Cutoff,]

#### PDX ####
PDX_compiled_promoter_accessibility <- compiled_promoter_accessibility[compiled_promoter_accessibility$Sample_Group == "PDX",]
PDX_freq_table <- table(PDX_compiled_promoter_accessibility$geneID, PDX_compiled_promoter_accessibility$Sample)
PDX_freq_table[PDX_freq_table > 1] <- 1
PDX_freq_table <- as.data.frame(PDX_freq_table)

PDX_freq_table <- PDX_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq
  )

PDX_freq_table$Var1 <- toupper(PDX_freq_table$Var1)
PDX_freq_table_samp <- PDX_freq_table
colnames(PDX_freq_table_samp)[!grepl("Var1", colnames(PDX_freq_table_samp))] <- paste0("PDX_", colnames(PDX_freq_table_samp)[!grepl("Var1", colnames(PDX_freq_table_samp))])

PDX_freq_table$row_sum <- rowSums(as.matrix(PDX_freq_table[ , -1, drop = FALSE]))
PDX_freq_table <- PDX_freq_table[PDX_freq_table$row_sum >= PDX_Cutoff,]

#### CellLines ####
CellLines_compiled_promoter_accessibility <- compiled_promoter_accessibility[compiled_promoter_accessibility$Sample_Group == "CellLines",]
CellLines_freq_table <- table(CellLines_compiled_promoter_accessibility$geneID, CellLines_compiled_promoter_accessibility$Sample)
CellLines_freq_table[CellLines_freq_table > 1] <- 1
CellLines_freq_table <- as.data.frame(CellLines_freq_table)

CellLines_freq_table <- CellLines_freq_table %>%
  tidyr::pivot_wider(
    names_from  = Var2,
    values_from = Freq
  )

CellLines_freq_table$Var1 <- toupper(CellLines_freq_table$Var1)
CellLines_freq_table_samp <- CellLines_freq_table
colnames(CellLines_freq_table_samp)[!grepl("Var1", colnames(CellLines_freq_table_samp))] <- paste0("CellLine_", colnames(CellLines_freq_table_samp)[!grepl("Var1", colnames(CellLines_freq_table_samp))])

CellLines_freq_table$row_sum <- rowSums(as.matrix(CellLines_freq_table[ , -1, drop = FALSE]))
CellLines_freq_table <- CellLines_freq_table[CellLines_freq_table$row_sum >= CellLine_Cutoff,]

# Compile to a matrix
TCGA_PDX_CellLine_freq_table <- merge(TCGA_freq_table_samp, PDX_freq_table_samp, by='Var1', all=TRUE)
TCGA_PDX_CellLine_freq_table <- merge(TCGA_PDX_CellLine_freq_table, CellLines_freq_table_samp, by='Var1', all=TRUE)
colnames(TCGA_PDX_CellLine_freq_table)[grepl("Var1", colnames(TCGA_PDX_CellLine_freq_table))] <- "TNBC_TR"

diff_TR <- setdiff(TCGA_METABRIC_TR$TR, TCGA_PDX_CellLine_freq_table$TNBC_TR)

new_rows <- as.data.frame(
  matrix(
    0,
    nrow = length(diff_TR),
    ncol = ncol(TCGA_PDX_CellLine_freq_table)
  )
)

colnames(new_rows) <- colnames(TCGA_PDX_CellLine_freq_table)

new_rows$TNBC_TR <- diff_TR

TCGA_PDX_CellLine_freq_table <- rbind(TCGA_PDX_CellLine_freq_table, new_rows)
TCGA_PDX_CellLine_freq_table <- TCGA_PDX_CellLine_freq_table %>%
  arrange(TNBC_TR)

#write.table(TCGA_PDX_CellLine_freq_table,
#  paste0("data/Figure2_TR_Chromatin/SupplementaryTables/Promoter_Accessibility_Across_Samples.tsv"),
#  col.names=TRUE,
#  row.names=FALSE,
#  sep='\t',
#  quote=FALSE)

# Now compare across TCGA, PDX and CellLines
TR_list <- list(TCGA = TCGA_freq_table$Var1,
    PDX = PDX_freq_table$Var1,
    CellLines = CellLines_freq_table$Var1)

TR_list_list_mat <- ComplexHeatmap::list_to_matrix(TR_list)
TR_list_list_mat <- TR_list_list_mat[,c("TCGA", "PDX", "CellLines")]

print(head(as.data.frame(TR_list_list_mat)))

# Plot the UpSet plot
TR_upset <- ComplexUpset::upset(
  data      = as.data.frame(TR_list_list_mat),
  intersect = c("TCGA", "PDX", "CellLines"),
  set_sizes = FALSE,  # Add this line to remove set size bars
  queries = list(
    ComplexUpset::upset_query("TCGA",     color = cohort_colours["TCGA"]),
    ComplexUpset::upset_query("PDX",      color = cohort_colours["PDX"]),
    ComplexUpset::upset_query("CellLines",color = cohort_colours["CellLines"])
  ),
  name   = "TRs with promoter accessibility\nin at least 50% of samples",
  #themes = upset_default_themes(text = element_text(size = 10, colour = "black")),
  themes = theme(text = element_text(size = 10, colour = "black")),
  base_annotations = list(
    "Number of TRs" = ComplexUpset::intersection_size() +
      labs(y = "Number of TRs") +
      theme_classic() +
      theme(
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  )
)

TR_upset <- TR_upset +
  theme(
    text           = element_text(colour = "black"),
    axis.text      = element_text(colour = "black"),
    axis.title     = element_text(colour = "black"),
    strip.text     = element_text(colour = "black"),
    legend.text    = element_text(colour = "black"),
    legend.title   = element_text(colour = "black"),
    plot.title     = element_text(colour = "black")
  )

#write.table(TR_list_list_mat, paste0(args$data_outdir, "/TR_PromoterAccessibility_Intersection_Sample_Groups.tsv"), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

cairo_pdf(paste0(args$visuals_outdir, "/Figure2B_TF_PromoterAccessibility_Intersection_Sample_Groups.pdf"))
print(TR_upset)
dev.off()

sessionInfo()
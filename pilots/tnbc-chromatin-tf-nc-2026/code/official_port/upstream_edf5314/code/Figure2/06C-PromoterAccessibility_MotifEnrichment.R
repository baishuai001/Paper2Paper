# Plot heatmap of promoter accessibility and motif enrichment for TNBC-specific TRs
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

## Inputs     
parser$add_argument("-compiled_motif_enrichment_JASPAR", "--compiled_motif_enrichment_JASPAR",
    default="/data/Figure2/Motif/Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
    help="Path to directory containing the combined motif enrichment for JASPAR")
parser$add_argument("-compiled_motif_enrichment_CISBP", "--compiled_motif_enrichment_CISBP",
    default="/data/Figure2/Motif/Combined_MotifEnrichment/CISBP_Combined_MotifEnrichment_TNBC.tsv",
    help="Path to directory containing the combined motif enrichment for CISBP")
parser$add_argument("-compiled_promoter_accessibility", "--compiled_promoter_accessibility",
    default="/data/Figure2/Promoter_Accessibility/Combined_Promoter_Accessibility_TNBC_TR_Genes.tsv",
    help="Path to the compiled promoter accessibility file")
parser$add_argument("-TCGA_RNA_NES", "--TCGA_RNA_NES",
    default="/data/Figure1/TCGA/TCGA_TR_activities_viper.tsv",
    help="Path to the NES for TCGA RNA-seq data")
parser$add_argument("-PDX_RNA_NES", "--PDX_RNA_NES",
    default="/data/Figure1/PDX/PDX73_TR_activities_viper_usingTCGAregulons.tsv",
    help='Path to the NES for PDX RNA-seq data')
parser$add_argument("-PDX_metadata", "--PDX_metadata",
    default="/data/Reference_Data/PDX/PDX_molecular_subtype.tsv",
    help='Path to the metadata for the PDXs')
parser$add_argument("-CellLines_RNA_NES", "--CellLines_RNA_NES",
    default="/data/Figure1/CellLines/CellLines82_TR_activities_viper_usingTCGAregulons.tsv",
    help='Path to the NES for CellLines RNA-seq data')
parser$add_argument("-CellLines_metadata", "--CellLines_metadata",
    default="/data/Reference_Data/CellLines_Marcotte_et_al_2016/CellLines_molecular_subtype.tsv",
    help='Path to the meta data for cell-lines')
parser$add_argument("-METABRIC_RNA_NES", "--METABRIC_RNA_NES",
    default="/data/Figure1/METABRIC/METABRIC_TR_activities_viper.tsv",
    help='Path to the NES for METABRIC RNA-seq data')
parser$add_argument("-METABRIC_metadata", "--METABRIC_metadata",
    default="/data/Reference_Data/METABRIC/METABRIC_molecular_subtype.tsv",
    help="Path to the molecular subtyping information for METABRIC")
parser$add_argument("-TR_file", "--TR_file", 
    default="/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
    help="Path to the file with the TNBC-specific TRs")
parser$add_argument("-HC_TRs", "--HC_TRs",
    default="/data/Figure2/High_RNA_NES.tsv",
    help="Path to the HC-TR file")
parser$add_argument("-UpSet_Cohort", "--UpSet_Cohort",
    default="/data/Figure2/Motif/JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv",
    help="Path to the UpSet plot comparision across TCGA, PDX and CellLine cohorts")

## Outputs
parser$add_argument("-visuals_outdir", "--visuals_outdir",
    default="/results/visuals/Figure2",
    help="Path to output plot")

args <- parser$parse_args()
print(args)
##################################################################################################
#### For Testing
#args <- list(compiled_motif_enrichment_JASPAR = "data/Figure2_TR_Chromatin/Motif/#Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
#    compiled_motif_enrichment_CISBP = "data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/#CISBP_Combined_MotifEnrichment_TNBC.tsv",
#    compiled_promoter_accessibility = "data/Figure2_TR_Chromatin/Promoter_Accessibility/#Combined_Promoter_Accessibility/Combined_Promoter_Accessibility_TNBC_TR_Genes.tsv",
#    TCGA_RNA_NES = "data/Fig1_TRNetwork/TCGA/TCGA_TR_activities_viper.tsv",
#    PDX_RNA_NES = "data/Fig1_TRNetwork/PDX/PDX73_TR_activities_viper_usingTCGAregulons.tsv",
#    PDX_metadata = "data/Reference_Data/PDX/PDX_molecular_subtype.tsv",
#    CellLines_RNA_NES = "data/Fig1_TRNetwork/CellLines/#CellLines82_TR_activities_viper_usingTCGAregulons.tsv",
#    CellLines_metadata = "data/Reference_Data/CellLines_Marcotte_et_al_2016/#CellLines_molecular_subtype.tsv",
#    METABRIC_RNA_NES = "data/Fig1_TxnalNetwork/METABRIC/METABRIC_TR_activities_viper.tsv",
#    METABRIC_metadata = "data/Reference_Data/METABRIC/METABRIC_molecular_subtype.tsv",
#    TR_file = "data/Fig1_TRNetwork/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
#    HC_TRs = "code/Brain_Dump/High_RNA_NES.tsv",
#    UpSet_Cohort = "data/Figure2_TR_Chromatin/Motif/#JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv",
#    visuals_outdir = "visuals/Figure2_TR_Chromatin")

#print(args)
##################################################################################################
dir.create(args$visuals_outdir, showWarnings=FALSE, recursive=TRUE)

# Read in the TNBC TR information
TR_file <- read.table(args$TR_file, header=TRUE, sep='\t')
TR_file$TR <- toupper(TR_file$TR)
TR_file <- TR_file[TR_file$Group %in% c("TNBC_TCGA_Only", "TNBC_METABRIC_TCGA_Shared"),]
colnames(TR_file) <- c("TR_Category", "geneID")
print(head(TR_file))
print(dim(TR_file))

# Read in UpSet plot motif enrichment
UpSet_Cohort <- read.table(args$UpSet_Cohort, header=TRUE, sep='\t')

# Read in the promoter accessibility
compiled_promoter_accessibility <- read.table(args$compiled_promoter_accessibility, header=TRUE, sep='\t') 
compiled_promoter_accessibility <- compiled_promoter_accessibility[,c("geneID", "Sample", "Sample_Group")]
compiled_promoter_accessibility <- compiled_promoter_accessibility[!duplicated(compiled_promoter_accessibility),]
compiled_promoter_accessibility$geneID <- toupper(compiled_promoter_accessibility$geneID)
compiled_promoter_accessibility$Promoter_Accessibility <- compiled_promoter_accessibility$Sample_Group

# Compiled motif enrichment
compiled_motif_enrichment_JASPAR <- read.table(args$compiled_motif_enrichment_JASPAR, header=TRUE, sep='\t')
compiled_motif_enrichment_JASPAR$Motif.Name <- toupper(compiled_motif_enrichment_JASPAR$Motif.Name)
compiled_motif_enrichment_JASPAR <- compiled_motif_enrichment_JASPAR[compiled_motif_enrichment_JASPAR$Motif.Name 
  %in% TR_file$geneID,]

compiled_motif_enrichment_CISBP <- read.table(args$compiled_motif_enrichment_CISBP, header=TRUE, sep='\t')
compiled_motif_enrichment_CISBP$Motif.Name <- toupper(compiled_motif_enrichment_CISBP$Motif.Name)
compiled_motif_enrichment_CISBP <- compiled_motif_enrichment_CISBP[compiled_motif_enrichment_CISBP$Motif.Name
  %in% TR_file$geneID,]

setdiff <- setdiff(TR_file$geneID, compiled_motif_enrichment_JASPAR$Motif.Name) 
compiled_motif_enrichment_CISBP_inclusion <- compiled_motif_enrichment_CISBP[compiled_motif_enrichment_CISBP$Motif.Name
  %in% setdiff,]

compiled_motif_enrichment <- rbind(compiled_motif_enrichment_JASPAR, compiled_motif_enrichment_CISBP_inclusion)
compiled_motif_enrichment <- compiled_motif_enrichment[,c("Motif.Name", "qvalue_Benjamini", "Sample", "Sample_Group")]
compiled_motif_enrichment$qvalue_Benjamini <- ifelse(compiled_motif_enrichment$qvalue_Benjamini <= 0.00001, "1", "0")

compiled_motif_enrichment <- compiled_motif_enrichment %>%
  mutate(
    Motif_enrichment_qvalue_Benjamini_Category =
      ifelse(qvalue_Benjamini == 1, Sample_Group, "NA")
  )

compiled_motif_enrichment$Motif.Name <- toupper(compiled_motif_enrichment$Motif.Name)
colnames(compiled_motif_enrichment) <- c("geneID", "Motif_enrichment_qvalue_Benjamini", "Sample", "Sample_Group", "Motif_enrichment_qvalue_Benjamini_Category")
summary(compiled_motif_enrichment)

diff <- setdiff(TR_file$geneID, compiled_motif_enrichment$geneID)

# Merge the information
merged_promoter_accessibility <- merge(
  TR_file,
  compiled_promoter_accessibility,
  by = c("geneID"),
  all = TRUE
)

merged_motif_enrichment <- merge(
  TR_file,
  compiled_motif_enrichment,
  by = c("geneID"),
  all = TRUE
)

merged_promoter_accessibility[is.na(merged_promoter_accessibility)] <- "NA"
merged_promoter_accessibility <- merged_promoter_accessibility[merged_promoter_accessibility$Sample_Group != "NA",]

merged_motif_enrichment[is.na(merged_motif_enrichment)] <- "NA"
merged_motif_enrichment <- merged_motif_enrichment[merged_motif_enrichment$Sample_Group != "NA",]

## NES RNAseq
#### TCGA ####
TCGA_RNA_NES <- read.table(args$TCGA_RNA_NES, header=TRUE, sep='\t')
TNBC_samples <- gsub("-", ".", unique(compiled_promoter_accessibility[compiled_promoter_accessibility$Sample_Group == "TCGA", "Sample"]))
TNBC_samples <- gsub("01A", "01", TNBC_samples)
TCGA_RNA_NES <- TCGA_RNA_NES[,TNBC_samples]
TCGA_RNA_NES$Avg_NES_RNA <- rowMeans(TCGA_RNA_NES)
TCGA_RNA_NES$geneID <- toupper(rownames(TCGA_RNA_NES))
TCGA_RNA_NES <- TCGA_RNA_NES[,c("geneID", "Avg_NES_RNA")]
TCGA_RNA_NES$Sample_Group <- "TCGA"
rownames(TCGA_RNA_NES) <- NULL

#### METABRIC ####
METABRIC_RNA_NES <- read.table(args$METABRIC_RNA_NES, header=TRUE, sep='\t', check.names=FALSE)
METABRIC_metadata <- read.table(args$METABRIC_metadata, header=TRUE, sep='\t')
METABRIC_metadata <- METABRIC_metadata[METABRIC_metadata$PAM50 == "Basal" & METABRIC_metadata$SCMOD2 == "Basal",]
METABRIC_RNA_NES <- METABRIC_RNA_NES[,METABRIC_metadata$Sample]
METABRIC_RNA_NES$Avg_NES_RNA <- rowMeans(METABRIC_RNA_NES)
METABRIC_RNA_NES$geneID <- toupper(rownames(METABRIC_RNA_NES))
METABRIC_RNA_NES <- METABRIC_RNA_NES[,c("geneID", "Avg_NES_RNA")]
METABRIC_RNA_NES$Sample_Group <- "METABRIC"
rownames(METABRIC_RNA_NES) <- NULL

#### PDX ####
PDX_RNA_NES <- read.table(args$PDX_RNA_NES, header=TRUE, sep='\t', check.names=FALSE)
PDX_metadata <- read.table(args$PDX_metadata, header=TRUE, sep='\t')
PDX_metadata <- PDX_metadata[PDX_metadata$PAM50 == "Basal" | PDX_metadata$SCMOD2 == "Basal",]
PDX_RNA_NES <- PDX_RNA_NES[,c(PDX_metadata$Sample)]
PDX_RNA_NES$Avg_NES_RNA <- rowMeans(PDX_RNA_NES)
PDX_RNA_NES$geneID <- toupper(rownames(PDX_RNA_NES))
PDX_RNA_NES <- PDX_RNA_NES[,c("geneID", "Avg_NES_RNA")]
PDX_RNA_NES$Sample_Group <- "PDX"
rownames(PDX_RNA_NES) <- NULL

#### CellLines ####
CellLines_RNA_NES <- read.table(args$CellLines_RNA_NES, header=TRUE, sep='\t', check.names=FALSE)
CellLines_metadata <- read.table(args$CellLines_metadata, header=TRUE, sep='\t')
CellLines_metadata <- CellLines_metadata[CellLines_metadata$PAM50 == "Basal" | CellLines_metadata$SCMOD2 == "Basal",]
CellLines_RNA_NES <- CellLines_RNA_NES[,c(CellLines_metadata$Sample)]
CellLines_RNA_NES$Avg_NES_RNA <- rowMeans(CellLines_RNA_NES)
CellLines_RNA_NES$geneID <- toupper(rownames(CellLines_RNA_NES))
CellLines_RNA_NES <- CellLines_RNA_NES[,c("geneID", "Avg_NES_RNA")]
CellLines_RNA_NES$Sample_Group <- "CellLines"
rownames(CellLines_RNA_NES) <- NULL

###### Samples
TCGA_samples <- unique(merged_promoter_accessibility[merged_promoter_accessibility$Sample_Group == "TCGA", "Sample"])
PDX_samples <- unique(merged_promoter_accessibility[merged_promoter_accessibility$Sample_Group == "PDX", "Sample"])
CellLines_samples <- unique(merged_promoter_accessibility[merged_promoter_accessibility$Sample_Group == "CellLines", "Sample"])

## Plot the tile matrix - RNA NES
merged_RNA_NES <- rbind(TCGA_RNA_NES, PDX_RNA_NES)
merged_RNA_NES <- rbind(merged_RNA_NES, CellLines_RNA_NES)
merged_RNA_NES <- rbind(merged_RNA_NES, METABRIC_RNA_NES)

merged_RNA_NES <- merged_RNA_NES[merged_RNA_NES$geneID %in% TR_file$geneID,]

merged_RNA_NES$Sample_Group <- factor(merged_RNA_NES$Sample_Group, levels = c("TCGA", "METABRIC", "PDX", "CellLines"))
merged_RNA_NES$NES_Category <- ifelse(merged_RNA_NES$Avg_NES_RNA < 0, "<0",
  ifelse(merged_RNA_NES$Avg_NES_RNA >= 0 & merged_RNA_NES$Avg_NES_RNA < 1, "0-1",
  ifelse(merged_RNA_NES$Avg_NES_RNA >= 1 & merged_RNA_NES$Avg_NES_RNA < 3, "1-3",
  ifelse(merged_RNA_NES$Avg_NES_RNA >= 3, ">3", "NA"))))

merged_RNA_NES$NES_Category <- factor(merged_RNA_NES$NES_Category, levels=c(">3", "1-3", "0-1", "<0"))

# Determine order of TRs
TR_RNA_NES <- merged_RNA_NES[,c("geneID", "Avg_NES_RNA", "Sample_Group")]
TR_RNA_NES <- TR_RNA_NES %>%
  tidyr::pivot_wider(
    id_cols     = geneID,
    names_from  = Sample_Group,
    values_from = Avg_NES_RNA
  )
TR_RNA_NES <- as.data.frame(TR_RNA_NES)
TR_RNA_NES$Avg_NES_RNA <- rowMeans(TR_RNA_NES[,c(2:4)])

TR_RNA_NES_NoMotif <- TR_RNA_NES[TR_RNA_NES$geneID %in% diff,]
TR_RNA_NES_WithMotif <- TR_RNA_NES[!(TR_RNA_NES$geneID %in% diff),]

TR_RNA_NES_NoMotif <- TR_RNA_NES_NoMotif[order(TR_RNA_NES_NoMotif$Avg_NES_RNA, decreasing=TRUE),]
TR_RNA_NES_WithMotif <- TR_RNA_NES_WithMotif[order(TR_RNA_NES_WithMotif$Avg_NES_RNA, decreasing=TRUE),]

merged_RNA_NES$geneID <- factor(merged_RNA_NES$geneID, levels = c(rev(TR_RNA_NES_NoMotif$geneID), rev(TR_RNA_NES_WithMotif$geneID)))

HC_TRs <- read.table(args$HC_TRs, header=FALSE, sep='\t')
#HC_TR_df <- merged_RNA_NES
#HC_TR_df$TR_Group <- ifelse(merged_RNA_NES$geneID %in% HC_TRs$V1, "HC-TR", "Not HC-TR")
merged_RNA_NES$TR_Group <- ifelse(merged_RNA_NES$geneID %in% HC_TRs$V1, "HC-TR", "Not HC-TR")
HC_TR_df <- merged_RNA_NES

HC_TR_df$geneID <- factor(merged_RNA_NES$geneID, levels = c(rev(TR_RNA_NES_NoMotif$geneID), rev(TR_RNA_NES_WithMotif$geneID)))

print(head(merged_RNA_NES))

HC_TR_df <- merged_RNA_NES[,c("geneID", "TR_Group")]
HC_TR_df <- HC_TR_df[!duplicated(HC_TR_df),]

# Add in annotation bar for HC-TR and Non HC-TR
HC_TR_anno_plt <- ggplot(HC_TR_df,
       aes(x = "TR Group", y = geneID, fill = TR_Group)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, fill="TR Group") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0,
                                colour = "black", size = 8),
    axis.text.y  = element_text(size = 3),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  ) +
  scale_fill_manual(values=c("HC-TR" = "#0CABA8", "Not HC-TR" = "grey"))

# NES - Normalized Enrichment Scores
RNA_NES_plt <- ggplot(merged_RNA_NES,
       aes(x = Sample_Group, y = geneID, fill = NES_Category)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, fill="NES") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0,
                                colour = "black", size = 8),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  ) +
  scale_fill_manual(values=NES_category_colours)

#write.table(merged_RNA_NES, "data/Figure2_TR_Chromatin/Tables_for_Plotting/Combined_RNA_NES_Heatmap.tsv", col.names=TRUE, row.names=FALSE, sep='\t')

## Plot the tile matrix - Promoter Accessibility
merged_promoter_accessibility$Promoter_Accessibility <- as.character(merged_promoter_accessibility$Promoter_Accessibility)
merged_promoter_accessibility$Sample_Group <- factor(merged_promoter_accessibility$Sample_Group, levels=c("TCGA", "PDX", "CellLines"))

merged_promoter_accessibility$Sample <- factor(merged_promoter_accessibility$Sample, levels = c(TCGA_samples, PDX_samples, CellLines_samples))

merged_promoter_accessibility$geneID <- factor(merged_promoter_accessibility$geneID, levels = c(rev(TR_RNA_NES_NoMotif$geneID), rev(TR_RNA_NES_WithMotif$geneID)))

Promoter_Accessibility_plt <- ggplot(merged_promoter_accessibility,
       aes(x = Sample, y = geneID, fill = Promoter_Accessibility)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    legend.position="none"
  ) +
  scale_fill_manual(values = cohort_colours)

## Plot the tile matrix - Motif Enrichment
merged_motif_enrichment$Motif_enrichment_qvalue_Benjamini <- as.character(merged_motif_enrichment$Motif_enrichment_qvalue_Benjamini)
merged_motif_enrichment$Sample_Group <- factor(merged_motif_enrichment$Sample_Group, levels=c("TCGA", "PDX", "CellLines"))

merged_motif_enrichment$Sample <- factor(merged_motif_enrichment$Sample, levels = c(TCGA_samples, PDX_samples, CellLines_samples))

merged_motif_enrichment <- merged_motif_enrichment[!(merged_motif_enrichment$geneID %in% diff),]
merged_motif_enrichment <- merged_motif_enrichment[,c("geneID", "Sample", "Sample_Group", "Motif_enrichment_qvalue_Benjamini_Category")]

no_motif_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(no_motif_df) <- c("geneID", "Sample", "Sample_Group", "Motif_enrichment_qvalue_Benjamini_Category")

for(group in c("TCGA", "PDX", "CellLines")){
  for(sample in unique(merged_promoter_accessibility[merged_promoter_accessibility$Sample_Group == group, "Sample"])){
    df <- data.frame(geneID = diff,
      Sample = sample,
      Sample_Group = group,
      Motif_enrichment_qvalue_Benjamini_Category = "No Motif")

    no_motif_df <- rbind(no_motif_df, df)
  }
}

merged_motif_enrichment <- rbind(merged_motif_enrichment, no_motif_df)
merged_motif_enrichment$geneID <- factor(merged_motif_enrichment$geneID, levels = c(rev(TR_RNA_NES_NoMotif$geneID), rev(TR_RNA_NES_WithMotif$geneID)))

head(merged_motif_enrichment)
table(merged_motif_enrichment$Motif_enrichment_qvalue_Benjamini_Category)

merged_motif_enrichment$Motif_enrichment_qvalue_Benjamini_Category <- factor(merged_motif_enrichment$Motif_enrichment_qvalue_Benjamini_Category,
  levels = c("TCGA", "PDX", "CellLines", "NA", "No Motif"))

Motif_Enrichment_plt <- ggplot(merged_motif_enrichment,
       aes(x = Sample, y = geneID, fill = Motif_enrichment_qvalue_Benjamini_Category)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, fill=NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  ) +
  scale_fill_manual(values = cohort_colours)

## TRs with generalizability
#UpSet_Cohort$geneID <- rownames(UpSet_Cohort)
#TR_file$General <- ifelse(TR_file$geneID %in% UpSet_Cohort$geneID, "1", "NA")
#TR_file$geneID <- factor(TR_file$geneID, levels=c(rev(TR_RNA_NES_NoMotif$geneID), rev(TR_RNA_NES_WithMotif$geneID)))

#General_plt <- ggplot(TR_file,
#       aes(x = "General", y = geneID, fill = General)) +
#  geom_tile() +
#  scale_x_discrete(position = "top") +
#  labs(x = NULL, y = NULL, fill=NULL) +
#  theme_minimal(base_size = 12) +
#  theme(
#    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = "black", size = 5),
#    axis.text.y = element_text(size = 8),
#    axis.title.x = element_blank(),
#    axis.title.y = element_blank(),
#    panel.grid   = element_blank()
#  ) +
#  scale_fill_manual(values = c("1" = "navy", "NA" = "white"), guide = "none")

## Merge all the plots together
#RNA_NES_plt <- RNA_NES_plt +
#  theme(
#    axis.text.y  = element_blank(),
#    axis.ticks.y = element_blank(),
#    axis.title.y = element_blank()
#  )

Promoter_Accessibility_plt <- Promoter_Accessibility_plt +
  scale_fill_manual(
    values = cohort_colours,
    guide  = "none"   # drop ONLY this legend, keep colors [web:120]
  ) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

Motif_Enrichment_plt <- Motif_Enrichment_plt +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

combined_plt <- (HC_TR_anno_plt + RNA_NES_plt + Promoter_Accessibility_plt + Motif_Enrichment_plt) +
  plot_layout(
    widths = c(0.05, 0.1, 1, 1),
    guides  = "collect"
  ) &
  theme(legend.position = "right")

cairo_pdf(paste0(args$visuals_outdir, "/Figure2C_JASPAR_CISBP_Combined_ChromatinAccessibility_MotifEnrichment_Heatmap.pdf"), 
  width = 12, height = 8)
print(combined_plt)
dev.off()

sessionInfo()
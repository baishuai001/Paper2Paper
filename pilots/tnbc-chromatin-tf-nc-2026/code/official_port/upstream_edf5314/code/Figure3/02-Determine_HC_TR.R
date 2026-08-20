# Determine the HC-TRs based on certain criteria
library(dplyr)
library(ComplexUpset)
library(ggplot2)

dir.create("/results/visuals/Figure2/", showWarnings=FALSE, recursive=TRUE)

## Read in the 150 TNBC-TRs
TNBC_TRs <- read.table("/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
    header=TRUE, sep='\t')
TNBC_TRs <- TNBC_TRs[TNBC_TRs$Group %in% c("TNBC_TCGA_Only", "TNBC_METABRIC_TCGA_Shared"),]
dim(TNBC_TRs)

## (1) Accessible promoter - TNBC tumours
### Read in the promoter accessibility information
promoter_accessibility <- read.table("/data/Figure2/Promoter_Accessibility/Combined_Promoter_Accessibility_TNBC_TR_Genes.tsv", header=TRUE, sep='\t')

promoter_accessibility_TNBC_tumours <- promoter_accessibility[promoter_accessibility$Sample_Group == "TCGA",]

promoter_accessibility_TNBC_tumours_freq <- table(promoter_accessibility_TNBC_tumours$Sample, promoter_accessibility_TNBC_tumours$geneID)
promoter_accessibility_TNBC_tumours_freq <- t(promoter_accessibility_TNBC_tumours_freq)
promoter_accessibility_TNBC_tumours_freq <- ifelse(promoter_accessibility_TNBC_tumours_freq >= 1, 1, 0)
promoter_accessibility_TNBC_tumours_freq <- as.data.frame(promoter_accessibility_TNBC_tumours_freq)
promoter_accessibility_TNBC_tumours_freq$sum <- rowSums(promoter_accessibility_TNBC_tumours_freq)
promoter_accessibility_TNBC_tumours_freq <- promoter_accessibility_TNBC_tumours_freq[promoter_accessibility_TNBC_tumours_freq$sum >= 4,]

length(intersect(TNBC_TRs$TR, rownames(promoter_accessibility_TNBC_tumours_freq)))

tcga_accessible <- intersect(TNBC_TRs$TR, rownames(promoter_accessibility_TNBC_tumours_freq))

## (2) Accessible promoter across TCGA tumours, PDXs and CellLines
promoter_accessibility_sub <- promoter_accessibility[promoter_accessibility$geneID %in% tcga_accessible,]

## PDX
promoter_accessibility_PDX <- promoter_accessibility_sub[promoter_accessibility_sub$Sample_Group == "PDX",]
promoter_accessibility_PDX_freq <- table(promoter_accessibility_PDX$Sample, promoter_accessibility_PDX$geneID)
promoter_accessibility_PDX_freq <- t(promoter_accessibility_PDX_freq)
promoter_accessibility_PDX_freq <- ifelse(promoter_accessibility_PDX_freq >= 1, 1, 0)
promoter_accessibility_PDX_freq <- as.data.frame(promoter_accessibility_PDX_freq)
promoter_accessibility_PDX_freq$sum <- rowSums(promoter_accessibility_PDX_freq)
promoter_accessibility_PDX_freq <- promoter_accessibility_PDX_freq[promoter_accessibility_PDX_freq$sum >= 38/2,]

pdx_accessible <- intersect(TNBC_TRs$TR, rownames(promoter_accessibility_PDX_freq))

## CellLines
promoter_accessibility_CellLines <- promoter_accessibility_sub[promoter_accessibility_sub$Sample_Group == "CellLines",]
promoter_accessibility_CellLines_freq <- table(promoter_accessibility_CellLines$Sample, promoter_accessibility_CellLines$geneID)
promoter_accessibility_CellLines_freq <- t(promoter_accessibility_CellLines_freq)
promoter_accessibility_CellLines_freq <- ifelse(promoter_accessibility_CellLines_freq >= 1, 1, 0)
promoter_accessibility_CellLines_freq <- as.data.frame(promoter_accessibility_CellLines_freq)
promoter_accessibility_CellLines_freq$sum <- rowSums(promoter_accessibility_CellLines_freq)
promoter_accessibility_CellLines_freq <- promoter_accessibility_CellLines_freq[promoter_accessibility_CellLines_freq$sum >= 26/2,]

CellLines_accessible <- intersect(TNBC_TRs$TR, rownames(promoter_accessibility_CellLines_freq))

## Accessible across all TCGA, PDX and CellLines
intersect_1 <- intersect(tcga_accessible, pdx_accessible)
intersect_2 <- intersect(intersect_1, CellLines_accessible)

length(intersect_2)

# (2) Remove TRs with low mean activity scores (NES < 0) across all patient tumours, PDXs and CellLines
RNA_NES <- read.table("/data/Figure2/Combined_RNA_NES_Heatmap.tsv", header=TRUE, sep='\t')

RNA_NES <- RNA_NES[RNA_NES$geneID %in% intersect_2,]

### UpSet plot
TR_list <- list("TCGA 1-3" = RNA_NES[RNA_NES$Sample_Group == "TCGA" & RNA_NES$NES_Category == "1-3", "geneID"],
    "TCGA >3" = RNA_NES[RNA_NES$Sample_Group == "TCGA" & RNA_NES$NES_Category == ">3", "geneID"],
    "TCGA 0-1" = RNA_NES[RNA_NES$Sample_Group == "TCGA" & RNA_NES$NES_Category == "0-1", "geneID"],
    "TCGA <0" = RNA_NES[RNA_NES$Sample_Group == "TCGA" & RNA_NES$NES_Category == "<0", "geneID"],
    "PDX 1-3" = RNA_NES[RNA_NES$Sample_Group == "PDX" & RNA_NES$NES_Category == "1-3", "geneID"],
    "PDX >3" = RNA_NES[RNA_NES$Sample_Group == "PDX" & RNA_NES$NES_Category == ">3", "geneID"],
    "PDX 0-1" = RNA_NES[RNA_NES$Sample_Group == "PDX" & RNA_NES$NES_Category == "0-1", "geneID"],
    "PDX <0" = RNA_NES[RNA_NES$Sample_Group == "PDX" & RNA_NES$NES_Category == "<0", "geneID"],
    "CellLines 1-3" = RNA_NES[RNA_NES$Sample_Group == "CellLines" & RNA_NES$NES_Category == "1-3", "geneID"],
    "CellLines >3" = RNA_NES[RNA_NES$Sample_Group == "CellLines" & RNA_NES$NES_Category == ">3", "geneID"],
    "CellLines 0-1" = RNA_NES[RNA_NES$Sample_Group == "CellLines" & RNA_NES$NES_Category == "0-1", "geneID"],
    "CellLines <0" = RNA_NES[RNA_NES$Sample_Group == "CellLines" & RNA_NES$NES_Category == "<0", "geneID"])

TR_list_list_mat <- ComplexHeatmap::list_to_matrix(TR_list)

# Plot the UpSet plot
TR_upset <- ComplexUpset::upset(
  data      = as.data.frame(TR_list_list_mat),
  intersect = colnames(TR_list_list_mat),
  name   = "NES category of TRs with promoter accessibility\nin at least 50% of samples",
  themes = theme(text = element_text(size = 10, colour = "black")),
  base_annotations = list(
    "Number of TRs" = ComplexUpset::intersection_size() +
      labs(y = "Number of TRs")
  ),
  set_sizes = FALSE,
  queries = list(
    ComplexUpset::upset_query(
      intersect = c("TCGA <0", "PDX <0", "CellLines <0"),
      color = "#BF0413",
      fill  = "#BF0413",
      only_components = c(
        "intersections_matrix",  # dots
        "Number of TRs"          # your bar annotation name
      )
    )
  )
)

cairo_pdf(paste0("/results/visuals/Figure2/SupplementaryFigure4B_NES_Average_Across_CellLines_PDX_TCGA.pdf"), width = 15)
print(TR_upset)
dev.off()

RNA_NES <- RNA_NES[RNA_NES$NES_Category != "<0",]

RNA_NES_freq <- table(RNA_NES$Sample_Group, RNA_NES$geneID)
RNA_NES_freq <- t(RNA_NES_freq)
RNA_NES_freq <- as.matrix(RNA_NES_freq)
RNA_NES_freq <- as.data.frame.matrix(RNA_NES_freq)
RNA_NES_freq <- RNA_NES_freq %>% select(-c("METABRIC"))
RNA_NES_freq$sum <- rowSums(RNA_NES_freq)
RNA_NES_freq <- RNA_NES_freq[RNA_NES_freq$sum >= 1,]

High_RNA_NES <- intersect(rownames(RNA_NES_freq), intersect_2)
length(unique(High_RNA_NES))

#write.table(as.data.frame(High_RNA_NES), "code/Brain_Dump/High_RNA_NES.tsv", col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)

sessionInfo()
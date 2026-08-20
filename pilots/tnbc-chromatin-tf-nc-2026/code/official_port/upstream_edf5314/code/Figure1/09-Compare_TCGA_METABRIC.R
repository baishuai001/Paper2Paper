# Compare the TRs identified as TNBC-specific or Non-TNBC specific
# Author: Shalini Bahl
# Description:
#   This script compares transcriptional regulators (TRs) identified as TNBC- or
#   Non-TNBC-specific in TCGA and METABRIC cohorts. It generates Venn diagrams
#   to visualize overlaps and extracts shared and unique TRs for each group.
#
# Useage: Rscript code/Figure1/08-Prepare_Annotation_ForPlotting.R
#
# Outputs:
#   - Venn diagrams (PDF) for TNBC and Non-TNBC TRs
#   - Tables listing shared and cohort-specific TRs (.tsv)

library(VennDiagram)
library(eulerr)

dir.create("/results/visuals/Figure1/", showWarnings=FALSE, recursive=TRUE)

# Read in the NES-Limma results from TCGA
TCGA_NES_Limma <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')
table(TCGA_NES_Limma$category)

# Read in the NES-Limma results from METABRIC
METABRIC_NES_Limma <- read.table("/data/Figure1/Tables_ForPlotting/METABRIC_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')
table(METABRIC_NES_Limma$category)

## Check TNBC
TCGA_TNBC <- TCGA_NES_Limma[TCGA_NES_Limma$category == "TNBC",]
METABRIC_TNBC <- METABRIC_NES_Limma[METABRIC_NES_Limma$category == "TNBC",]

TNBC_list <- list(TCGA = TCGA_TNBC$tf,
    METABRIC = METABRIC_TNBC$tf)

fit <- euler(TNBC_list)

venn_diagram_TNBC <- plot(
  fit,
  fills = list(fill = c("#435773", "#435773"), alpha = 0.7),
  labels = list(font = 2), # bold labels
  quantities = TRUE)

cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2B_VennDiagram_TNBC_TCGA_METABRIC.pdf", width = 3.5, height = 3)
print(venn_diagram_TNBC)
dev.off()

# Extract combinations and save to a table
TNBC_TCGA_METABRIC_comparision <- data.frame(matrix(nrow=0, ncol=2))
colnames(TNBC_TCGA_METABRIC_comparision) <- c("Group", "TR")

# Unique to TCGA
TNBC_only_TCGA <- setdiff(TNBC_list$TCGA, TNBC_list$METABRIC)
TNBC_TCGA_METABRIC_comparision <- rbind(TNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "TNBC_TCGA_Only", TR = TNBC_only_TCGA))

# Unique to METABRIC
TNBC_only_METABRIC <- setdiff(TNBC_list$METABRIC, TNBC_list$TCGA)
TNBC_TCGA_METABRIC_comparision <- rbind(TNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "TNBC_METABRIC_Only", TR = TNBC_only_METABRIC))

# Shared between both
TNBC_shared <- intersect(TNBC_list$TCGA, TNBC_list$METABRIC)
TNBC_TCGA_METABRIC_comparision <- rbind(TNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "TNBC_METABRIC_TCGA_Shared", TR = TNBC_shared))

#write.table(TNBC_TCGA_METABRIC_comparision,
#  "data/Fig1_TRNetwork/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv",
#  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

## Check Non-TNBC
TCGA_NonTNBC <- TCGA_NES_Limma[TCGA_NES_Limma$category == "Non-TNBC",]
METABRIC_NonTNBC <- METABRIC_NES_Limma[METABRIC_NES_Limma$category == "Non-TNBC",]

NonTNBC_list <- list(TCGA = TCGA_NonTNBC$tf,
    METABRIC = METABRIC_NonTNBC$tf)

fit <- euler(NonTNBC_list)

venn_diagram_NonTNBC <- plot(
  fit,
  fills = list(fill = c("#9bc1e5", "#9bc1e5"), alpha = 0.7),
  labels = list(font = 2), # bold labels
  quantities = TRUE)

cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2B_VennDiagram_NonTNBC_TCGA_METABRIC.pdf", width = 3.5, height = 3)
print(venn_diagram_NonTNBC)
dev.off()

# Extract combinations and save to a table
NonTNBC_TCGA_METABRIC_comparision <- data.frame(matrix(nrow=0, ncol=2))
colnames(NonTNBC_TCGA_METABRIC_comparision) <- c("Group", "TR")

# Unique to TCGA
NonTNBC_only_TCGA <- setdiff(NonTNBC_list$TCGA, NonTNBC_list$METABRIC)
NonTNBC_TCGA_METABRIC_comparision <- rbind(NonTNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "NonTNBC_TCGA_Only", TR = NonTNBC_only_TCGA))

# Unique to METABRIC
NonTNBC_only_METABRIC <- setdiff(NonTNBC_list$METABRIC, NonTNBC_list$TCGA)
NonTNBC_TCGA_METABRIC_comparision <- rbind(NonTNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "NonTNBC_METABRIC_Only", TR = NonTNBC_only_METABRIC))

# Shared between both
NonTNBC_shared <- intersect(NonTNBC_list$TCGA, NonTNBC_list$METABRIC)
NonTNBC_TCGA_METABRIC_comparision <- rbind(NonTNBC_TCGA_METABRIC_comparision, 
  data.frame(Group = "NonTNBC_METABRIC_TCGA_Shared", TR = NonTNBC_shared))

#write.table(NonTNBC_TCGA_METABRIC_comparision,
#  "data/Fig1_TRNetwork/Tables_ForPlotting/NonTNBC_TCGA_METABRIC_comparision.tsv",
#  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

sessionInfo()
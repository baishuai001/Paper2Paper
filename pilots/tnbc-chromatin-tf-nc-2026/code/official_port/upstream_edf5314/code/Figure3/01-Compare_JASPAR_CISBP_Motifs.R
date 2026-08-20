# Compare JASPAR and CISBP motif enrichments
library(VennDiagram)
library(eulerr)

# Read in the HC-TR
HC_TR <- read.table("code/Brain_Dump/High_RNA_NES.tsv", header=FALSE, sep='\t')
head(HC_TR)

# Compare CISBP and JASPAR results
CISBP <- read.table("data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/CISBP_Combined_MotifEnrichment_TNBC.tsv",
    header=TRUE, sep='\t')
CISBP$Motif.Name <- toupper(CISBP$Motif.Name)
CISBP <- CISBP[CISBP$Motif.Name %in% HC_TR$V1,]
length(unique(CISBP$Motif.Name))

JASPAR <- read.table("data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
    header=TRUE, sep='\t')
JASPAR$Motif.Name <- toupper(JASPAR$Motif.Name)
JASPAR <- JASPAR[JASPAR$Motif.Name %in% HC_TR$V1,]
length(unique(JASPAR$Motif.Name))

# Generate a venn diagram comparing JASPAR called and CISBP called
JASPAR_CISBP_TNBC_Motif <- list(JASPAR_HC_TR_TNBC = unique(JASPAR$Motif.Name),
    CISBP_HC_TR_TNBC = unique(CISBP$Motif.Name))

fit <- euler(JASPAR_CISBP_TNBC_Motif)

venn_diagram_JASPAR_CISBP_TNBC_Motif <- plot(
  fit,
  fills = list(fill = c("#ADD7D6", "#DA779A"), alpha = 0.7),
  labels = list(font = 2), # bold labels
  quantities = TRUE)

cairo_pdf("visuals/Figure2_TR_Chromatin/JASPAR_vs_CISBP_HC-TRs_VennDiagram.pdf",
    width = 6, height = 3)
print(venn_diagram_JASPAR_CISBP_TNBC_Motif)
dev.off()

setdiff(JASPAR_CISBP_TNBC_Motif$JASPAR_HC_TR_TNBC, JASPAR_CISBP_TNBC_Motif$CISBP_HC_TR_TNBC)
setdiff(JASPAR_CISBP_TNBC_Motif$CISBP_HC_TR_TNBC, JASPAR_CISBP_TNBC_Motif$JASPAR_HC_TR_TNBC)

# Save the JASPAR vs. CISBP comparision
jaspar_only <- setdiff(unique(JASPAR$Motif.Name), unique(CISBP$Motif.Name))
cisbp_only <- setdiff(unique(CISBP$Motif.Name), unique(JASPAR$Motif.Name))
both <- intersect(unique(JASPAR$Motif.Name), unique(CISBP$Motif.Name))

jaspar_only_df <- data.frame(TR = jaspar_only, Database = "JASPAR_only")
cisbp_only_df <- data.frame(TR = cisbp_only, Database = "CISBP_only")
both_df <- data.frame(TR = both, Database = "Both")

# Combine into long format
comparison_long <- rbind(jaspar_only_df, cisbp_only_df, both_df)

# Save long format version
write.csv(comparison_long, "data/Figure2_TR_Chromatin/Motif/JASPAR_CISBP_comparison_long.csv", row.names = FALSE)

sessionInfo()
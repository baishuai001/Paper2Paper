# Generate input files for Figure 1 and associated Supplementary Figures
# Author: Shalini Bahl
# Description:
#   This script generates input files and visualizations for Figure 1 and related
#   supplementary figures. It includes scatter plots and heatmaps of transcriptional
#   regulator (TR) activity across TCGA, PDX, and cell line cohorts, using VIPER
#   scores and molecular annotations.
#
# Outputs:
#   - Figure 1B: NES vs. log2FC scatter plot (PDF)
#   - Figure 1C–E: TR activity heatmaps for TCGA, PDX, and cell lines (PDF)
#   - Intermediate annotation and activity matrices (.tsv)

library(viper)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(ggrepel)
library(circlize)
library(dendextend)

source("/code/Static_Scripts/plotting_aesthetics.R")

dir.create("/results/visuals/Figure1", showWarnings=FALSE, recursive = TRUE)
dir.create("/results/Figure1/Basal_FisherTests", showWarnings=FALSE, recursive=TRUE)

## Figure 1B - VIPER TCGA TNBC vs. Non-TNBC
TCGA_NES <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')

# Add in text labels
TCGA_NES$label <- ifelse(TCGA_NES$category %in% c("TNBC", "Non-TNBC"), TCGA_NES$tf, NA)

# Plot the logFC results against the viper based NES for the TR regulators from PANGO-GO
TCGA_limma_logFC_viper_NES_plt <- TCGA_NES %>% ggplot(aes(x = TCGA_NES, y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    geom_point(aes(color = category), size = 0.75) +
    geom_text_repel(aes(label = label), size = 3, 
        min.segment.length = 0.1, max.overlaps = 40) +
    scale_color_manual(values = viper_cutoff_colours$VIPER_TR_Category) +
    geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
    stat_cor(method = "pearson", digits = 2, label.x = 2.5, label.y = -3.75) +
    theme_classic(base_size = 12) +
    labs(x = "NES", y = "log2FC Gene expression") +
    theme(legend.position = "none")

# Save the figure in PDF format
cairo_pdf("/results/visuals/Figure1/Figure1B_log2FC_Viper_TNBC.pdf", width = 3.5, height = 4)
print(TCGA_limma_logFC_viper_NES_plt)
dev.off()

## Figure 1C - Heatmap TCGA TR activity
TCGA_TF_activities <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_Activity_TNBC150_NonTNBC155.tsv",
    header=TRUE, sep='\t', check.names = FALSE)

TCGA_TF_activities <- t(scale(t(TCGA_TF_activities)))

TCGA_column_annotation <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_Activity_TNBC150_NonTNBC155_SampleAnnotation.tsv",
    header=TRUE, sep='\t')
rownames(TCGA_column_annotation) <- TCGA_column_annotation$Sample

TCGA_column_annotation$Cancer_Type <- ifelse(grepl("-06", TCGA_column_annotation$Sample) == TRUE,
    "Metastatic", "Primary")
TCGA_column_annotation$Cancer_Type <- factor(TCGA_column_annotation$Cancer_Type,
    levels = c("Primary", "Metastatic"))

TCGA_column_annotation <- TCGA_column_annotation[,-1]
TCGA_column_annotation <- TCGA_column_annotation[colnames(TCGA_TF_activities),]

TCGA_column_annotation$ER[is.na(TCGA_column_annotation$ER)] <- "NA"
TCGA_column_annotation$ER <- factor(TCGA_column_annotation$ER, 
    levels = c("Positive", "Negative", "NA"))

TCGA_column_annotation$PR[is.na(TCGA_column_annotation$PR)] <- "NA"
TCGA_column_annotation$PR <- factor(TCGA_column_annotation$PR, 
    levels = c("Positive", "Negative", "NA"))

TCGA_column_annotation$HER2[is.na(TCGA_column_annotation$HER2)] <- "NA"
TCGA_column_annotation$HER2 <- factor(TCGA_column_annotation$HER2, 
    levels = c("Positive", "Negative", "NA"))

TCGA_column_annotation$Age[is.na(TCGA_column_annotation$Age)] <- "NA"
TCGA_column_annotation$Age <- factor(TCGA_column_annotation$Age,
    levels = c("<=40", "(40-70)", ">=70", "NA"))

TCGA_column_annotation$VitalSt[is.na(TCGA_column_annotation$VitalSt)] <- "NA"
TCGA_column_annotation$VitalSt <- factor(TCGA_column_annotation$VitalSt,
    levels = c("Dead", "Alive", "NA"))

TCGA_column_annotation$IntClust <- factor(TCGA_column_annotation$IntClust,
    levels = paste0("iC", c(1:10)))

TCGA_row_annotation <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_Activity_TNBC150_NonTNBC155_TRAnnotation.tsv",
    header=TRUE, sep='\t')
TCGA_row_annotation$VIPER_TR_Category <- factor(TCGA_row_annotation$VIPER_TR_Category, levels = c("TNBC", "Non-TNBC"))

# Construct annotations
TCGA_row_anno <- HeatmapAnnotation(df = TCGA_row_annotation,
    col = viper_cutoff_colours,
    which = "row")

# Change column order
TCGA_column_annotation_main <- TCGA_column_annotation[,c("Cancer_Type","ER","PR","HER2","Age","VitalSt","PAM50","SCMOD2","IntClust")]
TCGA_col_anno <- HeatmapAnnotation(df = TCGA_column_annotation_main,
    col = clinical_annot_colours,
    annotation_name_side = "left",
    which = "column")

TCGA_TNBC_NonTNBC_ActivityHeatmap <- Heatmap(as.matrix(TCGA_TF_activities),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    left_annotation = TCGA_row_anno,
    top_annotation = TCGA_col_anno,
    heatmap_legend_param = list(title = "TR activity"),
    show_row_names = FALSE,
#    column_km=2,
    column_title_gp = gpar(fontsize = 10), # optional: smaller column title
    column_names_gp = gpar(fontsize = 8)   # optional: smaller column names
)

# Reduce the height of the column annotation
TCGA_col_anno@anno_size <- unit(1, "mm")  # set to a smaller value as needed

cairo_pdf("/results/visuals/Figure1/Figure1C_TCGA_TRActivity_TNBC150_NonTNBC155_test.pdf", width = 7.5, height = 6.5)
draw(TCGA_TNBC_NonTNBC_ActivityHeatmap)
dev.off()

TCGA_htmap <- draw(TCGA_TNBC_NonTNBC_ActivityHeatmap)

## NonTNBC_specific
dend_NonTNBC_Specific <- column_dend(TCGA_htmap)[[1]]
dend_df_NonTNBC <- as.data.frame(labels(dend_NonTNBC_Specific))
dend_df_NonTNBC$Category <- "NonTNBC_Specific_TR_Active_Cluster"
colnames(dend_df_NonTNBC) <- c("Sample", "TR_Cluster_Activity")

## TNBC_specific
dend_TNBC_Specific <- column_dend(TCGA_htmap)[[2]]
dend_df_TNBC <- as.data.frame(labels(dend_TNBC_Specific))
dend_df_TNBC$Category <- "TNBC_Specific_TR_Active_Cluster"
colnames(dend_df_TNBC) <- c("Sample", "TR_Cluster_Activity")

dend_df <- rbind(dend_df_NonTNBC, dend_df_TNBC)

#write.table(dend_df, paste0("/data/Figure1/Basal_FisherTests/TCGA_150TNBC_154NonTNBC_Dendrogram_Clusters.tsv"),
#    col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

## Supplementary Figure 1 - TCGA TR activity correlation plots
TCGA_column_annotation$HistologicalType <- factor(TCGA_column_annotation$HistologicalType)

TCGA_TF_activities_sub_pearson_corr <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_Activity_TNBC150_NonTNBC155_PearsonCor.tsv",
    header=TRUE, sep='\t')

TCGA_column_annotation_all <- TCGA_column_annotation[,c("Cancer_Type","HistologicalType", "ER","PR","HER2","Age","VitalSt","PAM50","SCMOD2","IntClust")]
TCGA_col_anno <- HeatmapAnnotation(df = TCGA_column_annotation_all,
    col = clinical_annot_colours,
    annotation_name_side = "left",
    which = "column")

col_fun = colorRamp2(c(0.5, 0.75, 1), hcl_palette = "Purple-Green", reverse = TRUE)

TCGA_TF_activities_sub_pearson_corr_heatmap <- Heatmap(as.matrix(TCGA_TF_activities_sub_pearson_corr),
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    top_annotation = TCGA_col_anno,
    heatmap_legend_param = list(title = "Pearson correlation"),
    show_row_names = FALSE)

cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2C_TCGA_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf", width = 7.5, height = 6)
draw(TCGA_TF_activities_sub_pearson_corr_heatmap)
dev.off()

## Figure 1D - Heatmap PDX TR activity scores
PDX_TF_activities <- read.table("/data/Figure1/PDX/PDX73_TR_activities_viper_usingTCGAregulons.tsv",
    header=TRUE, sep='\t', check.names = FALSE)

TRs_TCGA_Category <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')
TNBC <- TRs_TCGA_Category[TRs_TCGA_Category$category == "TNBC",]
NonTNBC <- TRs_TCGA_Category[TRs_TCGA_Category$category == "Non-TNBC",]

PDX_TF_activities <- PDX_TF_activities[c(TNBC$tf, NonTNBC$tf),]

PDX_TF_activities <-  t(scale(t(PDX_TF_activities)))

PDX_column_annotation <- read.table("/data/Reference_Data/PDX/PDX_molecular_subtype.tsv",
    header=TRUE, sep='\t')
rownames(PDX_column_annotation) <- PDX_column_annotation$Sample
PDX_column_annotation <- PDX_column_annotation[colnames(PDX_TF_activities),]

PDX_column_annotation <- PDX_column_annotation[,c("IntClust", "SCMOD2", "PAM50")]

PDX_col_anno <- HeatmapAnnotation(df = PDX_column_annotation,
    col = clinical_annot_colours,
    annotation_name_side = "left",
    which = "column")

PDX_TNBC_NonTNBC_ActivityHeatmap <- Heatmap(as.matrix(PDX_TF_activities),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    left_annotation = TCGA_row_anno,
    top_annotation = PDX_col_anno,
#    column_km=2,
    heatmap_legend_param = list(title = "TR activity"),
    show_row_names = FALSE,
    column_title_gp = gpar(fontsize = 10), # optional: smaller column title
    column_names_gp = gpar(fontsize = 8)   # optional: smaller column names
    )

cairo_pdf("/results/visuals/Figure1/Figure1D_PDX_TRActivity_TNBC150_NonTNBC154.pdf", width = 7, height = 6.5)
print(PDX_TNBC_NonTNBC_ActivityHeatmap)
dev.off()

PDX_htmap <- draw(PDX_TNBC_NonTNBC_ActivityHeatmap)

## NonTNBC_specific
dend_NonTNBC_Specific <- column_dend(PDX_htmap)[[1]]
dend_df_NonTNBC <- as.data.frame(labels(dend_NonTNBC_Specific))
dend_df_NonTNBC$Category <- "NonTNBC_Specific_TR_Active_Cluster"
colnames(dend_df_NonTNBC) <- c("Sample", "TR_Cluster_Activity")

## TNBC_specific
dend_TNBC_Specific <- column_dend(PDX_htmap)[[2]]
dend_df_TNBC <- as.data.frame(labels(dend_TNBC_Specific))
dend_df_TNBC$Category <- "TNBC_Specific_TR_Active_Cluster"
colnames(dend_df_TNBC) <- c("Sample", "TR_Cluster_Activity")

dend_df <- rbind(dend_df_NonTNBC, dend_df_TNBC)

write.table(dend_df, paste0("/data/Figure1/Basal_FisherTests/PDX_150TNBC_154NonTNBC_Dendrogram_Clusters.tsv"),
    col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

## Supplementary Figure 1B - PDX pearson corelation
head(PDX_column_annotation)

PDX_TF_activities_sub_pearson_corr <- read.table("/data/Figure1/Tables_ForPlotting/PDX_TR_Activity_TNBC150_NonTNBC154_PearsonCor.tsv",
    header=TRUE, sep='\t')

col_anno <- HeatmapAnnotation(df = PDX_column_annotation,
    col = clinical_annot_colours,
    annotation_name_side = "left",
    which = "column")

col_fun = colorRamp2(c(0.5, 0.75, 1), hcl_palette = "Purple-Green", reverse = TRUE)

PDX_TF_activities_sub_pearson_corr_heatmap <- Heatmap(as.matrix(PDX_TF_activities_sub_pearson_corr),
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    top_annotation = col_anno,
    heatmap_legend_param = list(title = "Pearson correlation"),
    show_row_names = FALSE)

cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2D_PDX_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf", width = 7.5, height = 6)
draw(PDX_TF_activities_sub_pearson_corr_heatmap)
dev.off()

## Figure 1E - Heatmap cell-lines TR activity scores
cell_line_TF_activities <- read.table("/data/Figure1/CellLines/CellLines82_TR_activities_viper_usingTCGAregulons.tsv",
    header=TRUE, sep='\t', check.names = FALSE)

TRs_TCGA_Category <- read.table("/data/Figure1/Tables_ForPlotting/TCGA_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')
TNBC <- TRs_TCGA_Category[TRs_TCGA_Category$category == "TNBC",]
NonTNBC <- TRs_TCGA_Category[TRs_TCGA_Category$category == "Non-TNBC",]

cell_line_TF_activities <- cell_line_TF_activities[c(TNBC$tf, NonTNBC$tf),]
cell_line_TF_activities <- t(scale(t(cell_line_TF_activities)))

### Read in annotations for column
cell_line_molecular_subtyping <- read.table("/data/Reference_Data/CellLines_Marcotte_et_al_2016/CellLines_molecular_subtype.tsv",
    header=TRUE, sep='\t')
rownames(cell_line_molecular_subtyping) <- cell_line_molecular_subtyping$Sample

## Read in the Lehmann annotations
Lehmann_annotations <- read.table("/data/Reference_Data/CellLines_Marcotte_et_al_2016/Lehmann_Annotation.tsv",
    header=TRUE, sep='\t')

cell_line_molecular_subtyping <- merge(cell_line_molecular_subtyping, Lehmann_annotations, by = "Sample")
colnames(cell_line_molecular_subtyping)[grep("Annotation", colnames(cell_line_molecular_subtyping))] <- "Lehmann"

cell_line_molecular_subtyping <- cell_line_molecular_subtyping[,-1]

cell_line_molecular_subtyping <- cell_line_molecular_subtyping[,c("Lehmann", "IntClust", "SCMOD2", "PAM50")] # Change the order of annotation bars
cell_line_molecular_subtyping$IntClust <- factor(cell_line_molecular_subtyping$IntClust,
    levels = paste0("iC", c(1:10)))

### Generate the column annotation bars
col_anno <- HeatmapAnnotation(df = cell_line_molecular_subtyping,
    col = clinical_annot_colours,
    annotation_name_side = "left",
    which = "column")

Cell_line_TNBC_NonTNBC_ActivityHeatmap <- Heatmap(as.matrix(cell_line_TF_activities),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    left_annotation = TCGA_row_anno,
    top_annotation = col_anno,
    heatmap_legend_param = list(title = "TR activity"),
    show_row_names = FALSE,
#    column_km = 2,
    column_title_gp = gpar(fontsize = 10), # optional: smaller column title
    column_names_gp = gpar(fontsize = 8)   # optional: smaller column names
    )

cairo_pdf("/results/visuals/Figure1/Figure1E_CellLine_TRActivity_TNBC150_NonTNBC155.pdf", width = 7, height = 6.5)
print(Cell_line_TNBC_NonTNBC_ActivityHeatmap)
dev.off()

CellLines_htmap <- draw(Cell_line_TNBC_NonTNBC_ActivityHeatmap)

## NonTNBC_specific
dend_NonTNBC_Specific <- column_dend(CellLines_htmap)[[1]]
dend_df_NonTNBC <- as.data.frame(labels(dend_NonTNBC_Specific))
dend_df_NonTNBC$Category <- "NonTNBC_Specific_TR_Active_Cluster"
colnames(dend_df_NonTNBC) <- c("Sample", "TR_Cluster_Activity")

## TNBC_specific
dend_TNBC_Specific <- column_dend(CellLines_htmap)[[2]]
dend_df_TNBC <- as.data.frame(labels(dend_TNBC_Specific))
dend_df_TNBC$Category <- "TNBC_Specific_TR_Active_Cluster"
colnames(dend_df_TNBC) <- c("Sample", "TR_Cluster_Activity")

dend_df <- rbind(dend_df_NonTNBC, dend_df_TNBC)

write.table(dend_df, paste0("/data/Figure1/Basal_FisherTests/CellLines_150TNBC_154NonTNBC_Dendrogram_Clusters.tsv"),
    col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

## Supplementary Figure 1C - Cell-line correlation plots
cell_lines_TF_activities_sub_pearson_corr <- read.table("/data/Figure1/Tables_ForPlotting/CellLines_TR_Activity_TNBC150_NonTNBC155_PearsonCor.tsv",
  header=TRUE, sep='\t')

col_fun = colorRamp2(c(0.5, 0.75, 1), hcl_palette = "Purple-Green", reverse = TRUE)

cell_lines_TF_activities_sub_pearson_corr_heatmap <- Heatmap(as.matrix(cell_lines_TF_activities_sub_pearson_corr),
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    top_annotation = col_anno,
    heatmap_legend_param = list(title = "Pearson correlation"),
    show_row_names = FALSE)

cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2E_Cell_Lines_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf", width = 7.5, height = 6)
draw(cell_lines_TF_activities_sub_pearson_corr_heatmap)
dev.off()

## Supplementary Figure 1D - METABRIC logFC and VIPER results
METABRIC_NES <- read.table("/data/Figure1/Tables_ForPlotting/METABRIC_TR_activities_NES_viper.tsv",
    header=TRUE, sep='\t')

# Add in text labels
METABRIC_NES$label <- ifelse(METABRIC_NES$category %in% c("TNBC", "Non-TNBC"), METABRIC_NES$tf, NA)

# Plot the logFC results against the viper based NES for the TR regulators from PANGO-GO
METABRIC_limma_logFC_viper_NES_plt <- METABRIC_NES %>% ggplot(aes(x = METABRIC_NES, y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    geom_point(aes(color = category), size = 0.75) +
    geom_text_repel(aes(label = label), size = 3, 
        min.segment.length = 0.1, max.overlaps = 40) +
    scale_color_manual(values = viper_cutoff_colours$VIPER_TR_Category) +
    geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
    stat_cor(method = "pearson", digits = 2, label.x = 2.5, label.y = -3.75) +
    theme_classic(base_size = 12) +
    labs(x = "NES", y = "log2FC Gene expression") +
    theme(legend.position = "none")

# Annotate the number of TRs in each category directly on the plot
# Find category centroids
centroids <- METABRIC_NES %>%
  group_by(category) %>%
  summarize(
    METABRIC_NES = min(METABRIC_NES), 
    logFC = max(logFC),
    n_TFs = n(), 
    .groups = 'drop'
  )
centroids <- centroids[centroids$category %in% c("TNBC", "Non-TNBC"),]

METABRIC_limma_logFC_viper_NES_plt <- METABRIC_limma_logFC_viper_NES_plt +
  geom_text(
    data = centroids,
    aes(x = METABRIC_NES, y = logFC, label = n_TFs, color = category)#,
    #size = 4, vjust = -0.5 # or adjust hjust for horizontal nudging
  )

# Save the figure in PDF format
cairo_pdf("/results/visuals/Figure1/SupplementaryFigure2A_METABRIC_log2FC_Viper.pdf", width = 3.5, height = 4)
print(METABRIC_limma_logFC_viper_NES_plt)
dev.off()

## Supplementary Figure 1F - TNBC-Basal TCGA tumors NES distribution across TNBC TRs
TCGA_Basal <- TCGA_column_annotation[TCGA_column_annotation$ER == "Negative" & TCGA_column_annotation$PR == "Negative" &
    TCGA_column_annotation$HER2 == "Negative",]
TCGA_Basal <- TCGA_Basal[TCGA_Basal$PAM50 == "Basal" | TCGA_Basal$SCMOD2 == "Basal",]

#TNBC_TCGA_METABRIC_TRs <- c(TNBC$tf, METABRIC_NES[METABRIC_NES$category == "TNBC", "tf"])
TNBC_TCGA_METABRIC_TRs <- TNBC$tf
TNBC_TCGA_METABRIC_TRs <- unique(TNBC_TCGA_METABRIC_TRs)

TCGA_TF_activities <- read.table("/data/Figure1/TCGA/TCGA_TR_activities_viper.tsv",
    header=TRUE, sep='\t')
colnames(TCGA_TF_activities) <- gsub("\\.", "-", colnames(TCGA_TF_activities))

TCGA_TF_activities_TNBC_Specific <- TCGA_TF_activities[rownames(TCGA_TF_activities) %in% TNBC_TCGA_METABRIC_TRs, rownames(TCGA_Basal)]
dim(TCGA_TF_activities_TNBC_Specific)

TCGA_TF_activities_TNBC_Specific$Median_NES <- apply(TCGA_TF_activities_TNBC_Specific, 1, median)
TCGA_TF_activities_TNBC_Specific$tf <- rownames(TCGA_TF_activities_TNBC_Specific)

library(tidyr)
TCGA_TF_activities_TNBC_Specific_long <- TCGA_TF_activities_TNBC_Specific %>% 
    select(-c("Median_NES")) %>%
    gather(key = "Patient", value = "NES", -"tf")

library(ggplot2)
library(dplyr)

# Sort TFs by median NES
sorted_tfs <- TCGA_TF_activities_TNBC_Specific_long %>%
  group_by(tf) %>%
  summarize(median_NES = median(NES, na.rm = TRUE)) %>%
  arrange(median_NES) %>%
  pull(tf)

# Add factor levels for sorted order
TCGA_TF_activities_TNBC_Specific_long <- TCGA_TF_activities_TNBC_Specific_long %>%
  mutate(tf = factor(tf, levels = sorted_tfs),
         color_group = ifelse(NES >= 0, "Positive", "Negative"))

# Plot
TCGA_TNBC_Basal <- ggplot(TCGA_TF_activities_TNBC_Specific_long, aes(x = NES, y = tf, color = color_group)) +
  geom_point(size = 0.10) +
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "navy")) +
  theme_classic() +
  labs(
       x = "NES",
       y = NULL) +
  theme(axis.text.y = element_text(size = 4.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none") +
  geom_vline(xintercept = 0, linetype="dashed", color = "black")

cairo_pdf(paste0("/results/visuals/Figure1/SupplementaryFigure2F_TCGA_Basal", length(unique(TCGA_TF_activities_TNBC_Specific_long$tf)), "TRs.pdf"),
     width = 3.75, height = 8)
print(TCGA_TNBC_Basal)
dev.off()

## Supplementary Figure 1F - TNBC-Basal PDX tumors NES distribution across TNBC TRs
PDX_column_annotation <- read.table("/data/Reference_Data/PDX/PDX_molecular_subtype.tsv",
    header=TRUE, sep='\t')

PDX_Basal <- PDX_column_annotation[PDX_column_annotation$PAM50 == "Basal" | PDX_column_annotation$SCMOD2 == "Basal",]

PDX_TF_activities <- read.table("/data/Figure1/PDX/PDX73_TR_activities_viper_usingTCGAregulons.tsv",
    header=TRUE, sep='\t', check.names = FALSE)

sorted_tfs_PDX <- unique(intersect(sorted_tfs, rownames(PDX_TF_activities)))

PDX_TF_activities_TNBC_Specific <- PDX_TF_activities[sorted_tfs_PDX, PDX_Basal$Sample]
dim(PDX_TF_activities_TNBC_Specific)

PDX_TF_activities_TNBC_Specific$Median_NES <- apply(PDX_TF_activities_TNBC_Specific, 1, median)
PDX_TF_activities_TNBC_Specific$tf <- rownames(PDX_TF_activities_TNBC_Specific)

PDX_TF_activities_TNBC_Specific_long <- TCGA_TF_activities_TNBC_Specific %>% 
    select(-c("Median_NES")) %>%
    gather(key = "Patient", value = "NES", -"tf")

PDX_TF_activities_TNBC_Specific_long <- PDX_TF_activities_TNBC_Specific_long %>%
  mutate(tf = factor(tf, levels = sorted_tfs_PDX),
         color_group = ifelse(NES >= 0, "Positive", "Negative"))

## Plot
PDX_TNBC_Basal <- ggplot(PDX_TF_activities_TNBC_Specific_long, aes(x = NES, y = tf, color = color_group)) +
  geom_point(size = 0.10) +
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "navy")) +
  theme_classic() +
  labs(
       x = "NES",
       y = NULL) +
  theme(axis.text.y = element_text(size = 4.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none") +
  geom_vline(xintercept = 0, linetype="dashed", color = "black")

cairo_pdf(paste0("/results/visuals/Figure1/SupplementaryFigure2F_PDX_Basal", length(unique(PDX_TF_activities_TNBC_Specific_long$tf)), "TRs.pdf"), 
    width = 3.75, height = 8)
print(PDX_TNBC_Basal)
dev.off()

## Supplementary Figure 1F - TNBC-Basal Cell lines tumors NES distribution across TNBC TRs
cell_line_molecular_subtyping <- read.table("/data/Reference_Data/CellLines_Marcotte_et_al_2016/CellLines_molecular_subtype.tsv",
    header=TRUE, sep='\t')
rownames(cell_line_molecular_subtyping) <- cell_line_molecular_subtyping$Sample

cell_line_Basal <- cell_line_molecular_subtyping[cell_line_molecular_subtyping$PAM50 == "Basal" | cell_line_molecular_subtyping$SCMOD2 == "Basal",]

cell_line_TF_activities <- read.table("/data/Figure1/CellLines/CellLines82_TR_activities_viper_usingTCGAregulons.tsv",
    header=TRUE, sep='\t', check.names = FALSE)

sorted_tfs_cell_lines <- unique(intersect(sorted_tfs, rownames(cell_line_TF_activities)))

CellLine_TF_activities_TNBC_Specific <- cell_line_TF_activities[sorted_tfs_cell_lines, cell_line_Basal$Sample]
dim(CellLine_TF_activities_TNBC_Specific)

CellLine_TF_activities_TNBC_Specific$Median_NES <- apply(CellLine_TF_activities_TNBC_Specific, 1, median)
CellLine_TF_activities_TNBC_Specific$tf <- rownames(CellLine_TF_activities_TNBC_Specific)

CellLine_TF_activities_TNBC_Specific_long <- CellLine_TF_activities_TNBC_Specific %>% 
    select(-c("Median_NES")) %>%
    gather(key = "Patient", value = "NES", -"tf")

CellLine_TF_activities_TNBC_Specific_long <- CellLine_TF_activities_TNBC_Specific_long %>%
  mutate(tf = factor(tf, levels = sorted_tfs_cell_lines),
         color_group = ifelse(NES >= 0, "Positive", "Negative"))

## Plot
CellLines_TNBC_Basal <- ggplot(CellLine_TF_activities_TNBC_Specific_long, aes(x = NES, y = tf, color = color_group)) +
  geom_point(size = 0.10) +
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "navy")) +
  theme_classic() +
  labs(
       x = "NES",
       y = NULL) +
  theme(axis.text.y = element_text(size = 4.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none") +
  geom_vline(xintercept = 0, linetype="dashed", color = "black")

cairo_pdf(paste0("/results/visuals/Figure1/SupplementaryFigure2F_CellLines_Basal",
    length(unique(CellLine_TF_activities_TNBC_Specific_long$tf)), "TRs.pdf"), width = 3.75, height = 8)
print(CellLines_TNBC_Basal)
dev.off()

## Supplementary Figure 1F - TNBC-Basal METABRIC tumors NES distribution across TNBC TRs
METABRIC_clinical <- read.table("/data/Reference_Data/METABRIC/METABRIC_TNBC233_clinical_info.tsv",
    header=TRUE, sep='\t')

METABRIC_molecular_subtype <- read.table("/data/Reference_Data/METABRIC/METABRIC_molecular_subtype.tsv",
    header=TRUE, sep='\t')
METABRIC_molecular_subtype <- METABRIC_molecular_subtype[METABRIC_molecular_subtype$Sample %in% METABRIC_clinical$SAMPLE_ID,]

METABRIC_molecular_subtype <- METABRIC_molecular_subtype[METABRIC_molecular_subtype$PAM50 == "Basal" | METABRIC_molecular_subtype$SCMOD2 == "Basal",]

#TNBC_TCGA_METABRIC_TRs <- c(TNBC$tf, METABRIC_NES[METABRIC_NES$category == "TNBC", "tf"])
TNBC_TCGA_METABRIC_TRs <- c(TNBC$tf)
TNBC_TCGA_METABRIC_TRs <- unique(TNBC_TCGA_METABRIC_TRs)

METABRIC_TF_activities <- read.table("/data/Figure1/METABRIC/METABRIC_TR_activities_viper.tsv",
    header=TRUE, sep='\t')
colnames(METABRIC_TF_activities) <- gsub("\\.", "-", colnames(METABRIC_TF_activities))

METABRIC_TF_activities_TNBC_Specific <- METABRIC_TF_activities[rownames(METABRIC_TF_activities) %in% TNBC_TCGA_METABRIC_TRs, METABRIC_molecular_subtype$Sample]
dim(METABRIC_TF_activities_TNBC_Specific)

METABRIC_TF_activities_TNBC_Specific$Median_NES <- apply(METABRIC_TF_activities_TNBC_Specific, 1, median)
METABRIC_TF_activities_TNBC_Specific$tf <- rownames(METABRIC_TF_activities_TNBC_Specific)

library(tidyr)
METABRIC_TF_activities_TNBC_Specific_long <- METABRIC_TF_activities_TNBC_Specific %>% 
    select(-c("Median_NES")) %>%
    gather(key = "Patient", value = "NES", -"tf")

library(ggplot2)
library(dplyr)

# Add factor levels for sorted order
sorted_tfs_metabric <- intersect(sorted_tfs, unique(METABRIC_TF_activities_TNBC_Specific_long$tf))

METABRIC_TF_activities_TNBC_Specific_long <- METABRIC_TF_activities_TNBC_Specific_long %>%
  mutate(tf = factor(METABRIC_TF_activities_TNBC_Specific_long$tf, levels = sorted_tfs_metabric),
    color_group = ifelse(NES >= 0, "Positive", "Negative"))

METABRIC_TF_activities_TNBC_Specific_long <- METABRIC_TF_activities_TNBC_Specific_long[!is.na(METABRIC_TF_activities_TNBC_Specific_long$tf),]

# Plot
METABRIC_TNBC_Basal <- ggplot(METABRIC_TF_activities_TNBC_Specific_long, aes(x = NES, y = tf, color = color_group)) +
  geom_point(size = 0.10) +
  scale_color_manual(values = c("Positive" = "darkred", "Negative" = "navy")) +
  theme_classic() +
  labs(
       x = "NES",
       y = NULL) +
  theme(axis.text.y = element_text(size = 2.6),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none") +
  geom_vline(xintercept = 0, linetype="dashed", color = "black")

cairo_pdf(paste0("/results/visuals/Figure1/SupplementaryFigure2F_METABRIC_Basal", length(unique(METABRIC_TF_activities_TNBC_Specific_long$tf)), 
    "TRs.pdf"), width = 3.75, height = 8)
print(METABRIC_TNBC_Basal)
dev.off()

sessionInfo()
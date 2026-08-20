# Plot pearson correlation of signal intensity across all samples
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

parser$add_argument("-cohort", "--cohort", required=TRUE,
                    help="Cohort name (e.g., TCGA, METABRIC, PDX, CellLines)")
parser$add_argument("-pearson_cor", "--pearson_cor", required=TRUE,
                    help="Path to the pearson correlation output file")
parser$add_argument("-meta_data", "--meta_data", required=FALSE,
                    help="Path to meta data file relevant for CellLines")
parser$add_argument("-file_prefix", "--file_prefix", required=TRUE,
                    help="file prefix - SupplementaryFigureX pattern")
parser$add_argument("-visual_outdir", "--visual_outdir", default = "/results/visuals/Figure2",
                    help="Output directory for plot of pearson correlation")

args <- parser$parse_args()
print(args)
##################################################################################################
## For testing
#args <- list(cohort = "TCGA",
#    pearson_cor = "data/Figure2_TR_Chromatin/BulkATAC_SignalMatrix/TCGA_SM/Stats/Cor_7.rds",
#    visual_outdir = "visuals/Figure2_TR_Chromatin/")

#args <- list(cohort = "PDX",
#    pearson_cor = "data/Figure2_TR_Chromatin/BulkATAC_SignalMatrix/PDX/Stats/Cor_38.rds",
#    visual_outdir = "visuals/Figure2_TR_Chromatin/")

#args <- list(cohort = "CellLines",
#    pearson_cor = "data/Figure2_TR_Chromatin/BulkATAC_SignalMatrix/CellLines/Stats/Cor_44.rds",
#    meta_data = "data/Reference_Data/UHN_PDX_Cell_Lines/config/samples.tsv",
#    visual_outdir = "visuals/Figure2_TR_Chromatin/")

#print(args)
##################################################################################################
dir.create(args$visual_outdir, showWarnings=FALSE, recursive=TRUE)

# Read in the pearson correlation file
pearson_cor <- readRDS(args$pearson_cor)
print(pearson_cor)

col_fun <- colorRamp2(
  breaks = seq(from = 0, to = 1, length.out = 10),
  colors = viridis(10, option = "magma")   # or more colors, e.g. viridis(10)
)

if(length(args$meta_data) != 0){
    meta_data <- read.table(args$meta_data, header=TRUE, sep='\t')

    meta_data <- meta_data[meta_data$sample_name %in% rownames(pearson_cor),]
    meta_data$SequencingType <- ifelse(meta_data$Cohort == "Nergiz_se", "SingleEnd", "PairedEnd")
    meta_data <- meta_data[,c("sample_name", "SequencingType")]
    rownames(meta_data) <- meta_data$sample_name
    meta_data <- meta_data[rownames(pearson_cor),]
    meta_data <- meta_data %>% select(-c("sample_name"))
    
    column_anno <- HeatmapAnnotation(
        df  = meta_data[, "SequencingType", drop = FALSE],
        col = list(SequencingType = SequencingType),
        which = "column")
}


if(length(args$meta_data) != 0){
    pearson_correlation_heatmap <- ComplexHeatmap::Heatmap(
        pearson_cor,
        name = "Pearson Correlation",
        col  = col_fun,
        top_annotation = column_anno,
        column_names_gp = gpar(fontsize=6),
        row_names_gp = gpar(fontsize=6)
    )
} else{

    pearson_correlation_heatmap <- ComplexHeatmap::Heatmap(
        pearson_cor,
        name = "Pearson Correlation",
        col  = col_fun,
        column_names_gp = gpar(fontsize=6),
        row_names_gp = gpar(fontsize=6)
    )

}

cairo_pdf(paste0(args$visual_outdir, "/", args$file_prefix, "_", args$cohort, "_Pearson_Correlation_Signal_ConsensusPeakSet.pdf"),
    width = 7, height = 6)
draw(pearson_correlation_heatmap, annotation_legend_side = "bottom")
dev.off()

sessionInfo()
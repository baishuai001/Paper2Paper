options(scipen=999)
# Plot the number of called peaks
## Pertains to the PDX and primary tumours
library(ggplot2)
library(dplyr)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

parser$add_argument("-cohort", "--cohort", required=TRUE,
                    help="Cohort name (TCGA, PDX)")
parser$add_argument("-Number_Called_Peaks_file", "--Number_Called_Peaks_file", required=TRUE,
                    help="Path to the directory containing the narrowpeak files")
parser$add_argument("-file_prefix", "--file_prefix", required=TRUE,
                    help="file prefix - SupplementaryFigureX pattern")
parser$add_argument("-visual_outdir", "--visual_outdir", default = "/results/visuals/Figure2/",
                    help="Output directory for plot of number of called peaks")

args <- parser$parse_args()
print(args)
##################################################################################################
## For testing
#args <- list(cohort = "PDX",
#    Number_Called_Peaks_file = "data/Figure2_TR_Chromatin/PDX_ATAC/Called_Peaks/PDX_Number_Peaks_PDX_TNBC_ATAC.tsv",
#    visual_outdir = "visuals/Figure2_TR_Chromatin/")

#print(args)
##################################################################################################
dir.create(args$visual_outdir, showWarnings=FALSE, recursive=TRUE)

# Read in the number of peaks data frame
number_peaks_df <- read.table(args$Number_Called_Peaks_file, header=TRUE, sep='\t')
head(number_peaks_df)
colnames(number_peaks_df) <- c("Number_Peaks", "Sample")

number_peaks_df <- number_peaks_df[order(number_peaks_df$Number_Peaks, decreasing=TRUE),]
number_peaks_df$Sample <- factor(number_peaks_df$Sample, levels = number_peaks_df$Sample)

number_peaks_df$Number_Peaks_1000 <- number_peaks_df$Number_Peaks/1000

# Plot the number of peaks
number_peaks_plt <- number_peaks_df %>% ggplot(aes(x = Sample, y = Number_Peaks_1000, fill = "IDR")) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = replicates_idr_colours) +
    scale_y_continuous(expand = c(0,0), limits=c(0, max(number_peaks_df$Number_Peaks_1000) + 50),
        n.breaks=10) +
    theme(axis.text.x = element_text(angle=270, hjust=0, vjust=0.5, colour = "black", size=6),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") +
    labs(x = "", y = "Number of Peaks (x1000)")

outfile <- paste0(args$visual_outdir, "/", args$file_prefix, "_", args$cohort, "_Number_Of_Called_Peaks.pdf")

cairo_pdf(outfile)
print(number_peaks_plt)
dev.off()

sessionInfo()
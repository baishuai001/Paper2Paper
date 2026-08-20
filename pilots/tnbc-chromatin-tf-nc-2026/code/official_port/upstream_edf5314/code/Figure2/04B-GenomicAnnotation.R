# Summarize the genomic annotation of peaks across the cohorts
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

## Inputs
parser$add_argument("-cohort", "--cohort", nargs='+',
    default=c("TCGA"),
    help="Cohort name (e.g., TCGA, PDX, CellLines)")              
parser$add_argument("-consensus_peakset_dir", "--consensus_peakset_dir",
    default = "data/Figure2_TR_Chromatin/Consensus_Cohort/",
    help="Path to consensus peak set parent directory")

## Reference
parser$add_argument("-genome_ref", "--genome_ref", 
    default="/data/Figure2/GenomicAnnotation/Genome_Wide_GenomicAnnotation_500bp.tsv.gz",
    help="Path to the genome reference")

## Outputs
parser$add_argument("-outfile_data", "--outfile_data",
    default="/data/Figure2/GenomicAnnotation/Compiled_GenomicAnnotation_Distribution.tsv",
    help="Path to output file for genomic annotation")
parser$add_argument("-outfile_plot", "--outfile_plot",
    default="/results/visuals/Figure2/Compiled_GenomicAnnotation.pdf",
    help="Path to output plot summarizing genomic annotation")

args <- parser$parse_args()
print(args)
##################################################################################################
## Testing
#args <- list(cohort = c("TCGA"),
#    consensus_peakset_dir = "data/Figure2_TR_Chromatin/Consensus_Cohort/",
#    genome_ref = "data/Figure2_TR_Chromatin/GenomicAnnotation/Genome_Wide_GenomicAnnotation_500bp.tsv",
#    outfile_data = "data/Figure2_TR_Chromatin/GenomicAnnotation/Compiled_GenomicAnnotation_Distribution.tsv",
#    outfile_plot = "visuals/Figure2_TR_Chromatin/TCGA_Compiled_GenomicAnnotation.pdf")

#args <- list(cohort = c("PDX"),
#    consensus_peakset_dir = "data/Figure2_TR_Chromatin/Consensus_Cohort/",
#    genome_ref = "data/Figure2_TR_Chromatin/GenomicAnnotation/Genome_Wide_GenomicAnnotation_500bp.tsv",
#    outfile_data = "data/Figure2_TR_Chromatin/GenomicAnnotation/Compiled_GenomicAnnotation_Distribution.tsv",
#    outfile_plot = "visuals/Figure2_TR_Chromatin/PDX_Compiled_GenomicAnnotation.pdf")

#args <- list(cohort = c("CellLines"),
#    consensus_peakset_dir = "data/Figure2_TR_Chromatin/Consensus_Cohort/",
#    genome_ref = "data/Figure2_TR_Chromatin/GenomicAnnotation/Genome_Wide_GenomicAnnotation_500bp.tsv",
#    outfile_data = "data/Figure2_TR_Chromatin/GenomicAnnotation/Compiled_GenomicAnnotation_Distribution.tsv",
#    outfile_plot = "visuals/Figure2_TR_Chromatin/CellLines_Compiled_GenomicAnnotation.pdf")

#print(args)
##################################################################################################
dir.create(dirname(args$outfile_data), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(args$outfile_plot), showWarnings=FALSE, recursive=TRUE)

genomic_annotation_df <- data.frame(matrix(nrow=0, ncol=6))
colnames(genomic_annotation_df) <- c("Annotation", "Number_of_Peaks", "Total_Peaks_Sample", "Proportion_Peaks", "Sample", "Cohort")

# Read in the genome reference file
genome_ref <- read.table(args$genome_ref, header=TRUE, sep='\t')
head(genome_ref)

## merge intergenic and distal together
genome_ref$type <- gsub("intergenic", "Distal", genome_ref$type)
genome_ref$type <- gsub("distal", "Distal", genome_ref$type)
genome_ref$type <- gsub("promoter", "Promoter", genome_ref$type)
genome_ref$type <- gsub("exonic", "Exonic", genome_ref$type)
genome_ref$type <- gsub("intronic", "Intronic", genome_ref$type)

genome_ref_freq <- as.data.frame(table(genome_ref$type))
genome_ref_total <- sum(genome_ref_freq$Freq)
colnames(genome_ref_freq) <- c("Annotation", "Number_of_Peaks")
genome_ref_freq$Proportion_Peaks <- genome_ref_freq$Number_of_Peaks/genome_ref_total
genome_ref_freq$Sample <- "Genome"
genome_ref_freq$Cohort <- "Genome"
genome_ref_freq$Total_Peaks_Sample <- genome_ref_total

genomic_annotation_df <- rbind(genomic_annotation_df, genome_ref_freq)

# Read in all the consensus peak files
consensus_files <- list.files(paste0(args$consensus_peakset, "/", args$cohort, "_results/"), 
    pattern="PeakCalls_unionPeakSet.bed.gz", full.names=TRUE)

for(file in consensus_files){
    cohort <- gsub("", "", basename(dirname(file)))
    cohort <- gsub("_results", "", cohort)

    file <- read.table(file, header=TRUE, sep='\t')
    file_df <- as.data.frame(table(file$GroupReplicate, file$peakType))
    colnames(file_df) <- c("Sample", "Annotation", "Number_of_Peaks")
    file_df$Cohort <- cohort

    file_df <- file_df %>% 
        group_by(Sample) %>% 
        mutate(Total_Peaks_Sample = sum(Number_of_Peaks)) %>% 
        ungroup() %>% 
        mutate(Proportion_Peaks = Number_of_Peaks / Total_Peaks_Sample)

    genomic_annotation_df <- rbind(genomic_annotation_df, file_df)
}

genomic_annotation_df_plt <- genomic_annotation_df[genomic_annotation_df$Cohort != "Genome",]
genomic_annotation_df_plt$Annotation <- factor(genomic_annotation_df_plt$Annotation, levels = c("Distal", "Intronic", "Exonic", "Promoter"))

# Calculate the yintercepts and colors for the dashed lines
genome_lines <- genomic_annotation_df %>%
  filter(Cohort == "Genome") %>%
  select(Proportion_Peaks, Annotation)

# Plot
# Count samples per Cohort
sample_counts <- genomic_annotation_df_plt %>%
  group_by(Cohort) %>%
  summarise(n = n_distinct(Sample))

# Calculate statistics separately for each cohort
comparisons <- list(c("Distal", "Intronic"),
    c("Distal", "Exonic"),
    c("Distal", "Promoter"),
    c("Intronic", "Exonic"),
    c("Exonic", "Promoter"),
    c("Exonic", "Promoter"))

stat_test <- genomic_annotation_df_plt %>%
  wilcox_test(Proportion_Peaks ~ Annotation, 
              comparisons = comparisons) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Annotation")

# Create plot with stat_pvalue_manual
genomic_annotation_plt <- genomic_annotation_df_plt %>%
  ggplot(aes(x = Annotation, y = Proportion_Peaks)) +
  geom_boxplot(aes(fill = Annotation), outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_bw(base_size = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_color_manual(values = genomic_anno) +
  scale_fill_manual(values = genomic_anno) +
  geom_hline(data = genome_lines,
             aes(yintercept = Proportion_Peaks, color = Annotation),
             linetype = "dashed", show.legend = FALSE) +
  labs(x = "", y = "Proportion of peaks") +
  facet_grid(. ~ Cohort) +
  geom_text(data = sample_counts,
            aes(x = 4, y = 0.6, label = paste0("n = ", n)),
            inherit.aes = FALSE, size = 3.5) +
  stat_pvalue_manual(stat_test, 
                     label = "p.adj.signif",
                     hide.ns = TRUE,
                     tip.length = 0.01) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 12))


cairo_pdf(args$outfile_plot, height = 4, width = 6)
print(genomic_annotation_plt)
dev.off()

# Save the data files
## Plotting info
#write.table(genomic_annotation_df_plt, paste0(dirname(args$outfile_data), "/", args$cohort, #"_GenomicAnnotation_SummarizedInfo.tsv"),
#    col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

## Statistics
stat_test <- as.data.frame(stat_test)
stat_test_clean <- stat_test %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = ", ")))

#write.table(stat_test_clean, paste0(dirname(args$outfile_data), "/", args$cohort, #"_GenomicAnnotation_Pairwise_Stats.tsv"),
#    col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

sessionInfo()
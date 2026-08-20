# Plot number of IDR peaks (if more replicates)
library(ggplot2)
library(dplyr)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")
source("/code/Static_Scripts/Subset_TNBC_Basal_CellLines_Helper.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

parser$add_argument("-cohort", "--cohort", default = "CellLines",
                    help="Cohort name with replicates - CellLines")
parser$add_argument("-compiled_peaks", "--compiled_peaks", default = "/data/Figure2/IDR/IDR_Results/Compiled_IDR_AnnotatedPeaks_hg38.tsv.gz", help="Path to compiled IDR peaks")
parser$add_argument("-sample", "--sample", default = "/data/Figure2/IDR/replicate_prefixes.tsv",
                    help="Path to compiled IDR peaks")
parser$add_argument("-cellLine_PAM50_assignment", "--cellLine_PAM50_assignment", default = "/data/Reference_Data/CellLines_UHN_ATAC/BCa_CellLines_classification.tsv", help="Path to compiled IDR peaks")
parser$add_argument("-sample_inclusion", "--sample_inclusion", default = "/data/Reference_Data/UHN_PDX_Cell_Lines/config/samples.tsv", help="Path to the file containing samples to include")
parser$add_argument("-visuals_outdir", "--visuals_outdir", default = "/results/visuals/Figure2", help="Output directory for Figure2 related plots")
parser$add_argument("-cell_line_inclusion_replicate_metadata", "--cell_line_inclusion_replicate_metadata", default = "/data/Figure2/IDR/Included_cellLines_IDR.csv", help="Path to the replicate metadata for cell-lines with bulk ATAC")

args <- parser$parse_args()
print(args)
##################################################################################################
dir.create(args$visuals_outdir, recursive=TRUE, showWarnings=FALSE)

# Read in the compiled peaks
compiled_peaks <- read.table(args$compiled_peaks, header=TRUE, sep='\t')
head(compiled_peaks)
number_IDR_replicates <- as.data.frame(table(compiled_peaks$IDR_Replicate))
number_IDR_replicates$Var1 <- gsub("_peaks.filtered.narrowPeak", "", number_IDR_replicates$Var1)

# Read in the replicate info stuff
sample <- read.table(args$sample, header=TRUE, sep='\t')

# build a lookup from all three sample columns to Replicate
lookup <- setNames(
  rep(sample$Replicate, 3),
  c(sample$IDR_prefix, sample$Sample1, sample$Sample2)
)

# then map Var1 through that lookup
number_IDR_replicates$CellLine <- lookup[number_IDR_replicates$Var1]

na_idx <- is.na(number_IDR_replicates$CellLine)

number_IDR_replicates$CellLine[na_idx] <- gsub("_.*", "", number_IDR_replicates$Var1)[na_idx]

number_IDR_replicates <- number_IDR_replicates[order(number_IDR_replicates$Freq, decreasing = TRUE),]
number_IDR_replicates$Var1 <- factor(number_IDR_replicates$Var1, levels = (number_IDR_replicates$Var1))
number_IDR_replicates$Freq <- number_IDR_replicates$Freq/1000

number_IDR_replicates$Category <- ifelse(grepl("IDR", number_IDR_replicates$Var1) == TRUE,
    "IDR", "Replicate")

number_IDR_replicates <- number_IDR_replicates %>%
  group_by(CellLine) %>%
  mutate(
    # 1. label replicates by decreasing Freq
    Replicate = case_when(
      Category == "Replicate" ~ paste0(
        "Rep",
        dense_rank(desc(Freq[Category == "Replicate"]))[
          match(row_number(), which(Category == "Replicate"))
        ]
      ),
      TRUE ~ "IDR"
    ),
    # 2. numeric order: within cell line, Rep1>Rep2>...>IDR
    .order = case_when(
      Replicate == "Rep1" ~ 1,
      Replicate == "Rep2" ~ 2,
      Replicate == "Rep3" ~ 3,
      Replicate == "Rep4" ~ 4,
      Replicate == "IDR"  ~ 5,
      TRUE                ~ 6
    )
  ) %>%
  ungroup() %>%
  mutate(
    Replicate = factor(Replicate, levels = c("Rep1","Rep2","Rep3","Rep4","IDR"))
  )

# 3. global ordering for x-axis: by CellLine, then .order, then Freq
number_IDR_replicates <- number_IDR_replicates %>%
  arrange(CellLine, .order, desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))

number_IDR_replicates$Replicate <- factor(number_IDR_replicates$Replicate,
    levels = c(paste0("Rep", c(1:4)), "IDR"))

# Subset for the TNBC (Basal) cell-lines
meta_data_sub <- meta_data[meta_data$sample_name %in% number_IDR_replicates$Var1,]
unique(meta_data_sub$CellLine)

meta_data_sub <- meta_data[meta_data$PAM50 == "Basal" & 
    meta_data$SCMOD2 == "Basal",]
unique(meta_data_sub$CellLine)

number_IDR_replicates_TNBC_Basal_Rep <- number_IDR_replicates[number_IDR_replicates$Var1 %in% meta_data_sub$sample_name,]
number_IDR_replicates_TNBC_Basal_IDR <- number_IDR_replicates[number_IDR_replicates$Category == "IDR" &
    number_IDR_replicates$CellLine %in% number_IDR_replicates_TNBC_Basal_Rep$CellLine,]

number_IDR_replicates_TNBC_Basal <- rbind(number_IDR_replicates_TNBC_Basal_Rep, number_IDR_replicates_TNBC_Basal_IDR)

## Inclusion ##
sample_inclusion <- read.table(args$sample_inclusion, header=TRUE, sep='\t')
sample_inclusion <- sample_inclusion[sample_inclusion$Inclusion == "Yes",]
unique(sample_inclusion$Replicate)

number_IDR_replicates_TNBC_Basal_included_Replicate <- number_IDR_replicates_TNBC_Basal[number_IDR_replicates_TNBC_Basal$Var1 %in% sample_inclusion$sample_name,]
number_IDR_replicates_TNBC_Basal_included_IDR <- number_IDR_replicates_TNBC_Basal[number_IDR_replicates_TNBC_Basal$Category == "IDR",]

# all Var1 values currently included
keep_var1 <- number_IDR_replicates_TNBC_Basal_included_Replicate$Var1

# rows of `sample` whose Sample1 and Sample2 are both in Var1
pairs_included <- sample %>%
  filter(Sample1 %in% keep_var1,
         Sample2 %in% keep_var1)

# now subset the IDR dataframe to those two replicates only
number_IDR_reps_subset <- number_IDR_replicates_TNBC_Basal_included_IDR[number_IDR_replicates_TNBC_Basal_included_IDR$Var1 %in% pairs_included$IDR_prefix,]

number_IDR_replicates_TNBC_Basal_included <- rbind(number_IDR_replicates_TNBC_Basal_included_Replicate, number_IDR_reps_subset)

singular_PE_bulkATAC <- c("HCC1806_2_S21", "HCC3153_1_S8", "Hs-578-T_2_S1", "MDA468_2_S11")

# Read in the re-assignment of replicate number
cell_line_inclusion_replicate_metadata <- read.table(args$cell_line_inclusion_replicate_metadata,
  header=TRUE, sep=',')

colnames(number_IDR_replicates_TNBC_Basal_included)[1] <- "SampleID"
number_IDR_replicates_TNBC_Basal_included <- merge(number_IDR_replicates_TNBC_Basal_included, cell_line_inclusion_replicate_metadata,
  by = "SampleID")

head(number_IDR_replicates_TNBC_Basal_included)

### Plot ###
IDR_plt <- number_IDR_replicates_TNBC_Basal_included %>% ggplot(aes(x = SampleID, y = Freq, fill = Replicate.x)) +
    scale_y_continuous(expand=c(0,0), n.breaks = 10, limits=c(0,ceiling(max(number_IDR_replicates$Freq) + 10))) +
    scale_fill_manual(values = replicates_idr_colours) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Number of peaks (x1000)", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=270, hjust = 0, vjust = 0.5, color = "black", size = 8),
        axis.text.y = element_text(color = "black"),
        legend.position="bottom") +
    geom_hline(yintercept=50, linetype="dashed")

cairo_pdf(paste0(args$visuals_outdir, "/SupplementaryFigure3C_IDR_ATAC_Peaks_Replicate_TNBC_Basal_CellLines_Inclusion.pdf"))
print(IDR_plt)
dev.off()

number_IDR_replicates_TNBC_Basal_included <- number_IDR_replicates_TNBC_Basal_included[,
  c("SampleID", "Freq", "CellLine.x", "Category", "Replicate.x", "Include")]
colnames(number_IDR_replicates_TNBC_Basal_included) <- c("SampleID", "Number_Of_Peaks_x1000",
  "CellLine", "Category", "Replicate", "Include")

#write.table(number_IDR_replicates_TNBC_Basal_included,
#  paste0("data/Figure2_TR_Chromatin/SupplementaryTables/CellLines_Number_Peaks_Replicate_IDR.tsv"),
#  col.names=TRUE,
#  row.names=FALSE,
#  sep='\t',
#  quote=FALSE)

sessionInfo()
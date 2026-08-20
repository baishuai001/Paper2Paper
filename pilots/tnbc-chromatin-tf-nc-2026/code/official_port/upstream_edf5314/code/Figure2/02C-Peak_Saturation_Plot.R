library(ggplot2)
library(data.table)
library(scales)
library(argparse)

set.seed(2025)

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

parser$add_argument("-cohort", "--cohort", required=TRUE,
                    help="Cohort name (e.g., TCGA, METABRIC, PDX, CellLines)")
parser$add_argument("-file_prefix", "--file_prefix", required=TRUE,
                    help="File name - Format SupplementaryFigureX")
parser$add_argument("-peak_saturation_dir", "--peak_saturation_dir", default = "/data/Figure2/peak_saturation_curves",
                    help="Output directory for saturation curve results")
parser$add_argument("-outdir_visuals", "--outdir_visuals", default = "/results/visuals/Figure2",
                    help="Output directory for saturation curve plots")

args <- parser$parse_args()
print(args)
##################################################################################################
dir.create(args$outdir_visuals, showWarnings=FALSE, recursive = TRUE)

# Read in the peak saturation file
peak_saturation <- list.files(args$peak_saturation_dir, pattern = ".tsv", full.names=TRUE)
peak_saturation <- peak_saturation[grepl(args$cohort, peak_saturation)]
print(peak_saturation)

bootstraps_long <- peak_saturation[grep("_peak-saturation.bootstraps", peak_saturation)]
bootstraps_long <- bootstraps_long[grepl(args$cohort, bootstraps_long)]
bootstraps_long <- as.data.frame(fread(bootstraps_long, header=TRUE, sep='\t'))

sample_ests <- peak_saturation[grep("_peak-saturation.model-estimates.tsv", peak_saturation)]
sample_ests <- sample_ests[grepl(args$cohort, sample_ests)]
sample_ests <- as.data.frame(fread(sample_ests, header=TRUE, sep='\t'))

saturation_ests <- peak_saturation[grep("_peak-saturation.saturation-estimates.tsv", peak_saturation)]
saturation_ests <- saturation_ests[grepl(args$cohort, saturation_ests)]
saturation_ests <- as.data.frame(fread(saturation_ests, header=TRUE, sep='\t'))

saturation_model <- readRDS(paste0(args$peak_saturation_dir, "/", args$cohort, "_saturation_model.rds"))

# Pre-compute summary values for annotations
label_x <- mean(saturation_ests$N_Samples[saturation_ests$Frac_Saturation >= 0.95])
label_y <- 0.9 * saturation_model$model$Asym
peaks_detected_y <- saturation_model$peaks$identified
asymptote_y <- saturation_model$model$Asym
peaks_label <- paste0("Peaks Detected: ", round(saturation_model$peaks$identified), 
                      " (", round(100 * saturation_model$Frac_Saturation, 1), "%)")
asym_label <- paste("Estimated Number of Peaks:", round(saturation_model$model$Asym))

gg <- ggplot() +
    geom_point(
        data = bootstraps_long,
        mapping = aes(x = N_Samples, y = N_Peaks, colour = N_Samples),
        position = position_jitter(width = 0.2, height = 0),
        alpha = 0.1
    ) +
    geom_boxplot(
        data = bootstraps_long,
        mapping = aes(x = N_Samples, y = N_Peaks, group = N_Samples, fill = N_Samples),
        outlier.shape = NA,
        alpha = 0.5
    ) +
    # asymptotic curve fit
    geom_path(
        data = sample_ests,
        aes(x = N_Samples, y = Est_N_Peaks),
        linetype = "dashed",
        size = 0.5,
        color = "darkgrey"
    ) +
    # annotate model fit label
    annotate("text",
        x = label_x,
        y = label_y,
        label = "Model Fit",
        vjust = 0,
        hjust = 0.5
    ) +
    # add asymptotic value
    geom_hline(
        yintercept = c(asymptote_y, peaks_detected_y),
        linetype = "dashed",
        colour = "black"
    ) +
    # labels for horizontal lines
    annotate("text",
        x = 0,
        y = peaks_detected_y,
        label = peaks_label,
        vjust = -1,
        hjust = 0,
        colour = "black"
    ) +
    annotate("text",
        x = 0,
        y = asymptote_y,
        label = asym_label,
        vjust = -1,
        hjust = 0,
        colour = "black"
    ) +
    # add bars for number of samples required to reach saturation
    geom_vline(
        data = saturation_ests,
        aes(xintercept = round(N_Samples)),
        linetype = "dashed",
        colour = "black"
    ) +
    # labels for vertical lines
    geom_text(
        data = saturation_ests,
        aes(
            x = round(N_Samples),
            y = 0,
            label = paste0(
                round(100 * Frac_Saturation),
                "%: ",
                round(N_Samples),
                " samples"
            )
        ),
        vjust = -0.5,
        hjust = 0,
        angle = 90,
        colour = "black"
    ) +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    scale_x_continuous(
        name = "Number of samples",
        limits = c(
            0,
            min(saturation_ests$N_Samples[saturation_ests$Frac_Saturation == 0.99] + 1, 100)
        )
    ) +
    scale_y_continuous(
        name = "Number of unique peaks",
        limits = c(0, 1.05 * saturation_model$model$Asym),
        labels = scales::comma
    ) +
    guides(fill = "none", colour = "none") +
    theme_minimal(base_size=8) +
    theme(
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, color = "black")
    ) +
    ggtitle(paste(args$cohort, " ATAC-seq"))

cairo_pdf(paste0(args$outdir_visuals, "/", args$file_prefix, "_", args$cohort, "_peak_saturation.pdf"), height = 4, width = 3.5)
print(gg)
dev.off()

sessionInfo()
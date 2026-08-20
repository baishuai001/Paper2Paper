# Plot key examples for top TRs
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(forcats)
library(stringr)
library(tidyr)
library(argparse)

source("/code/Static_Scripts/plotting_aesthetics.R")

#### Collect all the arguments ####################################################################
parser <- ArgumentParser()

## Inputs
parser$add_argument("-JASPAR_compiled_motif_enrichment", "--JASPAR_compiled_motif_enrichment",
    default = "/data/Figure2/Motif/Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
    help="Compiled JASPAR motif enrichment")
parser$add_argument("-CISBP_compiled_motif_enrichment", "--CISBP_compiled_motif_enrichment",
    default = "/data/Figure2/Motif/Combined_MotifEnrichment/CISBP_Combined_MotifEnrichment_TNBC.tsv",
    help="Compiled CISBP motif enrichment")
parser$add_argument("-JASPAR_TR_MotifEnrichment", "--JASPAR_TR_MotifEnrichment",
    default = "/data/Figure2/Motif/JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv",
    help="Sample groups motif enrichment")
parser$add_argument("-Motif_HC_TRs", "-Motif_HC_TRs",
    default = "/data/Figure2/Motif/JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv",
    help="Motif HC-TFs")

## Outputs
parser$add_argument("-visuals_outdir", "--visuals_outdir",
    default="/results/visuals/Figure2/",
    help="Path to output plot")

args <- parser$parse_args()
print(args)
##################################################################################################
dir.create(args$visuals_outdir, showWarnings=FALSE, recursive=TRUE)

#args <- list(JASPAR_compiled_motif_enrichment = "data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/JASPAR_Combined_MotifEnrichment_TNBC.tsv",
#    CISBP_compiled_motif_enrichment = "data/Figure2_TR_Chromatin/Motif/Combined_MotifEnrichment/CISBP_Combined_MotifEnrichment_TNBC.tsv",
#    Motif_HC_TRs = "data/Figure2_TR_Chromatin/Motif/JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv",
#    visuals_outdir = "visuals/Figure2_TR_Chromatin/")
#print(args)

# Read in HC-TRs
Motif_HC_TRs <- read.csv(args$Motif_HC_TRs, header=TRUE, sep='\t')
print(head(Motif_HC_TRs))

TR_list_list_mat <- read.table(args$JASPAR_TR_MotifEnrichment, header=TRUE, sep='\t')

# Read in the combined TR info for motif enrichment
compiled_motif_enrichment <- read.table(args$JASPAR_compiled_motif_enrichment, header=TRUE, sep='\t')
compiled_motif_enrichment$neglog10_padj <- -log10(compiled_motif_enrichment$qvalue_Benjamini)
compiled_motif_enrichment$TR <- toupper(compiled_motif_enrichment$Motif.Name)
compiled_motif_enrichment <- compiled_motif_enrichment[compiled_motif_enrichment$TR %in% rownames(Motif_HC_TRs),]

compiled_motif_enrichment_CISBP <- read.table(args$CISBP_compiled_motif_enrichment, header=TRUE, sep='\t')
compiled_motif_enrichment_CISBP$neglog10_padj <- -log10(compiled_motif_enrichment_CISBP$qvalue_Benjamini)
compiled_motif_enrichment_CISBP$TR <- toupper(compiled_motif_enrichment_CISBP$Motif.Name)
compiled_motif_enrichment_CISBP <- compiled_motif_enrichment_CISBP[compiled_motif_enrichment_CISBP$TR %in% rownames(Motif_HC_TRs[Motif_HC_TRs$hc_group == "CISBP only",]),]

compiled_motif_enrichment <- rbind(compiled_motif_enrichment, compiled_motif_enrichment_CISBP)

#compiled_motif_enrichment <- compiled_motif_enrichment[compiled_motif_enrichment$TR %in%
#    c("FOSL1", "ETV4", "TCF7L1", "FOXC1", "MYC", "NFIL3", "KLF5", "SOX10", "OTX1", "NFE2L3"),]

head(compiled_motif_enrichment)

df <- compiled_motif_enrichment

# Harmonize group labels (adjust mapping to your data)
df <- df %>%
  mutate(
    Sample_Group = recode(Sample_Group,
      "TCGA"="TCGA", "PDX"="PDX",
      "CellLine"="CellLines"
    ),
    Sample_Group = factor(Sample_Group, levels = c("TCGA","PDX","CellLines"))
  )

# --- Parse percentages if they are strings like "95.42%" ---
parse_pct <- function(x) {
  if (is.numeric(x)) return(x/100)
  as.numeric(str_replace(x, "%", "")) / 100
}

p_tgt <- parse_pct(df$Percent_TargetSequencesMotif)
p_bkg <- parse_pct(df$Percent_TargetSequencesMotif.1)

# JASPAR-like outputs often provide counts for hits (k1, k0) but not totals (n1, n0).
# We can back-calc totals using: n ≈ k / p
k1 <- df$Number_TargetSequencesMotif
k0 <- df$Number_BackgroundSequencesMotif

# Guard: if percentages missing, fall back to approximate totals using min denominator 1
n1 <- ifelse(!is.na(p_tgt) & p_tgt > 0, round(k1 / p_tgt), NA)
n0 <- ifelse(!is.na(p_bkg) & p_bkg > 0, round(k0 / p_bkg), NA)

# If any NA remain (no % available), you can provide known totals here, or infer from file metadata.
# As a safe fallback, assume totals are the same across motifs within a sample context
# by taking the max observed denominator among that context:
fill_totals <- function(n_vec, k_vec, grp) {
  d <- data.frame(n = n_vec, k = k_vec, g = grp)
  d %>%
    group_by(g) %>%
    mutate(n_filled = ifelse(is.na(n), max(k, na.rm = TRUE) / max(k / pmax(n, 1), na.rm = TRUE), n)) %>%
    pull(n_filled)
}
# If you are confident n1/n0 are recovered well, you can skip the filler above.

# --- Haldane–Anscombe correction (adds 0.5 to all cells of the 2x2 table) ---
# 2x2: [k1, n1-k1; k0, n0-k0]
eps <- 0.5

# Ensure denominators > 0
n1 <- ifelse(is.na(n1), k1 + 1, n1)
n0 <- ifelse(is.na(n0), k0 + 1, n0)

# Compute LOR and its standard error (Wald)
lor <- log2( ((k1 + eps) / (n1 - k1 + eps)) / ((k0 + eps) / (n0 - k0 + eps)) )
# Var(log(OR)) ~ 1/a + 1/b + 1/c + 1/d (with natural log); convert to log2 by dividing by ln(2).
var_log_or <- 1/(k1 + eps) + 1/(n1 - k1 + eps) + 1/(k0 + eps) + 1/(n0 - k0 + eps)
se_lor <- sqrt(var_log_or) / log(2)  # because lor = log2(OR) = ln(OR)/ln(2)

df_lor <- df %>%
  mutate(
    k1 = k1, n1 = n1, k0 = k0, n0 = n0,
    LOR = lor,
    SE_LOR = se_lor,
    LOR_low = LOR - 1.96 * SE_LOR,
    LOR_high = LOR + 1.96 * SE_LOR
  )

# Order TRs by overall median LOR (or any custom order you prefer)
order_df <- df_lor %>%
  group_by(Motif.Name) %>%
  summarise(med_LOR = median(LOR, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_LOR))

df_lor <- df_lor %>%
  mutate(Motif.Name = factor(Motif.Name, levels = order_df$Motif.Name))

#### Plot the log odds ratio and FDR
# Recompute mean_lor with total samples per group for accurate Pct
total_per_group <- df_lor %>%
  group_by(Sample_Group) %>%
  summarise(total = n_distinct(Sample), .groups = "drop")

mean_lor <- df_lor %>%
  group_by(TR, Sample_Group) %>%
  summarise(
    Mean_LOR  = mean(LOR, na.rm = TRUE),
    SD_LOR    = sd(LOR, na.rm = TRUE),
    n_samples = n_distinct(Sample),
    .groups   = "drop"
  ) %>%
  left_join(total_per_group, by = "Sample_Group") %>%
  mutate(Pct_samples = (n_samples / total) * 100)

# Reshape TR_list_list_mat to long format for joining
sig_long <- TR_list_list_mat %>%
  tibble::rownames_to_column("TR") %>%
  pivot_longer(cols = c("TCGA", "PDX", "CellLines"),
               names_to  = "Sample_Group",
               values_to = "sig_cohort") %>%
  mutate(sig_cohort = as.logical(sig_cohort))

# Join into mean_lor
mean_lor <- mean_lor %>%
  left_join(sig_long %>% select(TR, Sample_Group, sig_cohort),
            by = c("TR", "Sample_Group"))

# Plot with SD_LOR as size
Motif_logOdds_plt <- ggplot(mean_lor, aes(x = Mean_LOR, y = reorder(TR, Mean_LOR))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(
    aes(colour    = Sample_Group,
        size      = SD_LOR,
        shape     = sig_cohort,
        alpha     = sig_cohort,
        stroke    = sig_cohort)
  ) +
  scale_shape_manual(
    values = c("TRUE" = 19, "FALSE" = 1),
    name   = "Enriched in\n≥50% of samples",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 0.9, "FALSE" = 0.6),
    guide  = "none"
  ) +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  scale_discrete_manual(
    "stroke",
    values = c("TRUE" = 0.3, "FALSE" = 1.2),  # thicker border for hollow
    guide  = "none"
  ) +
  scale_size_continuous(name = "SD of LOR", range = c(6, 1)) +
  theme_bw() +
  labs(
    x      = "Mean Log2 Odds Ratio\n(target vs. background)",
    y      = "Transcription factor",
    colour = "Cohort"
  ) +
  theme(
    axis.text    = element_text(colour = "black", size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

cairo_pdf(paste0(args$visuals_outdir, "/Figure2E_Motif_Enrichment_HC-TRs.pdf"),
    height=8, width = 5)
print(Motif_logOdds_plt)
dev.off()

# Plot
pd <- position_dodge(width = 0.6)

# Make sure levels are consistent and define significance
df_lor <- df_lor %>%
  mutate(
    Sample_Group = fct_relevel(Sample_Group, "TCGA", "PDX", "CellLines"),
    Motif.Name   = fct_inorder(Motif.Name),
    sig          = qvalue_Benjamini < 0.05  # or (LOR_low > 0) | (LOR_high < 0)
  ) %>%
  filter(!is.na(Sample_Group))

col_grp <- cohort_colours
col_ns  <- "#F28E2B"  # orange for non‑sig

ns_lab <- "Non significant"

df_lor_sub <- df_lor[df_lor$TR %in% c("DMRTA2", "FOXN2", "E2F3", "VAX2", "ETV6", "FOXP4"),]

plt <- ggplot(df_lor_sub, aes(x = Sample_Group, y = LOR)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +

  geom_boxplot(aes(color = Sample_Group), width = 0.65,
               outlier.shape = NA, size = 0.6) +

  ggbeeswarm::geom_quasirandom(
    data = ~ dplyr::filter(.x, sig),
    aes(color = Sample_Group),
    width = 0.12, varwidth = FALSE, size = 1.8, alpha = 0.9
  ) +

  ggbeeswarm::geom_quasirandom(
    data = ~ dplyr::filter(.x, !sig),
    aes(color = ns_lab),
    width = 0.12, varwidth = FALSE, size = 1.8, alpha = 0.9
  ) +

  scale_color_manual(
    values = c(col_grp, setNames(col_ns, ns_lab)),
    breaks = c(names(col_grp), ns_lab)
  ) +
  labs(x = NULL, y = "Log2 Odds Ratio (LOR)\n(target vs background)", color = NULL,
       title = "") +
  facet_grid(Motif.Name ~ ., scales = "free_y", space = "free_y") +  # flipped: rows now
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),       # was major.x
    axis.text.y  = element_blank(),             # was axis.text.x
    axis.text.x  = element_text(colour = "black"),  # was axis.text.y
    axis.ticks.y = element_blank(),             # was axis.ticks.x
    strip.text.y = element_text(face = "bold", size = 12, angle = 0),  # horizontal strip labels
    strip.background = element_blank(),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    panel.spacing.y = unit(1.1, "lines")        # was panel.spacing.x
  )

cairo_pdf(paste0(args$visuals_outdir, "/Figure2E_Case_Examples_HC-TRs.pdf"), height=7.5, width = 5)
print(plt)
dev.off()

# Pairwise Wilcoxon per TR (post-hoc)
library(rstatix)

wilcox_results <- df_lor_sub %>%
  group_by(Motif.Name) %>%
  wilcox_test(LOR ~ Sample_Group, p.adjust.method = "BH") %>%
  add_significance("p.adj")

print(wilcox_results)

sessionInfo()
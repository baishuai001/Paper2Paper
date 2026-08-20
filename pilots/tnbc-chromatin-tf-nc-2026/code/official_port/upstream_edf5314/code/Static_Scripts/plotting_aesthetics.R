# Plotting aesthetics

## Aesthetics related to Figure 1
viper_cutoff_colours <- list(VIPER_TR_Category = 
    c("TNBC" = "#435773", 
    "Non-TNBC" = "#9bc1e5",
    "Non-Significant" = "black",
    "LogFC Non-Signifcant" = "#9d9b9b",
    "LogFC Significant" = "#9d9b9b",
    "VIPER_Significant" = "#dcdcdc"))

clinical_annot_colours <- list(
  PAM50 = c(
    Basal = "#49a12e",
    LumA = "#2836a0",
    LumB = "#85a6d7",
    HER2 = "#a02774",
    Normal = "#abafaf"
  ),
  ER = c(
    Positive = "black",
    Negative = "white",
    "NA" = "#981e39"
  ),
  PR = c(
    Positive = "black",
    Negative = "white",
    "NA" = "#981e39"
  ),
  HER2 = c(
    Positive = "black",
    Negative = "white",
    "NA" = "#981e39"
  ),
  VitalSt = c(
    Alive = "#D99AB1",
    Dead = "#8C3048",
    "NA" = "#e1e1dd"
  ),
  Age = c(
    "<=40" = "#d7afd2",
    "(40-70)" = "#a34d9d",
    ">=70" = "#560f55",
    "NA" = "#e1e1dd"
  ),
  SCMOD2 = c(
    Basal = "#49a12e",
    LumA = "#2836a0",
    LumB = "#85a6d7",
    HER2 = "#a02774"
  ),
  IntClust = c(
    iC1 = "#85a6d7",
    iC2 = "#325192",
    iC3 = "#2836a0",
    iC4 = "#f2a039",
    iC5 = "#a02774",
    iC6 = "#e1e555",
    iC7 = "#b9b037",
    iC8 = "#53bab0",
    iC9 = "#875591",
    iC10 = "#49a146"
  ),
  Cancer_Type = c(
    Primary = "#beb7ff",
    Metastatic = "#b20037"
  ),
  HistologicalType = c(
    IDC = "#fffab4",
    ILC = "#d4d34f",
    Other = "#919136",
    Mixed = "#b0b0b0",
    MucC = "#828282",
    MetC = "#787979",
    MedC = "#5f6060",
    NA_IC = "#0b0b0b"
  ),
  Lehmann = c(BL1 = "#2b6519",
  BL2 = "#49a146",
  IM = "#79f271",
  LAR = "#afd7a5",
  M = "#73a12f",
  MSL = "#a5d360",
  "Possible ER+" = "#85a6d7",
  UNS = "#e1e1e1")
)

# Sequencing Type Colours
SequencingType <- c(SingleEnd = "#D9BBA9", PairedEnd = "#3E498C")

# Replicate peaks colours
replicates_idr_colours <- c(Rep1 = "#DD9FA8",
  Rep2 = "#C16E8A",
  Rep3 = "#BBADD9",
  Rep4 = "#745570",
  IDR = "#341F36")

## Genomic Annotation colours
genomic_anno <- c("Distal" = "#798774",
  "Intronic" = "#D4947D",
  "Exonic" = "#96B9D9",
  "Promoter" = "#C23E34")

## Motif Enrichment Cohort Colours
HC_TR_Motif <- c("CISBP only" = "#DA779A",
  "JASPAR only" = "#ADD7D6",
  "Both" = "#B87EF2",
  "Not HC-TR" = "grey")

cohort_colours <- c("TCGA" = "#7b007c",
  "PDX" = "#6c8438",
  "CellLines" = "#147574",
  "NA" = "white",
  "No Motif" = "grey")

NES_category_colours <- c("<0" = "#c9cdd4",
  "0-1" = "#b36082",
  "1-3" = "#982c54",
  ">3" = "#731b42")

## Colours for activating and repressing interactions
activating_repressing_colours <- c("Activated Shared Targets" = "#981111",
  "Activated Non-Shared Targets" = "#efc8c7",
  "Repressed Shared Targets" = "#454ca0",
  "Repressed Non-Shared Targets" = "#c9c4e2")

## TR activity sample categories
sample_category_colours <- c("TNBC" = "#965DA6",
  "Non-TNBC" = "#99E2F2")

## Survival analysis
TR_cox_proportional_hazards <- c(Risk = "#A60311",
  Protective = "#04588C",
  "Not Significant" = "#BFBFBF")

TR_Kaplan_Meyer <- c("Low" = "#1b9e77",
  "High" = "#e7298a")

### Venn Diagrams
VennDiagram_Prognosis_colours <- list()

## Concordance Index
CI_colours <- c(Positive_Correlated_CI_Significant = "#A60311",
  Negative_Correlated_CI_Significant = "#04588C",
  Non_Significant = "#BFBFBF")

## Drug VennDiagrams
drug_dataset_colours <- c(GRAY = "#72779E",
  CTRPv2 = "#FFD1FF",
  GDSC1000 = "#82BD68")

## AAC Concordance Index colours
AAC_concordance_colours_Low <- c(GRAY = "#9bbcff", 
  GDSC1000 = "#2f6df6",
  CTRPv2   = "#1c3fb3")

AAC_concordance_colours_High <- c(GRAY = "#FCAE91",
    GDSC1000 = "#FB6A4A",
    CTRPv2   = "#CB181D")

sessionInfo()
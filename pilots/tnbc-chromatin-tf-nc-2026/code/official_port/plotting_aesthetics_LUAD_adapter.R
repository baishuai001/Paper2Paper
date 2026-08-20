# Cancer-label and annotation-key adapter for the official plotting palette.
# No colour value, font, theme, size, scale, or geometry is redesigned here.

official_aesthetics <- Sys.getenv("PAPER2PAPER_OFFICIAL_AESTHETICS", unset = "")
if (!nzchar(official_aesthetics) || !file.exists(official_aesthetics)) {
  stop("PAPER2PAPER_OFFICIAL_AESTHETICS must name the frozen edf5314 aesthetics file.")
}
source(official_aesthetics)

# Same two colours as TNBC/Non-TNBC in the official capsule; labels only.
viper_cutoff_colours$VIPER_TR_Category <- c(
  "LUAD" = "#435773",
  "LUSC" = "#9bc1e5",
  "Non-Significant" = "black",
  "LogFC Non-Signifcant" = "#9d9b9b",
  "LogFC Significant" = "#9d9b9b",
  "VIPER_Significant" = "#dcdcdc"
)
sample_category_colours <- c("LUAD" = "#965DA6", "LUSC" = "#99E2F2")

# LUAD-specific annotation keys reuse colours already present in the official
# breast-cancer palette.  Values are not newly selected or interpolated.
clinical_annot_colours$HistologicalType <- c(
  LUAD = "#fffab4", LUSC = "#d4d34f", Unavailable = "#e1e1dd"
)
clinical_annot_colours$Sex <- c(
  Female = "#D99AB1", Male = "#2836a0", Unavailable = "#e1e1dd"
)
clinical_annot_colours$Stage <- c(
  I = "#85a6d7", II = "#325192", III = "#2836a0", IV = "#a02774",
  Unavailable = "#e1e1dd"
)
clinical_annot_colours$Smoking <- c(
  Never = "#afd7a5", Former = "#73a12f", Current = "#2b6519",
  Unavailable = "#e1e1dd"
)
clinical_annot_colours$PublishedSubtype <- c(
  "LUAD: TRU" = "#49a12e", "LUAD: PI" = "#f2a039",
  "LUAD: PP" = "#875591", Unavailable = "#e1e1dd"
)
clinical_annot_colours$MolecularRecord <- c(
  Reported = "#49a146", Unavailable = "#e1e1dd"
)
clinical_annot_colours$Lehmann <- c(
  LUAD = "#2b6519", LUSC = "#49a146", Unavailable = "#e1e1dd"
)

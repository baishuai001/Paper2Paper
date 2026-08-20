# Shared visual grammar transcribed from the TNBC Code Ocean capsule.
#
# This file deliberately contains only presentation constants and small plot
# helpers.  LUAD wrappers must preserve the scientific object used by the
# corresponding TNBC panel; changing a palette is never a substitute for
# changing an incorrect statistic or graph construction rule.

suppressPackageStartupMessages({
  library(ggplot2)
})

# The capsule requested Helvetica.  The analysis server does not ship that
# proprietary face, so use the metrically compatible Liberation Sans unless a
# licensed Helvetica/Arial family is explicitly supplied at run time.  One
# family is used throughout to prevent Cairo from silently mixing Nimbus Roman
# and Nimbus Sans in the same figure.
anchor_font_family <- Sys.getenv("PAPER2PAPER_FONT_FAMILY", "Liberation Sans")

anchor_colours <- list(
  phenotype = c(LUAD = "#435773", LUSC = "#9bc1e5"),
  discovery = c(
    nonsignificant = "#000000",
    expression_only = "#9d9b9b",
    activity_only = "#dcdcdc",
    LUAD = "#435773",
    LUSC = "#9bc1e5"
  ),
  sample_category = c(LUAD = "#965DA6", LUSC = "#99E2F2"),
  atac_system = c(Patient = "#7b007c", PDX = "#6c8438", CellLine = "#147574"),
  motif_database = c(CISBP = "#DA779A", JASPAR = "#ADD7D6", Both = "#B87EF2"),
  nes_bin = c(`<0` = "#c9cdd4", `0-1` = "#b36082", `1-3` = "#982c54", `>3` = "#731b42"),
  regulon = c(
    activated_shared = "#981111",
    activated_private = "#efc8c7",
    repressed_shared = "#454ca0",
    repressed_private = "#c9c4e2"
  ),
  prognosis = c(risk = "#A60311", protective = "#04588C", nonsignificant = "#BFBFBF"),
  km = c(Low = "#1b9e77", High = "#e7298a"),
  drug_dataset = c(GRAY = "#72779E", CTRPv2 = "#FFD1FF", GDSC1000 = "#82BD68"),
  drug_low_ci = c(GRAY = "#9bbcff", GDSC1000 = "#2f6df6", CTRPv2 = "#1c3fb3"),
  drug_high_ci = c(GRAY = "#FCAE91", GDSC1000 = "#FB6A4A", CTRPv2 = "#CB181D"),
  replicate_idr = c(
    Rep1 = "#DD9FA8", Rep2 = "#C16E8A", Rep3 = "#BBADD9",
    Rep4 = "#745570", IDR = "#341F36"
  ),
  genomic_annotation = c(
    Distal = "#798774", Intronic = "#D4947D",
    Exonic = "#96B9D9", Promoter = "#C23E34"
  ),
  sequencing_type = c(SingleEnd = "#D9BBA9", PairedEnd = "#3E498C"),
  mki67 = c(
    significant = "#C0392B", nonsignificant = "#B3B3B3",
    low = "#FFFFFF", mid = "#FCBBA1", high = "#C0392B",
    missing = "#EBEBEB"
  )
)

anchor_theme_classic <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = anchor_font_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

anchor_theme_bw <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = anchor_font_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      panel.border = element_rect(linewidth = 0.45, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

anchor_theme_void <- function(base_size = 12) {
  theme_void(base_size = base_size, base_family = anchor_font_family) +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}

anchor_panel_tag_theme <- theme(
  plot.tag = element_text(
    family = anchor_font_family, face = "bold", size = 16,
    colour = "black", hjust = 0, vjust = 1
  ),
  plot.tag.position = c(0, 1)
)

anchor_cairo_pdf <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = anchor_font_family,
    onefile = TRUE,
    ...
  )
}

anchor_save_plot <- function(plot, filename, width, height, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    family = anchor_font_family,
    limitsize = FALSE,
    ...
  )
}

anchor_assert_columns <- function(x, required, object_name = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        object_name,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

anchor_not_supported_panel <- function(title, reason) {
  ggplot() +
    annotate("segment", x = 0.08, xend = 0.92, y = 0.70, yend = 0.70,
             linewidth = 0.5, colour = "#6F7478") +
    annotate("text", x = 0.08, y = 0.58, hjust = 0, vjust = 1,
             label = title, family = anchor_font_family, fontface = "bold",
             size = 4.2, colour = "#222222") +
    annotate("text", x = 0.08, y = 0.48, hjust = 0, vjust = 1,
             label = reason, family = anchor_font_family, size = 3.2,
             colour = "#55595C", lineheight = 1.05) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    anchor_theme_void(12)
}

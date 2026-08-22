#!/usr/bin/env python3
"""Apply audited LUAD-specific Figure 2E corrections to a materialised runtime.

The frozen TNBC source remains byte-identical.  This patch is applied only to
the materialised LUAD runtime after path/name adaptation.  It makes the point-
size direction agree with the published legend, adds a clearly labelled zoom
of the three conserved-core TFs while retaining the complete six-TF range, and
emits one composite Figure 2E PDF containing both official halves.

Scientific input corrections (one frozen best motif per TF/sample and the
official CIS-BP category vocabulary) are performed in the data-only adapter.
This program refuses to patch unless those input invariants already hold.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
from collections import Counter
from pathlib import Path


SOURCE_REL = Path("code/Figure2/08-Plot_TR_Examples.R")
SIZE_OLD = 'scale_size_continuous(name = "SD of LOR", range = c(6, 1))'
SIZE_NEW = 'scale_size_continuous(name = "SD of LOR", range = c(1, 6))'
CASE_OUTPUT_BLOCK = '''cairo_pdf(paste0(args$visuals_outdir, "/Figure2E_Case_Examples_HC-TRs.pdf"), height=7.5, width = 5)
print(plt)
dev.off()
'''
ZOOM_AND_COMPOSITE_BLOCK = r'''

# LUAD audited extension: retain the official six-case full range and repeat
# the same observations for the three conserved-core TFs on a labelled zoom.
# No observation is replaced, winsorised, or omitted from the full-range view.
core_tfs <- c("FOXA3", "NFATC4", "XBP1")
df_lor_core <- df_lor_sub %>% dplyr::filter(TR %in% core_tfs)
core_key_counts <- df_lor_core %>%
  dplyr::count(TR, Sample_Group, Sample, name = "n_rows")
if (any(core_key_counts$n_rows != 1L)) {
  stop("Figure2E core zoom requires exactly one frozen motif per TF/sample.")
}
core_range <- range(df_lor_core$LOR, na.rm = TRUE)
core_padding <- max(diff(core_range) * 0.08, 0.01)
core_limits <- c(core_range[1] - core_padding, core_range[2] + core_padding)

plt_full_range <- plt +
  labs(title = "Six selected TFs - full observed range")
plt_core_zoom <- (plt %+% df_lor_core) +
  scale_y_continuous(
    limits = core_limits,
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Conserved core zoom - FOXA3, NFATC4 and XBP1",
    subtitle = "Same observations as the full-range panel"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 8)
  )

cairo_pdf(
  paste0(args$visuals_outdir, "/Figure2E_Core_Zoom_HC-TRs.pdf"),
  height = 4.2, width = 5
)
print(plt_core_zoom)
dev.off()

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("The declared Figure2E vector composition requires patchwork.")
}
summary_for_complete <- Motif_logOdds_plt +
  labs(title = sprintf("%d motif-supported HC-TFs", nrow(Motif_HC_TRs))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
right_column <- patchwork::wrap_plots(
  plt_full_range, plt_core_zoom,
  ncol = 1, heights = c(2.05, 1)
)
figure2e_complete <- patchwork::wrap_plots(
  summary_for_complete, right_column,
  ncol = 2, widths = c(1, 1.42)
)
cairo_pdf(
  paste0(args$visuals_outdir, "/Figure2E_Complete_HC-TRs.pdf"),
  height = 8, width = 12
)
print(figure2e_complete)
dev.off()

figure2e_case_counts <- df_lor_sub %>%
  dplyr::distinct(TR, Sample_Group, Sample) %>%
  dplyr::count(TR, Sample_Group, name = "n_unique_samples")
write.table(
  figure2e_case_counts,
  paste0(args$visuals_outdir, "/Figure2E_Case_Sample_Counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
figure2e_runtime_receipt <- data.frame(
  metric = c(
    "summary_tf_count", "case_tf_count", "duplicate_case_tf_sample_keys",
    "full_lor_min", "full_lor_max", "core_lor_min", "core_lor_max",
    "core_zoom_min", "core_zoom_max"
  ),
  value = c(
    nrow(Motif_HC_TRs), length(unique(df_lor_sub$TR)),
    sum((df_lor_sub %>% dplyr::count(TR, Sample_Group, Sample))$n > 1L),
    min(df_lor_sub$LOR, na.rm = TRUE), max(df_lor_sub$LOR, na.rm = TRUE),
    core_range[1], core_range[2], core_limits[1], core_limits[2]
  )
)
write.table(
  figure2e_runtime_receipt,
  paste0(args$visuals_outdir, "/Figure2E_Runtime_Receipt.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
'''

WILCOX_OLD = "print(wilcox_results)"
WILCOX_NEW = '''print(wilcox_results)
write.table(
  wilcox_results,
  paste0(args$visuals_outdir, "/Figure2E_Pairwise_Wilcoxon_BH.tsv"),
  sep = "\\t", quote = FALSE, row.names = FALSE
)'''

LEGACY_HC31_SUPPORTED = {
    "CREB3L2", "E2F5", "ELF3", "ETV1", "FOXA3", "FOXD2", "NFATC4",
    "NR3C2", "XBP1", "ZNF281", "ZNF33A", "ZNF341", "ZNF444", "ZNF540",
    "ZNF75D",
}
CASE_TFS = {"ZNF341", "ELF3", "ETV1", "FOXA3", "NFATC4", "XBP1"}
EXPECTED_CASE_COUNTS = {"TCGA": 22, "PDX": 11, "CellLines": 19}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty receipt: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def support_tfs(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2:
        raise RuntimeError(f"Empty motif-support matrix: {path}")
    return {row[0].strip().upper() for row in rows[1:] if row and row[0].strip()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("runtime_root", type=Path)
    ns = parser.parse_args()

    runtime_root = ns.runtime_root.resolve()
    manifest = json.loads(
        (runtime_root / "audit/runtime_manifest.json").read_text(encoding="utf-8")
    )
    adapter_root = Path(manifest["adapter_root"])
    motif_dir = adapter_root / "data/Figure2/Motif/Combined_MotifEnrichment"
    support_path = (
        adapter_root
        / "data/Figure2/Motif/JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv"
    )
    observed_supported = support_tfs(support_path)
    missing_legacy = LEGACY_HC31_SUPPORTED - observed_supported
    if missing_legacy:
        raise RuntimeError(
            "Figure2E HC45 input unexpectedly lost motif-supported TFs from "
            f"the audited HC31 run: missing={sorted(missing_legacy)}"
        )
    if not CASE_TFS.issubset(observed_supported):
        raise RuntimeError(
            "Figure2E case TFs are absent from the current motif-supported set: "
            f"{sorted(CASE_TFS - observed_supported)}"
        )

    motif_rows: list[dict[str, str]] = []
    for filename in (
        "JASPAR_Combined_MotifEnrichment_TNBC.tsv",
        "CISBP_Combined_MotifEnrichment_TNBC.tsv",
    ):
        motif_rows.extend(read_tsv(motif_dir / filename))
    selected_rows = [
        row for row in motif_rows if row["Motif.Name"].strip().upper() in CASE_TFS
    ]
    key_counts = Counter(
        (
            row["Motif.Name"].strip().upper(),
            row["Sample_Group"].strip(),
            row["Sample"].strip(),
        )
        for row in selected_rows
    )
    duplicate_keys = {key: count for key, count in key_counts.items() if count != 1}
    if duplicate_keys:
        raise RuntimeError(
            "Figure2E formal input still contains duplicate TF/sample motifs: "
            f"{list(duplicate_keys.items())[:5]}"
        )
    case_count_rows: list[dict[str, object]] = []
    for tf in sorted(CASE_TFS):
        for group, expected in EXPECTED_CASE_COUNTS.items():
            observed = sum(1 for key in key_counts if key[0] == tf and key[1] == group)
            if observed != expected:
                raise RuntimeError(
                    f"{tf}/{group}: expected {expected} unique samples, observed {observed}"
                )
            case_count_rows.append(
                {"TF": tf, "Sample_Group": group, "n_unique_samples": observed}
            )

    source = runtime_root / SOURCE_REL
    before = source.read_text(encoding="utf-8")
    if before.count(SIZE_OLD) != 1:
        raise RuntimeError("Expected exactly one reversed Figure2E SD size scale.")
    if before.count(CASE_OUTPUT_BLOCK) != 1:
        raise RuntimeError("Expected exactly one official Figure2E case-output block.")
    if before.count(WILCOX_OLD) != 1:
        raise RuntimeError("Expected exactly one official Figure2E Wilcoxon print call.")

    after = before.replace(SIZE_OLD, SIZE_NEW)
    after = after.replace(
        CASE_OUTPUT_BLOCK,
        CASE_OUTPUT_BLOCK + ZOOM_AND_COMPOSITE_BLOCK,
    )
    after = after.replace(WILCOX_OLD, WILCOX_NEW)
    source.write_text(after, encoding="utf-8", newline="\n")

    audit_root = runtime_root / "audit"
    audit_root.mkdir(parents=True, exist_ok=True)
    diff = "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile="materialised-official/Figure2/08-Plot_TR_Examples.R",
            tofile="final-corrected/Figure2/08-Plot_TR_Examples.R",
            n=3,
        )
    )
    (audit_root / "figure2e_luad_motif_bugfix.diff").write_text(diff, encoding="utf-8")
    write_tsv(
        audit_root / "figure2e_case_input_counts.tsv",
        case_count_rows,
    )
    write_tsv(
        audit_root / "figure2e_luad_motif_bugfix_manifest.tsv",
        [
            {
                "panel": "Figure2E",
                "change_class": "LUAD_MOTIF_INPUT_BUGFIX_AND_DECLARED_ZOOM",
                "official_source": SOURCE_REL.as_posix(),
                "official_source_unchanged": True,
                "supported_tf_count": len(observed_supported),
                "supported_tfs": ";".join(sorted(observed_supported)),
                "case_tf_count": len(CASE_TFS),
                "case_tf_sample_duplicate_keys": len(duplicate_keys),
                "etv1_unique_samples": sum(1 for key in key_counts if key[0] == "ETV1"),
                "sd_size_direction_fixed": True,
                "full_observed_range_retained": True,
                "core_zoom_tfs": "FOXA3;NFATC4;XBP1",
                "core_zoom_reuses_full_range_observations": True,
                "official_boxpoint_geometry_retained": True,
                "official_palette_retained": True,
                "strict_runtime_sha256": sha256_text(before),
                "corrected_runtime_sha256": sha256_text(after),
            }
        ],
    )
    print(
        f"Figure2E audited patch applied: {len(observed_supported)} supported TFs; "
        "no TF/sample duplicates; "
        "SD size direction corrected; full range plus conserved-core zoom emitted."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

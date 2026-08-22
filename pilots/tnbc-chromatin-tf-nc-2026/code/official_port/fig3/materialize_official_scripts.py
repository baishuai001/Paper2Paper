#!/usr/bin/env python3
"""Validate and minimally materialize the official CodeOcean Figure 3 scripts.

The plots are still built by the official R scripts. This program performs only
explicitly enumerated portability/adaptation edits and writes a unified diff for
every materialized script. Any upstream source checksum or expected text change
causes a hard failure.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
from pathlib import Path


SOURCE_SPECS = {
    "03-Pearson_Correlation_TRActivity.R": {
        "sha256": "61f46a08820e68a0b3ae5fcc6aa73f540039505b28424a08e29664ce943fd2a8",
        "panels": "Figure3C",
        "lines": "44-156",
    },
    "04-Summarize_Regulons_HC_TR.R": {
        "sha256": "c846a6375fca87914cbef46837c4b31d273937fdc2d40d9604f5894046d73e44",
        "panels": "Figure3A",
        "lines": "50-150;157-222",
    },
    "05-Generate_Regulon_Network.R": {
        "sha256": "2f241e60d1c21dc381942abb7ee92170b8a2343c85d38951cc680fb1d646e800",
        "panels": "Figure3B;Figure3D;Figure3E",
        "lines": "44-104;110-227;239-375",
    },
    "06-Activity_Score_Heterogeneity_Across_Patients.R": {
        "sha256": "ea0ad4d093d1af0e65d2e81e25181d9f02e4b4614229acb37b4da44a8bf0df63",
        "panels": "Figure3F;Figure3G;SupplementaryFigure5A",
        "lines": "63-189;191-269",
    },
    "07B-Plot_MKi67_HC-TR.R": {
        "sha256": "7df80b0ae7e17ecee9d039f5ae039ca48ed5c4fae7e072cf5e8ca98d438d853d",
        "panels": "SupplementaryFigure4D_volcanoes",
        "lines": "31-75",
    },
    "07C-Plot_Heatmap_MKI67_HC-TR.R": {
        "sha256": "5b44ad2058e2dd556e851966a1d24d4942eb9a838df8269fe6bbe18853f11214",
        "panels": "SupplementaryFigure4D_heatmap",
        "lines": "16-123",
    },
    "08D-Plot_ESTIMATE_Volcano.R": {
        "sha256": "758248868b1375bbab20970e414294c0f085efae0104619a70019bee71a160ff",
        "panels": "SupplementaryFigure5B_volcanoes",
        "lines": "45-196",
    },
    "08E-Plot_Heatmap_ESTIMATE.R": {
        "sha256": "43569ca10503416e48abace611445148c5e706b6d926cf1c7da58a025d0e8c38",
        "panels": "SupplementaryFigure5B_heatmap",
        "lines": "53-249",
    },
}

AESTHETICS_SHA256 = "52f1fc590c6cad927ee0d1b030ddc215e99c506d15b465f3c805a9cc0ca11c1d"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {count}")
    return text.replace(old, new, 1)


def replace_exact_count(text: str, old: str, new: str, expected: int, label: str) -> str:
    count = text.count(old)
    if count != expected:
        raise RuntimeError(f"{label}: expected {expected} matches, found {count}")
    return text.replace(old, new)


def patch_script(name: str, source: str, figure3b_min_shared_tfs: int) -> str:
    patched = source
    source_old = 'source("/code/Static_Scripts/plotting_aesthetics.R")'
    source_new = (
        'source(file.path(Sys.getenv("PAPER2PAPER_OFFICIAL_PORT_WORK"), '
        '"Static_Scripts", "plotting_aesthetics.R"))'
    )

    if name in {
        "04-Summarize_Regulons_HC_TR.R",
        "05-Generate_Regulon_Network.R",
        "06-Activity_Score_Heterogeneity_Across_Patients.R",
    }:
        patched = replace_once(patched, source_old, source_new, f"{name}: source path")

    if name == "04-Summarize_Regulons_HC_TR.R":
        patched = replace_once(
            patched,
            'labs(x = "Number of interactions in TNBC regulatory network", \n'
            '        y = "TNBC-specific high-confidence TRs", fill = "")',
            'labs(x = "Number of interactions in LUAD regulatory network", \n'
            '        y = "LUAD-specific high-confidence TFs", fill = "")',
            f"{name}: LUAD axis labels",
        )
        patched = replace_once(
            patched,
            '#write.table(VIPER_interactions, paste0(args$data_outdir, "/HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv"),\n'
            "#    col.names=TRUE, row.names=FALSE, sep='\\t', quote=FALSE)",
            'write.table(VIPER_interactions, paste0(args$data_outdir, "/HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv"),\n'
            "    col.names=TRUE, row.names=FALSE, sep='\\t', quote=FALSE)",
            f"{name}: expose official downstream table",
        )

    if name == "05-Generate_Regulon_Network.R":
        patched = replace_once(
            patched,
            '#write.table(edge_df, paste0(args$data_outdir, "/HC_TR_60_Plotted_Network.tsv"),\n'
            '#  col.names=TRUE, row.names=FALSE, sep=\'\\t\', quote=FALSE)',
            'write.table(edge_df, paste0(args$data_outdir, "/HC_TR_60_Plotted_Network.tsv"),\n'
            '  col.names=TRUE, row.names=FALSE, sep=\'\\t\', quote=FALSE)',
            f"{name}: expose official Figure3D edge table",
        )
        patched = replace_once(
            patched,
            "scatter_df$percent_HC_TF_overlap <- scatter_df$n_partners/94",
            "scatter_df$percent_HC_TF_overlap <- scatter_df$n_partners/nrow(HC_TRs)",
            f"{name}: dynamic HC-TF denominator",
        )
        if figure3b_min_shared_tfs != 3:
            patched = replace_once(
                patched,
                "min_shared_tfs <- 3  # Adjust threshold - targets must share at least 3 TFs",
                (
                    f"min_shared_tfs <- {figure3b_min_shared_tfs}  # Paper2Paper frozen LUAD parameter: "
                    f"targets share at least {figure3b_min_shared_tfs} HC-TFs"
                ),
                f"{name}: explicit Figure 3B overlap parameter",
            )

    if name == "06-Activity_Score_Heterogeneity_Across_Patients.R":
        patched = replace_once(
            patched,
            'source(file.path(Sys.getenv("PAPER2PAPER_OFFICIAL_PORT_WORK"), '
            '"Static_Scripts", "plotting_aesthetics.R"))',
            'source(file.path(Sys.getenv("PAPER2PAPER_OFFICIAL_PORT_WORK"), '
            '"Static_Scripts", "plotting_aesthetics.R"))\n'
            'sample_category_colours <- setNames(unname(sample_category_colours), '
            'c("LUAD", "Non-LUAD"))',
            f"{name}: semantic palette labels",
        )
        patched = replace_once(
            patched,
            'TR_connections_top5 <- c("FOXC1", "MYBL2", "PML", "RCOR2")',
            'TR_connections_top5 <- readLines(Sys.getenv("PAPER2PAPER_REPRESENTATIVE_TFS_FILE"), warn = FALSE)',
            f"{name}: LUAD representatives",
        )
        patched = replace_exact_count(
            patched,
            '"TNBC", "Non-TNBC"',
            '"LUAD", "Non-LUAD"',
            3,
            f"{name}: subtype labels",
        )
        patched = replace_once(
            patched,
            'subtype == "TNBC"',
            'subtype == "LUAD"',
            f"{name}: skewness target group",
        )
        patched = patched.replace("Skewness (TNBC samples)", "Skewness (LUAD samples)")

    if name == "07C-Plot_Heatmap_MKI67_HC-TR.R":
        patched = replace_once(
            patched,
            'args <- list(visuals_out = "/results/visuals/Figure3/")',
            'args <- list(visuals_out = Sys.getenv("PAPER2PAPER_OFFICIAL_VISUALS_OUT"))',
            f"{name}: output path",
        )
        patched = replace_once(
            patched,
            'results_files <- list(\n'
            '    TCGA     = "/data/Figure3/MKi67_Correlation/TCGA_MKi67_Correlation_HC-TR_results.tsv",\n'
            '    METABRIC = "/data/Figure3/MKi67_Correlation/METABRIC_MKi67_Correlation_HC-TR_results.tsv",\n'
            '    PDX      = "/data/Figure3/MKi67_Correlation/PDX_MKi67_Correlation_HC-TR_results.tsv",\n'
            '    CellLines = "/data/Figure3/MKi67_Correlation/CellLines_MKi67_Correlation_HC-TR_results.tsv"\n'
            ')',
            'mki67_dir <- file.path(Sys.getenv("PAPER2PAPER_OFFICIAL_ADAPTER_INPUTS"), "MKI67")\n'
            'results_files <- list(\n'
            '    TCGA      = file.path(mki67_dir, "TCGA_MKi67_Correlation_HC-TR_results.tsv"),\n'
            '    GSE41271  = file.path(mki67_dir, "GSE41271_MKi67_Correlation_HC-TR_results.tsv"),\n'
            '    PDX       = file.path(mki67_dir, "PDX_MKi67_Correlation_HC-TR_results.tsv"),\n'
            '    CellLines = file.path(mki67_dir, "CellLines_MKi67_Correlation_HC-TR_results.tsv")\n'
            ')',
            f"{name}: four-system LUAD input paths",
        )
        patched = replace_once(
            patched,
            'cohort_order <- c("TCGA", "METABRIC", "PDX", "CellLines")',
            'cohort_order <- c("TCGA", "GSE41271", "PDX", "CellLines")',
            f"{name}: independent patient cohort label",
        )

    return patched


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--official-code-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--figure3b-min-shared-tfs", type=int, default=3)
    args = parser.parse_args()
    if args.figure3b_min_shared_tfs < 1:
        raise SystemExit("--figure3b-min-shared-tfs must be a positive integer")

    official_root = args.official_code_root.resolve()
    figure3_root = official_root / "Figure3"
    aesthetics_path = official_root / "Static_Scripts" / "plotting_aesthetics.R"
    work_root = args.output_root.resolve() / "work"
    script_out = work_root / "Figure3"
    static_out = work_root / "Static_Scripts"
    audit_root = args.output_root.resolve() / "audit"
    diff_root = audit_root / "official_diffs"
    for path in (script_out, static_out, diff_root):
        path.mkdir(parents=True, exist_ok=True)

    aesthetics_payload = aesthetics_path.read_bytes()
    aesthetics_hash = sha256_bytes(aesthetics_payload)
    if aesthetics_hash != AESTHETICS_SHA256:
        raise RuntimeError(
            "Official plotting_aesthetics.R checksum mismatch: "
            f"expected {AESTHETICS_SHA256}, found {aesthetics_hash}"
        )
    (static_out / "plotting_aesthetics.R").write_bytes(aesthetics_payload)

    manifest_rows: list[dict[str, str]] = []
    patch_rows: list[dict[str, str | int]] = []
    for name, spec in SOURCE_SPECS.items():
        source_path = figure3_root / name
        source_payload = source_path.read_bytes()
        source_hash = sha256_bytes(source_payload)
        if source_hash != spec["sha256"]:
            raise RuntimeError(
                f"Official source checksum mismatch for {name}: "
                f"expected {spec['sha256']}, found {source_hash}"
            )
        source_text = source_payload.decode("utf-8")
        patched_text = patch_script(name, source_text, args.figure3b_min_shared_tfs)
        patched_path = script_out / name
        patched_path.write_text(patched_text, encoding="utf-8", newline="")

        diff = "".join(
            difflib.unified_diff(
                source_text.splitlines(keepends=True),
                patched_text.splitlines(keepends=True),
                fromfile=f"official/{name}",
                tofile=f"official_port/{name}",
            )
        )
        (diff_root / f"{name}.diff").write_text(
            diff if diff else "# Byte-identical: no adaptation required.\n",
            encoding="utf-8",
        )
        changed_lines = sum(
            1
            for line in diff.splitlines()
            if (line.startswith("+") or line.startswith("-"))
            and not line.startswith("+++")
            and not line.startswith("---")
        )
        manifest_rows.append(
            {
                "source_file": name,
                "official_absolute_path": str(source_path),
                "official_sha256": source_hash,
                "official_line_count": str(len(source_text.splitlines())),
                "panels": spec["panels"],
                "constructor_dataflow_lines": spec["lines"],
                "materialized_sha256": sha256_bytes(patched_text.encode("utf-8")),
                "diff_file": str(diff_root / f"{name}.diff"),
            }
        )
        patch_rows.append(
            {
                "source_file": name,
                "changed_plus_minus_lines": changed_lines,
                "plot_constructor_changed": "NO",
                "analysis_parameter_changed": (
                    "YES" if name == "05-Generate_Regulon_Network.R" and args.figure3b_min_shared_tfs != 3 else "NO"
                ),
                "allowed_change_classes": (
                    "none"
                    if not diff
                    else (
                        "path;LUAD semantic labels;expose commented official tables;LUAD representative IDs;"
                        "explicit Figure3B shared-TF threshold"
                    )
                ),
            }
        )

    with (audit_root / "official_source_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    with (audit_root / "official_patch_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(patch_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(patch_rows)

    print(f"Materialized {len(SOURCE_SPECS)} checksum-locked official Figure 3 scripts in {script_out}")


if __name__ == "__main__":
    main()

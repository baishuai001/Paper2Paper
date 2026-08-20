#!/usr/bin/env python3
"""Independent structural verification for the strict official Figure 3 port."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def components(nodes: set[str], edges: set[tuple[str, str]]) -> int:
    parent = {node: node for node in nodes}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: str, right: str) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left, right in edges:
        union(left, right)
    return len({find(node) for node in nodes})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()

    root = args.output_root.resolve()
    data = root / "data"
    inputs = data / "inputs"
    official_tables = data / "official_tables"
    visuals = root / "results" / "visuals" / "Figure3"
    audit = root / "audit"
    audit.mkdir(parents=True, exist_ok=True)

    hc_tfs = {
        line.strip()
        for line in (inputs / "High_RNA_NES.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if len(hc_tfs) != 31:
        raise RuntimeError(f"Expected 31 HC-TFs; found {len(hc_tfs)}")

    regulon_rows = read_tsv(official_tables / "HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv")
    required_regulon = {"TR_Source", "Target_Gene", "VIPER_Weight", "Label"}
    if not regulon_rows or not required_regulon.issubset(regulon_rows[0]):
        raise RuntimeError("Official Figure 3A regulon table is absent or malformed")
    # Official script 04 merges the motif annotation with all=TRUE after the
    # regulon rows have been built. Motif-only merge rows therefore carry R's
    # literal NA in TR_Source. They are official annotation placeholders, not
    # a 32nd regulator, and must be excluded from the source-TF identity check.
    missing_tokens = {"", "NA", "NaN", "nan", "None", "NULL"}
    regulon_source_tfs = {
        row["TR_Source"].strip()
        for row in regulon_rows
        if row.get("TR_Source", "").strip() not in missing_tokens
    }
    if regulon_source_tfs != hc_tfs:
        raise RuntimeError("Official Figure 3A regulon table does not contain exactly 31 HC-TFs")

    activated_shared: dict[str, set[str]] = defaultdict(set)
    for row in regulon_rows:
        if row["Label"] == "Activated Shared Targets":
            activated_shared[row["TR_Source"]].add(row["Target_Gene"])

    # Figure 3D: every HC-TF pair sharing >=1 activated shared target.
    expected_tf_edges: dict[tuple[str, str], int] = {}
    for left, right in itertools.combinations(sorted(hc_tfs), 2):
        weight = len(activated_shared[left] & activated_shared[right])
        if weight >= 1:
            expected_tf_edges[(left, right)] = weight

    official_edge_rows = read_tsv(official_tables / "HC_TR_60_Plotted_Network.tsv")
    observed_tf_edges = {
        tuple(sorted((row["from"], row["to"]))): int(float(row["weight"]))
        for row in official_edge_rows
    }
    if observed_tf_edges != expected_tf_edges:
        raise RuntimeError("Official Figure 3D edge table is not the complete >=1 overlap graph")
    figure3d_nodes = {node for edge in observed_tf_edges for node in edge}
    figure3d_components = components(figure3d_nodes, set(observed_tf_edges))
    if (len(figure3d_nodes), len(observed_tf_edges), figure3d_components) != (31, 134, 1):
        raise RuntimeError(
            "Figure 3D expected 31 nodes / 134 edges / 1 component; found "
            f"{len(figure3d_nodes)} / {len(observed_tf_edges)} / {figure3d_components}"
        )

    # Figure 3B: exact official transposed incidence, >=3 shared TFs, and no
    # explicit vertices argument (only qualifying endpoints are nodes).
    target_to_tfs: dict[str, set[str]] = defaultdict(set)
    for tf, targets in activated_shared.items():
        for target in targets:
            target_to_tfs[target].add(tf)
    target_edges: dict[tuple[str, str], int] = {}
    for left, right in itertools.combinations(sorted(target_to_tfs), 2):
        weight = len(target_to_tfs[left] & target_to_tfs[right])
        if weight >= 3:
            target_edges[(left, right)] = weight
    figure3b_nodes = {node for edge in target_edges for node in edge}
    if (len(figure3b_nodes), len(target_edges)) != (6, 3):
        raise RuntimeError(
            "Figure 3B expected the frozen official topology 6 nodes / 3 edges; found "
            f"{len(figure3b_nodes)} / {len(target_edges)}"
        )

    representatives = read_tsv(audit / "figure3G_representative_tfs.tsv")
    representative_tfs = [row["TF"] for row in representatives]
    if len(representative_tfs) != 4 or len(set(representative_tfs)) != 4:
        raise RuntimeError("Figure 3G does not have four unique LUAD representative TFs")

    expected_pdfs = [
        "Figure3A_RegulonInteractions_TCGA_VIPER.pdf",
        "Figure3B_Target_Network_Shared_TFs.pdf",
        "Figure3C_TCGA_PearsonCorrelation_TR_Activity.pdf",
        "Figure3C_Top_Bottom_TF_PearsonCorrelation.pdf",
        "Figure3D_TR_Network_Shared_ActivatingTargetGenes.pdf",
        "Figure3E_Scatter_RegulonSize_vs_Partners.pdf",
        "Figure3F_Skewness_TCGA.pdf",
        "SupplementaryFigure5A_Skewness_PDX_CellLines.pdf",
        "SupplementaryFigure4D_TCGA_MKi67_Correlation_HC-TR.pdf",
        "SupplementaryFigure4D_GSE41271_MKi67_Correlation_HC-TR.pdf",
        "SupplementaryFigure4D_PDX_MKi67_Correlation_HC-TR.pdf",
        "SupplementaryFigure4D_CellLines_MKi67_Correlation_HC-TR.pdf",
        "SupplementaryFigure4D_All_MKi67_Correlation_HC-TR_Heatmap.pdf",
        "SupplementaryFigure5B_All_ESTIMATE_Volcano_HC-TR.pdf",
        "SupplementaryFigure5B_ESTIMATE_Combined_Heatmap_HC-TR.pdf",
    ]
    expected_pdfs.extend(
        f"Figure3G_{tf}_TR_Activity_Distribution.pdf" for tf in representative_tfs
    )

    output_rows: list[dict[str, str | int]] = []
    for name in expected_pdfs:
        path = visuals / name
        if not path.is_file():
            raise RuntimeError(f"Missing official PDF: {path}")
        size = path.stat().st_size
        if size < 1000 or path.read_bytes()[:5] != b"%PDF-":
            raise RuntimeError(f"Invalid or suspiciously small PDF: {path}")
        output_rows.append(
            {"artifact": name, "bytes": size, "sha256": sha256(path), "status": "PASS"}
        )

    source_manifest = read_tsv(audit / "official_source_manifest.tsv")
    patch_manifest = read_tsv(audit / "official_patch_manifest.tsv")
    if len(source_manifest) != 8 or len(patch_manifest) != 8:
        raise RuntimeError("Official source/patch manifest is incomplete")
    if any(row["plot_constructor_changed"] != "NO" for row in patch_manifest):
        raise RuntimeError("At least one official plot constructor was marked as changed")

    with (audit / "official_output_sha256.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    receipt = {
        "status": "PASS",
        "implementation": "checksum-locked official CodeOcean scripts with data-only LUAD adapters",
        "hc_tf_count": len(hc_tfs),
        "regulon_interaction_rows": len(regulon_rows),
        "figure3b": {
            "nodes": len(figure3b_nodes),
            "edges": len(target_edges),
            "threshold": "at least 3 shared activating HC-TFs",
            "vertices_argument": False,
        },
        "figure3d": {
            "nodes": len(figure3d_nodes),
            "edges": len(observed_tf_edges),
            "components": figure3d_components,
            "isolates": len(hc_tfs - figure3d_nodes),
            "threshold": "at least 1 shared activated target",
            "display_edge_filter": None,
        },
        "figure3g_representatives": representative_tfs,
        "pdf_count": len(output_rows),
        "official_source_files": len(source_manifest),
        "plot_constructors_changed": False,
    }
    (audit / "official_port_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

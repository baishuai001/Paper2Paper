#!/usr/bin/env python3
"""Independently verify the official-code Figure 3B HC45/overlap>=2 run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def component_sizes(nodes: set[str], edges: dict[tuple[str, str], int]) -> list[int]:
    neighbours: dict[str, set[str]] = {node: set() for node in nodes}
    for left, right in edges:
        neighbours[left].add(right)
        neighbours[right].add(left)
    unseen = set(nodes)
    sizes: list[int] = []
    while unseen:
        seed = unseen.pop()
        stack = [seed]
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            new = neighbours[node] & unseen
            unseen -= new
            stack.extend(new)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    inputs = root / "data/inputs"
    tables = root / "data/official_tables"
    visuals = root / "results/visuals/Figure3"
    audit = root / "audit"
    audit.mkdir(parents=True, exist_ok=True)

    hc = {
        line.strip()
        for line in (inputs / "High_RNA_NES.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if len(hc) != 45:
        raise RuntimeError(f"Expected 45 HC-TFs; found {len(hc)}")
    rows = read_tsv(tables / "HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv")
    activated_shared: dict[str, set[str]] = defaultdict(set)
    missing = {"", "NA", "NaN", "NULL"}
    source_tfs = {
        row["TR_Source"] for row in rows if row.get("TR_Source", "") not in missing
    }
    if source_tfs != hc:
        raise RuntimeError("Official intermediate does not contain exactly the HC45 set")
    for row in rows:
        if row.get("Label") == "Activated Shared Targets":
            activated_shared[row["TR_Source"]].add(row["Target_Gene"])

    target_to_tfs: dict[str, set[str]] = defaultdict(set)
    for tf, targets in activated_shared.items():
        for target in targets:
            target_to_tfs[target].add(tf)
    edges: dict[tuple[str, str], int] = {}
    for left, right in itertools.combinations(sorted(target_to_tfs), 2):
        weight = len(target_to_tfs[left] & target_to_tfs[right])
        if weight >= 2:
            edges[(left, right)] = weight
    nodes = {node for pair in edges for node in pair}
    if not edges:
        raise RuntimeError("The user-frozen >=2 Figure 3B graph is empty")
    sizes = component_sizes(nodes, edges)
    with (audit / "Figure3B_HC45_shared_ge2_edges.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["from", "to", "shared_activating_HC_TFs"], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {"from": pair[0], "to": pair[1], "shared_activating_HC_TFs": weight}
            for pair, weight in sorted(edges.items())
        )

    pdf = visuals / "Figure3B_Target_Network_Shared_TFs.pdf"
    if not pdf.is_file() or pdf.stat().st_size < 1000 or pdf.read_bytes()[:5] != b"%PDF-":
        raise RuntimeError(f"Official Figure 3B PDF is missing or invalid: {pdf}")
    receipt = {
        "status": "PASS",
        "implementation": "checksum-locked TNBC Figure3/04 and Figure3/05; LUAD data adapter; user-frozen overlap threshold",
        "HC_TFs": len(hc),
        "activated_shared_targets_before_pair_filter": len(target_to_tfs),
        "minimum_shared_activating_HC_TFs": 2,
        "nodes": len(nodes),
        "edges": len(edges),
        "components": len(sizes),
        "largest_component": sizes[0],
        "component_sizes": sizes,
        "maximum_shared_HC_TFs": max(edges.values()),
        "inference_note": "The >=2 network is descriptive. The prior degree-preserving permutation test found zero FDR-significant target pairs.",
        "pdf": str(pdf),
        "pdf_sha256": sha256(pdf),
    }
    (audit / "figure3b_official_hc45_ge2_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

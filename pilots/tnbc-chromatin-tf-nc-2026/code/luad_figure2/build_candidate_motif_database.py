#!/usr/bin/env python3
"""Build frozen combined and source-split candidate motif databases."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from parse_harmonized_motifs import (
    LINEAGE_CONTROLS,
    build_token_map,
    read_tsv,
    sha256,
    write_tsv,
)


def read_homer_records(path: Path) -> list[tuple[str, list[str]]]:
    records: list[tuple[str, list[str]]] = []
    token = ""
    lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if token:
                    records.append((token, lines))
                fields = line[1:].strip().split()
                if len(fields) < 2:
                    raise RuntimeError(f"Malformed HOMER header: {line.rstrip()}")
                token = fields[1]
                lines = [line]
            elif token:
                lines.append(line)
    if token:
        records.append((token, lines))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    full_database = root / "reference/motifs/full_reference_both_databases.homer"
    inventory = read_tsv(root / "results/motif_gate/hc_tf_motif_inventory.tsv")
    hc_tfs = {row["TF"] for row in inventory if row["motif_testable"] == "TRUE"}
    token_map, control_catalog = build_token_map(root, hc_tfs, full_database)
    wanted = set(token_map)
    records = read_homer_records(full_database)
    observed = {token for token, _ in records}
    missing = sorted(wanted - observed)
    if missing:
        raise RuntimeError(f"Candidate motif tokens absent from full database: {missing}")

    output = root / "reference/motifs/candidate_luad_hc_and_lineage_controls.homer"
    source_outputs = {
        source: root / f"reference/motifs/candidate_luad_hc_and_lineage_controls.{source}.homer"
        for source in ("JASPAR2024", "CIS-BP2.00")
    }
    selected_records = [(token, lines) for token, lines in records if token in wanted]
    with output.open("w", encoding="utf-8", newline="") as handle:
        for _, lines in selected_records:
            handle.writelines(lines)
    for source, source_output in source_outputs.items():
        with source_output.open("w", encoding="utf-8", newline="") as handle:
            for token, lines in selected_records:
                if token.startswith(f"{source}|"):
                    handle.writelines(lines)

    mapping_rows = []
    for token in sorted(wanted):
        for tf in sorted(token_map[token]):
            mapping_rows.append(
                {
                    "motif_token": token,
                    "TF": tf,
                    "is_luad_hc_tf": str(tf in hc_tfs).upper(),
                    "is_lineage_positive_control": str(tf in LINEAGE_CONTROLS).upper(),
                }
            )
    mapping_path = root / "audit/motif_equivalence/candidate_motif_token_map.tsv"
    write_tsv(
        mapping_path,
        mapping_rows,
        ["motif_token", "TF", "is_luad_hc_tf", "is_lineage_positive_control"],
    )
    source_counts = Counter(token.split("|", 1)[0] for token in wanted)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "candidate motif family frozen once and emitted both combined and source-split; source-split files reproduce the TNBC author's independent JASPAR/CIS-BP correction strategy",
        "full_database_sha256": sha256(full_database),
        "candidate_database": str(output),
        "candidate_database_sha256": sha256(output),
        "source_split_databases": {
            source: {
                "path": str(source_output),
                "sha256": sha256(source_output),
                "motif_model_count": source_counts[source],
            }
            for source, source_output in source_outputs.items()
        },
        "candidate_motif_model_count": len(wanted),
        "candidate_motif_models_by_source": dict(sorted(source_counts.items())),
        "motif_testable_hc_tf_count": len(hc_tfs),
        "lineage_positive_controls": list(LINEAGE_CONTROLS),
        "lineage_controls_with_catalogued_motifs": sorted(tf for tf in LINEAGE_CONTROLS if control_catalog.get(tf)),
        "token_mapping_sha256": sha256(mapping_path),
    }
    receipt_path = root / "audit/motif_equivalence/candidate_motif_database_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Inventory JASPAR/CIS-BP motifs and freeze the primary HC-TF motif set."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PRIMARY_DEFINITION = "PRIMARY_anchor_promoter_tcga_cpm1_both"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_meme(path: Path) -> tuple[list[str], list[tuple[str, str, list[str]]]]:
    header: list[str] = []
    motifs: list[tuple[str, str, list[str]]] = []
    current: list[str] = []
    seen_motif = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MOTIF "):
                if current:
                    first = current[0].strip().split(maxsplit=2)
                    motifs.append((first[1], first[2] if len(first) > 2 else "", current))
                current = [line]
                seen_motif = True
            elif seen_motif:
                current.append(line)
            else:
                header.append(line)
    if current:
        first = current[0].strip().split(maxsplit=2)
        motifs.append((first[1], first[2] if len(first) > 2 else "", current))
    if not header or not motifs:
        raise RuntimeError(f"Could not parse MEME motif database: {path}")
    return header, motifs


def matched_tfs(motif_id: str, alternate_name: str, hc_tfs: list[str]) -> list[str]:
    text = f" {motif_id} {alternate_name} ".upper()
    output: list[str] = []
    for tf in hc_tfs:
        # Require symbol boundaries so e.g. E2F3 does not match E2F3A.
        if re.search(rf"(?<![A-Z0-9]){re.escape(tf.upper())}(?![A-Z0-9])", text):
            output.append(tf)
    return output


def meme_matrix(block: list[str]) -> list[list[float]]:
    """Extract an A/C/G/T probability matrix from one MEME motif block."""
    matrix_start = next(
        (index for index, line in enumerate(block) if line.startswith("letter-probability matrix:")),
        None,
    )
    if matrix_start is None:
        raise RuntimeError(f"MEME motif lacks a probability matrix: {block[0].strip()}")
    width_match = re.search(r"\bw=\s*(\d+)", block[matrix_start])
    if not width_match:
        raise RuntimeError(f"MEME motif lacks matrix width: {block[0].strip()}")
    width = int(width_match.group(1))
    matrix: list[list[float]] = []
    for line in block[matrix_start + 1 :]:
        fields = line.strip().split()
        if len(fields) != 4:
            if matrix:
                break
            continue
        try:
            row = [float(value) for value in fields]
        except ValueError:
            if matrix:
                break
            continue
        matrix.append(row)
        if len(matrix) == width:
            break
    if len(matrix) != width:
        raise RuntimeError(
            f"MEME motif matrix width mismatch for {block[0].strip()}: expected {width}, got {len(matrix)}"
        )
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    motif_root = root / "reference/motifs"
    out = root / "results/motif_gate"
    audit = root / "audit/motif_gate"
    out.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)

    gate_rows = read_tsv(root / "results/promoter_gate/tf_promoter_gate_summary.tsv")
    hc_tfs = sorted(
        row["TF"] for row in gate_rows
        if row["analysis_definition"] == PRIMARY_DEFINITION and row["HC_TF_promoter_activity_definition"] == "TRUE"
    )
    if not hc_tfs:
        raise RuntimeError("No HC-TFs passed the promoter/activity gate; motif analysis cannot proceed")

    databases = {
        "JASPAR2024": motif_root / "downloads/JASPAR2024_CORE_vertebrates_non-redundant.meme",
        "CIS-BP2.00": motif_root / "CIS-BP_2.00_Homo_sapiens.meme",
    }
    parsed: dict[str, tuple[list[str], list[tuple[str, str, list[str]]]]] = {
        name: parse_meme(path) for name, path in databases.items()
    }
    matches: dict[str, dict[str, list[tuple[str, str, str, list[str]]]]] = {
        tf: defaultdict(list) for tf in hc_tfs
    }
    full_entries: dict[str, list[tuple[str, str, str, list[str]]]] = defaultdict(list)
    for database, (_header, motifs) in parsed.items():
        for index, (motif_id, alternate_name, block) in enumerate(motifs, start=1):
            safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", motif_id)
            homer_token = f"{database}|{index:06d}|{safe_id}"
            full_entries[database].append((homer_token, motif_id, alternate_name, block))
            for tf in matched_tfs(motif_id, alternate_name, hc_tfs):
                matches[tf][database].append((homer_token, motif_id, alternate_name, block))

    inventory_rows: list[dict[str, object]] = []
    selected: list[tuple[str, str, str, str, str, list[str]]] = []
    for tf in hc_tfs:
        jaspar = matches[tf].get("JASPAR2024", [])
        cisbp = matches[tf].get("CIS-BP2.00", [])
        if jaspar and cisbp:
            category = "both"
        elif jaspar:
            category = "JASPAR_only"
        elif cisbp:
            category = "CIS-BP_only"
        else:
            category = "none"
        primary_source = "JASPAR2024" if jaspar else ("CIS-BP2.00" if cisbp else "none")
        inventory_rows.append(
            {
                "TF": tf,
                "JASPAR_motif_count": len(jaspar),
                "CISBP_motif_count": len(cisbp),
                "total_tested_motif_count": len(jaspar) + len(cisbp),
                "motif_database_category": category,
                "primary_source": primary_source,
                "primary_motif_count": len(jaspar if jaspar else cisbp),
                "motif_testable": str(primary_source != "none").upper(),
            }
        )
        # The anchor ran JASPAR and CIS-BP independently; a TF represented in
        # both databases therefore contributes motifs to both frozen tests.
        for homer_token, motif_id, alternate_name, block in jaspar:
            selected.append(("JASPAR2024", tf, homer_token, motif_id, alternate_name, block))
        for homer_token, motif_id, alternate_name, block in cisbp:
            selected.append(("CIS-BP2.00", tf, homer_token, motif_id, alternate_name, block))

    with (out / "hc_tf_motif_inventory.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(inventory_rows)

    selected_meme = motif_root / "selected_hc_tf_primary_motifs.meme"
    # Both source files are DNA MEME v4; retain a single canonical header.
    header = parsed["JASPAR2024"][0]
    with selected_meme.open("w", encoding="utf-8") as handle:
        handle.writelines(header)
        if header and not header[-1].endswith("\n"):
            handle.write("\n")
        for source, tf, _homer_token, motif_id, _alternate_name, block in selected:
            safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", motif_id)
            handle.write(f"MOTIF {source}|{tf}|{safe_id} {tf}\n")
            handle.writelines(block[1:])
            if block and not block[-1].endswith("\n"):
                handle.write("\n")

    # HOMER requires its own PWM format and a detection threshold.  HOMER's
    # bundled parseJasparMatrix.pl uses threshold 0 for imported JASPAR PWMs;
    # use that documented importer convention for both frozen reference sets.
    full_homer = motif_root / "full_reference_both_databases.homer"
    homer_by_database = {
        "JASPAR2024": motif_root / "full_reference_JASPAR2024.homer",
        "CIS-BP2.00": motif_root / "full_reference_CIS-BP2.00.homer",
    }
    bases = "ACGT"
    homer_handles = {
        database: path.open("w", encoding="utf-8") for database, path in homer_by_database.items()
    }
    try:
        with full_homer.open("w", encoding="utf-8") as combined_handle:
            for source in ("JASPAR2024", "CIS-BP2.00"):
                for homer_token, _motif_id, _alternate_name, block in full_entries[source]:
                    matrix = meme_matrix(block)
                    consensus = "".join(bases[max(range(4), key=lambda index: row[index])] for row in matrix)
                    rendered = [f">{consensus}\t{homer_token}\t0\n"]
                    rendered.extend("\t".join(f"{value:.9g}" for value in row) + "\n" for row in matrix)
                    combined_handle.writelines(rendered)
                    homer_handles[source].writelines(rendered)
    finally:
        for handle in homer_handles.values():
            handle.close()

    mapping_rows = [
        {
            "database": source,
            "TF": tf,
            "homer_motif_token": homer_token,
            "original_motif_id": motif_id,
            "original_alternate_name": alternate_name,
            "selected_meme_id": f"{source}|{tf}|{re.sub(r'[^A-Za-z0-9_.-]+', '_', motif_id)}",
        }
        for source, tf, homer_token, motif_id, alternate_name, _block in selected
    ]
    mapping_columns = ["database", "TF", "homer_motif_token", "original_motif_id", "original_alternate_name", "selected_meme_id"]
    with (out / "selected_motif_mapping.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mapping_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(mapping_rows)

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "HC_TFs": len(hc_tfs),
        "motif_testable_HC_TFs": sum(row["motif_testable"] == "TRUE" for row in inventory_rows),
        "selected_motif_records": len(selected),
        "selected_motif_records_by_database": {
            database: sum(source == database for source, _tf, _token, _motif_id, _alternate, _block in selected)
            for database in ("JASPAR2024", "CIS-BP2.00")
        },
        "full_reference_motif_records_by_database": {
            database: len(full_entries[database]) for database in ("JASPAR2024", "CIS-BP2.00")
        },
        "selection_rule": "HC-TF support is mapped after HOMER tests each complete reference database independently",
        "homer_import_threshold": 0,
        "homer_import_threshold_basis": "HOMER bundled parseJasparMatrix.pl convention for external JASPAR matrices",
        "categories": {
            category: sum(row["motif_database_category"] == category for row in inventory_rows)
            for category in ("JASPAR_only", "CIS-BP_only", "both", "none")
        },
        "selected_meme": str(selected_meme),
        "full_homer_combined_audit_copy": str(full_homer),
        "full_homer_by_database": {database: str(path) for database, path in homer_by_database.items()},
    }
    (audit / "hc_tf_motif_catalog_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

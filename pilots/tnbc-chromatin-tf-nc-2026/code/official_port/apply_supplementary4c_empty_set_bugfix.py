#!/usr/bin/env python3
"""Make the official Supplementary Figure 4C set table empty-set safe.

The frozen upstream source remains untouched. The capsule constructs three
two-column data.frames by pairing a possibly empty TF vector with one scalar
database label. R rejects the 0-versus-1 row count when a set is empty. This
patch changes only those three label expressions to length-matched ``rep``
calls; it does not add a dummy TF or change the eulerr constructor, palette, or
set definitions.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
from pathlib import Path


REPLACEMENTS = (
    (
        'jaspar_only_df <- data.frame(TR = jaspar_only, Database = "JASPAR_only")',
        'jaspar_only_df <- data.frame(TR = jaspar_only, Database = rep("JASPAR_only", length(jaspar_only)))',
        47,
    ),
    (
        'cisbp_only_df <- data.frame(TR = cisbp_only, Database = "CISBP_only")',
        'cisbp_only_df <- data.frame(TR = cisbp_only, Database = rep("CISBP_only", length(cisbp_only)))',
        48,
    ),
    (
        'both_df <- data.frame(TR = both, Database = "Both")',
        'both_df <- data.frame(TR = both, Database = rep("Both", length(both)))',
        49,
    ),
)


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def motif_names(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Motif.Name" not in reader.fieldnames:
            raise RuntimeError(f"Motif.Name column missing from {path}")
        return {
            row["Motif.Name"].strip().upper()
            for row in reader
            if row["Motif.Name"].strip()
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("runtime_root", type=Path)
    ns = parser.parse_args()

    runtime_root = ns.runtime_root.resolve()
    source = runtime_root / "code/Figure3/01-Compare_JASPAR_CISBP_Motifs.R"
    before = source.read_text(encoding="utf-8")
    before_lines = before.splitlines()
    after = before
    for old, new, expected_line in REPLACEMENTS:
        observed = after.count(old)
        if observed != 1:
            raise RuntimeError(
                f"Expected one empty-set-unsafe expression at line {expected_line}; "
                f"observed {observed}"
            )
        after = after.replace(old, new)

    after_lines = after.splitlines()
    changed_lines = [
        index + 1
        for index, (left, right) in enumerate(zip(before_lines, after_lines))
        if left != right
    ]
    if changed_lines != [47, 48, 49]:
        raise RuntimeError(
            f"ANCHOR_EMPTY_SET_BUGFIX changed unexpected lines: {changed_lines}"
        )

    runtime_manifest = json.loads(
        (runtime_root / "audit/runtime_manifest.json").read_text(encoding="utf-8")
    )
    adapter_root = Path(runtime_manifest["adapter_root"])
    motif_root = (
        adapter_root
        / "data/Figure2/Motif/Combined_MotifEnrichment"
    )
    jaspar_path = (
        motif_root / "JASPAR_Combined_MotifEnrichment_DatabaseAvailability.tsv"
    )
    cisbp_path = (
        motif_root / "CISBP_Combined_MotifEnrichment_DatabaseAvailability.tsv"
    )
    jaspar = motif_names(jaspar_path)
    cisbp = motif_names(cisbp_path)
    jaspar_only = sorted(jaspar - cisbp)
    cisbp_only = sorted(cisbp - jaspar)
    both = sorted(jaspar & cisbp)
    if jaspar_only:
        raise RuntimeError(
            "Frozen LUAD JASPAR-only set is no longer empty; re-audit before applying bugfix."
        )

    source.write_text(after, encoding="utf-8", newline="\n")
    audit_root = runtime_root / "audit"
    audit_root.mkdir(parents=True, exist_ok=True)
    diff = "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile="strict-replay/code/Figure3/01-Compare_JASPAR_CISBP_Motifs.R",
            tofile="final-corrected/code/Figure3/01-Compare_JASPAR_CISBP_Motifs.R",
            n=2,
        )
    )
    (audit_root / "supplementary4c_empty_set_bugfix.diff").write_text(
        diff, encoding="utf-8"
    )
    receipt = {
        "panel": "SupplementaryFigure4C",
        "change_class": "ANCHOR_EMPTY_SET_BUGFIX",
        "official_source": "Figure3/01-Compare_JASPAR_CISBP_Motifs.R",
        "official_source_lines": "47;48;49",
        "changed_line_count": 3,
        "jaspar_only_count": len(jaspar_only),
        "cisbp_only_count": len(cisbp_only),
        "both_count": len(both),
        "jaspar_only_members": ";".join(jaspar_only),
        "cisbp_only_members": ";".join(cisbp_only),
        "both_members": ";".join(both),
        "dummy_tf_added": False,
        "eulerr_constructor_changed": False,
        "palette_changed": False,
        "strict_runtime_sha256": sha256(before),
        "corrected_runtime_sha256": sha256(after),
        "scientific_reason": (
            "Represent a genuine empty JASPAR-only set as a zero-row data.frame "
            "without inventing a TF."
        ),
        "jaspar_input": jaspar_path.as_posix(),
        "cisbp_input": cisbp_path.as_posix(),
    }
    with (audit_root / "supplementary4c_empty_set_bugfix_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(receipt), delimiter="\t")
        writer.writeheader()
        writer.writerow(receipt)
    print(
        "ANCHOR_EMPTY_SET_BUGFIX applied to official source lines 47-49 only; "
        f"sets JASPAR_only={len(jaspar_only)}, CISBP_only={len(cisbp_only)}, "
        f"Both={len(both)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

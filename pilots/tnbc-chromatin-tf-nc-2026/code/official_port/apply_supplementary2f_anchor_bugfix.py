#!/usr/bin/env python3
"""Apply the one acknowledged TNBC capsule bugfix to a materialised runtime.

The frozen upstream source remains untouched.  This changes exactly one object
name in Figure1/10 so the PDX long table is built from the PDX matrix rather
than the TCGA matrix, and records a machine-readable receipt plus a one-line
unified diff.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
from pathlib import Path


OLD = (
    "PDX_TF_activities_TNBC_Specific_long <- "
    "TCGA_TF_activities_TNBC_Specific %>%"
)
NEW = (
    "PDX_TF_activities_TNBC_Specific_long <- "
    "PDX_TF_activities_TNBC_Specific %>%"
)


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("runtime_root", type=Path)
    ns = parser.parse_args()

    runtime_root = ns.runtime_root.resolve()
    source = runtime_root / "code/Figure1/10-Plot_Figure1.R"
    before = source.read_text(encoding="utf-8")
    observed = before.count(OLD)
    if observed != 1:
        raise RuntimeError(
            f"Expected exactly one Supplementary Figure 2F capsule bug; observed {observed}"
        )
    after = before.replace(OLD, NEW)
    before_lines = before.splitlines()
    after_lines = after.splitlines()
    changed_lines = [
        index + 1
        for index, (left, right) in enumerate(zip(before_lines, after_lines))
        if left != right
    ]
    if changed_lines != [475]:
        raise RuntimeError(f"ANCHOR_BUGFIX changed unexpected lines: {changed_lines}")
    source.write_text(after, encoding="utf-8", newline="\n")

    audit_root = runtime_root / "audit"
    audit_root.mkdir(parents=True, exist_ok=True)
    diff = "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile="strict-replay/code/Figure1/10-Plot_Figure1.R",
            tofile="final-corrected/code/Figure1/10-Plot_Figure1.R",
            n=2,
        )
    )
    (audit_root / "supplementary2f_anchor_bugfix.diff").write_text(
        diff, encoding="utf-8"
    )
    receipt = {
        "panel": "SupplementaryFigure2F_PDX",
        "change_class": "ANCHOR_BUGFIX",
        "official_source": "Figure1/10-Plot_Figure1.R",
        "official_source_line": 475,
        "changed_line_count": 1,
        "old_text": OLD,
        "new_text": NEW,
        "strict_runtime_sha256": sha256(before),
        "corrected_runtime_sha256": sha256(after),
        "scientific_reason": "Build the PDX long table from the PDX activity object.",
    }
    with (audit_root / "supplementary2f_anchor_bugfix_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(receipt), delimiter="\t")
        writer.writeheader()
        writer.writerow(receipt)
    print("ANCHOR_BUGFIX applied to official source line 475 only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

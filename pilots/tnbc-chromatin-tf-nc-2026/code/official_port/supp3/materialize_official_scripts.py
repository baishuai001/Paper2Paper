#!/usr/bin/env python3
"""Lock and materialise the authors' Supplementary Figure 3 scripts.

The plotting code is not rewritten.  The only runtime substitutions are:
1. the absolute path used to source the authors' plotting-aesthetics file;
2. correction of the 04B parser/path variable typo
   ``args$consensus_peakset`` -> ``args$consensus_peakset_dir``.
Both substitutions are outside the plotting constructors and are recorded as
machine-readable receipts plus a unified diff.
"""

from __future__ import annotations

import csv
import difflib
import hashlib
import json
import shutil
import sys
from pathlib import Path


LOCKS = {
    "Figure2/02J2-Plot_Number_Peaks.R": "e3f975ffb877ba689155fc6a7b42c9a8d4d5a6e0ed8840e9c6c2028d8d9cd19f",
    "Figure2/02C-Peak_Saturation_Plot.R": "e5d7984de178318d6e5982d91cd65ac96647fa178f0eed2b7e905ebd2191dafa",
    "Figure2/03D-Pearson_Correlation_Accessibility_ConsensusPeakSet.R": "44745a0df1650246a50da41f65ad6f27a399034e835f7ad151ebec5ce6a3818f",
    "Figure2/04B-GenomicAnnotation.R": "411150cc1d6e3b7c5545f69e0b8e541ff91bc23c353aa0eed93a13841a9bf9f6",
    "Static_Scripts/plotting_aesthetics.R": "52f1fc590c6cad927ee0d1b030ddc215e99c506d15b465f3c805a9cc0ca11c1d",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def replace_exact(text: str, old: str, new: str, expected: int, source: str, kind: str):
    observed = text.count(old)
    if observed != expected:
        raise RuntimeError(
            f"{source}: expected {expected} occurrence(s) for {kind}, observed {observed}"
        )
    return text.replace(old, new), {
        "source_file": source,
        "change_kind": kind,
        "old": old,
        "new": new,
        "expected_occurrences": expected,
        "observed_occurrences": observed,
    }


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: materialize_official_scripts.py OFFICIAL_CODE_ROOT WORK_ROOT")
    official_root = Path(sys.argv[1]).resolve()
    work_root = Path(sys.argv[2]).resolve()
    code_out = work_root / "runtime" / "code"
    audit_out = work_root / "runtime" / "audit"
    code_out.mkdir(parents=True, exist_ok=True)
    audit_out.mkdir(parents=True, exist_ok=True)

    source_rows = []
    originals: dict[str, str] = {}
    for rel, expected_hash in LOCKS.items():
        source = official_root / rel
        if not source.is_file():
            raise FileNotFoundError(source)
        observed_hash = sha256(source)
        if observed_hash.lower() != expected_hash.lower():
            raise RuntimeError(
                f"official source hash mismatch for {rel}: {observed_hash} != {expected_hash}"
            )
        source_rows.append(
            {
                "source_file": rel,
                "official_path": str(source),
                "sha256": observed_hash,
                "bytes": source.stat().st_size,
            }
        )
        if rel.endswith(".R"):
            originals[rel] = source.read_text(encoding="utf-8")

    aesthetic_source = official_root / "Static_Scripts" / "plotting_aesthetics.R"
    aesthetic_target = code_out / "Static_Scripts" / "plotting_aesthetics.R"
    aesthetic_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(aesthetic_source, aesthetic_target)
    aesthetic_runtime = str(aesthetic_target)

    substitutions = []
    diffs = []
    for rel in [
        "Figure2/02J2-Plot_Number_Peaks.R",
        "Figure2/02C-Peak_Saturation_Plot.R",
        "Figure2/03D-Pearson_Correlation_Accessibility_ConsensusPeakSet.R",
        "Figure2/04B-GenomicAnnotation.R",
    ]:
        source_text = originals[rel]
        runtime_text = source_text
        if 'source("/code/Static_Scripts/plotting_aesthetics.R")' in runtime_text:
            runtime_text, receipt = replace_exact(
                runtime_text,
                'source("/code/Static_Scripts/plotting_aesthetics.R")',
                f'source("{aesthetic_runtime}")',
                1,
                rel,
                "PATH_ONLY",
            )
            substitutions.append(receipt)
        if rel == "Figure2/04B-GenomicAnnotation.R":
            runtime_text, receipt = replace_exact(
                runtime_text,
                "args$consensus_peakset",
                "args$consensus_peakset_dir",
                1,
                rel,
                "PATH_ARGUMENT_TYPO_ONLY",
            )
            substitutions.append(receipt)

        target = code_out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(runtime_text, encoding="utf-8", newline="")
        if runtime_text != source_text:
            diffs.extend(
                difflib.unified_diff(
                    source_text.splitlines(keepends=True),
                    runtime_text.splitlines(keepends=True),
                    fromfile=f"official/{rel}",
                    tofile=f"runtime/{rel}",
                )
            )

    with (audit_out / "official_source_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(source_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(source_rows)
    with (audit_out / "runtime_substitution_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(substitutions[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(substitutions)
    (audit_out / "runtime_unified.diff").write_text("".join(diffs), encoding="utf-8")
    (audit_out / "materialization_receipt.json").write_text(
        json.dumps(
            {
                "official_code_root": str(official_root),
                "work_root": str(work_root),
                "official_sources_verified": len(source_rows),
                "runtime_substitutions": len(substitutions),
                "plotting_constructor_rewrites": 0,
                "allowed_change_kinds": ["PATH_ONLY", "PATH_ARGUMENT_TYPO_ONLY"],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

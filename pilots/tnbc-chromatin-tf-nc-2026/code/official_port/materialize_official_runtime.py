#!/usr/bin/env python3
"""Materialise an auditable LUAD runtime from frozen TNBC plot sources.

The upstream files are never edited.  This program applies only three generic
path substitutions plus the explicitly whitelisted content substitutions in
``allowed_content_substitutions.json``.  It writes hashes, changed-line spans,
and a unified diff beside the runtime copy.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import re
from pathlib import Path
from typing import Any


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def line_numbers_with(text: str, token: str) -> list[int]:
    return [i for i, line in enumerate(text.splitlines(), 1) if token in line]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("adapter_root", type=Path)
    parser.add_argument("runtime_root", type=Path)
    parser.add_argument(
        "--official-port-root",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    parser.add_argument(
        "--validation-label",
        default="GSE41271",
        help="Dataset label inserted into label-only adapter substitutions.",
    )
    ns = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9._]*", ns.validation_label):
        raise ValueError("--validation-label must be a safe R name and display token")

    port_root = ns.official_port_root.resolve()
    upstream_code = port_root / "upstream_edf5314" / "code"
    adapter_root = ns.adapter_root.resolve()
    runtime_root = ns.runtime_root.resolve()
    audit_root = runtime_root / "audit"
    runtime_code = runtime_root / "code"
    runtime_code.mkdir(parents=True, exist_ok=True)
    audit_root.mkdir(parents=True, exist_ok=True)

    rules: list[dict[str, Any]] = json.loads(
        (port_root / "allowed_content_substitutions.json").read_text(encoding="utf-8")
    )
    rules_by_source: dict[str, list[dict[str, Any]]] = {}
    for rule in rules:
        rules_by_source.setdefault(rule["source"], []).append(rule)

    source_files = sorted(
        p for p in upstream_code.rglob("*") if p.is_file()
    )
    lock_path = port_root / "OFFICIAL_SOURCE_LOCK.tsv"
    with lock_path.open("r", encoding="utf-8", newline="") as handle:
        source_lock = {
            row["relative_path"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    observed_relpaths = {
        path.relative_to(upstream_code).as_posix() for path in source_files
    }
    if observed_relpaths != set(source_lock):
        raise RuntimeError(
            "Frozen source inventory differs from OFFICIAL_SOURCE_LOCK.tsv: "
            f"extra={sorted(observed_relpaths - set(source_lock))}, "
            f"missing={sorted(set(source_lock) - observed_relpaths)}"
        )
    source_manifest: list[dict[str, Any]] = []
    substitution_manifest: list[dict[str, Any]] = []
    all_diffs: list[str] = []

    adapter_data = (adapter_root / "data").as_posix().rstrip("/") + "/"
    adapter_results = (adapter_root / "results").as_posix().rstrip("/") + "/"
    aesthetics_adapter = (port_root / "plotting_aesthetics_LUAD_adapter.R").as_posix()

    for source_path in source_files:
        rel = source_path.relative_to(upstream_code).as_posix()
        original = source_path.read_text(encoding="utf-8")
        locked = source_lock[rel]
        if sha256_text(original) != locked["sha256"]:
            raise RuntimeError(f"Frozen upstream hash mismatch: {rel}")
        if len(original.splitlines()) != int(locked["line_count"]):
            raise RuntimeError(f"Frozen upstream line-count mismatch: {rel}")
        changed = original

        generic_rules = [
            {
                "rule_id": "PATH_DATA_ROOT",
                "old": "/data/",
                "new": adapter_data,
                "class": "PATH_ONLY",
                "reason": "Point the official reader at the isolated LUAD adapter tree.",
            },
            {
                "rule_id": "PATH_RESULTS_ROOT",
                "old": "/results/",
                "new": adapter_results,
                "class": "PATH_ONLY",
                "reason": "Point the official writer at the isolated official-port result tree.",
            },
            {
                "rule_id": "PATH_AESTHETICS_ADAPTER",
                "old": 'source("/code/Static_Scripts/plotting_aesthetics.R")',
                "new": f'source("{aesthetics_adapter}")',
                "class": "PATH_ONLY",
                "reason": "Load the frozen palette through the LUAD label-key adapter.",
            },
        ]

        for rule in generic_rules:
            count = changed.count(rule["old"])
            if count:
                before_lines = line_numbers_with(changed, rule["old"])
                changed = changed.replace(rule["old"], rule["new"])
                substitution_manifest.append(
                    {
                        "source": rel,
                        "rule_id": rule["rule_id"],
                        "class": rule["class"],
                        "expected_count": count,
                        "observed_count": count,
                        "source_lines": ",".join(map(str, before_lines)),
                        "old_text": rule["old"],
                        "new_text": rule["new"],
                        "reason": rule["reason"],
                    }
                )

        for rule in rules_by_source.get(rel, []):
            old_text = (
                rule["old"]
                .replace("@@ADAPTER_DATA_ROOT@@", adapter_data.rstrip("/"))
                .replace("@@ADAPTER_RESULTS_ROOT@@", adapter_results.rstrip("/"))
                .replace("@@VALIDATION_LABEL@@", ns.validation_label)
            )
            observed = changed.count(old_text)
            expected = int(rule["expected_count"])
            if observed != expected:
                raise RuntimeError(
                    f"{rel} {rule['rule_id']}: expected {expected} occurrences, "
                    f"observed {observed}"
                )
            before_lines = line_numbers_with(changed, old_text.splitlines()[0])
            new_text = (
                rule["new"]
                .replace("@@ADAPTER_DATA_ROOT@@", adapter_data.rstrip("/"))
                .replace("@@ADAPTER_RESULTS_ROOT@@", adapter_results.rstrip("/"))
                .replace("@@VALIDATION_LABEL@@", ns.validation_label)
            )
            changed = changed.replace(old_text, new_text)
            substitution_manifest.append(
                {
                    "source": rel,
                    "rule_id": rule["rule_id"],
                    "class": rule["class"],
                    "expected_count": expected,
                    "observed_count": observed,
                    "source_lines": ",".join(map(str, before_lines)),
                    "old_text": old_text.replace("\n", "\\n"),
                    "new_text": new_text.replace("\n", "\\n"),
                    "reason": rule["reason"],
                }
            )

        destination = runtime_code / rel
        if "@@VALIDATION_LABEL@@" in changed:
            raise RuntimeError(f"Unresolved validation-label placeholder: {rel}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(changed, encoding="utf-8", newline="\n")

        diff = list(
            difflib.unified_diff(
                original.splitlines(keepends=True),
                changed.splitlines(keepends=True),
                fromfile=f"upstream_edf5314/code/{rel}",
                tofile=f"runtime/code/{rel}",
            )
        )
        all_diffs.extend(diff)
        source_manifest.append(
            {
                "source": rel,
                "upstream_sha256": sha256_text(original),
                "runtime_sha256": sha256_text(changed),
                "upstream_lines": len(original.splitlines()),
                "runtime_lines": len(changed.splitlines()),
                "changed": original != changed,
                "content_rule_count": len(rules_by_source.get(rel, [])),
            }
        )

    unused_sources = sorted(set(rules_by_source) - {p["source"] for p in source_manifest})
    if unused_sources:
        raise RuntimeError(f"Rules name missing upstream sources: {unused_sources}")

    def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(audit_root / "runtime_source_manifest.tsv", source_manifest)
    write_tsv(audit_root / "runtime_substitution_manifest.tsv", substitution_manifest)
    (audit_root / "runtime_unified.diff").write_text("".join(all_diffs), encoding="utf-8")
    (audit_root / "runtime_manifest.json").write_text(
        json.dumps(
            {
                "official_commit": "edf5314ce5b9ee0e2f88b2310e7c2df5619ad888",
                "upstream_code_root": upstream_code.as_posix(),
                "adapter_root": adapter_root.as_posix(),
                "runtime_root": runtime_root.as_posix(),
                "validation_label": ns.validation_label,
                "sources": source_manifest,
                "substitutions": substitution_manifest,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Materialised {len(source_manifest)} official source files in {runtime_code}")
    print(f"Recorded {len(substitution_manifest)} substitutions in {audit_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

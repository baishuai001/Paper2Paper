#!/usr/bin/env python3
"""Verify and materialize line-exact TNBC CodeOcean plotting blocks.

This program is intentionally mechanical.  It never interprets R and it never
generates plotting code.  Every output is a source-line slice from the frozen
CodeOcean commit, optionally followed by literal substitutions declared in
``source_manifest.tsv``.  A unified diff and cryptographic receipt make those
substitutions auditable.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def replacements(spec: str) -> list[tuple[str, str]]:
    answer: list[tuple[str, str]] = []
    if not spec.strip():
        return answer
    for item in spec.split(";"):
        if "=>" not in item:
            raise ValueError(f"Malformed replacement: {item!r}")
        old, new = item.split("=>", 1)
        if not old:
            raise ValueError("Empty replacement source is forbidden")
        answer.append((old, new))
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--official-code-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=HERE / "source_manifest.tsv")
    args = parser.parse_args()

    code_root = args.official_code_root.resolve()
    out = args.output_dir.resolve()
    exact_dir = out / "exact"
    patched_dir = out / "patched"
    diff_dir = out / "diff"
    for path in (exact_dir, patched_dir, diff_dir):
        path.mkdir(parents=True, exist_ok=True)

    receipts: list[dict[str, object]] = []
    with args.manifest.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {
        "block_id", "panel", "upstream_file", "sha256", "line_start",
        "line_end", "mode", "allowed_replacements", "note",
    }
    if not rows or set(rows[0]) != required:
        raise SystemExit("source_manifest.tsv has an unexpected schema")

    verified_files: dict[Path, str] = {}
    for row in rows:
        source = (code_root / row["upstream_file"]).resolve()
        try:
            source.relative_to(code_root)
        except ValueError as exc:
            raise SystemExit(f"Source escapes official root: {source}") from exc
        if not source.is_file():
            raise SystemExit(f"Missing official source: {source}")
        observed = verified_files.setdefault(source, sha256(source))
        expected = row["sha256"].lower()
        if observed.lower() != expected:
            raise SystemExit(
                f"Checksum mismatch for {row['upstream_file']}: "
                f"expected {expected}, observed {observed}"
            )

        raw = source.read_bytes()
        lines = raw.splitlines(keepends=True)
        start, end = int(row["line_start"]), int(row["line_end"])
        if start < 1 or end < start or end > len(lines):
            raise SystemExit(f"Invalid line range for {row['block_id']}: {start}-{end}")
        exact_bytes = b"".join(lines[start - 1 : end])
        exact_path = exact_dir / f"{row['block_id']}.R"
        exact_path.write_bytes(exact_bytes)

        exact_text = exact_bytes.decode("utf-8")
        patched_text = exact_text
        applied: list[dict[str, object]] = []
        for old, new in replacements(row["allowed_replacements"]):
            count = patched_text.count(old)
            if count == 0:
                raise SystemExit(
                    f"Declared replacement {old!r} is absent from {row['block_id']}"
                )
            patched_text = patched_text.replace(old, new)
            applied.append({"from": old, "to": new, "count": count})

        patched_path = patched_dir / f"{row['block_id']}.R"
        patched_path.write_text(patched_text, encoding="utf-8", newline="")
        diff = "".join(
            difflib.unified_diff(
                exact_text.splitlines(keepends=True),
                patched_text.splitlines(keepends=True),
                fromfile=f"upstream/{row['upstream_file']}:{start}-{end}",
                tofile=f"patched/{row['block_id']}.R",
            )
        )
        diff_path = diff_dir / f"{row['block_id']}.diff"
        diff_path.write_text(diff, encoding="utf-8", newline="")

        receipts.append(
            {
                **row,
                "source_sha256_observed": observed,
                "exact_block_sha256": sha256(exact_path),
                "patched_block_sha256": sha256(patched_path),
                "replacement_receipt": applied,
                "exact_path": str(exact_path),
                "patched_path": str(patched_path),
                "diff_path": str(diff_path),
            }
        )

    receipt_path = out / "materialization_receipt.json"
    receipt_path.write_text(
        json.dumps(
            {
                "capsule": "7227095/v1",
                "commit": "edf5314",
                "official_code_root": str(code_root),
                "manifest_sha256": sha256(args.manifest),
                "blocks": receipts,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(receipt_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

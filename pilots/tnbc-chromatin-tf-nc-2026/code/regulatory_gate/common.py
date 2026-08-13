#!/usr/bin/env python3
"""Shared, dependency-light helpers for the CRC regulatory gate."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest().upper()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="\n", dir=path.parent, delete=False
    ) as stream:
        stream.write(text)
        temporary = Path(stream.name)
    os.replace(temporary, path)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    p = np.clip(values[finite], 0.0, 1.0)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[finite] = restored
    return result


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "<na>"} else text


def read_h5_column(group: h5py.Group, name: str) -> np.ndarray:
    """Read a plain or AnnData categorical HDF5 column."""
    item = group[name]
    if isinstance(item, h5py.Dataset):
        return np.asarray([clean_text(value) for value in item[:]], dtype=object)
    if not isinstance(item, h5py.Group) or "codes" not in item or "categories" not in item:
        raise ValueError(f"Unsupported H5AD column representation: {name}")
    codes = np.asarray(item["codes"][:], dtype=int)
    categories = np.asarray([clean_text(value) for value in item["categories"][:]], dtype=object)
    values = np.full(len(codes), "", dtype=object)
    valid = (codes >= 0) & (codes < len(categories))
    values[valid] = categories[codes[valid]]
    return values


def output_manifest(paths: list[Path]) -> dict[str, dict[str, object]]:
    return {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in sorted(paths)
    }

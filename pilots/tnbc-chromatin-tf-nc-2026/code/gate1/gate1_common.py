#!/usr/bin/env python3
"""Shared fail-closed helpers for the CRC TF gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


MISSING_TEXT = {"", "nan", "none", "na", "n/a", "unknown", "not available"}


def decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in "SUO":
        return np.asarray(
            [value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value) for value in values],
            dtype=object,
        )
    return values.astype(object)


def read_obs_column(obs: h5py.Group, column: str) -> np.ndarray:
    if column not in obs:
        raise ValueError(f"H5AD obs missing required column: {column}")
    item = obs[column]
    if isinstance(item, h5py.Dataset):
        return decode(item[:])
    if not isinstance(item, h5py.Group) or "codes" not in item or "categories" not in item:
        raise ValueError(f"unsupported H5AD obs encoding for {column}: {type(item).__name__}")
    codes = np.asarray(item["codes"][:], dtype=np.int64)
    categories = decode(item["categories"][:])
    result = np.empty(codes.shape[0], dtype=object)
    result[:] = ""
    valid = codes >= 0
    if valid.any():
        if int(codes[valid].max()) >= len(categories):
            raise ValueError(f"categorical code exceeds categories for {column}")
        result[valid] = categories[codes[valid]]
    return result


def clean_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in MISSING_TEXT else text


def join_unique(values: pd.Series) -> str:
    return "|".join(sorted({clean_text(value) for value in values if clean_text(value)}))


def unique_count(values: pd.Series) -> int:
    return len({clean_text(value) for value in values if clean_text(value)})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def bh_fdr(pvalues: np.ndarray | pd.Series) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    selected = values[finite]
    order = np.argsort(selected)
    ranks = np.arange(1, len(selected) + 1, dtype=float)
    adjusted = selected[order] * len(selected) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(selected)
    restored[order] = np.minimum(adjusted, 1.0)
    result[finite] = restored
    return result


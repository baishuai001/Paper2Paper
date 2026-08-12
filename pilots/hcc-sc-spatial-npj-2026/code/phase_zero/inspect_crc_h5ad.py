#!/usr/bin/env python3
"""Inspect a CRC-atlas H5AD without loading its full expression matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def json_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.ndarray):
        return [json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def dataframe_columns(group: h5py.Group) -> list[str]:
    order = group.attrs.get("column-order")
    if order is not None:
        return [str(json_value(value)) for value in order]
    return sorted(key for key in group.keys() if key != "_index")


def index_length(group: h5py.Group) -> int:
    index_name = json_value(group.attrs.get("_index", "_index"))
    if index_name in group and isinstance(group[index_name], h5py.Dataset):
        return int(group[index_name].shape[0])
    if "_index" in group and isinstance(group["_index"], h5py.Dataset):
        return int(group["_index"].shape[0])
    raise ValueError(f"cannot determine dataframe index length for {group.name}")


def encoding_type(item: h5py.Group | h5py.Dataset) -> str:
    return str(json_value(item.attrs.get("encoding-type", "dataset" if isinstance(item, h5py.Dataset) else "group")))


def column_receipt(group: h5py.Group, column: str) -> dict[str, object]:
    item = group[column]
    receipt: dict[str, object] = {
        "column": column,
        "storage": "dataset" if isinstance(item, h5py.Dataset) else "group",
        "encoding": encoding_type(item),
    }
    if isinstance(item, h5py.Dataset):
        receipt.update({"shape": "x".join(map(str, item.shape)), "dtype": str(item.dtype), "categories": ""})
    else:
        codes = item.get("codes")
        categories = item.get("categories")
        receipt.update(
            {
                "shape": "x".join(map(str, codes.shape)) if isinstance(codes, h5py.Dataset) else "",
                "dtype": str(codes.dtype) if isinstance(codes, h5py.Dataset) else "",
                "categories": int(categories.shape[0]) if isinstance(categories, h5py.Dataset) else "",
            }
        )
    return receipt


def sample_numeric(dataset: h5py.Dataset, limit: int = 100_000) -> dict[str, object]:
    if dataset.size == 0 or dataset.dtype.kind not in "iufb":
        return {"sample_n": 0, "sample_min": "", "sample_max": "", "nonnegative": "", "fraction_noninteger": ""}
    if dataset.ndim == 1:
        values = dataset[: min(limit, dataset.shape[0])]
    elif dataset.ndim == 2:
        rows = max(1, min(dataset.shape[0], limit // max(1, min(dataset.shape[1], 1000))))
        cols = min(dataset.shape[1], max(1, limit // rows))
        values = dataset[:rows, :cols]
    else:
        values = dataset[tuple(slice(0, min(size, 10)) for size in dataset.shape)]
    array = np.asarray(values).ravel()
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return {"sample_n": 0, "sample_min": "", "sample_max": "", "nonnegative": "", "fraction_noninteger": ""}
    return {
        "sample_n": int(finite.size),
        "sample_min": float(np.min(finite)),
        "sample_max": float(np.max(finite)),
        "nonnegative": bool(np.min(finite) >= 0),
        "fraction_noninteger": float(np.mean(np.abs(finite - np.rint(finite)) > 1e-6)),
    }


def matrix_receipt(name: str, item: h5py.Group | h5py.Dataset) -> dict[str, object]:
    row: dict[str, object] = {
        "matrix": name,
        "storage": "dataset" if isinstance(item, h5py.Dataset) else "group",
        "encoding": encoding_type(item),
        "shape": "",
        "dtype": "",
        "sample_n": 0,
        "sample_min": "",
        "sample_max": "",
        "nonnegative": "",
        "fraction_noninteger": "",
    }
    if isinstance(item, h5py.Dataset):
        row["shape"] = "x".join(map(str, item.shape))
        row["dtype"] = str(item.dtype)
        row.update(sample_numeric(item))
    elif "data" in item and isinstance(item["data"], h5py.Dataset):
        shape = item.attrs.get("shape", "")
        if isinstance(shape, np.ndarray):
            row["shape"] = "x".join(map(str, shape.tolist()))
        else:
            row["shape"] = str(json_value(shape))
        row["dtype"] = str(item["data"].dtype)
        row.update(sample_numeric(item["data"]))
    return row


def write_tsv(rows: list[dict[str, object]], fields: list[str], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inspect_h5ad(h5ad: Path, output_dir: Path, expected_bytes: int | None = None) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD is missing or empty: {h5ad}")
    actual_bytes = h5ad.stat().st_size
    if expected_bytes is not None and actual_bytes != expected_bytes:
        raise ValueError(f"H5AD byte mismatch: expected {expected_bytes}, got {actual_bytes}")
    output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5ad, "r") as handle:
        for required in ("obs", "var"):
            if required not in handle or not isinstance(handle[required], h5py.Group):
                raise ValueError(f"H5AD missing required dataframe group: {required}")
        if "X" not in handle and "layers" not in handle and "raw" not in handle:
            raise ValueError("H5AD contains no X, layers or raw matrix")
        obs = handle["obs"]
        var = handle["var"]
        n_obs = index_length(obs)
        n_vars = index_length(var)
        obs_rows = [column_receipt(obs, column) for column in dataframe_columns(obs)]
        var_rows = [column_receipt(var, column) for column in dataframe_columns(var)]
        matrix_rows: list[dict[str, object]] = []
        if "X" in handle:
            matrix_rows.append(matrix_receipt("X", handle["X"]))
        if "layers" in handle and isinstance(handle["layers"], h5py.Group):
            for key in sorted(handle["layers"].keys()):
                matrix_rows.append(matrix_receipt(f"layers/{key}", handle["layers"][key]))
        if "raw" in handle and isinstance(handle["raw"], h5py.Group) and "X" in handle["raw"]:
            matrix_rows.append(matrix_receipt("raw/X", handle["raw"]["X"]))
        root_keys = sorted(handle.keys())
        root_attrs = {str(key): json_value(value) for key, value in handle.attrs.items()}

    obs_path = output_dir / "P0_h5ad_obs_columns.tsv"
    var_path = output_dir / "P0_h5ad_var_columns.tsv"
    matrix_path = output_dir / "P0_h5ad_matrices.tsv"
    receipt_path = output_dir / "P0_h5ad_receipt.json"
    fields = ["column", "storage", "encoding", "shape", "dtype", "categories"]
    write_tsv(obs_rows, fields, obs_path)
    write_tsv(var_rows, fields, var_path)
    matrix_fields = [
        "matrix",
        "storage",
        "encoding",
        "shape",
        "dtype",
        "sample_n",
        "sample_min",
        "sample_max",
        "nonnegative",
        "fraction_noninteger",
    ]
    write_tsv(matrix_rows, matrix_fields, matrix_path)

    receipt: dict[str, object] = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "bytes": actual_bytes,
        "sha256": sha256_file(h5ad),
        "n_obs": n_obs,
        "n_vars": n_vars,
        "root_keys": root_keys,
        "root_attrs": root_attrs,
        "obs_columns": len(obs_rows),
        "var_columns": len(var_rows),
        "matrix_locations": [row["matrix"] for row in matrix_rows],
        "outputs": {},
    }
    for path in (obs_path, var_path, matrix_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-bytes", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = inspect_h5ad(args.h5ad, args.output_dir, args.expected_bytes)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

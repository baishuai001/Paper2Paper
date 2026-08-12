#!/usr/bin/env python3
"""Run CopyKAT and SCEVAN over a frozen real-data CRC sample panel."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from prepare_crc_cna_inputs import prepare_cna_input
from summarize_crc_cna_concordance import summarize as summarize_concordance
from summarize_crc_copykat import summarize as summarize_copykat


REQUIRED_PANEL_COLUMNS = {
    "panel_order",
    "dataset",
    "sample_id",
    "analysis_patient_id",
    "sample_type",
    "tissue",
    "eligible",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def write_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_panel(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"CNA panel missing or empty: {path}")
    panel = pd.read_csv(path, sep="\t")
    missing = REQUIRED_PANEL_COLUMNS - set(panel.columns)
    if missing:
        raise ValueError(f"CNA panel fields missing: {sorted(missing)}")
    if panel.empty:
        raise ValueError("CNA panel is empty")
    if panel["panel_order"].duplicated().any() or panel["sample_id"].isna().any():
        raise ValueError("CNA panel order or sample IDs are missing or duplicated")
    if panel.duplicated(["sample_id", "analysis_patient_id"]).any():
        raise ValueError("CNA panel contains duplicate sample-patient analysis units")
    eligible = panel["eligible"].astype(str).str.lower().isin({"true", "1", "yes"})
    if not eligible.all():
        raise ValueError("CNA panel contains ineligible rows")
    return panel.sort_values("panel_order").reset_index(drop=True)


def run_command(command: list[str], log_path: Path, environment: dict[str, str] | None = None) -> dict[str, object]:
    started = utc_now()
    with log_path.open("wb") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    return {
        "command": command,
        "started_at": started,
        "finished_at": utc_now(),
        "exit_code": result.returncode,
        "log": str(log_path.resolve()),
        "log_bytes": log_path.stat().st_size,
        "log_sha256": sha256_file(log_path),
    }


def run_panel(
    h5ad: Path,
    panel_path: Path,
    output_root: Path,
    cna_runner: Path,
    scevan_r_lib: Path,
    rscript: str,
    cores: int,
    seed: int,
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    if not cna_runner.is_file() or cna_runner.stat().st_size == 0:
        raise FileNotFoundError(f"R CNA runner missing or empty: {cna_runner}")
    if not scevan_r_lib.is_dir():
        raise FileNotFoundError(f"SCEVAN R library missing: {scevan_r_lib}")
    if cores < 1:
        raise ValueError("cores must be positive")
    panel = load_panel(panel_path)
    if output_root.exists():
        raise FileExistsError(f"refusing to mix with existing CNA panel output: {output_root}")
    output_root.mkdir(parents=True)
    state_path = output_root / "P0_cna_panel_run_state.json"
    receipt_path = output_root / "P0_cna_panel_run_receipt.json"
    state: dict[str, object] = {
        "status": "running",
        "started_at": utc_now(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "h5ad": str(h5ad.resolve()),
        "panel": str(panel_path.resolve()),
        "panel_sha256": sha256_file(panel_path),
        "cna_runner": str(cna_runner.resolve()),
        "cna_runner_sha256": sha256_file(cna_runner),
        "scevan_r_lib": str(scevan_r_lib.resolve()),
        "rscript": rscript,
        "cores": cores,
        "seed": seed,
        "units": [],
    }
    write_json(state_path, state)
    failed_units = 0
    for row in panel.itertuples(index=False):
        order = int(row.panel_order)
        unit_dir = output_root / f"unit_{order:02d}"
        input_dir = unit_dir / "input"
        copykat_dir = unit_dir / "copykat"
        scevan_dir = unit_dir / "scevan"
        summary_dir = unit_dir / "summary"
        unit_dir.mkdir()
        unit: dict[str, object] = {
            "panel_order": order,
            "dataset": str(row.dataset),
            "sample_id": str(row.sample_id),
            "patient_id": str(row.analysis_patient_id),
            "sample_type": str(row.sample_type),
            "tissue": str(row.tissue),
            "status": "running",
            "started_at": utc_now(),
            "unit_dir": str(unit_dir.resolve()),
        }
        state["units"].append(unit)
        write_json(state_path, state)
        try:
            prepare_receipt = prepare_cna_input(
                h5ad,
                input_dir,
                str(row.sample_id),
                str(row.analysis_patient_id),
                max_cancer=300,
                max_epithelial=100,
                max_reference=300,
                min_cancer=50,
                min_reference=50,
                seed=seed + order,
            )
            unit["prepare"] = prepare_receipt
        except Exception as exc:  # preserve the unit failure and continue independent samples
            unit["status"] = "prepare_failed"
            unit["error"] = f"{type(exc).__name__}: {exc}"
            unit["finished_at"] = utc_now()
            failed_units += 1
            write_json(state_path, state)
            continue

        sample_name = f"panel_{order:02d}"
        copykat_command = [
            rscript,
            str(cna_runner),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(copykat_dir),
            "--method",
            "copykat",
            "--sample-name",
            sample_name,
            "--cores",
            str(cores),
            "--seed",
            str(seed + order),
        ]
        copykat_run = run_command(copykat_command, unit_dir / "copykat.log")
        unit["copykat_run"] = copykat_run
        if copykat_run["exit_code"] != 0:
            unit["status"] = "copykat_failed"
            unit["finished_at"] = utc_now()
            failed_units += 1
            write_json(state_path, state)
            continue

        scevan_environment = os.environ.copy()
        scevan_environment["R_LIBS_USER"] = str(scevan_r_lib)
        scevan_command = [
            rscript,
            str(cna_runner),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(scevan_dir),
            "--method",
            "scevan",
            "--sample-name",
            sample_name,
            "--cores",
            str(cores),
            "--seed",
            str(seed + order),
        ]
        scevan_run = run_command(scevan_command, unit_dir / "scevan.log", scevan_environment)
        unit["scevan_run"] = scevan_run
        if scevan_run["exit_code"] != 0:
            unit["status"] = "scevan_failed"
            unit["finished_at"] = utc_now()
            failed_units += 1
            write_json(state_path, state)
            continue

        try:
            copykat_receipt = summarize_copykat(
                input_dir / "cell_metadata.tsv",
                copykat_dir / "copykat_prediction.tsv",
                summary_dir / "copykat",
            )
            concordance_receipt = summarize_concordance(
                input_dir / "cell_metadata.tsv",
                copykat_dir / "copykat_prediction.tsv",
                scevan_dir / "scevan_prediction.tsv",
                summary_dir / "concordance",
            )
            unit["copykat_summary"] = copykat_receipt
            unit["concordance"] = concordance_receipt
            unit["status"] = "completed"
        except Exception as exc:
            unit["status"] = "summary_failed"
            unit["error"] = f"{type(exc).__name__}: {exc}"
            failed_units += 1
        unit["finished_at"] = utc_now()
        write_json(state_path, state)

    state["status"] = "completed" if failed_units == 0 else "completed_with_failures"
    state["finished_at"] = utc_now()
    state["units_total"] = len(panel)
    state["units_completed"] = sum(unit["status"] == "completed" for unit in state["units"])
    state["units_failed"] = failed_units
    state["claim_boundary"] = (
        "Panel execution qualifies the selected diagnostic samples only. Calls use author Cancer "
        "cell candidates and predefined normal references; neither algorithm is DNA truth, and "
        "state comparisons remain observational and potentially dataset-confounded."
    )
    write_json(state_path, state)
    write_json(receipt_path, state)
    return state


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--cna-runner", required=True, type=Path)
    parser.add_argument("--scevan-r-lib", required=True, type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = run_panel(
        args.h5ad,
        args.panel,
        args.output_root,
        args.cna_runner,
        args.scevan_r_lib,
        args.rscript,
        args.cores,
        args.seed,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if receipt["units_failed"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

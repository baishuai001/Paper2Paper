#!/usr/bin/env python3
"""Download and freeze a signed human CollecTRI network from OmniPath."""

from __future__ import annotations

import argparse
import io
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import numpy as np

from gate1_common import sha256_file, write_json


URL = (
    "https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&format=tsv&genesymbols=yes&"
    "fields=sources,references,curation_effort&license=academic"
)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)


def run(output_dir: Path) -> dict[str, object]:
    request = urllib.request.Request(URL, headers={"User-Agent": "Paper2Paper-CRC-gate1/1.0", "Accept": "text/tab-separated-values"})
    with urllib.request.urlopen(request, timeout=180) as response:
        raw = response.read()
        status = int(response.status)
        headers = dict(response.headers.items())
    frame = pd.read_csv(io.BytesIO(raw), sep="\t", low_memory=False)
    required = {"source_genesymbol", "target_genesymbol", "is_directed", "consensus_stimulation", "consensus_inhibition"}
    if required - set(frame):
        raise ValueError(f"CollecTRI missing columns: {sorted(required-set(frame))}")
    frame = frame[numeric(frame["is_directed"]).eq(1)].copy()
    frame["source"] = frame["source_genesymbol"].astype(str).str.strip()
    frame["target"] = frame["target_genesymbol"].astype(str).str.strip()
    stimulation = numeric(frame["consensus_stimulation"])
    inhibition = numeric(frame["consensus_inhibition"])
    frame["sign"] = np.where((stimulation == 1) & (inhibition == 0), 1, np.where((stimulation == 0) & (inhibition == 1), -1, 0))
    signed = frame[(frame["source"].ne("")) & (frame["target"].ne("")) & (frame["sign"].ne(0))][["source", "target", "sign"]]
    signed = signed.groupby(["source", "target"], as_index=False)["sign"].sum()
    signed["sign"] = signed["sign"].clip(-1, 1)
    signed = signed[signed["sign"].ne(0)].sort_values(["source", "target"]).reset_index(drop=True)
    if len(signed) < 10000:
        raise ValueError(f"Implausibly few signed CollecTRI edges: {len(signed)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / "G4_collectri_raw.tsv"
    network_path = output_dir / "G4_collectri_signed.tsv"
    raw_path.write_bytes(raw)
    signed.to_csv(network_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "passed", "generated_at": datetime.now(timezone.utc).isoformat(), "url": URL,
        "http_status": status, "http_headers": headers, "license_query": "academic", "raw_rows": int(len(frame)),
        "signed_edges": int(len(signed)), "TFs": int(signed["source"].nunique()), "targets": int(signed["target"].nunique()),
        "outputs": {p.name: {"bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in (raw_path, network_path)},
    }
    write_json(output_dir / "G4_collectri_receipt.json", receipt)
    print(receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--output-dir", required=True, type=Path); args = parser.parse_args()
    run(args.output_dir); return 0


if __name__ == "__main__":
    raise SystemExit(main())

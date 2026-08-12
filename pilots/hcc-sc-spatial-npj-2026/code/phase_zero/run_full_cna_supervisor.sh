#!/usr/bin/env bash
set -euo pipefail

# Reusable wrapper. The exact server-bound script that produced the recorded
# run is retained separately under evidence/phase0-paired-cna-review/.
ROOT=${PAPER2PAPER_SERVER_ROOT:?Set PAPER2PAPER_SERVER_ROOT to the repository root}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PHASE=${P2P_PHASE_ZERO_ROOT:-"$ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero"}
CODE=${P2P_PHASE_ZERO_CODE:-"$SCRIPT_DIR"}
SCREEN="$PHASE/outputs/cna_paired_screen_v1"
MAN="$PHASE/manifests/paired_cna_full_v1"
OUT="$PHASE/outputs/cna_paired_full_v1"
PY=${P2P_BOOTSTRAP_PYTHON:-"$PHASE/env/bootstrap/bin/python3"}
H5AD=${P2P_CRC_H5AD:-"$PHASE/raw/crc_atlas_cellxgene.h5ad"}
RUNNER="$CODE/run_crc_cna_panel_configurable.py"
COMBINER="$CODE/combine_crc_paired_cna_summaries.py"
CNA_R="$CODE/run_crc_cna.R"
SCEVAN_LIB=${P2P_SCEVAN_LIB:-"$PHASE/env/R-scevan-lib"}

test -s "$SCREEN/POST_SCREEN_COMPLETE"
test -s "$MAN/P0_paired_cna_full_shards_receipt.json"
test ! -e "$OUT"
test -s "$H5AD"
test -x "$PY"
test -s "$RUNNER"
test -s "$COMBINER"
test -s "$CNA_R"
mkdir -p "$OUT/logs"

run_shard() {
  local shard="$1"
  local seed=$((20260811 + shard * 1000))
  PYTHONPATH="$CODE" OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
    nice -n 19 "$PY" "$RUNNER" \
      --h5ad "$H5AD" \
      --panel "$MAN/P0_paired_cna_full_shard_$(printf '%02d' "$shard").tsv" \
      --output-root "$OUT/shard_$(printf '%02d' "$shard")" \
      --cna-runner "$CNA_R" \
      --scevan-r-lib "$SCEVAN_LIB" \
      --rscript Rscript \
      --cores 4 \
      --seed "$seed" \
      --max-cancer 100000 \
      --max-epithelial 100 \
      --max-reference 300 \
      --min-cancer 50 \
      --min-reference 50 \
      > "$OUT/logs/shard_$(printf '%02d' "$shard").stdout.log" 2>&1
}

failure=0
for first in 1 3; do
  second=$((first + 1))
  run_shard "$first" &
  first_pid=$!
  echo "$first_pid" > "$OUT/shard_$(printf '%02d' "$first").pid"
  run_shard "$second" &
  second_pid=$!
  echo "$second_pid" > "$OUT/shard_$(printf '%02d' "$second").pid"
  wait "$first_pid" || failure=1
  wait "$second_pid" || failure=1
  if [[ "$failure" -ne 0 ]]; then
    date -u +'%Y-%m-%dT%H:%M:%SZ' > "$OUT/FULL_CNA_FAILED"
    exit 2
  fi
done

for shard in 1 2 3 4; do
  receipt="$OUT/shard_$(printf '%02d' "$shard")/P0_cna_panel_run_receipt.json"
  "$PY" - "$receipt" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.is_file() or path.stat().st_size == 0:
    raise SystemExit(f"missing full-CNA shard receipt: {path}")
receipt = json.loads(path.read_text())
if receipt.get("status") != "completed" or receipt.get("units_failed") != 0:
    raise SystemExit(f"full-CNA shard did not complete cleanly: {path}")
PY
  PYTHONPATH="$CODE" "$PY" "$CODE/summarize_crc_cna_panel.py" \
    --panel "$MAN/P0_paired_cna_full_shard_$(printf '%02d' "$shard").tsv" \
    --run-root "$OUT/shard_$(printf '%02d' "$shard")" \
    --output-dir "$OUT/summary_$(printf '%02d' "$shard")" \
    > "$OUT/logs/summary_$(printf '%02d' "$shard").log" 2>&1
done

PYTHONPATH="$CODE" "$PY" "$COMBINER" \
  --panel "$MAN/P0_paired_cna_full_allowed.tsv" \
  --summary-dir "$OUT/summary_01" \
  --summary-dir "$OUT/summary_02" \
  --summary-dir "$OUT/summary_03" \
  --summary-dir "$OUT/summary_04" \
  --output-dir "$OUT/review" \
  --stage full \
  > "$OUT/logs/combine_full_review.log" 2>&1

date -u +'%Y-%m-%dT%H:%M:%SZ' > "$OUT/FULL_CNA_COMPLETE"

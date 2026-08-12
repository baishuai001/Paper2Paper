#!/usr/bin/env bash
set -euo pipefail

ROOT=/media/desk16/iy13202/projects/Paper2Paper
PHASE="$ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero"
CODE="$ROOT/pilots/hcc-sc-spatial-npj-2026/code/phase_zero"
MAN="$PHASE/manifests/paired_cna_v1"
OUT="$PHASE/outputs/cna_paired_screen_v1"
PY="$PHASE/env/bootstrap/bin/python3"

for shard in liu other; do
  test -s "$OUT/$shard.pid"
done

while :; do
  live=0
  for shard in liu other; do
    pid=$(cat "$OUT/$shard.pid")
    if kill -0 "$pid" 2>/dev/null; then
      live=1
    fi
  done
  if [[ "$live" -eq 0 ]]; then
    break
  fi
  date -u +'%Y-%m-%dT%H:%M:%SZ screen still running'
  sleep 60
done

"$PY" - "$OUT" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
for shard, expected in (("liu", 25), ("other", 18)):
    path = root / shard / "P0_cna_panel_run_receipt.json"
    if not path.is_file():
        raise SystemExit(f"missing run receipt: {path}")
    state = json.loads(path.read_text())
    if state.get("status") != "completed" or state.get("units_completed") != expected:
        raise SystemExit(
            f"screen shard failed: {shard} status={state.get('status')} "
            f"completed={state.get('units_completed')} failed={state.get('units_failed')}"
        )
print("screen receipts complete")
PY

PYTHONPATH="$CODE" "$PY" "$CODE/summarize_crc_cna_panel.py" \
  --panel "$MAN/P0_paired_cna_panel_liu.tsv" \
  --run-root "$OUT/liu" \
  --output-dir "$OUT/summary_liu"
PYTHONPATH="$CODE" "$PY" "$CODE/summarize_crc_cna_panel.py" \
  --panel "$MAN/P0_paired_cna_panel_other.tsv" \
  --run-root "$OUT/other" \
  --output-dir "$OUT/summary_other"
PYTHONPATH="$CODE" "$PY" "$MAN/combine_crc_paired_cna_summaries.py" \
  --panel "$MAN/P0_paired_cna_panel.tsv" \
  --summary-dir "$OUT/summary_liu" \
  --summary-dir "$OUT/summary_other" \
  --output-dir "$OUT/review" \
  --stage screen

date -u +'%Y-%m-%dT%H:%M:%SZ' > "$OUT/POST_SCREEN_COMPLETE"

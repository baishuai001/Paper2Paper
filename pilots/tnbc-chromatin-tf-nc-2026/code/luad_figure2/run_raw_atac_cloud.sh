#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_raw_atac_cloud.sh <figure2_run_root>" >&2
  exit 2
fi
RUN_ROOT=$(realpath "$1")
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
MANIFEST="$RUN_ROOT/audit/manifests/raw_atac/figure2_raw_atac_runs.tsv"
WORKER="$PROJECT_ROOT/pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2/process_raw_atac_sample.sh"
PARALLEL_JOBS=${PARALLEL_JOBS:-4}
ONLY_SYSTEM=${ONLY_SYSTEM:-}
[[ -s "$MANIFEST" ]] || { echo "Missing raw ATAC manifest: $MANIFEST" >&2; exit 2; }

running=0
failed=0
while IFS=$'\t' read -r system sample_id sample_slug run layout host_depletion source unit include; do
  [[ "$include" == "TRUE" ]] || continue
  [[ -z "$ONLY_SYSTEM" || "$system" == "$ONLY_SYSTEM" ]] || continue
  bash "$WORKER" "$RUN_ROOT" "$system" "$sample_id" "$sample_slug" "$run" "$layout" &
  running=$((running + 1))
  if (( running >= PARALLEL_JOBS )); then
    if ! wait -n; then failed=$((failed + 1)); fi
    running=$((running - 1))
  fi
done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')

while (( running > 0 )); do
  if ! wait -n; then failed=$((failed + 1)); fi
  running=$((running - 1))
done

if (( failed > 0 )); then
  echo "RAW_ATAC_WORKERS_FAILED=$failed" >&2
  exit 1
fi
mkdir -p "$RUN_ROOT/audit/raw_atac_qc"
scope=${ONLY_SYSTEM:-all}
printf 'completed_utc=%s\nworkers_failed=0\n' "$(date -u +%FT%TZ)" \
  > "$RUN_ROOT/audit/raw_atac_qc/raw_atac_driver.${scope}.complete"
if [[ -z "$ONLY_SYSTEM" ]]; then
  cp "$RUN_ROOT/audit/raw_atac_qc/raw_atac_driver.${scope}.complete" \
    "$RUN_ROOT/audit/raw_atac_qc/raw_atac_driver.complete"
fi
echo "RAW_ATAC_WORKERS_COMPLETE scope=$scope"

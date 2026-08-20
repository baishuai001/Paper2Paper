#!/usr/bin/env bash
set -uo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_anchor_consensus_motif_primary_driver.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOG_ROOT="$RUN_ROOT/logs/motif_primary/anchor_consensus_author"
EXIT_FILE="$LOG_ROOT/driver.exit"
mkdir -p "$LOG_ROOT"
rm -f "$EXIT_FILE"

status=99
finish() {
  printf '%s\n' "$status" > "$EXIT_FILE"
}
trap finish EXIT

PARALLEL_JOBS=${PARALLEL_JOBS:-6} \
THREADS_PER_JOB=${THREADS_PER_JOB:-2} \
bash "$CODE_DIR/run_anchor_consensus_motif_primary.sh" "$RUN_ROOT"
status=$?
exit "$status"

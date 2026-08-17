#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_motif_equivalence_audit.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$RUN_ROOT/logs/motif_equivalence"

python3 "$CODE_DIR/harmonize_motif_inputs.py" "$RUN_ROOT" \
  2>&1 | tee "$RUN_ROOT/logs/motif_equivalence/01_harmonize_inputs.log"

bash "$CODE_DIR/run_harmonized_homer.sh" "$RUN_ROOT" \
  2>&1 | tee "$RUN_ROOT/logs/motif_equivalence/02_homer.log"

python3 "$CODE_DIR/parse_harmonized_motifs.py" "$RUN_ROOT" \
  2>&1 | tee "$RUN_ROOT/logs/motif_equivalence/03_parse.log"

echo "MOTIF_EQUIVALENCE_AUDIT_COMPLETE"

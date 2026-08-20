#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 OFFICIAL_CODE_ROOT FIGURE1_ROOT FIGURE4_ROOT FIGURE5_ROOT PUBLICATION_ROOT WORK_ROOT RESULT_ROOT" >&2
  exit 64
fi

OFFICIAL_CODE_ROOT=$1
FIGURE1_ROOT=$2
FIGURE4_ROOT=$3
FIGURE5_ROOT=$4
PUBLICATION_ROOT=$5
WORK_ROOT=$6
RESULT_ROOT=$7
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if [[ -e "$WORK_ROOT" || -e "$RESULT_ROOT" ]]; then
  echo "Refusing to overwrite an existing audited run: $WORK_ROOT or $RESULT_ROOT" >&2
  exit 73
fi

mkdir -p "$WORK_ROOT/generated" "$WORK_ROOT/adapted" "$RESULT_ROOT"

python3 "$SCRIPT_DIR/materialize_official_blocks.py" \
  --official-code-root "$OFFICIAL_CODE_ROOT" \
  --output-dir "$WORK_ROOT/generated"

Rscript "$SCRIPT_DIR/prepare_luad_inputs.R" \
  "$FIGURE1_ROOT" "$FIGURE4_ROOT" "$FIGURE5_ROOT" \
  "$PUBLICATION_ROOT" "$WORK_ROOT/adapted"

Rscript "$SCRIPT_DIR/render_official_fig45.R" \
  "$WORK_ROOT/generated" "$WORK_ROOT/adapted" "$RESULT_ROOT"

python3 "$SCRIPT_DIR/verify_official_port.py" \
  --code-dir "$SCRIPT_DIR" \
  --generated-dir "$WORK_ROOT/generated" \
  --results-dir "$RESULT_ROOT" \
  --rscript "$(command -v Rscript)"

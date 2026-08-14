#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 CLOUD_RUN_ROOT" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
AUDIT_DIR="$RUN_ROOT/audit/figure1_annotations"
LOG_DIR="$RUN_ROOT/logs"
mkdir -p "$AUDIT_DIR" "$LOG_DIR"

find "$RUN_ROOT/results/tables" -maxdepth 1 -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum > "$AUDIT_DIR/frozen_numeric_sha256.before.tsv"

python "$SCRIPT_DIR/fetch_pdmr_annotations.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_DIR/fetch_pdmr_annotations.log"
Rscript "$SCRIPT_DIR/build_figure1_annotations.R" "$RUN_ROOT" \
  2>&1 | tee "$LOG_DIR/build_figure1_annotations.log"
Rscript "$SCRIPT_DIR/plot_figure1_annotations_and_supplements.R" "$RUN_ROOT" \
  2>&1 | tee "$LOG_DIR/plot_figure1_annotations_and_supplements.log"

find "$RUN_ROOT/results/tables" -maxdepth 1 -type f \
  ! -name 'SupplementaryFigure*' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "$AUDIT_DIR/frozen_numeric_sha256.after.tsv"

# The before list predates new supplementary tables.  Compare only pre-existing
# paths so a new plot-support table is not mistaken for a numerical mutation.
awk '{print $2}' "$AUDIT_DIR/frozen_numeric_sha256.before.tsv" > "$AUDIT_DIR/frozen_numeric_paths.before.txt"
while IFS= read -r path; do sha256sum "$path"; done < "$AUDIT_DIR/frozen_numeric_paths.before.txt" \
  > "$AUDIT_DIR/frozen_numeric_sha256.after_comparable.tsv"
if ! diff -u "$AUDIT_DIR/frozen_numeric_sha256.before.tsv" "$AUDIT_DIR/frozen_numeric_sha256.after_comparable.tsv" \
  > "$AUDIT_DIR/frozen_numeric_sha256.diff"; then
  echo "Frozen Figure 1 numerical outputs changed; refusing to accept annotation run." >&2
fi

python "$SCRIPT_DIR/verify_figure1_annotations.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_DIR/verify_figure1_annotations.log"

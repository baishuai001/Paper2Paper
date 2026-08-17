#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: wait_and_compile_motif_equivalence.sh <luad_root> <tnbc_root> <anchor_candidate_pid> <combined_candidate_pid> <full_pid>" >&2
  exit 2
fi

LUAD_ROOT=$(realpath "$1")
TNBC_ROOT=$(realpath "$2")
ANCHOR_CANDIDATE_PID=$3
COMBINED_CANDIDATE_PID=$4
FULL_PID=$5
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANCHOR_CANDIDATE_RECEIPT="$LUAD_ROOT/audit/motif_equivalence/candidate_anchor_equivalent_motif_receipt.json"
COMBINED_CANDIDATE_RECEIPT="$LUAD_ROOT/audit/motif_equivalence/candidate_database_motif_receipt.json"
FULL_RECEIPT="$LUAD_ROOT/audit/motif_equivalence/harmonized_motif_receipt.json"
MAX_POLLS=720  # six hours at 30 seconds per poll

for ((poll=1; poll<=MAX_POLLS; poll++)); do
  if [[ -s "$ANCHOR_CANDIDATE_RECEIPT" && -s "$COMBINED_CANDIDATE_RECEIPT" && -s "$FULL_RECEIPT" ]]; then
    python3 "$CODE_DIR/compile_motif_equivalence_audit.py" "$LUAD_ROOT" "$TNBC_ROOT"
    python3 "$CODE_DIR/validate_motif_equivalence_audit.py" "$LUAD_ROOT" "$TNBC_ROOT"
    echo "MOTIF_EQUIVALENCE_FINAL_COMPILE_COMPLETE"
    exit 0
  fi
  if [[ ! -s "$ANCHOR_CANDIDATE_RECEIPT" ]] && ! kill -0 "$ANCHOR_CANDIDATE_PID" 2>/dev/null; then
    echo "Anchor-equivalent candidate run stopped before producing its receipt" >&2
    exit 1
  fi
  if [[ ! -s "$COMBINED_CANDIDATE_RECEIPT" ]] && ! kill -0 "$COMBINED_CANDIDATE_PID" 2>/dev/null; then
    echo "Combined candidate-database run stopped before producing its receipt" >&2
    exit 1
  fi
  if [[ ! -s "$FULL_RECEIPT" ]] && ! kill -0 "$FULL_PID" 2>/dev/null; then
    echo "Full-database run stopped before producing its receipt" >&2
    exit 1
  fi
  sleep 30
done

echo "Timed out waiting for motif equivalence receipts" >&2
exit 1

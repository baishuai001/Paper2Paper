#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 STATE_ROOT FIGURE2_ROOT CODE_ROOT" >&2
  exit 64
fi

state_root="$(realpath "$1")"
figure2_root="$(realpath "$2")"
code_root="$(realpath "$3")"
cell_receipt="${state_root}/audit/state_classifier/state_classifier_receipt.json"
motif_receipt="${figure2_root}/audit/motif_background_sensitivity/final_receipt.json"
output_receipt="${state_root}/audit/tru_state_diagnostic/final_receipt.json"
log_root="${state_root}/logs"
exit_file="${log_root}/tru_state_diagnostic.exit"

mkdir -p "${log_root}"
rm -f "${exit_file}"
echo "$(date --iso-8601=seconds) waiting for cell-state and motif-background receipts"
while [[ ! -s "${cell_receipt}" || ! -s "${motif_receipt}" ]]; do
  sleep 300
done

if [[ -s "${output_receipt}" ]]; then
  echo "$(date --iso-8601=seconds) TRU diagnostic receipt already exists"
  printf '0\n' > "${exit_file}"
  exit 0
fi

echo "$(date --iso-8601=seconds) upstream receipts detected; starting TRU diagnostic"
if python3 "${code_root}/summarize_tru_state_diagnostic.py" "${state_root}" "${figure2_root}"; then
  printf '0\n' > "${exit_file}"
  echo "$(date --iso-8601=seconds) TRU diagnostic completed"
  exit 0
else
  status="$?"
  printf '%s\n' "${status}" > "${exit_file}"
  echo "$(date --iso-8601=seconds) TRU diagnostic failed with exit ${status}" >&2
  exit "${status}"
fi

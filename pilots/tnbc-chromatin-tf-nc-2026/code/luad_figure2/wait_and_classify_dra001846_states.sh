#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 STATE_ROOT FIGURE1_ROOT FIGURE2_ROOT" >&2
  exit 64
fi

state_root="$1"
figure1_root="$2"
figure2_root="$3"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
receipt="${state_root}/audit/dra001846_quant/dra001846_quantification_receipt.json"
audit_dir="${state_root}/audit/state_classifier"
log_dir="${state_root}/logs"

mkdir -p "${audit_dir}" "${log_dir}"
rm -f "${log_dir}/dra001846_state_classifier.exit"

echo "$(date --iso-8601=seconds) waiting for ${receipt}"
while [[ ! -s "${receipt}" ]]; do
  if ! pgrep -f "quantify_dra001846_baseline.py ${state_root}" >/dev/null; then
    echo "$(date --iso-8601=seconds) quantification process is absent before receipt was written" >&2
    printf '1\n' > "${log_dir}/dra001846_state_classifier.exit"
    exit 1
  fi
  sleep 60
done

echo "$(date --iso-8601=seconds) quantification receipt detected; starting frozen classifier"
if Rscript "${script_dir}/classify_luad_states.R" "${state_root}" "${figure1_root}" "${figure2_root}"; then
  printf '0\n' > "${log_dir}/dra001846_state_classifier.exit"
  echo "$(date --iso-8601=seconds) frozen cell-line classifier completed"
  exit 0
else
  status="$?"
  printf '%s\n' "${status}" > "${log_dir}/dra001846_state_classifier.exit"
  echo "$(date --iso-8601=seconds) frozen cell-line classifier failed with exit ${status}" >&2
  exit "${status}"
fi

#!/usr/bin/env bash
set -uo pipefail

# cron has a minimal PATH; Rscript is installed in /usr/local/bin on this host.
export PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 REPOSITORY_ROOT" >&2
  exit 64
fi

repo_root="$(realpath "$1")"
code_root="${repo_root}/pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2"
run_base="${repo_root}/tmp/tnbc-chromatin-tf-nc-2026"
state_root="${run_base}/luad-state-audit"
figure2_root="${run_base}/luad-figure2"
log_root="${state_root}/logs"

mkdir -p "${log_root}"

if [[ ! -s "${state_root}/audit/cloud_keeper_complete.timestamp" ]] \
  && ! pgrep -f "[m]onitor_luad_cloud_jobs.sh ${repo_root}" >/dev/null; then
  nohup "${code_root}/monitor_luad_cloud_jobs.sh" "${repo_root}" \
    >> "${log_root}/cloud_keeper.log" 2>&1 &
fi

if [[ ! -s "${state_root}/audit/tru_state_diagnostic/final_receipt.json" ]] \
  && ! pgrep -f "[w]ait_and_run_tru_state_diagnostic.sh ${state_root}" >/dev/null; then
  nohup "${code_root}/wait_and_run_tru_state_diagnostic.sh" \
    "${state_root}" "${figure2_root}" "${code_root}" \
    >> "${log_root}/tru_state_diagnostic_driver.log" 2>&1 &
fi

"${code_root}/cloud_status_luad.sh" "${repo_root}"

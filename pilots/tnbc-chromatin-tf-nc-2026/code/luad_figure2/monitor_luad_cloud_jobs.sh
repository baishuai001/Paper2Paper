#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 REPOSITORY_ROOT" >&2
  exit 64
fi

repo_root="$(realpath "$1")"
pilot_root="${repo_root}/pilots/tnbc-chromatin-tf-nc-2026"
code_root="${pilot_root}/code/luad_figure2"
run_base="${repo_root}/tmp/tnbc-chromatin-tf-nc-2026"
state_root="${run_base}/luad-state-audit"
figure1_root="${run_base}/luad-figure1"
figure2_root="${run_base}/luad-figure2"
log_root="${state_root}/logs"
lock_path="${state_root}/audit/cloud_keeper.lock"

mkdir -p "${log_root}" "$(dirname "${lock_path}")"
exec 9>"${lock_path}"
if ! flock -n 9; then
  echo "Another LUAD cloud monitor already holds ${lock_path}" >&2
  exit 73
fi

quant_receipt="${state_root}/audit/dra001846_quant/dra001846_quantification_receipt.json"
cell_receipt="${state_root}/audit/state_classifier/state_classifier_receipt.json"
motif_receipt="${figure2_root}/audit/motif_background_sensitivity/final_receipt.json"

log() {
  printf '%s %s\n' "$(date --iso-8601=seconds)" "$*"
}

quant_running() {
  pgrep -f "[q]uantify_dra001846_baseline.py ${state_root}" >/dev/null
}

classifier_waiter_running() {
  pgrep -f "[w]ait_and_classify_dra001846_states.sh ${state_root}" >/dev/null
}

motif_running() {
  pgrep -f "[r]un_motif_background_sensitivity_driver.sh ${figure2_root}" >/dev/null
}

start_quantification() {
  log "Restarting resumable DRA001846 download and Salmon quantification"
  nohup python3 "${code_root}/quantify_dra001846_baseline.py" "${state_root}" \
    --manifest "${state_root}/audit/manifests/dra_rna/ddbj_dra001846_baseline_libraries.tsv" \
    --salmon-index "${state_root}/reference/gencode_v36/salmon_index" \
    --workers 3 --threads-per-sample 6 \
    >> "${log_root}/dra001846_quant_driver.log" 2>&1 &
}

start_classifier_waiter() {
  log "Restarting automatic frozen cell-line classifier waiter"
  nohup "${code_root}/wait_and_classify_dra001846_states.sh" \
    "${state_root}" "${figure1_root}" "${figure2_root}" \
    >> "${log_root}/dra001846_state_classifier_driver.log" 2>&1 &
}

start_motif_sensitivity() {
  log "Restarting resumable B1-B3 motif-background sensitivity driver"
  nohup "${code_root}/run_motif_background_sensitivity_driver.sh" "${figure2_root}" \
    >> "${figure2_root}/logs/motif_background_sensitivity/driver.nohup.log" 2>&1 &
}

log "LUAD cloud monitor started"
while true; do
  if [[ ! -s "${quant_receipt}" ]] && ! quant_running; then
    start_quantification
    sleep 5
  fi

  if [[ ! -s "${cell_receipt}" ]] && ! classifier_waiter_running; then
    start_classifier_waiter
  fi

  if [[ ! -s "${motif_receipt}" ]] && ! motif_running; then
    start_motif_sensitivity
  fi

  if [[ -s "${cell_receipt}" && -s "${motif_receipt}" ]]; then
    log "Cell-line state classifier and B1-B3 motif sensitivity both completed"
    printf '%s\n' "$(date --iso-8601=seconds)" \
      > "${state_root}/audit/cloud_keeper_complete.timestamp"
    exit 0
  fi

  sleep 300
done

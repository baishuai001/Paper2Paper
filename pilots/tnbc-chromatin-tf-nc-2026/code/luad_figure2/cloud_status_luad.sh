#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 REPOSITORY_ROOT" >&2
  exit 64
fi

repo_root="$(realpath "$1")"
run_base="${repo_root}/tmp/tnbc-chromatin-tf-nc-2026"
state_root="${run_base}/luad-state-audit"
figure2_root="${run_base}/luad-figure2"

count_files() {
  local root="$1"
  local name="$2"
  find "$root" -type f -name "$name" 2>/dev/null | wc -l
}

process_count() {
  local pattern="$1"
  pgrep -fc "$pattern" 2>/dev/null || true
}

printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
printf 'quant_driver_processes\t%s\n' "$(process_count '[q]uantify_dra001846_baseline.py')"
printf 'quantified_cell_lines\t%s/19\n' "$(count_files "${state_root}/data/processed/dra001846_quant" quant.sf)"
printf 'quant_receipts\t%s/19\n' "$(find "${state_root}/audit/dra001846_quant" -maxdepth 1 -type f -name '*.json' ! -name 'dra001846_quantification_receipt.json' 2>/dev/null | wc -l)"
printf 'classifier_waiter_processes\t%s\n' "$(process_count '[w]ait_and_classify_dra001846_states.sh')"
printf 'cell_classifier_complete\t%s\n' "$([[ -s "${state_root}/audit/state_classifier/state_classifier_receipt.json" ]] && echo TRUE || echo FALSE)"
printf 'motif_driver_processes\t%s\n' "$(process_count '[r]un_motif_background_sensitivity_driver.sh')"
for variant in B1 B2 B3; do
  for source in JASPAR2024 CIS-BP2.00; do
    root="${figure2_root}/results/motif_background_sensitivity/${variant}/homer/${source}"
    printf '%s_%s_completed_samples\t%s\n' "$variant" "$source" "$(count_files "$root" .complete)"
  done
done
printf 'motif_sensitivity_complete\t%s\n' "$([[ -s "${figure2_root}/audit/motif_background_sensitivity/final_receipt.json" ]] && echo TRUE || echo FALSE)"
printf 'motif_driver_exit\t%s\n' "$(cat "${figure2_root}/logs/motif_background_sensitivity/driver.exit" 2>/dev/null || echo RUNNING_OR_NOT_WRITTEN)"
printf 'cloud_keeper_processes\t%s\n' "$(process_count '[m]onitor_luad_cloud_jobs.sh')"
printf 'TRU_diagnostic_waiter_processes\t%s\n' "$(process_count '[w]ait_and_run_tru_state_diagnostic.sh')"
printf 'TRU_diagnostic_complete\t%s\n' "$([[ -s "${state_root}/audit/tru_state_diagnostic/final_receipt.json" ]] && echo TRUE || echo FALSE)"
printf 'TRU_diagnostic_exit\t%s\n' "$(cat "${state_root}/logs/tru_state_diagnostic.exit" 2>/dev/null || echo RUNNING_OR_NOT_WRITTEN)"
printf 'available_disk\t%s\n' "$(df -h --output=avail "${repo_root}" | tail -n 1 | xargs)"

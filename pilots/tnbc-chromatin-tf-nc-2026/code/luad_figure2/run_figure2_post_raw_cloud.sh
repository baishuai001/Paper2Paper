#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_figure2_post_raw_cloud.sh <figure2_run_root> <project_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
PROJECT_ROOT=$(realpath "$2")
CODE="$PROJECT_ROOT/pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2"
LOG="$RUN_ROOT/logs/figure2_atomic"
mkdir -p "$LOG"

run_step() {
  local label=$1
  shift
  echo "STEP_START $label $(date -u +%FT%TZ)"
  "$@" > >(tee "$LOG/${label}.log") 2>&1
  echo "STEP_COMPLETE $label $(date -u +%FT%TZ)"
}

aggregate_args=("$RUN_ROOT")
if [[ -n "${EARLY_STOP_LOCKED_SYSTEM:-}" ]]; then
  case "$EARLY_STOP_LOCKED_SYSTEM" in
    PDX|cell_line) ;;
    *) echo "Invalid EARLY_STOP_LOCKED_SYSTEM: $EARLY_STOP_LOCKED_SYSTEM" >&2; exit 2 ;;
  esac
  aggregate_args+=(--early-stop-locked-system "$EARLY_STOP_LOCKED_SYSTEM")
fi
run_step 01_aggregate_raw_qc python3 "$CODE/aggregate_raw_atac_qc.py" "${aggregate_args[@]}"
raw_status=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["raw_data_gate"])' \
  "$RUN_ROOT/audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json")
if [[ "$raw_status" != "PASS" ]]; then
  run_step 02_finalize_data_failure python3 "$CODE/finalize_figure2_gate.py" "$RUN_ROOT"
  if [[ "${EARLY_STOP_LOCKED_SYSTEM:-}" == "cell_line" ]]; then
    run_step 03_verify_data_failure python3 "$CODE/verify_figure2_data_failure.py" "$RUN_ROOT"
  fi
  printf 'completed_utc=%s\nverdict=FAIL_DATA\n' "$(date -u +%FT%TZ)" \
    > "$RUN_ROOT/audit/final_gate/figure2_data_failure_run.complete"
  echo "FIGURE2_STOP_FAIL_DATA"
  exit 3
fi

run_step 02_build_atomic_manifest python3 "$CODE/build_figure2_sample_manifest.py" "$RUN_ROOT"
run_step 03_build_supp3_inputs bash "$CODE/build_supp3_accessibility_inputs.sh" "$RUN_ROOT"
run_step 04_compute_supp3 python3 "$CODE/compute_supp3_saturation_correlations.py" "$RUN_ROOT"
run_step 05_compute_promoter_gate python3 "$CODE/compute_promoter_gate_and_annotations.py" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 06_promoter_matched_permutation python3 "$CODE/run_promoter_matched_permutation.py" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 07_prepare_motif_catalog python3 "$CODE/prepare_hc_motif_catalog.py" "$RUN_ROOT"
run_step 08_homer_motif_gate bash "$CODE/run_homer_motif_gate.sh" "$RUN_ROOT"
run_step 09_parse_motif_gate python3 "$CODE/parse_homer_motif_gate.py" "$RUN_ROOT"
run_step 10_ataqv_qc bash "$CODE/run_ataqv_raw_qc.sh" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 10a_tn5_tss_qc python3 "$CODE/compute_tn5_tss_enrichment.py" "$RUN_ROOT"
run_step 10b_parse_ataqv_qc python3 "$CODE/parse_ataqv_metrics.py" "$RUN_ROOT"
run_step 11_plot_atomic Rscript "$CODE/plot_figure2_and_supplementary.R" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 12_finalize_gate python3 "$CODE/finalize_figure2_gate.py" "$RUN_ROOT"
run_step 13_verify_atomic python3 "$CODE/verify_figure2_atomic.py" "$RUN_ROOT" "$PROJECT_ROOT"

printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" > "$RUN_ROOT/audit/final_gate/figure2_atomic_run.complete"
echo "FIGURE2_ATOMIC_GATE_COMPLETE"

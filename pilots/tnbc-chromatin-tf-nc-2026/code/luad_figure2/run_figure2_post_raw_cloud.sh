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

run_step 01_aggregate_raw_qc python3 "$CODE/aggregate_raw_atac_qc.py" "$RUN_ROOT"
run_step 02_build_atomic_manifest python3 "$CODE/build_figure2_sample_manifest.py" "$RUN_ROOT"
run_step 03_build_supp3_inputs bash "$CODE/build_supp3_accessibility_inputs.sh" "$RUN_ROOT"
run_step 04_compute_supp3 python3 "$CODE/compute_supp3_saturation_correlations.py" "$RUN_ROOT"
run_step 05_compute_promoter_gate python3 "$CODE/compute_promoter_gate_and_annotations.py" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 06_prepare_motif_catalog python3 "$CODE/prepare_hc_motif_catalog.py" "$RUN_ROOT"
run_step 07_homer_motif_gate bash "$CODE/run_homer_motif_gate.sh" "$RUN_ROOT"
run_step 08_parse_motif_gate python3 "$CODE/parse_homer_motif_gate.py" "$RUN_ROOT"
run_step 09_ataqv_qc bash "$CODE/run_ataqv_raw_qc.sh" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 09a_tn5_tss_qc python3 "$CODE/compute_tn5_tss_enrichment.py" "$RUN_ROOT"
run_step 09b_parse_ataqv_qc python3 "$CODE/parse_ataqv_metrics.py" "$RUN_ROOT"
run_step 10_plot_atomic Rscript "$CODE/plot_figure2_and_supplementary.R" "$RUN_ROOT" "$PROJECT_ROOT"
run_step 11_finalize_report python3 "$CODE/finalize_figure2_gate.py" "$RUN_ROOT"
run_step 12_verify_atomic python3 "$CODE/verify_figure2_atomic.py" "$RUN_ROOT" "$PROJECT_ROOT"

printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" > "$RUN_ROOT/audit/final_gate/figure2_atomic_run.complete"
echo "FIGURE2_ORIGINAL_STYLE_ANALYSIS_COMPLETE"

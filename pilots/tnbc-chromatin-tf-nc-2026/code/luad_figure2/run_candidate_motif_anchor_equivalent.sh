#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_candidate_motif_anchor_equivalent.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOG_ROOT="$RUN_ROOT/logs/motif_equivalence/candidate_anchor_equivalent"
mkdir -p "$LOG_ROOT"

python3 "$CODE_DIR/build_candidate_motif_database.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_ROOT/01_build_database.log"

for source in JASPAR2024 CIS-BP2.00; do
  MOTIF_DB="$RUN_ROOT/reference/motifs/candidate_luad_hc_and_lineage_controls.${source}.homer" \
  RESULTS_DIR="$RUN_ROOT/results/motif_equivalence/candidate_anchor_homer/${source}" \
  LOG_DIR="$LOG_ROOT/per_sample/${source}" \
  PARALLEL_JOBS=${PARALLEL_JOBS:-2} \
  THREADS_PER_JOB=${THREADS_PER_JOB:-1} \
  bash "$CODE_DIR/run_harmonized_homer.sh" "$RUN_ROOT" \
    2>&1 | tee "$LOG_ROOT/02_homer_${source}.log"
done

python3 "$CODE_DIR/parse_harmonized_motifs.py" "$RUN_ROOT" \
  --homer-results-relative results/motif_equivalence/candidate_anchor_homer/JASPAR2024 \
  --homer-results-relative results/motif_equivalence/candidate_anchor_homer/CIS-BP2.00 \
  --output-relative results/motif_equivalence/candidate_anchor_equivalent \
  --motif-database-relative reference/motifs/candidate_luad_hc_and_lineage_controls.homer \
  --receipt-relative audit/motif_equivalence/candidate_anchor_equivalent_motif_receipt.json \
  --analysis-label "anchor-equivalent candidate motifs with independent JASPAR2024/CIS-BP2.00 correction and per-TF JASPAR preference with CIS-BP fallback" \
  --prefer-jaspar-then-cisbp \
  2>&1 | tee "$LOG_ROOT/03_parse.log"

echo "CANDIDATE_MOTIF_ANCHOR_EQUIVALENT_COMPLETE"

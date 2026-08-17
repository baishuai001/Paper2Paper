#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_candidate_motif_equivalence.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOG_ROOT="$RUN_ROOT/logs/motif_equivalence/candidate_database"
mkdir -p "$LOG_ROOT"

python3 "$CODE_DIR/build_candidate_motif_database.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_ROOT/01_build_database.log"

MOTIF_DB="$RUN_ROOT/reference/motifs/candidate_luad_hc_and_lineage_controls.homer" \
RESULTS_DIR="$RUN_ROOT/results/motif_equivalence/candidate_homer" \
LOG_DIR="$LOG_ROOT/per_sample" \
PARALLEL_JOBS=${PARALLEL_JOBS:-2} \
THREADS_PER_JOB=${THREADS_PER_JOB:-1} \
bash "$CODE_DIR/run_harmonized_homer.sh" "$RUN_ROOT" \
  2>&1 | tee "$LOG_ROOT/02_homer.log"

python3 "$CODE_DIR/parse_harmonized_motifs.py" "$RUN_ROOT" \
  --homer-results-relative results/motif_equivalence/candidate_homer \
  --output-relative results/motif_equivalence/candidate_database \
  --motif-database-relative reference/motifs/candidate_luad_hc_and_lineage_controls.homer \
  --receipt-relative audit/motif_equivalence/candidate_database_motif_receipt.json \
  --analysis-label "combined 47-model candidate sensitivity analysis with joint JASPAR2024/CIS-BP2.00 correction" \
  2>&1 | tee "$LOG_ROOT/03_parse.log"

echo "CANDIDATE_MOTIF_EQUIVALENCE_COMPLETE"

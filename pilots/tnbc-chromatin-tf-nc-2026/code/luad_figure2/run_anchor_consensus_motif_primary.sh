#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_anchor_consensus_motif_primary.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INPUT="$RUN_ROOT/data/processed/motif_primary/anchor_consensus_200bp"
MANIFEST="$INPUT/anchor_consensus_motif_manifest.tsv"
BACKGROUND="$INPUT/lung_accessible_consensus_background.bed"
RESULTS="$RUN_ROOT/results/motif_primary/anchor_consensus_author"
HOMER_RESULTS="$RUN_ROOT/results/motif_primary/anchor_consensus_homer"
LOG_ROOT="$RUN_ROOT/logs/motif_primary/anchor_consensus_author"
mkdir -p "$LOG_ROOT" "$RESULTS" "$HOMER_RESULTS"

python3 "$CODE_DIR/build_anchor_consensus_motif_inputs.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_ROOT/01_build_inputs.log"

for source in JASPAR2024 CIS-BP2.00; do
  motif_db="$RUN_ROOT/reference/motifs/full_reference_${source}.homer"
  [[ -s "$motif_db" ]] || { echo "Missing full motif database: $motif_db" >&2; exit 2; }
  INPUT_DIR="$INPUT" \
  MANIFEST_PATH="$MANIFEST" \
  BACKGROUND_PATH="$BACKGROUND" \
  PREPARSED_DIR="$RUN_ROOT/cache/motif_primary/anchor_consensus/$source" \
  NOT_TESTABLE_STATUS=NOT_TESTABLE_NO_USABLE_CANONICAL_Q0.01_PEAKS \
  MOTIF_DB="$motif_db" \
  RESULTS_DIR="$HOMER_RESULTS/$source" \
  LOG_DIR="$LOG_ROOT/per_sample/$source" \
  PARALLEL_JOBS=${PARALLEL_JOBS:-6} \
  THREADS_PER_JOB=${THREADS_PER_JOB:-2} \
  bash "$CODE_DIR/run_harmonized_homer.sh" "$RUN_ROOT" \
    2>&1 | tee "$LOG_ROOT/02_homer_${source}.log"
done

python3 "$CODE_DIR/parse_harmonized_motifs.py" "$RUN_ROOT" \
  --homer-results-relative results/motif_primary/anchor_consensus_homer/JASPAR2024 \
  --homer-results-relative results/motif_primary/anchor_consensus_homer/CIS-BP2.00 \
  --output-relative results/motif_primary/anchor_consensus_author \
  --motif-database-relative reference/motifs/full_reference_both_databases.homer \
  --receipt-relative audit/motif_primary/anchor_consensus_author_receipt.json \
  --manifest-relative data/processed/motif_primary/anchor_consensus_200bp/anchor_consensus_motif_manifest.tsv \
  --background-relative data/processed/motif_primary/anchor_consensus_200bp/lung_accessible_consensus_background.bed \
  --output-prefix anchor_consensus \
  --analysis-label "primary TNBC-author-analog LUAD motif analysis: all usable 200-bp peaks, one shared lung accessible consensus background, full JASPAR2024 and CIS-BP2.00 run separately, per-TF JASPAR preference with CIS-BP fallback" \
  --prefer-jaspar-then-cisbp \
  --not-testable-status NOT_TESTABLE_NO_USABLE_CANONICAL_Q0.01_PEAKS \
  2>&1 | tee "$LOG_ROOT/03_parse.log"

python3 "$CODE_DIR/validate_anchor_consensus_motif_primary.py" "$RUN_ROOT" \
  2>&1 | tee "$LOG_ROOT/04_validate.log"

echo "ANCHOR_CONSENSUS_MOTIF_PRIMARY_COMPLETE"

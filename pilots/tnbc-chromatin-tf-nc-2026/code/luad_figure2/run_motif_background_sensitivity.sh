#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_motif_background_sensitivity.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
B0_INPUT="$RUN_ROOT/data/processed/motif_primary/anchor_consensus_200bp"
MANIFEST="$B0_INPUT/anchor_consensus_motif_manifest.tsv"
LOG_ROOT="$RUN_ROOT/logs/motif_background_sensitivity"
mkdir -p "$LOG_ROOT"

python3 "$CODE_DIR/build_motif_background_sensitivity.py" "$RUN_ROOT" \
  > "$LOG_ROOT/01_build_backgrounds.log" 2>&1

for variant in B1 B2 B3; do
  BACKGROUND="$RUN_ROOT/data/processed/motif_background_sensitivity/$variant/${variant}_shared_accessible_background_200bp.bed"
  for source in JASPAR2024 CIS-BP2.00; do
    motif_db="$RUN_ROOT/reference/motifs/full_reference_${source}.homer"
    [[ -s "$motif_db" ]] || { echo "Missing motif database: $motif_db" >&2; exit 2; }
    INPUT_DIR="$B0_INPUT" \
    MANIFEST_PATH="$MANIFEST" \
    BACKGROUND_PATH="$BACKGROUND" \
    PREPARSED_DIR="$RUN_ROOT/cache/motif_background_sensitivity/$variant/$source" \
    NOT_TESTABLE_STATUS=NOT_TESTABLE_NO_USABLE_CANONICAL_Q0.01_PEAKS \
    MOTIF_DB="$motif_db" \
    RESULTS_DIR="$RUN_ROOT/results/motif_background_sensitivity/$variant/homer/$source" \
    LOG_DIR="$LOG_ROOT/$variant/per_sample/$source" \
    PARALLEL_JOBS=${PARALLEL_JOBS:-6} \
    THREADS_PER_JOB=${THREADS_PER_JOB:-2} \
    bash "$CODE_DIR/run_harmonized_homer.sh" "$RUN_ROOT" \
      > "$LOG_ROOT/${variant}_homer_${source}.log" 2>&1
  done

  python3 "$CODE_DIR/parse_harmonized_motifs.py" "$RUN_ROOT" \
    --homer-results-relative "results/motif_background_sensitivity/$variant/homer/JASPAR2024" \
    --homer-results-relative "results/motif_background_sensitivity/$variant/homer/CIS-BP2.00" \
    --output-relative "results/motif_background_sensitivity/$variant/parsed" \
    --motif-database-relative reference/motifs/full_reference_both_databases.homer \
    --receipt-relative "audit/motif_background_sensitivity/${variant}_motif_receipt.json" \
    --manifest-relative data/processed/motif_primary/anchor_consensus_200bp/anchor_consensus_motif_manifest.tsv \
    --background-relative "data/processed/motif_background_sensitivity/$variant/${variant}_shared_accessible_background_200bp.bed" \
    --output-prefix "$variant" \
    --analysis-label "prespecified $variant shared-background sensitivity; immutable B0 targets" \
    --prefer-jaspar-then-cisbp \
    --not-testable-status NOT_TESTABLE_NO_USABLE_CANONICAL_Q0.01_PEAKS \
    > "$LOG_ROOT/${variant}_parse.log" 2>&1
done

python3 "$CODE_DIR/summarize_motif_background_sensitivity.py" "$RUN_ROOT" \
  > "$LOG_ROOT/99_summarize.log" 2>&1

echo "MOTIF_BACKGROUND_SENSITIVITY_COMPLETE"


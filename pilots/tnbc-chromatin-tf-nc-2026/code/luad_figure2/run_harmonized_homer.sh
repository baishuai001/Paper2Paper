#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_harmonized_homer.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
HOMER="$RUN_ROOT/tools/homer"
REF="$RUN_ROOT/reference"
INPUT="$RUN_ROOT/data/processed/motif_equivalence/harmonized_200bp"
MANIFEST="$INPUT/harmonized_motif_manifest.tsv"
BACKGROUND="$INPUT/LUSC_accessible_common_GC_matched_background.bed"
MOTIF_DB=${MOTIF_DB:-"$REF/motifs/full_reference_both_databases.homer"}
RESULTS=${RESULTS_DIR:-"$RUN_ROOT/results/motif_equivalence/harmonized_homer"}
LOG=${LOG_DIR:-"$RUN_ROOT/logs/motif_equivalence"}
PARALLEL_JOBS=${PARALLEL_JOBS:-4}
THREADS_PER_JOB=${THREADS_PER_JOB:-2}

mkdir -p "$RESULTS" "$LOG"
for required in "$HOMER/bin/findMotifsGenome.pl" "$REF/hg38.fa" "$MANIFEST" "$BACKGROUND" "$MOTIF_DB"; do
  [[ -s "$required" ]] || { echo "Missing required input: $required" >&2; exit 2; }
done
export PATH="$HOMER/bin:$PATH"

manifest_sha=$(sha256sum "$MANIFEST" | cut -d' ' -f1)
background_sha=$(sha256sum "$BACKGROUND" | cut -d' ' -f1)
motif_sha=$(sha256sum "$MOTIF_DB" | cut -d' ' -f1)

run_one() {
  local system=$1 sample=$2 slug=$3 matchable=$4 target=$5 target_sha=$6
  local out="$RESULTS/$system/$slug"
  mkdir -p "$out"
  if [[ "$matchable" != "TRUE" ]]; then
    printf 'status=NOT_TESTABLE_PEAK_COUNT_MATCH\nmanifest_sha256=%s\n' "$manifest_sha" > "$out/.not_testable"
    printf 'status=NOT_TESTABLE_PEAK_COUNT_MATCH\nmanifest_sha256=%s\nbackground_sha256=%s\nmotif_sha256=%s\n' \
      "$manifest_sha" "$background_sha" "$motif_sha" > "$out/.complete"
    echo "HARMONIZED_HOMER_NOT_TESTABLE $system $sample"
    return 0
  fi
  if [[ -s "$out/.complete" && -s "$out/knownResults.txt" ]] \
    && grep -qx "manifest_sha256=$manifest_sha" "$out/.complete" \
    && grep -qx "background_sha256=$background_sha" "$out/.complete" \
    && grep -qx "motif_sha256=$motif_sha" "$out/.complete" \
    && grep -qx "target_sha256=$target_sha" "$out/.complete"; then
    echo "HARMONIZED_HOMER_SKIP_COMPLETE $system $sample"
    return 0
  fi
  rm -f "$out/.complete" "$out/.not_testable"
  findMotifsGenome.pl "$target" "$REF/hg38.fa" "$out" \
    -bg "$BACKGROUND" -size 200 -nomotif -mknown "$MOTIF_DB" -p "$THREADS_PER_JOB" \
    > "$LOG/${system}.${slug}.log" 2>&1
  [[ -s "$out/knownResults.txt" ]] || { echo "HOMER failed: $system/$sample" >&2; return 1; }
  printf 'status=COMPLETE\nmanifest_sha256=%s\nbackground_sha256=%s\nmotif_sha256=%s\ntarget_sha256=%s\n' \
    "$manifest_sha" "$background_sha" "$motif_sha" "$target_sha" > "$out/.complete"
  echo "HARMONIZED_HOMER_COMPLETE $system $sample"
}

running=0
failed=0
while IFS=$'\t' read -r system sample slug source source_rows invalid_removed blacklist_removed non_acgt_removed harmonized_count floor matchable reason matched_count original_gc_mean original_gc_sd matched_gc_mean matched_gc_sd target target_sha; do
  run_one "$system" "$sample" "$slug" "$matchable" "$target" "$target_sha" &
  running=$((running + 1))
  if (( running >= PARALLEL_JOBS )); then
    if ! wait -n; then failed=$((failed + 1)); fi
    running=$((running - 1))
  fi
done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')
while (( running > 0 )); do
  if ! wait -n; then failed=$((failed + 1)); fi
  running=$((running - 1))
done
if (( failed > 0 )); then
  echo "HARMONIZED_HOMER_FAILURES=$failed" >&2
  exit 1
fi
echo "HARMONIZED_HOMER_ALL_COMPLETE"

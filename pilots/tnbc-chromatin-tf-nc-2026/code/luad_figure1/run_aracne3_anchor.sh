#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: run_aracne3_anchor.sh ARACNE_REPO EXPRESSION_TSV REGULATORS OUTPUT_DIR RUN_ID THREADS" >&2
  exit 2
fi

repo=$1
expression=$2
regulators=$3
output=$4
run_id=$5
threads=$6
expected_commit=3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14
seed=2025

actual_commit=$(git -C "$repo" rev-parse HEAD)
if [[ "$actual_commit" != "$expected_commit" ]]; then
  echo "ARACNe3 commit $actual_commit does not match $expected_commit" >&2
  exit 3
fi
if [[ ! -s "$expression" || ! -s "$regulators" ]]; then
  echo "Missing expression matrix or regulator list" >&2
  exit 4
fi

header_fields=$(awk -F'\t' 'NR==1 {print NF; exit}' "$expression")
row_fields=$(awk -F'\t' 'NR==2 {print NF; exit}' "$expression")
if [[ "$header_fields" -ne "$row_fields" ]]; then
  echo "ARACNe3 input has $header_fields header fields and $row_fields row fields" >&2
  exit 5
fi

binary="$repo/build-system/ARACNe3_app_release"
if [[ ! -x "$binary" ]]; then
  mkdir -p "$repo/build-system"
  g++ -O3 -std=c++17 -fopenmp -I"$repo/include/ARACNe3" \
    "$repo/src/app/cmdline_parser.cpp" "$repo/src/app/stopwatch.cpp" \
    "$repo/src/app/io.cpp" "$repo/src/app/algorithms.cpp" \
    "$repo/src/app/apmi_nullmodel.cpp" "$repo/src/app/subnet_operations.cpp" \
    "$repo/src/app/ARACNe3.cpp" -lstdc++fs -o "$binary"
fi

mkdir -p "$output"
subnet="$output/subnets/subnet1_${run_id}.tsv"
log="$output/log_${run_id}.txt"
if [[ -s "$subnet" ]] && grep -q 'SUCCESS!' "$log" 2>/dev/null; then
  echo "ARACNe3 output already completed; preserving $subnet"
  exit 0
fi

# The TNBC capsule invokes ARACNe3 with only -e/-r/-o, hence one subnet,
# 1-exp(-1) subsampling, alpha=0.05, FDR pruning and MaxEnt/DPI pruning.  We
# expose only the thread count and freeze the otherwise wall-clock seed.
"$binary" \
  -e "$expression" \
  -r "$regulators" \
  -o "$output" \
  --seed "$seed" \
  --threads "$threads" \
  --runid "$run_id"

grep -q 'SUCCESS!' "$log"
test -s "$subnet"
sha256sum "$binary" "$expression" "$regulators" "$subnet" "$log" \
  > "$output/aracne3_sha256.tsv"

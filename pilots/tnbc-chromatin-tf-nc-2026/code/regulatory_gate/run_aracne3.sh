#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: run_aracne3.sh ARACNE_REPO EXPRESSION_TSV REGULATORS OUTPUT_DIR THREADS SUBNETWORKS" >&2
  exit 2
fi

repo=$1
expression=$2
regulators=$3
output=$4
threads=$5
subnetworks=$6
expected_commit=3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14
seed=1729

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
network="$output/consolidated-net_crc.tsv"
if [[ -s "$network" ]] && grep -q 'SUCCESS!' "$output/log_crc.txt" 2>/dev/null; then
  echo "ARACNe3 completed output already exists; preserving it."
  exit 0
fi

"$binary" \
  -e "$expression" \
  -r "$regulators" \
  -o "$output" \
  -x "$subnetworks" \
  --subsample 0.63212 \
  --alpha 0.05 \
  --seed "$seed" \
  --threads "$threads" \
  --runid crc

grep -q 'SUCCESS!' "$output/log_crc.txt"
actual_subnets=$(find "$output/subnets" -maxdepth 1 -type f -name 'subnet*_crc.tsv' | wc -l)
if [[ "$actual_subnets" -ne "$subnetworks" ]]; then
  echo "Expected $subnetworks subnetworks, found $actual_subnets" >&2
  exit 6
fi
sha256sum "$binary" "$network" "$output/log_crc.txt"

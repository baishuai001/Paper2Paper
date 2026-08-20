#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 FIGURE1_ROOT FIGURE2_ROOT OUTPUT_ROOT" >&2
  exit 2
fi

fig1=$1
fig2=$2
out=$3
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out/logs" "$out/audit" "$out/results/tables" "$out/results/figures"

deadline=$(( $(date +%s) + 86400 ))
while [[ ! -s "$fig1/audit/figure1_independent_verification.json" ]]; do
  if (( $(date +%s) >= deadline )); then
    echo "Timed out waiting for corrected Figure 1 verification" >&2
    exit 9
  fi
  sleep 30
done
grep -q '"verification_status": "passed"' "$fig1/audit/figure1_independent_verification.json"

required=(
  "$fig1/audit/figure1_independent_verification.json"
  "$fig1/results/tcga/TCGA_regulon.rds"
  "$fig1/results/tcga/TCGA_viper_activity.tsv.gz"
  "$fig2/results/promoter_gate/tf_promoter_gate_summary.tsv"
  "$fig2/results/motif_threshold_sensitivity/motif_threshold_sensitivity_triple_system_tfs.tsv"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Required Figure 3 input missing: $path" >&2; exit 10; }
done

sha256sum "${required[@]}" > "$out/audit/figure3_input_sha256.tsv"
Rscript "$script_dir/run_luad_figure3.R" "$fig1" "$fig2" "$out" \
  2>&1 | tee "$out/logs/run_luad_figure3.log"
python3 "$script_dir/assemble_figure3.py" "$out" \
  2>&1 | tee "$out/logs/assemble_figure3.log"
python3 "$script_dir/verify_luad_figure3.py" "$out" \
  2>&1 | tee "$out/logs/verify_luad_figure3.log"
find "$out/results" "$out/audit" -type f -print0 | sort -z | xargs -0 sha256sum \
  > "$out/audit/figure3_output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" > "$out/audit/run_completion.tsv"

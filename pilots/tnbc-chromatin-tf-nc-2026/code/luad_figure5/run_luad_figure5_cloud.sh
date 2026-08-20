#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 4 ]]; then echo "Usage: $0 FIG1_ROOT FIG4_ROOT PDXE_TSV_DIR OUTPUT_ROOT" >&2; exit 2; fi
fig1=$1; fig4=$2; pdxe=$3; out=$4
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out/logs" "$out/audit" "$out/results/tables" "$out/results/figures"
export R_LIBS_USER="${R_LIBS_USER:-$PWD/$out/r-lib}"
required=(
  "$fig1/audit/figure1_independent_verification.json"
  "$fig4/audit/figure4_independent_verification.json"
  "$out/audit/figure5_cellline_input_receipt.json"
  "$out/audit/pdxe_lung_classifier_receipt.json"
  "$out/data/processed/figure5_cellline_inputs.rds"
  "$out/results/pdxe_projection/PDXE_V2_viper_activity.tsv.gz"
)
for path in "${required[@]}"; do [[ -s "$path" ]] || { echo "Missing required input: $path" >&2; exit 10; }; done
grep -q '"verification_status": "passed"' "$fig1/audit/figure1_independent_verification.json"
grep -q '"verification_status": "passed"' "$fig4/audit/figure4_independent_verification.json"
sha256sum "${required[@]}" > "$out/audit/figure5_input_sha256.tsv"
Rscript "$script_dir/analyze_luad_figure5.R" "$fig4" "$pdxe" "$out" 2>&1 | tee "$out/logs/analyze_luad_figure5.log"
python3 "$script_dir/assemble_figure5.py" "$out" 2>&1 | tee "$out/logs/assemble_figure5.log"
python3 "$script_dir/verify_luad_figure5.py" "$out" 2>&1 | tee "$out/logs/verify_luad_figure5.log"
find "$out/results" "$out/audit" -type f -print0 | sort -z | xargs -0 sha256sum > "$out/audit/figure5_output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" > "$out/audit/run_completion.tsv"


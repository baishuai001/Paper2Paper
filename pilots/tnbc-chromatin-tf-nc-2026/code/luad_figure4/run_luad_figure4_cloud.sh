#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 FIGURE1_ROOT FIGURE2_ROOT FIGURE3_ROOT OUTPUT_ROOT THREADS" >&2
  exit 2
fi
fig1=$1
fig2=$2
fig3=$3
out=$4
threads=$5
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out/logs" "$out/audit" "$out/results/tables" "$out/results/figures"

required=(
  "$fig1/audit/figure1_complete/figure1_complete_independent_verification.json"
  "$fig1/results/gse41271/GSE41271_viper_activity.tsv.gz"
  "$fig1/data/processed/gse41271_manifest.tsv"
  "$fig1/data/processed/annotations/tcga_figure1_annotations.tsv"
  "$fig2/results/promoter_gate/tf_promoter_gate_summary.tsv"
  "$fig3/audit/figure3_independent_verification.json"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Required Figure 4 input missing: $path" >&2; exit 10; }
done
grep -q '"verification_status": "passed"' "${required[0]}"
grep -q '"verification_status": "passed"' "${required[5]}"
sha256sum "${required[@]}" > "$out/audit/figure4_input_sha256.tsv"

python3 "$script_dir/fetch_tcga_luad_smoking.py" "$out" \
  2>&1 | tee "$out/logs/fetch_tcga_luad_smoking.log"

mkdir -p "$out/results/gse41271_tcga_projection"
Rscript "$script_dir/../luad_figure1/project_viper.R" \
  "$fig1/results/tcga/TCGA_regulon.rds" \
  "$fig1/data/processed/gse41271_aracne_expression.tsv" \
  "$fig1/data/processed/gse41271_manifest.tsv" \
  GSE41271_TCGA_PROJECTED "$out/results/gse41271_tcga_projection" "$threads" \
  2>&1 | tee "$out/logs/project_tcga_regulon_to_gse41271.log"

Rscript "$script_dir/run_luad_figure4.R" "$fig1" "$fig2" "$out" \
  2>&1 | tee "$out/logs/run_luad_figure4.log"
Rscript "$script_dir/run_figure4_permutations.R" "$out" 5000 "$threads" \
  2>&1 | tee "$out/logs/run_figure4_permutations.log"
python3 "$script_dir/assemble_figure4.py" "$out" \
  2>&1 | tee "$out/logs/assemble_figure4.log"
python3 "$script_dir/verify_luad_figure4.py" "$out" \
  2>&1 | tee "$out/logs/verify_luad_figure4.log"
find "$out/results" "$out/audit" -type f -print0 | sort -z | xargs -0 sha256sum \
  > "$out/audit/figure4_output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" > "$out/audit/run_completion.tsv"

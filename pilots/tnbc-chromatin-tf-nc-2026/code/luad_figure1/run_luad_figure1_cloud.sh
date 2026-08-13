#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: run_luad_figure1_cloud.sh CLOUD_RUN_ROOT ARACNE3_REPO THREADS" >&2
  exit 2
fi

root=$1
aracne_repo=$2
threads=$3
code="$root/code"
processed="$root/data/processed"
results="$root/results"
audit="$root/audit"
logs="$root/logs"
mkdir -p "$results/tcga" "$results/gse81089" "$results/pdmr" "$results/depmap" "$audit" "$logs"

if [[ ! -s "$processed/pdmr_tpm_symbols.tsv" ]]; then
  echo "PDMR matrix is incomplete" >&2
  exit 10
fi
if [[ ! -s "$root/data/raw/depmap/depmap_TPM_22Q2.tsv.gz" ]]; then
  echo "DepMap matrix is incomplete" >&2
  exit 11
fi
if [[ ! -s "$root/data/raw/tcga/TCGA-LUAD.star_tpm.tsv.gz" || ! -s "$root/data/raw/tcga/TCGA-LUSC.star_tpm.tsv.gz" ]]; then
  echo "TCGA TPM matrices are incomplete" >&2
  exit 12
fi

required_processed=(
  tcga_aracne_expression.tsv tcga_counts_symbols.tsv.gz tcga_manifest.tsv
  gse81089_aracne_expression.tsv gse81089_counts_symbols.tsv.gz gse81089_manifest.tsv
  pdmr_viper_expression.tsv.gz pdmr_manifest.tsv
  depmap_viper_expression.tsv.gz depmap_manifest.tsv aracne_regulators.txt
)
processed_checkpoint_ok=true
for required_file in "${required_processed[@]}"; do
  [[ -s "$processed/$required_file" ]] || processed_checkpoint_ok=false
done
if [[ "$processed_checkpoint_ok" == true && -s "$audit/expression_preparation_receipt.json" && -s "$audit/expression_scale_receipt.tsv" ]]; then
  echo "Preserving verified processed-expression checkpoint" | tee "$logs/prepare_expression.log"
else
  Rscript "$code/prepare_expression.R" "$root" 2>&1 | tee "$logs/prepare_expression.log"
fi

# Both networks use exactly the anchor ARACNe3 settings. They are independent
# cohorts, so running them concurrently does not alter either result.
bash "$code/run_aracne3_anchor.sh" \
  "$aracne_repo" "$processed/tcga_aracne_expression.tsv" "$processed/aracne_regulators.txt" \
  "$results/tcga/aracne3" TCGA_LUAD_LUSC "$threads" > "$logs/aracne_tcga.log" 2>&1 &
tcga_aracne_pid=$!
bash "$code/run_aracne3_anchor.sh" \
  "$aracne_repo" "$processed/gse81089_aracne_expression.tsv" "$processed/aracne_regulators.txt" \
  "$results/gse81089/aracne3" GSE81089_LUAD_LUSC "$threads" > "$logs/aracne_gse81089.log" 2>&1 &
gse_aracne_pid=$!
wait "$tcga_aracne_pid"
wait "$gse_aracne_pid"

Rscript "$code/run_viper_anchor.R" \
  "$results/tcga/aracne3/subnets/subnet1_TCGA_LUAD_LUSC.tsv" \
  "$processed/tcga_aracne_expression.tsv" "$processed/tcga_manifest.tsv" \
  TCGA "$results/tcga" discover "$threads" 1000 2>&1 | tee "$logs/viper_tcga.log"
Rscript "$code/run_viper_anchor.R" \
  "$results/gse81089/aracne3/subnets/subnet1_GSE81089_LUAD_LUSC.tsv" \
  "$processed/gse81089_aracne_expression.tsv" "$processed/gse81089_manifest.tsv" \
  GSE81089 "$results/gse81089" discover "$threads" 1000 2>&1 | tee "$logs/viper_gse81089.log"

Rscript "$code/project_viper.R" \
  "$results/tcga/TCGA_regulon.rds" "$processed/pdmr_viper_expression.tsv.gz" \
  "$processed/pdmr_manifest.tsv" PDMR "$results/pdmr" "$threads" 2>&1 | tee "$logs/viper_pdmr.log"
Rscript "$code/project_viper.R" \
  "$results/tcga/TCGA_regulon.rds" "$processed/depmap_viper_expression.tsv.gz" \
  "$processed/depmap_manifest.tsv" DEPMAP "$results/depmap" "$threads" 2>&1 | tee "$logs/viper_depmap.log"

Rscript "$code/select_and_plot_figure1.R" "$root" 2>&1 | tee "$logs/select_and_plot.log"
python3 "$code/verify_figure1.py" "$root" 2>&1 | tee "$logs/independent_verification.log"

find "$results" "$audit" -type f -print0 | sort -z | xargs -0 sha256sum > "$audit/luad_figure1_output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" > "$audit/run_completion.tsv"
printf 'aracne3_commit\t%s\n' "$(git -C "$aracne_repo" rev-parse HEAD)" >> "$audit/run_completion.tsv"
printf 'threads\t%s\n' "$threads" >> "$audit/run_completion.tsv"

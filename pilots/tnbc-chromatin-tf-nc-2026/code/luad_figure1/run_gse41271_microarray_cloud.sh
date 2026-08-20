#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: run_gse41271_microarray_cloud.sh CLOUD_RUN_ROOT ARACNE3_REPO THREADS" >&2
  exit 2
fi

root=$1
aracne_repo=$2
threads=$3
code="$root/code"
raw="$root/data/raw/gse41271"
processed="$root/data/processed"
results="$root/results"
logs="$root/logs"
audit="$root/audit/gse41271"
mkdir -p "$raw" "$processed" "$results/gse41271" "$logs" "$audit"

download() {
  local url=$1
  local output=$2
  if [[ -s "$output" ]]; then
    echo "Preserving $output"
    return
  fi
  curl -L --fail --retry 5 --retry-delay 5 -o "$output.part" "$url"
  mv "$output.part" "$output"
}

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41271/soft/GSE41271_family.soft.gz" \
  "$raw/GSE41271_family.soft.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41271/matrix/GSE41271_series_matrix.txt.gz" \
  "$raw/GSE41271_series_matrix.txt.gz"

Rscript "$code/prepare_gse41271_microarray.R" "$root" \
  2>&1 | tee "$logs/prepare_gse41271.log"

bash "$code/run_aracne3_anchor.sh" \
  "$aracne_repo" "$processed/gse41271_aracne_expression.tsv" \
  "$processed/gse41271_aracne_regulators.txt" \
  "$results/gse41271/aracne3" GSE41271_LUAD_LUSC "$threads" \
  > "$logs/aracne_gse41271.log" 2>&1

Rscript "$code/run_viper_anchor.R" \
  "$results/gse41271/aracne3/subnets/subnet1_GSE41271_LUAD_LUSC.tsv" \
  "$processed/gse41271_aracne_expression.tsv" \
  "$processed/gse41271_manifest.tsv" \
  GSE41271 "$results/gse41271" discover "$threads" 1000 \
  2>&1 | tee "$logs/viper_gse41271.log"

# The cross-platform summary also uses the corrected GSE81089 tables. Wait for
# the independent Figure 1 revision to finish instead of reading a partial run.
deadline=$(( $(date +%s) + 86400 ))
while [[ ! -s "$root/audit/figure1_independent_verification.json" ]] || \
      ! grep -q '"verification_status": "passed"' "$root/audit/figure1_independent_verification.json"; do
  if (( $(date +%s) >= deadline )); then
    echo "Timed out waiting for corrected Figure 1 verification" >&2
    exit 20
  fi
  sleep 60
done

Rscript "$code/analyze_gse41271_microarray.R" "$root" \
  2>&1 | tee "$logs/analyze_gse41271.log"
python3 "$code/verify_gse41271_microarray.py" "$root" \
  2>&1 | tee "$logs/verify_gse41271.log"

find "$results/gse41271" "$audit" -type f -print0 | sort -z | xargs -0 sha256sum \
  > "$audit/gse41271_output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" > "$audit/run_completion.tsv"
printf 'aracne3_commit\t%s\n' "$(git -C "$aracne_repo" rev-parse HEAD)" >> "$audit/run_completion.tsv"
printf 'threads\t%s\n' "$threads" >> "$audit/run_completion.tsv"


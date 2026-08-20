#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 CLOUD_RUN_ROOT" >&2
  exit 2
fi
root=$1
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$root/logs" "$root/audit/figure1_complete"

deadline=$(( $(date +%s) + 86400 ))
while [[ ! -s "$root/audit/figure1_independent_verification.json" ]] || \
      [[ ! -s "$root/audit/gse41271/gse41271_independent_verification.json" ]]; do
  if (( $(date +%s) >= deadline )); then
    echo "Timed out waiting for corrected Figure 1 base and GSE41271 verification" >&2
    exit 20
  fi
  sleep 30
done
grep -q '"verification_status": "passed"' "$root/audit/figure1_independent_verification.json"
grep -q '"verification_status": "passed"' "$root/audit/gse41271/gse41271_independent_verification.json"

bash "$script_dir/run_figure1_annotations_supp_cloud.sh" "$root" \
  2>&1 | tee "$root/logs/finalize_annotations_supplementary.log"
Rscript "$script_dir/render_figure1_revision.R" "$root" \
  2>&1 | tee "$root/logs/render_figure1_revision.log"
python3 "$script_dir/assemble_figure1_complete.py" "$root" \
  2>&1 | tee "$root/logs/assemble_figure1_complete.log"
python3 "$script_dir/verify_figure1_complete.py" "$root" \
  2>&1 | tee "$root/logs/verify_figure1_complete.log"

find "$root/results/figures" "$root/results/tables" "$root/audit/figure1_complete" \
  -type f -print0 | sort -z | xargs -0 sha256sum \
  > "$root/audit/figure1_complete/output_sha256.tsv"
printf 'completed_at_utc\t%s\n' "$(date -u +%FT%TZ)" \
  > "$root/audit/figure1_complete/run_completion.tsv"


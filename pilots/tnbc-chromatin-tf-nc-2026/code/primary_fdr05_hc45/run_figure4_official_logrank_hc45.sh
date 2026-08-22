#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 OFFICIAL_CODE_ROOT FIGURE1_ROOT FIGURE2_ROOT FIGURE4_SOURCE_ROOT HC45_AUDIT_TSV OUTPUT_ROOT R_LIBRARY" >&2
  exit 64
fi

official_code_root=$(realpath "$1")
figure1_root=$(realpath "$2")
figure2_root=$(realpath "$3")
figure4_source=$(realpath "$4")
hc45_audit=$(realpath "$5")
output_root=$6
r_library=$(realpath "$7")
script_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fig45_port=$(realpath "$script_root/../official_port/fig45")

if [[ -e "$output_root" ]]; then
  echo "Refusing to overwrite audited Figure 4 run: $output_root" >&2
  exit 73
fi
mkdir -p "$output_root"/{analysis,generated,logs,results}
output_root=$(realpath "$output_root")
export R_LIBS_USER="$r_library"

Rscript "$script_root/run_figure4_nominal_logrank_hc45.R" \
  "$figure1_root" "$figure2_root" "$figure4_source" "$hc45_audit" \
  "$output_root/analysis" \
  >"$output_root/logs/01_analysis.log" 2>&1

python3 "$fig45_port/materialize_official_blocks.py" \
  --official-code-root "$official_code_root" \
  --output-dir "$output_root/generated" \
  >"$output_root/logs/02_materialize_official_blocks.log" 2>&1

Rscript "$script_root/render_figure4_official_logrank.R" \
  "$output_root/generated" "$output_root/analysis" "$output_root/results" \
  >"$output_root/logs/03_render.log" 2>&1

pdf_count=$(find "$output_root/results/atomic" -maxdepth 1 -type f -name '*.pdf' | wc -l)
if [[ "$pdf_count" -lt 20 ]]; then
  echo "Expected at least 20 official Figure 4 PDF atoms; found $pdf_count" >&2
  exit 1
fi
printf 'status\tvalue\nrun_status\tPASS\npdf_atoms\t%s\n' "$pdf_count" \
  >"$output_root/results/audit/run_completion.tsv"
echo "Official-code Figure 4 HC45/log-rank run complete: $output_root"

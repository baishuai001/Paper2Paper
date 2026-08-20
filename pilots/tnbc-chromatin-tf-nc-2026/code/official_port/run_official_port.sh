#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: run_official_port.sh FIGURE1_ROOT FIGURE2_ROOT PUBLICATION_ROOT WORK_ROOT" >&2
  exit 2
fi

figure1_root="$1"
figure2_root="$2"
publication_root="$3"
work_root="$4"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
strict_adapter_root="${work_root}/adapter_strict_replay"
final_adapter_root="${work_root}/adapter_final_corrected"
strict_runtime_root="${work_root}/runtime_strict_replay"
final_runtime_root="${work_root}/runtime_final_corrected"
secondary_adapter_root="${final_adapter_root}/secondary_validation/GSE81089"
secondary_runtime_root="${work_root}/runtime_secondary_GSE81089_corrected"

mkdir -p "${work_root}"

Rscript "${script_dir}/adapt_luad_to_official_schema.R" \
  "${figure1_root}" \
  "${figure2_root}" \
  "${publication_root}" \
  "${strict_adapter_root}"

mkdir -p "${final_adapter_root}"
cp -a "${strict_adapter_root}/." "${final_adapter_root}/"

python3 "${script_dir}/materialize_official_runtime.py" \
  "${strict_adapter_root}" \
  "${strict_runtime_root}" \
  --official-port-root "${script_dir}" \
  --validation-label GSE41271

Rscript "${script_dir}/run_figure1_official_blocks.R" \
  "${strict_runtime_root}" \
  "${script_dir}" \
  PRIMARY_METABRIC_ANALOGUE \
  GSE41271 \
  CAPSULE_BUG_RETAINED

python3 "${script_dir}/materialize_official_runtime.py" \
  "${final_adapter_root}" \
  "${final_runtime_root}" \
  --official-port-root "${script_dir}" \
  --validation-label GSE41271

python3 "${script_dir}/apply_supplementary2f_anchor_bugfix.py" \
  "${final_runtime_root}"

python3 "${script_dir}/apply_supplementary4c_empty_set_bugfix.py" \
  "${final_runtime_root}"

Rscript "${script_dir}/run_figure1_official_blocks.R" \
  "${final_runtime_root}" \
  "${script_dir}" \
  PRIMARY_METABRIC_ANALOGUE \
  GSE41271 \
  ANCHOR_BUGFIX

Rscript "${script_dir}/run_figure2_official_scripts.R" \
  "${final_runtime_root}" \
  "${script_dir}"

python3 "${script_dir}/materialize_official_runtime.py" \
  "${secondary_adapter_root}" \
  "${secondary_runtime_root}" \
  --official-port-root "${script_dir}" \
  --validation-label GSE81089

python3 "${script_dir}/apply_supplementary2f_anchor_bugfix.py" \
  "${secondary_runtime_root}"

Rscript "${script_dir}/run_figure1_official_blocks.R" \
  "${secondary_runtime_root}" \
  "${script_dir}" \
  SECONDARY_INDEPENDENT_VALIDATION \
  GSE81089 \
  ANCHOR_BUGFIX

echo "Strict replay complete: ${strict_runtime_root}"
echo "Final corrected Figure 1-2 port complete: ${final_runtime_root}"
echo "Secondary GSE81089 official-block replay complete: ${secondary_runtime_root}"

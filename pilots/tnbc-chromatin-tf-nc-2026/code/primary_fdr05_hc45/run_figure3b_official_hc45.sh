#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 OFFICIAL_CODE_ROOT FIGURE2_ROOT HC45_EDGES_TSV OUTPUT_ROOT R_LIBRARY" >&2
  exit 64
fi

official_code_root=$(realpath "$1")
figure2_root=$(realpath "$2")
edge_path=$(realpath "$3")
output_root=$4
r_library=$(realpath "$5")
script_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fig3_port=$(realpath "$script_root/../official_port/fig3")

if [[ -e "$output_root" ]]; then
  echo "Refusing to overwrite audited Figure 3B run: $output_root" >&2
  exit 73
fi
mkdir -p "$output_root"/{audit,data,logs,results/visuals/Figure3,work}
output_root=$(realpath "$output_root")
export R_LIBS_USER="$r_library"
export PAPER2PAPER_OFFICIAL_PORT_WORK="$output_root/work"

python3 "$fig3_port/materialize_official_scripts.py" \
  --official-code-root "$official_code_root" \
  --output-root "$output_root" \
  --figure3b-min-shared-tfs 2 \
  >"$output_root/logs/00_materialize.log" 2>&1

Rscript "$script_root/prepare_figure3b_official_inputs.R" \
  "$edge_path" "$figure2_root" "$output_root/data" \
  >"$output_root/logs/01_prepare_inputs.log" 2>&1

inputs="$output_root/data/inputs"
tables="$output_root/data/official_tables"
visuals="$output_root/results/visuals/Figure3"
mkdir -p "$tables"

Rscript "$output_root/work/Figure3/04-Summarize_Regulons_HC_TR.R" \
  -cohort TCGA \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -JASPAR_CISBP_compare "$inputs/JASPAR_CISBP_HC_TF_motif.tsv" \
  -VIPER_results "$inputs/TCGA_Regulon_TF_Target_Weights_Viper.tsv" \
  -visuals_outdir "$visuals" \
  -data_outdir "$tables" \
  >"$output_root/logs/02_official_figure3A_dataflow.log" 2>&1

Rscript "$output_root/work/Figure3/05-Generate_Regulon_Network.R" \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -HC_TR_VIPER_results "$tables/HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv" \
  -data_outdir "$tables" \
  -visuals_outdir "$visuals" \
  >"$output_root/logs/03_official_figure3B.log" 2>&1

python3 "$script_root/verify_figure3b_official_hc45.py" \
  --run-root "$output_root" \
  >"$output_root/logs/04_verify.log" 2>&1
printf 'status\tvalue\nrun_status\tPASS\n' >"$output_root/audit/run_completion.tsv"
echo "Official-code Figure 3B HC45/shared>=2 complete: $output_root"

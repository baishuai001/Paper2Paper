#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  cat >&2 <<'USAGE'
Usage: run_official_fig3_port_cloud.sh \
  OFFICIAL_CODE_ROOT FIGURE1_ROOT FIGURE2_ROOT FIGURE3_ROOT \
  PUBLICATION_ROOT OUTPUT_ROOT R_LIBRARY
USAGE
  exit 2
fi

official_code_root=$(realpath "$1")
figure1_root=$(realpath "$2")
figure2_root=$(realpath "$3")
figure3_root=$(realpath "$4")
publication_root=$(realpath "$5")
output_root=$6
r_library=$(realpath "$7")

mkdir -p "$output_root"/{audit,data,logs,results/visuals/Figure3,work}
output_root=$(realpath "$output_root")

script_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export R_LIBS_USER="$r_library"
export PAPER2PAPER_OFFICIAL_PORT_WORK="$output_root/work"
export PAPER2PAPER_OFFICIAL_ADAPTER_INPUTS="$output_root/data/inputs"
export PAPER2PAPER_OFFICIAL_VISUALS_OUT="$output_root/results/visuals/Figure3"
export PAPER2PAPER_REPRESENTATIVE_TFS_FILE="$output_root/data/representative_tfs.txt"

python3 "$script_root/materialize_official_scripts.py" \
  --official-code-root "$official_code_root" \
  --output-root "$output_root" \
  --figure3b-min-shared-tfs "${FIGURE3B_MIN_SHARED_TFS:-3}" \
  >"$output_root/logs/00_materialize.log" 2>&1

Rscript "$script_root/prepare_luad_official_inputs.R" \
  --fig1-root "$figure1_root" \
  --fig2-root "$figure2_root" \
  --fig3-root "$figure3_root" \
  --publication-root "$publication_root" \
  --output-root "$output_root/data" \
  >"$output_root/logs/01_prepare_inputs.log" 2>&1

official_scripts="$output_root/work/Figure3"
inputs="$output_root/data/inputs"
visuals="$output_root/results/visuals/Figure3"
official_tables="$output_root/data/official_tables"
mkdir -p "$official_tables"

# Figure 3A: exact official constructor and official Shared/Non-Shared data flow.
Rscript "$official_scripts/04-Summarize_Regulons_HC_TR.R" \
  -cohort TCGA \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -JASPAR_CISBP_compare "$inputs/JASPAR_CISBP_HC_TF_motif.tsv" \
  -VIPER_results "$inputs/TCGA_Regulon_TF_Target_Weights_Viper.tsv" \
  -visuals_outdir "$visuals" \
  -data_outdir "$official_tables" \
  >"$output_root/logs/02_figure3A_official04.log" 2>&1

# Figures 3B, 3D and 3E: one execution of the official network script. The
# Figure 3D edge list is complete at min_shared=1; Figure 3B remains at the
# author's min_shared_tfs=3 and graph_from_data_frame endpoint behavior.
Rscript "$official_scripts/05-Generate_Regulon_Network.R" \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -HC_TR_VIPER_results "$official_tables/HC_TR_94_TCGA_VIPER_RegulonNetwork.tsv" \
  -data_outdir "$official_tables" \
  -visuals_outdir "$visuals" \
  >"$output_root/logs/03_figure3BDE_official05.log" 2>&1

Rscript "$script_root/select_luad_representatives.R" \
  "$official_tables/HC_TR_60_Plotted_Network.tsv" \
  "$PAPER2PAPER_REPRESENTATIVE_TFS_FILE" \
  "$output_root/audit/figure3G_representative_tfs.tsv" \
  >"$output_root/logs/04_select_representatives.log" 2>&1

# Figure 3C: exact official Pearson/ComplexHeatmap/bar constructors.
Rscript "$official_scripts/03-Pearson_Correlation_TRActivity.R" \
  -cohort TCGA \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -TR_activity_score "$inputs/TCGA_TR_activities_viper.tsv" \
  -TCGA_subtype_info "$inputs/TCGA_LUAD.tsv" \
  -visuals_outdir "$visuals" \
  >"$output_root/logs/05_figure3C_official03.log" 2>&1

# Figures 3F/G and Supplementary Figure 5A: exact official density, skewness
# and cross-system constructors; only the four LUAD representative IDs and
# subtype names are adapted.
Rscript "$official_scripts/06-Activity_Score_Heterogeneity_Across_Patients.R" \
  -HC_TRs "$inputs/High_RNA_NES.tsv" \
  -HC_TR_Network_Ranking "$official_tables/HC_TR_60_Plotted_Network.tsv" \
  -TCGA_TR_Activity_Score "$inputs/TCGA_TR_activities_viper.tsv" \
  -TCGA_subtype_info "$inputs/TCGA_LUAD.tsv" \
  -PDX_TR_Activity_Score "$inputs/PDX_TR_activities_viper.tsv" \
  -PDX_subtype_info "$inputs/PDX_molecular_subtype.tsv" \
  -CellLines_TR_Activity_Score "$inputs/CellLines_TR_activities_viper.tsv" \
  -CellLines_subtype_info "$inputs/CellLines_molecular_subtype.tsv" \
  -data_outdir "$official_tables" \
  -visuals_outdir "$visuals" \
  >"$output_root/logs/06_figure3FG_supp5A_official06.log" 2>&1

# Linked Supplementary Figure 4D: official volcano and heatmap constructors.
Rscript "$official_scripts/07B-Plot_MKi67_HC-TR.R" \
  -results \
  "$inputs/MKI67/TCGA_MKi67_Correlation_HC-TR_results.tsv" \
  "$inputs/MKI67/GSE41271_MKi67_Correlation_HC-TR_results.tsv" \
  "$inputs/MKI67/PDX_MKi67_Correlation_HC-TR_results.tsv" \
  "$inputs/MKI67/CellLines_MKi67_Correlation_HC-TR_results.tsv" \
  -visuals_out "$visuals" \
  >"$output_root/logs/07_supp4D_volcano_official07B.log" 2>&1

Rscript "$official_scripts/07C-Plot_Heatmap_MKI67_HC-TR.R" \
  >"$output_root/logs/08_supp4D_heatmap_official07C.log" 2>&1

# Linked Supplementary Figure 5B: official volcano and heatmap constructors.
estimate_tcga="$inputs/ESTIMATE/TCGA_ESTIMATE_all_results.tsv"
estimate_gse="$inputs/ESTIMATE/GSE41271_ESTIMATE_all_results.tsv"
Rscript "$official_scripts/08D-Plot_ESTIMATE_Volcano.R" \
  -results "$estimate_tcga" "$estimate_gse" \
  -cohort_order "TCGA,GSE41271" \
  -score_order "StromalScore,ImmuneScore,TumorPurity" \
  -visuals_out "$visuals" \
  >"$output_root/logs/09_supp5B_volcano_official08D.log" 2>&1

Rscript "$official_scripts/08E-Plot_Heatmap_ESTIMATE.R" \
  -results "$estimate_tcga,$estimate_gse" \
  -cohort_order "TCGA,GSE41271" \
  -score_order "StromalScore,ImmuneScore,TumorPurity" \
  -plot_title "LUAD HC-TF activity correlation" \
  -out_name "SupplementaryFigure5B_ESTIMATE_Combined_Heatmap_HC-TR" \
  -visuals_out "$visuals" \
  >"$output_root/logs/10_supp5B_heatmap_official08E.log" 2>&1

python3 "$script_root/verify_official_fig3_port.py" \
  --output-root "$output_root" \
  >"$output_root/logs/11_verify.log" 2>&1

printf 'status\tvalue\nrun_status\tPASS\n' >"$output_root/audit/run_completion.tsv"
echo "Strict official Figure 3 port completed: $output_root"

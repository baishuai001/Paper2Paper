#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -lt 4 || "$#" -gt 5 ]]; then
  echo "usage: run_official_supp3_cloud.sh OFFICIAL_CODE_ROOT LUAD_FIGURE2_ROOT PORT_CODE_DIR WORK_ROOT [R_LIB]" >&2
  exit 2
fi

OFFICIAL_CODE_ROOT="$(readlink -f "$1")"
LUAD_FIGURE2_ROOT="$(readlink -f "$2")"
PORT_CODE_DIR="$(readlink -f "$3")"
WORK_ROOT="$(readlink -m "$4")"
R_LIB="${5:-/media/desk16/iy13202/R/x86_64-pc-linux-gnu-library/4.4}"

mkdir -p "$WORK_ROOT" "$WORK_ROOT/logs" "$WORK_ROOT/results/panels" "$WORK_ROOT/results/rendered_png"
export R_LIBS_USER="$R_LIB"

python3 "$PORT_CODE_DIR/materialize_official_scripts.py" "$OFFICIAL_CODE_ROOT" "$WORK_ROOT"
Rscript "$PORT_CODE_DIR/prepare_luad_inputs.R" "$LUAD_FIGURE2_ROOT" "$WORK_ROOT/adapter" \
  >"$WORK_ROOT/logs/00_prepare_luad_inputs.log" 2>&1

RUNTIME_F2="$WORK_ROOT/runtime/code/Figure2"
PANELS="$WORK_ROOT/results/panels"
ADAPTER="$WORK_ROOT/adapter"

declare -a COHORTS=("TCGA" "PDX" "CellLines")
declare -a PEAK_PANELS=("SupplementaryFigure3A" "SupplementaryFigure3B" "SupplementaryFigure3C")
declare -a SAT_PANELS=("SupplementaryFigure3D" "SupplementaryFigure3E" "SupplementaryFigure3F")
declare -a COR_PANELS=("SupplementaryFigure3G" "SupplementaryFigure3H" "SupplementaryFigure3I")
declare -a GENOMIC_PANELS=("SupplementaryFigure3J" "SupplementaryFigure3K" "SupplementaryFigure3L")

for i in 0 1 2; do
  cohort="${COHORTS[$i]}"
  panel="${PEAK_PANELS[$i]}"
  Rscript "$RUNTIME_F2/02J2-Plot_Number_Peaks.R" \
    -cohort "$cohort" \
    -Number_Called_Peaks_file "$ADAPTER/number_peaks/${cohort}_Number_Peaks.tsv" \
    -file_prefix "$panel" \
    -visual_outdir "$PANELS" \
    >"$WORK_ROOT/logs/${panel}_official_02J2.log" 2>&1
done

for i in 0 1 2; do
  cohort="${COHORTS[$i]}"
  panel="${SAT_PANELS[$i]}"
  Rscript "$RUNTIME_F2/02C-Peak_Saturation_Plot.R" \
    -cohort "$cohort" \
    -file_prefix "$panel" \
    -peak_saturation_dir "$ADAPTER/saturation" \
    -outdir_visuals "$PANELS" \
    >"$WORK_ROOT/logs/${panel}_official_02C.log" 2>&1
done

for i in 0 1 2; do
  cohort="${COHORTS[$i]}"
  panel="${COR_PANELS[$i]}"
  extra=()
  if [[ "$cohort" == "CellLines" ]]; then
    extra=(-meta_data "$ADAPTER/correlation/CellLines_metadata.tsv")
  fi
  Rscript "$RUNTIME_F2/03D-Pearson_Correlation_Accessibility_ConsensusPeakSet.R" \
    -cohort "$cohort" \
    -pearson_cor "$ADAPTER/correlation/${cohort}_pearson_cor.rds" \
    -file_prefix "$panel" \
    -visual_outdir "$PANELS" \
    "${extra[@]}" \
    >"$WORK_ROOT/logs/${panel}_official_03D.log" 2>&1
done

GENOME_BASELINE="$(cat "$ADAPTER/GENOME_BASELINE_PATH.txt")"
for i in 0 1 2; do
  cohort="${COHORTS[$i]}"
  panel="${GENOMIC_PANELS[$i]}"
  Rscript "$RUNTIME_F2/04B-GenomicAnnotation.R" \
    -cohort "$cohort" \
    -consensus_peakset_dir "$ADAPTER/consensus" \
    -genome_ref "$GENOME_BASELINE" \
    -outfile_data "$WORK_ROOT/results/tables/${panel}_${cohort}_GenomicAnnotation.tsv" \
    -outfile_plot "$PANELS/${panel}_${cohort}_GenomicAnnotation.pdf" \
    >"$WORK_ROOT/logs/${panel}_official_04B.log" 2>&1
done

while IFS= read -r pdf; do
  stem="$(basename "${pdf%.pdf}")"
  pdftoppm -png -r 150 -singlefile "$pdf" "$WORK_ROOT/results/rendered_png/$stem" >/dev/null 2>&1
done < <(find "$PANELS" -maxdepth 1 -type f -name '*.pdf' | sort)

python3 "$PORT_CODE_DIR/verify_official_supp3.py" "$WORK_ROOT" "$PORT_CODE_DIR/PANEL_STATUS.tsv"
echo "Official Supplementary Figure 3 port complete: $WORK_ROOT"

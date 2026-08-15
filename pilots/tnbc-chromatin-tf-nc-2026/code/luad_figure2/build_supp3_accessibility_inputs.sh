#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: build_supp3_accessibility_inputs.sh <figure2_run_root>" >&2
  exit 2
fi
RUN_ROOT=$(realpath "$1")
MANIFEST="$RUN_ROOT/audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv"
OUT="$RUN_ROOT/data/processed/supplementary_figure3"
LOG="$RUN_ROOT/logs/supplementary_figure3"
mkdir -p "$OUT" "$LOG"
[[ -s "$MANIFEST" ]] || { echo "Missing Figure 2 sample manifest" >&2; exit 2; }

build_system() {
  local system=$1
  local manifest="$OUT/${system}_sample_files.tsv"
  awk -F '\t' -v target="$system" 'BEGIN{OFS="\t"; print "sample_id","sample_slug","peak_path","bam_path"}
    NR>1 && $1==target {print $2,$3,$4,$5}' "$MANIFEST" > "$manifest"
  local n
  n=$(( $(wc -l < "$manifest") - 1 ))
  if (( n < 1 )); then echo "No analysis samples for $system" >&2; exit 1; fi

  local combined="$OUT/${system}_all_sample_peaks.bed"
  : > "$combined"
  while IFS=$'\t' read -r sample slug peak bam; do
    [[ -s "$peak" ]] || { echo "Missing $system peak file: $peak" >&2; exit 1; }
    if [[ "$peak" == *.gz ]]; then gzip -dc "$peak"; else cat "$peak"; fi \
      | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ {print $1,$2,$3}' >> "$combined"
  done < <(tail -n +2 "$manifest" | sed 's/\r$//')
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$combined" \
    | bedtools merge -i - > "$OUT/${system}_consensus_peaks.bed"
  rm -f "$combined"

  local binary="$OUT/${system}_consensus_peak_presence.tsv"
  {
    printf 'chrom\tstart\tend'
    while IFS=$'\t' read -r sample slug peak bam; do printf '\t%s' "$sample"; done < <(tail -n +2 "$manifest" | sed 's/\r$//')
    printf '\n'
    presence_files=()
    while IFS=$'\t' read -r sample slug peak bam; do
      presence_file="$OUT/.${system}.${slug}.presence"
      presence_files+=("$presence_file")
      if [[ "$peak" == *.gz ]]; then
        bedtools intersect -a "$OUT/${system}_consensus_peaks.bed" -b <(gzip -dc "$peak") -c | awk '{print ($NF>0)?1:0}' > "$presence_file"
      else
        bedtools intersect -a "$OUT/${system}_consensus_peaks.bed" -b "$peak" -c | awk '{print ($NF>0)?1:0}' > "$presence_file"
      fi
    done < <(tail -n +2 "$manifest" | sed 's/\r$//')
    paste "$OUT/${system}_consensus_peaks.bed" <(paste "${presence_files[@]}")
  } > "$binary"
  rm -f "$OUT"/."${system}".*.presence

  if [[ "$system" != "patient" ]]; then
    mapfile -t bams < <(tail -n +2 "$manifest" | sed 's/\r$//' | cut -f4)
    bedtools multicov -bed "$OUT/${system}_consensus_peaks.bed" -bams "${bams[@]}" \
      > "$OUT/${system}_consensus_peak_read_counts.tsv"
  fi
  echo "SUPP3_INPUT_COMPLETE system=$system samples=$n peaks=$(wc -l < "$OUT/${system}_consensus_peaks.bed")"
}

build_system patient
build_system PDX
build_system cell_line
echo "SUPP3_ALL_ACCESSIBILITY_INPUTS_COMPLETE"

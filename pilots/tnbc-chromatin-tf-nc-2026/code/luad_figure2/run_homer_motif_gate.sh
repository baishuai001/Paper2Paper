#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: run_homer_motif_gate.sh <figure2_run_root>" >&2
  exit 2
fi
RUN_ROOT=$(realpath "$1")
HOMER="$RUN_ROOT/tools/homer"
REF="$RUN_ROOT/reference"
MOTIF_ROOT="$REF/motifs"
MANIFEST="$RUN_ROOT/audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv"
TCGA_ZIP="$RUN_ROOT/data/raw/tcga_atac_gdc/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip"
INPUT="$RUN_ROOT/data/processed/motif_inputs"
RESULTS="$RUN_ROOT/results/motif_gate/homer_per_sample"
LOG="$RUN_ROOT/logs/motif_gate"
JASPAR_HOMER="$MOTIF_ROOT/selected_hc_tf_JASPAR2024.homer"
CISBP_HOMER="$MOTIF_ROOT/selected_hc_tf_CIS-BP2.00.homer"
PARALLEL_JOBS=${PARALLEL_JOBS:-4}
THREADS_PER_JOB=${THREADS_PER_JOB:-2}
mkdir -p "$INPUT/targets" "$RESULTS" "$LOG"
for required in "$HOMER/bin/findMotifsGenome.pl" "$MANIFEST" "$TCGA_ZIP" "$JASPAR_HOMER" "$CISBP_HOMER" "$REF/hg38.fa"; do
  [[ -s "$required" ]] || { echo "Missing motif input/tool: $required" >&2; exit 2; }
done
export PATH="$HOMER/bin:$PATH"

BACKGROUND_BED="$INPUT/lung_accessible_CRE_background.bed"
if [[ ! -s "$BACKGROUND_BED" ]]; then
  temporary="$INPUT/.background_unsorted.bed"
  : > "$temporary"
  for cancer_type in LUAD LUSC; do
    unzip -p "$TCGA_ZIP" "${cancer_type}_raw_counts.txt" \
      | awk 'BEGIN{OFS="\t"} NR>1 && $1 ~ /^chr([0-9]+|X)$/ {print $1,$2,$3}' >> "$temporary"
  done
  while IFS=$'\t' read -r system sample sample_slug peak bam peak_definition idr; do
    if [[ "$peak" == *.gz ]]; then gzip -dc "$peak"; else cat "$peak"; fi \
      | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ {print $1,$2,$3}' >> "$temporary"
  done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$temporary" | bedtools merge -i - > "$BACKGROUND_BED"
  rm -f "$temporary"
fi
BACKGROUND_POS="$INPUT/lung_accessible_CRE_background.homer.pos"
awk 'BEGIN{OFS="\t"} {print "background_"NR,$1,$2,$3,"+"}' "$BACKGROUND_BED" > "$BACKGROUND_POS"

prepare_target() {
  local system=$1 sample=$2 slug=$3 peak=$4
  local pos="$INPUT/targets/${system}.${slug}.homer.pos"
  if [[ ! -s "$pos" ]]; then
    if [[ "$peak" == *.gz ]]; then gzip -dc "$peak"; else cat "$peak"; fi \
      | awk -v prefix="${system}_${slug}_" 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ {print prefix NR,$1,$2,$3,"+"}' > "$pos"
  fi
  [[ -s "$pos" ]] || { echo "Empty HOMER target for $system/$sample" >&2; return 1; }
}

run_sample() {
  local system=$1 sample=$2 slug=$3 database=$4 motif_file=$5
  local pos="$INPUT/targets/${system}.${slug}.homer.pos"
  local out="$RESULTS/${database}/${system}/${slug}"
  mkdir -p "$out"
  if [[ -s "$out/knownResults.txt" && -s "$out/.complete" ]]; then
    echo "SKIP_HOMER_COMPLETE $database $system $sample"
    return 0
  fi
  findMotifsGenome.pl "$pos" "$REF/hg38.fa" "$out" \
    -size given -bg "$BACKGROUND_POS" -nomotif -mknown "$motif_file" -p "$THREADS_PER_JOB" \
    > "$LOG/${database}.${system}.${slug}.log" 2>&1
  [[ -s "$out/knownResults.txt" ]] || { echo "HOMER failed for $database/$system/$sample" >&2; return 1; }
  printf 'completed_utc=%s\ndatabase=%s\n' "$(date -u +%FT%TZ)" "$database" > "$out/.complete"
  echo "HOMER_COMPLETE $database $system $sample"
}

while IFS=$'\t' read -r system sample sample_slug peak bam peak_definition idr; do
  prepare_target "$system" "$sample" "$sample_slug" "$peak"
done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')

running=0
failed=0
while IFS=$'\t' read -r system sample sample_slug peak bam peak_definition idr; do
  for database in JASPAR2024 CIS-BP2.00; do
    if [[ "$database" == "JASPAR2024" ]]; then
      motif_file="$JASPAR_HOMER"
    else
      motif_file="$CISBP_HOMER"
    fi
    run_sample "$system" "$sample" "$sample_slug" "$database" "$motif_file" &
    running=$((running + 1))
    if (( running >= PARALLEL_JOBS )); then
      if ! wait -n; then failed=$((failed + 1)); fi
      running=$((running - 1))
    fi
  done
done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')
while (( running > 0 )); do
  if ! wait -n; then failed=$((failed + 1)); fi
  running=$((running - 1))
done
if (( failed > 0 )); then echo "HOMER_SAMPLE_FAILURES=$failed" >&2; exit 1; fi
echo "HOMER_MOTIF_GATE_ALL_SAMPLES_COMPLETE"

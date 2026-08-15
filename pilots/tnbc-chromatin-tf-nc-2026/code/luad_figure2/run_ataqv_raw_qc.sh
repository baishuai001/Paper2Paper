#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_ataqv_raw_qc.sh <figure2_run_root> <project_root>" >&2
  exit 2
fi
RUN_ROOT=$(realpath "$1")
PROJECT_ROOT=$(realpath "$2")
TOOLS="$RUN_ROOT/tools/ataqv"
DOWNLOADS="$TOOLS/downloads"
MANIFEST="$RUN_ROOT/audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv"
TSS="$RUN_ROOT/reference/figure2_features/gencode.v47.protein_coding.gene_tss.bed"
BLACKLIST="$RUN_ROOT/reference/hg38-blacklist.v2.bed"
OUT="$RUN_ROOT/audit/ataqv"
LOG="$RUN_ROOT/logs/ataqv"
PARALLEL_JOBS=${PARALLEL_JOBS:-4}
THREADS_PER_JOB=${THREADS_PER_JOB:-4}
mkdir -p "$DOWNLOADS" "$OUT" "$LOG"
[[ -s "$MANIFEST" && -s "$BLACKLIST" ]] || { echo "Missing ataqv manifest/reference" >&2; exit 2; }

if [[ ! -x "$TOOLS/root/usr/bin/ataqv" ]]; then
  (cd "$DOWNLOADS" && apt-get download ataqv=1.3.0+ds-1) > "$LOG/ataqv_apt_download.log" 2>&1
  deb=$(find "$DOWNLOADS" -maxdepth 1 -type f -name 'ataqv_1.3.0+ds-1_amd64.deb' -print -quit)
  [[ -s "$deb" ]] || { echo "ataqv .deb download failed" >&2; exit 1; }
  mkdir -p "$TOOLS/root"
  dpkg-deb -x "$deb" "$TOOLS/root"
fi
ATAQV="$TOOLS/root/usr/bin/ataqv"
[[ -x "$ATAQV" ]] || { echo "ataqv binary unavailable" >&2; exit 1; }

if [[ ! -s "$TSS" ]]; then
  python3 "$PROJECT_ROOT/pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2/make_gencode_gene_tss_bed.py" \
    "$RUN_ROOT/reference/downloads/gencode.v47.basic.annotation.gtf.gz" "$TSS"
fi

run_one() {
  local system=$1 sample=$2 slug=$3 bam=$4 peak=$5
  local sample_out="$OUT/$system/$slug"
  mkdir -p "$sample_out"
  if [[ -s "$sample_out/metrics.ataqv.json.gz" && -s "$sample_out/.complete" ]]; then
    echo "SKIP_ATAQV_COMPLETE $system $sample"
    return 0
  fi
  "$ATAQV" --threads "$THREADS_PER_JOB" --peak-file "$peak" --tss-file "$TSS" \
    --excluded-region-file "$BLACKLIST" --name "$sample" --ignore-read-groups \
    --metrics-file "$sample_out/metrics.ataqv.json.gz" human "$bam" \
    > "$sample_out/summary.txt" 2> "$LOG/${system}.${slug}.log"
  [[ -s "$sample_out/metrics.ataqv.json.gz" ]] || { echo "ataqv failed: $system/$sample" >&2; return 1; }
  printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" > "$sample_out/.complete"
  echo "ATAQV_COMPLETE $system $sample"
}

running=0
failed=0
while IFS=$'\t' read -r system sample slug peak bam peak_definition idr; do
  [[ "$system" != "patient" ]] || continue
  run_one "$system" "$sample" "$slug" "$bam" "$peak" &
  running=$((running + 1))
  if (( running >= PARALLEL_JOBS )); then
    if ! wait -n; then failed=$((failed + 1)); fi
    running=$((running - 1))
  fi
done < <(tail -n +2 "$MANIFEST" | sed 's/\r$//')
while (( running > 0 )); do
  if ! wait -n; then failed=$((failed + 1)); fi
  running=$((running - 1))
done
if (( failed > 0 )); then echo "ATAQV_FAILURES=$failed" >&2; exit 1; fi
{
  echo "ataqv=$($ATAQV --version 2>&1 | head -1)"
  sha256sum "$ATAQV" "$TSS"
} > "$OUT/ataqv_tool_reference_receipt.txt"
echo "ATAQV_ALL_RAW_SAMPLES_COMPLETE"

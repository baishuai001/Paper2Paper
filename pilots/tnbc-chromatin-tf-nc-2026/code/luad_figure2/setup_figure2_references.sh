#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: setup_figure2_references.sh <figure2_run_root>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
REF="$RUN_ROOT/reference"
LOG="$RUN_ROOT/logs"
TOOLS="$RUN_ROOT/tools"
mkdir -p "$REF/downloads" "$REF/hg38_bowtie2" "$REF/mm10_bowtie2" "$LOG" "$TOOLS"

download() {
  local url=$1
  local out=$2
  if [[ ! -s "$out" ]]; then
    curl -fL --retry 8 --retry-delay 3 "$url" -o "${out}.part"
    mv "${out}.part" "$out"
  fi
}

download \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
  "$REF/downloads/hg38.fa.gz"
download \
  "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" \
  "$REF/downloads/mm10.fa.gz"
download \
  "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz" \
  "$REF/downloads/hg38-blacklist.v2.bed.gz"
download \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz" \
  "$REF/downloads/gencode.v47.basic.annotation.gtf.gz"

if [[ ! -s "$REF/hg38.fa" ]]; then gzip -dc "$REF/downloads/hg38.fa.gz" > "$REF/hg38.fa"; fi
if [[ ! -s "$REF/mm10.fa" ]]; then gzip -dc "$REF/downloads/mm10.fa.gz" > "$REF/mm10.fa"; fi
if [[ ! -s "$REF/hg38-blacklist.v2.bed" ]]; then
  gzip -dc "$REF/downloads/hg38-blacklist.v2.bed.gz" > "$REF/hg38-blacklist.v2.bed"
fi
samtools faidx "$REF/hg38.fa"
samtools faidx "$REF/mm10.fa"

index_complete() {
  local prefix=$1
  local extension
  for extension in bt2 bt2l; do
    if [[ -s "${prefix}.1.${extension}" && -s "${prefix}.2.${extension}" && \
          -s "${prefix}.3.${extension}" && -s "${prefix}.4.${extension}" && \
          -s "${prefix}.rev.1.${extension}" && -s "${prefix}.rev.2.${extension}" ]]; then
      return 0
    fi
  done
  return 1
}

if ! index_complete "$REF/hg38_bowtie2/hg38"; then
  bowtie2-build --threads 24 "$REF/hg38.fa" "$REF/hg38_bowtie2/hg38" \
    > "$LOG/bowtie2_build_hg38.log" 2>&1
fi
if ! index_complete "$REF/mm10_bowtie2/mm10"; then
  bowtie2-build --threads 24 "$REF/mm10.fa" "$REF/mm10_bowtie2/mm10" \
    > "$LOG/bowtie2_build_mm10.log" 2>&1
fi

if [[ ! -x "$TOOLS/cutadapt-venv/bin/cutadapt" ]]; then
  python3 -m venv "$TOOLS/cutadapt-venv"
  "$TOOLS/cutadapt-venv/bin/pip" install --disable-pip-version-check "cutadapt==4.9" \
    > "$LOG/install_cutadapt.log" 2>&1
fi
if [[ ! -x "$TOOLS/macs2-venv/bin/macs2" ]]; then
  python3 -m venv "$TOOLS/macs2-venv"
  "$TOOLS/macs2-venv/bin/pip" install --disable-pip-version-check \
    "numpy==1.26.4" "MACS2==2.2.9.1" > "$LOG/install_macs2.log" 2>&1
fi

{
  printf 'resource\tsha256\n'
  sha256sum \
    "$REF/downloads/hg38.fa.gz" \
    "$REF/downloads/mm10.fa.gz" \
    "$REF/downloads/hg38-blacklist.v2.bed.gz" \
    "$REF/downloads/gencode.v47.basic.annotation.gtf.gz" \
    | awk '{print $2"\t"$1}'
} > "$REF/reference_sha256.tsv"

{
  echo "bowtie2=$(bowtie2 --version | head -1)"
  echo "samtools=$(samtools --version | head -1)"
  echo "bedtools=$(bedtools --version)"
  echo "macs2=$($TOOLS/macs2-venv/bin/macs2 --version)"
  echo "cutadapt=$($TOOLS/cutadapt-venv/bin/cutadapt --version)"
  echo "gencode=47 basic annotation"
  echo "reference_downloaded_utc=$(date -u +%FT%TZ)"
} > "$REF/reference_tool_versions.txt"

echo "REFERENCE_SETUP_COMPLETE"
cat "$REF/reference_sha256.tsv"

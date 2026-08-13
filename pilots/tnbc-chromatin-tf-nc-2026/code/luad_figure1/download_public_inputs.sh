#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: download_public_inputs.sh CLOUD_RUN_ROOT" >&2
  exit 2
fi

run_root=$1
raw="$run_root/data/raw"
audit="$run_root/audit/downloads"
mkdir -p "$raw/tcga" "$raw/gse81089" "$raw/depmap" "$raw/pdmr" "$audit"

download() {
  local url=$1
  local destination=$2
  local remote_size local_size
  remote_size=$(curl -fsSIL --retry 3 "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" {size=$2} END {print size+0}')
  local_size=0
  if [[ -f "$destination" ]]; then
    local_size=$(stat -c '%s' "$destination")
  fi
  if [[ "$remote_size" -le 0 || "$local_size" -ne "$remote_size" ]]; then
    aria2c \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --continue=true \
      --max-connection-per-server=8 \
      --min-split-size=4M \
      --split=8 \
      --max-tries=5 \
      --retry-wait=3 \
      --dir="$(dirname "$destination")" \
      --out="$(basename "$destination")" \
      "$url"
  fi
  if [[ "$remote_size" -gt 0 && $(stat -c '%s' "$destination") -ne "$remote_size" ]]; then
    echo "Downloaded size mismatch for $destination" >&2
    exit 6
  fi
}

download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_counts.tsv.gz" \
  "$raw/tcga/TCGA-LUAD.star_counts.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.star_counts.tsv.gz" \
  "$raw/tcga/TCGA-LUSC.star_counts.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz" \
  "$raw/tcga/TCGA-LUAD.star_tpm.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.star_tpm.tsv.gz" \
  "$raw/tcga/TCGA-LUSC.star_tpm.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.clinical.tsv.gz" \
  "$raw/tcga/TCGA-LUAD.clinical.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.clinical.tsv.gz" \
  "$raw/tcga/TCGA-LUSC.clinical.tsv.gz"
download \
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap" \
  "$raw/tcga/gencode.v36.annotation.gtf.gene.probemap"

geo_base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089"
download \
  "$geo_base/suppl/GSE81089_FPKM_cufflinks.tsv.gz" \
  "$raw/gse81089/GSE81089_FPKM_cufflinks.tsv.gz"
download \
  "$geo_base/suppl/GSE81089_readcounts_featurecounts.tsv.gz" \
  "$raw/gse81089/GSE81089_readcounts_featurecounts.tsv.gz"
download \
  "$geo_base/soft/GSE81089_family.soft.gz" \
  "$raw/gse81089/GSE81089_family.soft.gz"

# DepMap 22Q2 is retrieved through Bioconductor ExperimentHub so that the
# resource version and provider metadata are frozen in the receipt.
Rscript "$(dirname "$0")/fetch_depmap_22q2.R" "$raw/depmap" "$audit"

find "$raw" -type f -print0 | sort -z | xargs -0 sha256sum > "$audit/input_sha256.tsv"
find "$raw" -type f -printf '%p\t%s\n' | sort > "$audit/input_sizes.tsv"

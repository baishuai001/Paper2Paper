#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: setup_figure2_motifs.sh <figure2_run_root>" >&2
  exit 2
fi
RUN_ROOT=$(realpath "$1")
TOOLS="$RUN_ROOT/tools"
MOTIFS="$RUN_ROOT/reference/motifs"
LOG="$RUN_ROOT/logs"
HOMER="$TOOLS/homer"
mkdir -p "$HOMER" "$MOTIFS/downloads" "$LOG"

download() {
  local url=$1
  local out=$2
  if [[ ! -s "$out" ]]; then
    curl -fL --retry 8 --retry-delay 5 "$url" -o "${out}.part"
    mv "${out}.part" "$out"
  fi
}

download "http://homer.ucsd.edu/homer/configureHomer.pl" "$HOMER/configureHomer.pl"
if [[ ! -x "$HOMER/bin/findMotifsGenome.pl" ]]; then
  (cd "$HOMER" && perl configureHomer.pl -install homer) > "$LOG/install_homer.log" 2>&1
fi
[[ -x "$HOMER/bin/findMotifsGenome.pl" ]] || { echo "HOMER installation failed" >&2; exit 1; }

download \
  "https://jaspar2024.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt" \
  "$MOTIFS/downloads/JASPAR2024_CORE_vertebrates_non-redundant.meme"
download \
  "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.27.tgz" \
  "$MOTIFS/downloads/motif_databases.12.27.tgz"

if [[ ! -s "$MOTIFS/CIS-BP_2.00_Homo_sapiens.meme" ]]; then
  # Do not exit awk early under pipefail: that sends SIGPIPE to tar and makes
  # an otherwise valid archive lookup terminate the setup script.
  member=$(tar -tzf "$MOTIFS/downloads/motif_databases.12.27.tgz" \
    | awk '/CIS-BP_2\.00\/Homo_sapiens\.meme$/ {found=$0} END {print found}')
  [[ -n "$member" ]] || { echo "CIS-BP Homo_sapiens.meme not found in MEME database snapshot" >&2; exit 1; }
  tar -xOzf "$MOTIFS/downloads/motif_databases.12.27.tgz" "$member" \
    > "$MOTIFS/CIS-BP_2.00_Homo_sapiens.meme"
fi

{
  printf 'resource\tsha256\n'
  sha256sum \
    "$HOMER/configureHomer.pl" \
    "$MOTIFS/downloads/JASPAR2024_CORE_vertebrates_non-redundant.meme" \
    "$MOTIFS/downloads/motif_databases.12.27.tgz" \
    "$MOTIFS/CIS-BP_2.00_Homo_sapiens.meme" \
    | awk '{print $2"\t"$1}'
} > "$MOTIFS/motif_resource_sha256.tsv"

{
  echo "homer_findMotifsGenome=$($HOMER/bin/findMotifsGenome.pl 2>&1 | head -1 || true)"
  echo "homer_config=$(grep -m1 '^homer' "$HOMER/config.txt" 2>/dev/null || true)"
  echo "jaspar=JASPAR 2024 CORE vertebrates non-redundant"
  echo "cisbp=MEME motif_databases.12.27 CIS-BP 2.00 Homo sapiens"
  echo "completed_utc=$(date -u +%FT%TZ)"
} > "$MOTIFS/motif_tool_versions.txt"
echo "MOTIF_SETUP_COMPLETE"
cat "$MOTIFS/motif_resource_sha256.tsv"

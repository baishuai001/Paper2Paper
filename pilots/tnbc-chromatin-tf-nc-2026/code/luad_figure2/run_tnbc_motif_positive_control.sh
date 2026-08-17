#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_tnbc_motif_positive_control.sh <tnbc_replay_root> <author_code_root>" >&2
  exit 2
fi

TNBC_ROOT=$(realpath "$1")
AUTHOR_ROOT=$(realpath "$2")
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_DIR="$TNBC_ROOT/data/codeocean"
HOMER_DIR="$SOURCE_DIR/HOMER"
REPLAY_DIR="$TNBC_ROOT/results/author_code_replay"

mkdir -p "$HOMER_DIR" "$REPLAY_DIR/combined" "$TNBC_ROOT/audit"

if [[ ! -f "$SOURCE_DIR/capsule-7227095-HOMER.zip" ]]; then
  echo "Missing Code Ocean HOMER archive: $SOURCE_DIR/capsule-7227095-HOMER.zip" >&2
  exit 1
fi

unzip -q -o "$SOURCE_DIR/capsule-7227095-HOMER.zip" -d "$HOMER_DIR"

Rscript "$AUTHOR_ROOT/code/Figure2/06B-Compare_TR_Motif_Site_Enrichment.R" \
  --category JASPAR \
  --motif_enrich_dir "$HOMER_DIR/JASPAR" \
  --output_dir "$REPLAY_DIR/combined" \
  > "$REPLAY_DIR/06B_JASPAR.log" 2>&1

Rscript "$AUTHOR_ROOT/code/Figure2/06B-Compare_TR_Motif_Site_Enrichment.R" \
  --category CISBP \
  --motif_enrich_dir "$HOMER_DIR/CISBP" \
  --output_dir "$REPLAY_DIR/combined" \
  > "$REPLAY_DIR/06B_CISBP.log" 2>&1

python3 "$CODE_DIR/replay_tnbc_author_motifs.py" "$TNBC_ROOT" \
  2>&1 | tee "$REPLAY_DIR/replay.log"

echo "TNBC_MOTIF_POSITIVE_CONTROL_COMPLETE"

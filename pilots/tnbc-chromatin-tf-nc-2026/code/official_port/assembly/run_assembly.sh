#!/usr/bin/env bash
set -euo pipefail

dry_run=0
require_publication_ready=0
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --dry-run) dry_run=1 ;;
    --require-publication-ready) require_publication_ready=1 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

if [[ $# -lt 8 || $# -gt 9 ]]; then
  cat >&2 <<'EOF'
Usage:
  run_assembly.sh [--dry-run] [--require-publication-ready] \
    FIG12_ROOT FIG3_ROOT FIG45_ROOT SUPP3_ROOT \
    TNBC_MAIN_PDF TNBC_SUPPLEMENT_PDF TNBC_REFERENCE_DIR PILOT_ROOT [UNICODE_TTF]

Use - for an unavailable SUPP3_ROOT, TNBC PDF, reference directory, or font.
Final PDFs are fixed to PILOT_ROOT/execution/luad-official-code-port/final/.
Default mode generates review PDFs and records the upstream release status.
--require-publication-ready activates the non-zero release gate.
EOF
  exit 2
fi

fig12_root=$1
fig3_root=$2
fig45_root=$3
supp3_root=$4
tnbc_main_pdf=$5
tnbc_supp_pdf=$6
tnbc_reference_dir=$7
pilot_root=$8
font_path=${9:--}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
execution_root="$pilot_root/execution/luad-official-code-port"
output_dir="$execution_root/final"
audit_dir="$execution_root/assembly-audit"
qa_dir="$execution_root/assembly-qa"
mkdir -p "$output_dir" "$audit_dir" "$qa_dir"

assembly_args=(
  --fig12-root "$fig12_root"
  --fig3-root "$fig3_root"
  --fig45-root "$fig45_root"
  --output-dir "$output_dir"
  --audit-dir "$audit_dir"
)

if [[ "$supp3_root" != "-" ]]; then
  assembly_args+=(--supp3-root "$supp3_root")
fi
if [[ "$tnbc_main_pdf" != "-" ]]; then
  assembly_args+=(--tnbc-main-pdf "$tnbc_main_pdf")
fi
if [[ "$tnbc_supp_pdf" != "-" ]]; then
  assembly_args+=(--tnbc-supplement-pdf "$tnbc_supp_pdf")
fi
if [[ "$tnbc_reference_dir" != "-" ]]; then
  assembly_args+=(--tnbc-reference-dir "$tnbc_reference_dir")
fi
if [[ "$font_path" != "-" ]]; then
  assembly_args+=(--font-path "$font_path")
fi

if [[ "$dry_run" -eq 1 ]]; then
  assembly_args+=(--dry-run)
fi
if [[ "$require_publication_ready" -eq 1 ]]; then
  assembly_args+=(--require-publication-ready)
fi

python3 "$script_dir/assemble_official_pdfs.py" "${assembly_args[@]}"

if [[ "$dry_run" -eq 1 ]]; then
  printf 'Input-contract dry-run complete: %s\n' "$audit_dir"
  exit 0
fi

qa_args=(
  --output-dir "$output_dir"
  --audit-dir "$audit_dir"
  --qa-dir "$qa_dir"
  --strict
)
if [[ "$require_publication_ready" -eq 1 ]]; then
  qa_args+=(--require-publication-ready)
fi
python3 "$script_dir/verify_assembled_pdfs.py" "${qa_args[@]}"

printf 'Review-mode layout-only assembly complete: %s\n' "$output_dir"
printf 'Authoritative release state: %s\n' "$audit_dir/release_status.json"

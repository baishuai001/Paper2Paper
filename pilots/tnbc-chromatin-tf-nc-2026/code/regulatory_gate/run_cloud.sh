#!/usr/bin/env bash
set -euo pipefail

repo=${PAPER2PAPER_REPO:-/media/desk16/iy13202/projects/Paper2Paper}
work=${REGULATORY_GATE_WORK:-$repo/tmp/tnbc-chromatin-tf-nc-2026/regulatory-gate-M-vs-rest}
h5ad=${CRC_ATLAS_H5AD:-$repo/tmp/hcc-sc-spatial-npj-2026/phase-zero/raw/crc_atlas_cellxgene.h5ad}
supplement=${TNBC_SUPPLEMENT_DATA2:-$repo/tmp/tnbc-chromatin-tf-nc-2026/external/source/41467_2026_76385_MOESM4_ESM.xlsx}
aracne_repo=${ARACNE3_REPO:-$repo/tmp/tnbc-chromatin-tf-nc-2026/external/code/ARACNe3}
anchor_code_repo=${TNBC_ANCHOR_CODE_REPO:-$repo/tmp/tnbc-chromatin-tf-nc-2026/external/code/TNBC_CodeOcean_7227095_v1}
code=$repo/pilots/tnbc-chromatin-tf-nc-2026/code/regulatory_gate
manifest=${REGULATORY_GATE_MANIFEST:-$repo/pilots/tnbc-chromatin-tf-nc-2026/analysis/phenotypes/m-vs-rest.json}
outputs=$work/outputs
network_outputs=${REGULATORY_GATE_SHARED_NETWORK_OUTPUTS:-$outputs}
logs=$work/logs
threads=${REGULATORY_GATE_THREADS:-24}
python_bin=${REGULATORY_GATE_PYTHON:-$repo/tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap/bin/python}
mkdir -p "$outputs" "$network_outputs" "$logs"
lock_file=$work/.regulatory-gate.lock
exec 9>"$lock_file"
if ! flock -n 9; then
  echo "Another regulatory-gate pipeline holds $lock_file; refusing concurrent output writes." >&2
  exit 73
fi
printf 'pid=%s started=%s\n' "$$" "$(date --iso-8601=seconds)" 1>&9
network_lock_file=$network_outputs/.regulatory-network.lock
exec 8>"$network_lock_file"
if ! flock -n 8; then
  echo "Another regulatory-gate pipeline holds $network_lock_file; refusing concurrent network reads/writes." >&2
  exit 73
fi
export PYTHONPATH=$code
if [[ ! -x "$python_bin" ]]; then
  echo "Configured analysis Python is not executable: $python_bin" >&2
  exit 2
fi

run_logged() {
  local name=$1
  shift
  local log=$logs/${name}.log
  echo "[$(date --iso-8601=seconds)] START $name" | tee -a "$log"
  "$@" 2>&1 | tee -a "$log"
  echo "[$(date --iso-8601=seconds)] END $name" | tee -a "$log"
}

audit_hash_args=()
if [[ -n "${REGULATORY_GATE_VERIFIED_H5AD_RECEIPT:-}" ]]; then
  audit_hash_args=(--verified-h5ad-receipt "$REGULATORY_GATE_VERIFIED_H5AD_RECEIPT")
fi

run_logged 00_input_audit "$python_bin" "$code/audit_inputs.py" \
  --h5ad "$h5ad" \
  --expected-h5ad-bytes 30875155333 \
  --expected-h5ad-sha256 718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6 \
  --aracne-repo "$aracne_repo" \
  --expected-aracne-commit 3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14 \
  --anchor-code-repo "$anchor_code_repo" \
  --expected-anchor-code-commit edf5314ce5b9ee0e2f88b2310e7c2df5619ad888 \
  --supplement "$supplement" \
  --manifest "$manifest" \
  --pipeline-code-dir "$code" \
  "${audit_hash_args[@]}" \
  --output "$outputs/input_audit.json"

run_logged 01_pango "$python_bin" "$code/extract_pango_regulators.py" \
  --xlsx "$supplement" --output-dir "$network_outputs/pango"

if [[ ! -s "$network_outputs/tcga/tcga_receipt.json" ]]; then
  run_logged 02_tcga Rscript "$code/prepare_tcga_crc.R" "$network_outputs/tcga"
fi

run_logged 03_aracne_inputs "$python_bin" "$code/prepare_aracne_inputs.py" \
  --expression "$network_outputs/tcga/tcga_crc_tpm.tsv" \
  --pango "$network_outputs/pango/pango_regulators.txt" \
  --output-dir "$network_outputs/aracne_inputs"

if [[ ! -s "$outputs/pseudobulk/pseudobulk_receipt.json" ]]; then
  run_logged 04_pseudobulk "$python_bin" "$code/build_pseudobulk.py" \
    --h5ad "$h5ad" --manifest "$manifest" --output-dir "$outputs/pseudobulk"
fi

if [[ ! -s "$network_outputs/aracne3_receipt.json" || ! -s "$network_outputs/aracne3/subnets/subnet1_crc.tsv" ]]; then
  run_logged 05_aracne bash "$code/run_aracne3.sh" \
    "$aracne_repo" "$network_outputs/tcga/tcga_crc_tpm.tsv" \
    "$network_outputs/aracne_inputs/aracne_regulators.txt" "$network_outputs/aracne3" "$threads" 1
fi

run_logged 06_aracne_audit "$python_bin" "$code/audit_aracne_run.py" \
  --repo "$aracne_repo" --expression "$network_outputs/tcga/tcga_crc_tpm.tsv" \
  --regulators "$network_outputs/aracne_inputs/aracne_regulators.txt" \
  --run-dir "$network_outputs/aracne3" --expected-subnetworks 1 \
  --threads "$threads" \
  --output "$network_outputs/aracne3_receipt.json"

run_logged 07_viper Rscript "$code/viper_msviper.R" \
  "$network_outputs/tcga/tcga_crc_tpm.rds" "$network_outputs/aracne3/subnets/subnet1_crc.tsv" \
  "$outputs/pseudobulk/log2cpm.tsv.gz" "$outputs/pseudobulk/patient_metadata.tsv" \
  "$outputs/viper" "$threads" 1000 "$manifest"

run_logged 08_statistics "$python_bin" "$code/statistics.py" \
  --activity "$outputs/viper/viper_activity.tsv.gz" \
  --metadata "$outputs/pseudobulk/patient_metadata.tsv" \
  --msviper "$outputs/viper/msviper_results.tsv" \
  --manifest "$manifest" --output-dir "$outputs/statistics" \
  --permutations 500 --bootstraps 1000

run_logged 09_decision "$python_bin" "$code/decide_gate.py" \
  --manifest "$manifest" \
  --input-audit "$outputs/input_audit.json" \
  --pango "$network_outputs/pango/pango_receipt.json" \
  --tcga "$network_outputs/tcga/tcga_receipt.json" \
  --aracne-input "$network_outputs/aracne_inputs/aracne_input_receipt.json" \
  --aracne-run "$network_outputs/aracne3_receipt.json" \
  --pseudobulk "$outputs/pseudobulk/pseudobulk_receipt.json" \
  --viper "$outputs/viper/viper_receipt.json" \
  --statistics "$outputs/statistics/statistics_receipt.json" \
  --output "$outputs/gate_decision.json"

echo "Decision written to $outputs/gate_decision.json"

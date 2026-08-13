#!/usr/bin/env bash
set -euo pipefail

repo=${PAPER2PAPER_REPO:-/media/desk16/iy13202/projects/Paper2Paper}
work=${REGULATORY_GATE_WORK:-$repo/tmp/tnbc-chromatin-tf-nc-2026/regulatory-gate-M-vs-rest}
h5ad=${CRC_ATLAS_H5AD:-$repo/tmp/hcc-sc-spatial-npj-2026/phase-zero/raw/crc_atlas_cellxgene.h5ad}
supplement=${TNBC_SUPPLEMENT_DATA2:-$repo/tmp/tnbc-chromatin-tf-nc-2026/external/source/41467_2026_76385_MOESM4_ESM.xlsx}
aracne_repo=${ARACNE3_REPO:-$repo/tmp/tnbc-chromatin-tf-nc-2026/external/code/ARACNe3}
code=$repo/pilots/tnbc-chromatin-tf-nc-2026/code/regulatory_gate
manifest=$repo/pilots/tnbc-chromatin-tf-nc-2026/analysis/phenotypes/m-vs-rest.json
outputs=$work/outputs
logs=$work/logs
threads=${REGULATORY_GATE_THREADS:-24}
python_bin=${REGULATORY_GATE_PYTHON:-$repo/tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap/bin/python}
mkdir -p "$outputs" "$logs"
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

run_logged 00_input_audit "$python_bin" "$code/audit_inputs.py" \
  --h5ad "$h5ad" \
  --expected-h5ad-bytes 30875155333 \
  --expected-h5ad-sha256 718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6 \
  --aracne-repo "$aracne_repo" \
  --expected-aracne-commit 3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14 \
  --supplement "$supplement" \
  --output "$outputs/input_audit.json"

run_logged 01_pango "$python_bin" "$code/extract_pango_regulators.py" \
  --xlsx "$supplement" --output-dir "$outputs/pango"

if [[ ! -s "$outputs/tcga/tcga_receipt.json" ]]; then
  run_logged 02_tcga Rscript "$code/prepare_tcga_crc.R" "$outputs/tcga"
fi

run_logged 03_aracne_inputs "$python_bin" "$code/prepare_aracne_inputs.py" \
  --expression "$outputs/tcga/tcga_crc_tpm.tsv" \
  --pango "$outputs/pango/pango_regulators.txt" \
  --output-dir "$outputs/aracne_inputs"

if [[ ! -s "$outputs/pseudobulk/pseudobulk_receipt.json" ]]; then
  run_logged 04_pseudobulk "$python_bin" "$code/build_pseudobulk.py" \
    --h5ad "$h5ad" --manifest "$manifest" --output-dir "$outputs/pseudobulk"
fi

run_logged 05_aracne "$code/run_aracne3.sh" \
  "$aracne_repo" "$outputs/tcga/tcga_crc_tpm.tsv" \
  "$outputs/aracne_inputs/aracne_regulators.txt" "$outputs/aracne3" "$threads" 100

run_logged 06_aracne_audit "$python_bin" "$code/audit_aracne_run.py" \
  --repo "$aracne_repo" --expression "$outputs/tcga/tcga_crc_tpm.tsv" \
  --regulators "$outputs/aracne_inputs/aracne_regulators.txt" \
  --run-dir "$outputs/aracne3" --expected-subnetworks 100 \
  --output "$outputs/aracne3_receipt.json"

run_logged 07_viper Rscript "$code/viper_msviper.R" \
  "$outputs/tcga/tcga_crc_tpm.rds" "$outputs/aracne3/consolidated-net_crc.tsv" \
  "$outputs/pseudobulk/log2cpm.tsv.gz" "$outputs/pseudobulk/patient_metadata.tsv" \
  "$outputs/viper" "$threads" 500

run_logged 08_statistics "$python_bin" "$code/statistics.py" \
  --activity "$outputs/viper/viper_activity.tsv.gz" \
  --metadata "$outputs/pseudobulk/patient_metadata.tsv" \
  --msviper "$outputs/viper/msviper_results.tsv" \
  --manifest "$manifest" --output-dir "$outputs/statistics" \
  --permutations 500 --bootstraps 1000

run_logged 09_decision "$python_bin" "$code/decide_gate.py" \
  --input-audit "$outputs/input_audit.json" \
  --pango "$outputs/pango/pango_receipt.json" \
  --tcga "$outputs/tcga/tcga_receipt.json" \
  --aracne-input "$outputs/aracne_inputs/aracne_input_receipt.json" \
  --aracne-run "$outputs/aracne3_receipt.json" \
  --pseudobulk "$outputs/pseudobulk/pseudobulk_receipt.json" \
  --viper "$outputs/viper/viper_receipt.json" \
  --statistics "$outputs/statistics/statistics_receipt.json" \
  --output "$outputs/gate_decision.json"

echo "Decision written to $outputs/gate_decision.json"

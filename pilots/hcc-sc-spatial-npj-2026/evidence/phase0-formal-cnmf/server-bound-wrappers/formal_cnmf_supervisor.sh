#!/usr/bin/env bash
set -euo pipefail

ROOT=/media/desk16/iy13202/projects/Paper2Paper/tmp/hcc-sc-spatial-npj-2026/phase-zero
MANIFEST=$ROOT/manifests/formal_cnmf_v1
CODE=$MANIFEST/code
BASE_CODE=/media/desk16/iy13202/projects/Paper2Paper/pilots/hcc-sc-spatial-npj-2026/code/phase_zero
OUT=$ROOT/outputs/formal_cnmf_v1
LOG=$ROOT/logs/formal_cnmf_v1
H5AD=$ROOT/raw/crc_atlas_cellxgene.h5ad
CNMF_PY=$ROOT/env/cnmf/bin/python
EXPECTED_H5AD_SHA256=718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6
RUN_NAME=crc_liu_formal_cnmf

export PYTHONPATH=$CODE:$BASE_CODE
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

if [[ -e "$OUT" || -e "$LOG" ]]; then
  printf 'refusing to mix with an existing formal cNMF run\n' >&2
  exit 2
fi
mkdir -p "$OUT" "$LOG"
date --iso-8601=seconds > "$OUT/STARTED"

(
  cd "$MANIFEST"
  sha256sum -c code_and_contract.sha256
) > "$LOG/code_hash_check.log" 2>&1
"$CNMF_PY" --version > "$MANIFEST/cnmf_python_version.txt" 2>&1
"$CNMF_PY" -m pip freeze > "$MANIFEST/cnmf_pip_freeze_runtime.txt"

"$CNMF_PY" -m unittest discover -s "$CODE" -p 'test_formal_cnmf_*.py' -v \
  > "$LOG/preflight_tests.log" 2>&1

actual_h5ad_sha256=$(sha256sum "$H5AD" | awk '{print toupper($1)}')
printf '%s  %s\n' "$actual_h5ad_sha256" "$H5AD" > "$MANIFEST/h5ad.sha256"
if [[ "$actual_h5ad_sha256" != "$EXPECTED_H5AD_SHA256" ]]; then
  printf 'H5AD SHA256 changed: %s\n' "$actual_h5ad_sha256" >&2
  exit 3
fi

evidence_args=()
for shard in 01 02 03 04; do
  evidence_args+=(
    --cna-evidence
    "$ROOT/outputs/cna_paired_full_v1/summary_${shard}/P0_cna_panel_cell_evidence.tsv"
  )
done

prepare_input() {
  local seed=$1
  local destination=$2
  "$CNMF_PY" "$CODE/prepare_crc_formal_cnmf_input.py" \
    --h5ad "$H5AD" \
    "${evidence_args[@]}" \
    --output-dir "$destination" \
    --h5ad-sha256 "$actual_h5ad_sha256" \
    --cohort-dataset Liu_2024_mixCD45PosCD45Neg \
    --analysis-role liu_primary_discovery \
    --expected-patients 11 \
    --expected-samples 25 \
    --expected-dual-cells 23766 \
    --max-cells-per-patient-state 500 \
    --seed "$seed"
}

prepare_input 20260812 "$OUT/input_main" > "$LOG/prepare_input_main.log" 2>&1
prepare_input 20260813 "$OUT/input_resample" > "$LOG/prepare_input_resample.log" 2>&1

"$CNMF_PY" - "$OUT/input_main/P0_formal_cnmf_input_receipt.json" \
  "$OUT/input_resample/P0_formal_cnmf_input_receipt.json" <<'PY'
import json
import sys
from pathlib import Path

expected_high = {
    "Liu_2024_mixCD45PosCD45Neg.Pt05_T1",
    "Liu_2024_mixCD45PosCD45Neg.Pt07_T2",
    "Liu_2024_mixCD45PosCD45Neg.Pt15_T1",
}
receipts = [json.loads(Path(value).read_text(encoding="utf-8")) for value in sys.argv[1:]]
for receipt in receipts:
    assert receipt["status"] == "formal_liu_discovery_input_frozen_pending_cnmf"
    assert receipt["available_dual_supported_cells"] == 23766
    assert receipt["selected_cells"] == 9370
    assert receipt["patients"] == 11
    assert receipt["samples"] == 25
    assert receipt["states"] == ["metastasis", "tumor"]
    assert set(receipt["high_reference_unresolved_samples_gt_0_5"]) == expected_high
assert receipts[0]["cell_sampling_seed"] == 20260812
assert receipts[1]["cell_sampling_seed"] == 20260813
assert (
    receipts[0]["outputs"]["P0_formal_cnmf_selected_cells.tsv"]["sha256"]
    != receipts[1]["outputs"]["P0_formal_cnmf_selected_cells.tsv"]["sha256"]
)
print("formal input receipts passed")
PY

"$CNMF_PY" "$CODE/split_crc_formal_patient_holdouts.py" \
  --input-h5ad "$OUT/input_main/crc_liu_dual_cna_patient_state_balanced_counts.h5ad" \
  --output-dir "$OUT/patient_holdouts" \
  --seed 20260812 > "$LOG/split_patient_holdouts.log" 2>&1

"$CNMF_PY" "$CODE/run_crc_cnmf.py" \
  --counts-h5ad "$ROOT/outputs/cnmf_cna_input_v1/crc_two_method_cna_supported_counts.h5ad" \
  --output-dir "$OUT/parallel_smoke" \
  --run-name crc_parallel_smoke \
  --components 5 \
  --n-iter 5 \
  --seed 20260812 \
  --total-workers 2 \
  --run-role diagnostic > "$LOG/parallel_smoke.log" 2>&1

run_one() {
  local role=$1
  local input=$2
  local destination=$3
  local nmf_seed=$4
  "$CNMF_PY" "$CODE/run_crc_cnmf.py" \
    --counts-h5ad "$input" \
    --output-dir "$destination" \
    --run-name "$RUN_NAME" \
    --components 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 \
    --n-iter 100 \
    --seed "$nmf_seed" \
    --density-threshold 0.1 \
    --num-highvar-genes 2000 \
    --max-nmf-iter 1000 \
    --total-workers 8 \
    --run-role "$role"
}

declare -a names=(
  seed_20260812
  seed_20260813
  cell_resample_20260813
  patient_holdout_A
  patient_holdout_B
)
declare -a roles=(
  liu_primary_seed
  liu_primary_seed
  liu_cell_resample
  liu_patient_holdout
  liu_patient_holdout
)
declare -a inputs=(
  "$OUT/input_main/crc_liu_dual_cna_patient_state_balanced_counts.h5ad"
  "$OUT/input_main/crc_liu_dual_cna_patient_state_balanced_counts.h5ad"
  "$OUT/input_resample/crc_liu_dual_cna_patient_state_balanced_counts.h5ad"
  "$OUT/patient_holdouts/crc_liu_patient_holdout_A_counts.h5ad"
  "$OUT/patient_holdouts/crc_liu_patient_holdout_B_counts.h5ad"
)
declare -a seeds=(20260812 20260813 20260812 20260812 20260812)
declare -a pids=()

mkdir -p "$OUT/runs"
for index in "${!names[@]}"; do
  name=${names[$index]}
  run_one "${roles[$index]}" "${inputs[$index]}" "$OUT/runs/$name" "${seeds[$index]}" \
    > "$LOG/$name.log" 2>&1 &
  pid=$!
  pids+=("$pid")
  printf '%s\n' "$pid" > "$LOG/$name.pid"
done

failed=0
for index in "${!names[@]}"; do
  name=${names[$index]}
  if wait "${pids[$index]}"; then
    printf '0\n' > "$LOG/$name.exit_code"
  else
    code=$?
    printf '%s\n' "$code" > "$LOG/$name.exit_code"
    failed=1
  fi
done
if [[ "$failed" -ne 0 ]]; then
  date --iso-8601=seconds > "$OUT/FAILED"
  exit 4
fi

"$CNMF_PY" "$CODE/evaluate_crc_cnmf_stability.py" \
  --run "seed_20260812=$OUT/runs/seed_20260812/$RUN_NAME" \
  --run "seed_20260813=$OUT/runs/seed_20260813/$RUN_NAME" \
  --run "cell_resample_20260813=$OUT/runs/cell_resample_20260813/$RUN_NAME" \
  --run "patient_holdout_A=$OUT/runs/patient_holdout_A/$RUN_NAME" \
  --run "patient_holdout_B=$OUT/runs/patient_holdout_B/$RUN_NAME" \
  --output-dir "$OUT/stability" \
  --top-n 50 \
  --permutations 100 \
  --seed 20260812 > "$LOG/evaluate_stability.log" 2>&1

"$CNMF_PY" - "$OUT" "$MANIFEST" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

out = Path(sys.argv[1])
manifest = Path(sys.argv[2])

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()

run_receipts = sorted(out.glob("runs/*/P0_cnmf_run_receipt.json"))
if len(run_receipts) != 5:
    raise RuntimeError(f"expected five run receipts, observed {len(run_receipts)}")
for path in run_receipts:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if receipt["status"] != "formal_liu_candidate_run_completed_pending_k_and_replication_review":
        raise RuntimeError(f"non-formal run receipt: {path}")
    if receipt["components"] != list(range(5, 21)) or receipt["n_iter"] != 100:
        raise RuntimeError(f"run contract mismatch: {path}")

stability = out / "stability" / "P0_cnmf_stability_receipt.json"
stability_receipt = json.loads(stability.read_text(encoding="utf-8"))
if stability_receipt["status"] != "candidate_seed_cell_resample_and_patient_holdout_comparison":
    raise RuntimeError("formal stability receipt lacks the required perturbations")

receipt = {
    "status": "formal_liu_cnmf_ready_for_k_review",
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "formal_cnmf_started": True,
    "formal_cnmf_completed": True,
    "program_names_frozen": False,
    "che_external_projection_started": False,
    "contract_sha256": sha256(manifest / "formal-cnmf-execution-contract.md"),
    "run_receipts": {str(path.relative_to(out)): sha256(path) for path in run_receipts},
    "stability_receipt": {
        "path": str(stability.relative_to(out)),
        "sha256": sha256(stability),
    },
    "claim_boundary": (
        "Five Liu candidate runs and their stability comparison completed. K, program identities, "
        "biological names and Che external replication remain pending human review."
    ),
}
(out / "formal_cnmf_execution_receipt.json").write_text(
    json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
PY

date --iso-8601=seconds > "$OUT/READY_FOR_K_REVIEW"

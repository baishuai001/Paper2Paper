#!/usr/bin/env bash
set -euo pipefail

ROOT=/media/desk16/iy13202/projects/Paper2Paper/tmp/hcc-sc-spatial-npj-2026/phase-zero
MANIFEST=$ROOT/manifests/formal_cnmf_sensitivity_v1
CODE=$MANIFEST/code
BASE_CODE=/media/desk16/iy13202/projects/Paper2Paper/pilots/hcc-sc-spatial-npj-2026/code/phase_zero
PRIMARY=$ROOT/outputs/formal_cnmf_v1
OUT=$ROOT/outputs/formal_cnmf_sensitivity_v1
LOG=$ROOT/logs/formal_cnmf_sensitivity_v1
H5AD=$ROOT/raw/crc_atlas_cellxgene.h5ad
CNMF_PY=$ROOT/env/cnmf/bin/python
EXPECTED_H5AD_SHA256=718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6
PRIMARY_RUN=crc_liu_formal_cnmf
SENSITIVITY_RUN=crc_liu_sensitivity_cnmf

export PYTHONPATH=$CODE:$BASE_CODE
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

if [[ -e "$OUT" || -e "$LOG" ]]; then
  printf 'refusing to mix with an existing cNMF sensitivity run\n' >&2
  exit 2
fi
if [[ ! -s "$PRIMARY/READY_FOR_K_REVIEW" || ! -s "$PRIMARY/formal_cnmf_execution_receipt.json" ]]; then
  printf 'primary cNMF run has not reached the frozen K-review point\n' >&2
  exit 3
fi
mkdir -p "$OUT" "$LOG"
date --iso-8601=seconds > "$OUT/STARTED"

(
  cd "$MANIFEST"
  sha256sum -c code_and_contract.sha256
) > "$LOG/code_hash_check.log" 2>&1
"$CNMF_PY" -m unittest discover -s "$CODE" -p 'test_formal_cnmf_*.py' -v \
  > "$LOG/preflight_tests.log" 2>&1

actual_h5ad_sha256=$(sha256sum "$H5AD" | awk '{print toupper($1)}')
printf '%s  %s\n' "$actual_h5ad_sha256" "$H5AD" > "$MANIFEST/h5ad.sha256"
if [[ "$actual_h5ad_sha256" != "$EXPECTED_H5AD_SHA256" ]]; then
  printf 'H5AD SHA256 changed: %s\n' "$actual_h5ad_sha256" >&2
  exit 4
fi

evidence_args=()
for shard in 01 02 03 04; do
  evidence_args+=(
    --cna-evidence
    "$ROOT/outputs/cna_paired_full_v1/summary_${shard}/P0_cna_panel_cell_evidence.tsv"
  )
done

"$CNMF_PY" "$CODE/prepare_crc_formal_cnmf_input.py" \
  --h5ad "$H5AD" \
  "${evidence_args[@]}" \
  --output-dir "$OUT/input_exclude_high_unresolved" \
  --h5ad-sha256 "$actual_h5ad_sha256" \
  --cohort-dataset Liu_2024_mixCD45PosCD45Neg \
  --analysis-role liu_reference_qc_exclusion_sensitivity \
  --expected-patients 11 \
  --expected-samples 22 \
  --expected-dual-cells 22338 \
  --max-cells-per-patient-state 500 \
  --seed 20260812 \
  --exclude-sample-id Liu_2024_mixCD45PosCD45Neg.Pt05_T1 \
  --exclude-sample-id Liu_2024_mixCD45PosCD45Neg.Pt07_T2 \
  --exclude-sample-id Liu_2024_mixCD45PosCD45Neg.Pt15_T1 \
  --allow-incomplete-patient-states \
  --expected-patient-state-units 19 \
  --input-filename crc_liu_exclude_high_unresolved_counts.h5ad \
  > "$LOG/prepare_exclude_high_unresolved.log" 2>&1

"$CNMF_PY" "$CODE/prepare_crc_formal_cnmf_input.py" \
  --h5ad "$H5AD" \
  "${evidence_args[@]}" \
  --output-dir "$OUT/input_include_one_method" \
  --h5ad-sha256 "$actual_h5ad_sha256" \
  --cohort-dataset Liu_2024_mixCD45PosCD45Neg \
  --analysis-role liu_one_method_expansion_sensitivity \
  --expected-patients 11 \
  --expected-samples 25 \
  --expected-dual-cells 27705 \
  --max-cells-per-patient-state 500 \
  --seed 20260812 \
  --selection-mode dual_plus_one_method \
  --expected-patient-state-units 22 \
  --input-filename crc_liu_dual_plus_one_method_counts.h5ad \
  > "$LOG/prepare_include_one_method.log" 2>&1

"$CNMF_PY" - "$OUT" <<'PY'
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
excluded = json.loads((out / "input_exclude_high_unresolved/P0_formal_cnmf_input_receipt.json").read_text())
expanded = json.loads((out / "input_include_one_method/P0_formal_cnmf_input_receipt.json").read_text())
expected_exclusions = {
    "Liu_2024_mixCD45PosCD45Neg.Pt05_T1",
    "Liu_2024_mixCD45PosCD45Neg.Pt07_T2",
    "Liu_2024_mixCD45PosCD45Neg.Pt15_T1",
}
assert excluded["status"] == "formal_liu_sensitivity_input_frozen_pending_cnmf"
assert excluded["available_selected_cells"] == 22338
assert excluded["selected_cells"] == 8067
assert excluded["patients"] == 11 and excluded["samples"] == 22
assert excluded["patient_state_units"] == 19
assert set(excluded["excluded_sample_ids"]) == expected_exclusions
assert excluded["retained_high_reference_unresolved_samples_gt_0_5"] == []
assert expanded["status"] == "formal_liu_sensitivity_input_frozen_pending_cnmf"
assert expanded["available_selected_cells"] == 27705
assert expanded["available_dual_supported_cells"] == 23766
assert expanded["available_one_method_supported_cells"] == 3939
assert expanded["selected_cells"] == 9893
assert expanded["patients"] == 11 and expanded["samples"] == 25
assert expanded["patient_state_units"] == 22
selected = out / "input_include_one_method/P0_formal_cnmf_selected_cells.tsv"
import pandas as pd
classes = pd.read_csv(selected, sep="\t")["cna_selection_class"].value_counts().to_dict()
assert classes == {
    "author_cancer_two_method_malignancy_support": 8308,
    "author_cancer_one_method_malignancy_support": 1585,
}
print("two frozen sensitivity inputs passed")
PY

run_one() {
  local role=$1
  local input=$2
  local destination=$3
  "$CNMF_PY" "$CODE/run_crc_cnmf.py" \
    --counts-h5ad "$input" \
    --output-dir "$destination" \
    --run-name "$SENSITIVITY_RUN" \
    --components 9 10 11 \
    --n-iter 100 \
    --seed 20260812 \
    --density-threshold 0.1 \
    --num-highvar-genes 2000 \
    --max-nmf-iter 1000 \
    --total-workers 8 \
    --run-role "$role"
}

mkdir -p "$OUT/runs"
run_one \
  liu_reference_qc_exclusion_sensitivity \
  "$OUT/input_exclude_high_unresolved/crc_liu_exclude_high_unresolved_counts.h5ad" \
  "$OUT/runs/exclude_high_unresolved" \
  > "$LOG/exclude_high_unresolved.log" 2>&1 &
pid_exclude=$!
printf '%s\n' "$pid_exclude" > "$LOG/exclude_high_unresolved.pid"

run_one \
  liu_one_method_expansion_sensitivity \
  "$OUT/input_include_one_method/crc_liu_dual_plus_one_method_counts.h5ad" \
  "$OUT/runs/include_one_method" \
  > "$LOG/include_one_method.log" 2>&1 &
pid_expand=$!
printf '%s\n' "$pid_expand" > "$LOG/include_one_method.pid"

failed=0
if wait "$pid_exclude"; then printf '0\n' > "$LOG/exclude_high_unresolved.exit_code"; else printf '%s\n' "$?" > "$LOG/exclude_high_unresolved.exit_code"; failed=1; fi
if wait "$pid_expand"; then printf '0\n' > "$LOG/include_one_method.exit_code"; else printf '%s\n' "$?" > "$LOG/include_one_method.exit_code"; failed=1; fi
if [[ "$failed" -ne 0 ]]; then
  date --iso-8601=seconds > "$OUT/FAILED"
  exit 5
fi

"$CNMF_PY" "$CODE/evaluate_crc_cnmf_stability.py" \
  --run "primary_seed_20260812=$PRIMARY/runs/seed_20260812/$PRIMARY_RUN" \
  --run "exclude_high_unresolved=$OUT/runs/exclude_high_unresolved/$SENSITIVITY_RUN" \
  --run "include_one_method=$OUT/runs/include_one_method/$SENSITIVITY_RUN" \
  --output-dir "$OUT/comparison" \
  --top-n 50 \
  --permutations 100 \
  --seed 20260812 > "$LOG/evaluate_sensitivities.log" 2>&1

"$CNMF_PY" - "$OUT" "$MANIFEST" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

out = Path(sys.argv[1])
manifest = Path(sys.argv[2])
matches = pd.read_csv(out / "comparison/P0_cnmf_stability.tsv", sep="\t")
direct = matches.loc[
    matches["comparison_type"].eq("same_k_cross_input")
    & (matches["left_k"].astype(int) == matches["right_k"].astype(int))
    & (
        matches["left_run"].eq("primary_seed_20260812")
        | matches["right_run"].eq("primary_seed_20260812")
    )
].copy()
direct["sensitivity"] = direct["left_run"].where(
    ~direct["left_run"].eq("primary_seed_20260812"), direct["right_run"]
)
direct["k"] = direct["left_k"].astype(int)
summary = direct.groupby(["sensitivity", "k"], as_index=False).agg(
    matched_pairs=("exceeds_both_nulls", "size"),
    pairs_exceeding_both_nulls=("exceeds_both_nulls", "sum"),
    median_cosine=("cosine", "median"),
    median_top_gene_jaccard=("top_gene_jaccard", "median"),
    minimum_cosine=("cosine", "min"),
    minimum_top_gene_jaccard=("top_gene_jaccard", "min"),
)
summary["pairs_exceeding_both_nulls"] = summary["pairs_exceeding_both_nulls"].astype(int)
summary["pass_fraction"] = summary["pairs_exceeding_both_nulls"] / summary["matched_pairs"]
if len(summary) != 6 or not (summary["matched_pairs"] == summary["k"]).all():
    raise RuntimeError("same-K primary-versus-sensitivity comparison is incomplete")
direct.to_csv(out / "P0_cnmf_sensitivity_same_k_pairs.tsv", sep="\t", index=False)
summary.to_csv(out / "P0_cnmf_sensitivity_by_k.tsv", sep="\t", index=False)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()

run_receipts = sorted(out.glob("runs/*/P0_cnmf_run_receipt.json"))
if len(run_receipts) != 2:
    raise RuntimeError("expected exactly two sensitivity run receipts")
for path in run_receipts:
    receipt = json.loads(path.read_text())
    assert receipt["status"] == "formal_liu_sensitivity_run_completed_pending_k_review"
    assert receipt["components"] == [9, 10, 11] and receipt["n_iter"] == 100
receipt = {
    "status": "formal_liu_cnmf_sensitivities_ready_for_human_review",
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "contract_sha256": sha256(manifest / "formal-cnmf-sensitivity-contract.md"),
    "program_names_frozen": False,
    "che_projection_started": False,
    "run_receipts": {str(p.relative_to(out)): sha256(p) for p in run_receipts},
    "summary_sha256": sha256(out / "P0_cnmf_sensitivity_by_k.tsv"),
    "claim_boundary": (
        "Both predeclared malignant-set sensitivities completed at K=9-11. Their direct same-K "
        "matches are ready for human review; K and program names remain unfrozen."
    ),
}
(out / "formal_cnmf_sensitivity_execution_receipt.json").write_text(
    json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
PY

date --iso-8601=seconds > "$OUT/READY_FOR_SENSITIVITY_REVIEW"

# Gastric NRRS real-data pilot

This pilot audits and executes a bounded public-data reproduction of the
eight-gene nerve-related risk score (NRRS) from the 2025 BMC Cancer anchor.
It is deliberately marked `route_role=training`: it is a completed worked
example for fixed bulk-signature reconstruction. It does not validate
Paper2Paper as a whole, unrelated modalities, or another anchor paper, and it
does not by itself define a non-duplicate paper.

## What is reproduced

- Cohort: GSE62254 / ACRG, 300 patients and 152 deaths.
- Expression: GEO-distributed RMA Series Matrix.
- Outcome: ACRG sheet in Supplementary Data 2 of DOI
  `10.1038/s41467-018-04179-8`.
- Probe mapping: official GPL570 annotation.
- Score: published eight coefficients, with a fully stated reconstructed
  transform, probe rule, missing-gene rule and median cutoff.

The paper does not disclose the original probe-collapsing rule, exact z-score
convention or cutoff. Therefore this is a transparent reconstruction, not an
exact reproduction of the authors' figure. GSE62254 also influenced the
anchor's feature-screening stage, so it is not untouched external validation.

## Run from the repository root

```powershell
python -m pip install -r pilots/gastric-nrrs/code/requirements-gse62254.txt

python pilots/gastric-nrrs/code/gse62254_nrrs_spike.py download `
  --manifest pilots/gastric-nrrs/config/gse62254_resources.tsv `
  --output-dir pilots/gastric-nrrs/raw_data/GSE62254

python -m unittest `
  pilots/gastric-nrrs/code/test_gse62254_nrrs_spike.py -v

python pilots/gastric-nrrs/code/gse62254_nrrs_spike.py analyze `
  --manifest pilots/gastric-nrrs/config/gse62254_resources.tsv `
  --data-dir pilots/gastric-nrrs/raw_data/GSE62254 `
  --output-dir pilots/gastric-nrrs/outputs/gse62254_spike
```

Downloads are excluded from Git. The manifest pins expected byte counts and
SHA256 checksums; a changed public file fails closed and requires re-audit.
The analysis also compares key results with
`config/gse62254_expected_results.json` and fails if the locked result changes.

## Auditable outputs

- `GSE62254_nrrs_source_table.tsv`: generated locally with one row per patient;
  it is intentionally not committed until source-data redistribution terms are
  fully audited. Its row count, columns, byte count and SHA256 are retained in
  `GSE62254_local_artifact_receipt.tsv`.
- `GSE62254_probe_mapping.tsv`: retained probes and sensitivity selection.
- `GSE62254_nrrs_summary.tsv`: effect estimates and discrimination.
- `GSE62254_nrrs_kaplan_meier.png`: reconstructed median-split survival plot.
- `GSE62254_run_metadata.json`: code hash, package versions, resource hashes,
  identity checks, model contract, results and claim boundary.
- `GSE62254_local_artifact_receipt.tsv`: repository-visible proof of the
  local-only patient table without redistributing patient-level rows.

See `reports/pilot-run-summary.md` for interpretation and
`reports/pilot-findings.md` for issues classified locally and, where justified,
linked to separately reviewed module/core promotions.

The MIT license covers Paper2Paper-authored code only. GEO files, the ACRG
supplement, GPL annotation, and derived patient-level records retain their
source-specific terms and attribution requirements.

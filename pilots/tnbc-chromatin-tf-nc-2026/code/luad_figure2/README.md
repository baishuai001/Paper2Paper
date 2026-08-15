# LUAD Figure 2 — TNBC-anchor-aligned execution

This directory transfers the anchor Figure 2 chain to all 158 LUAD-specific
TFs discovered in Figure 1. Figure 2 and Supplementary Figures 3 and 4A-C are
one atomic deliverable and share the same 22 patient / 13 PDX / 19 cell-line
sample manifest.

The primary rules are the anchor rules:

- a TF promoter is system-supported when it is accessible in at least
  `ceiling(n/2)` samples;
- a triple-promoter-supported TF is excluded from the HC-TF set only when its
  mean LUAD NES is below zero in patient, PDX, **and** cell-line activity
  cohorts;
- JASPAR and CIS-BP are tested independently with HOMER, using each complete
  motif reference rather than an HC-TF-only reduced database;
- adjusted `p < 1e-5` in at least `ceiling(n/2)` samples defines system-level
  motif support.

Read depth, peak count, FRiP, ataqv, Tn5-TSS enrichment and saturation are
reported as diagnostics. The former 1-million-read/10,000-peak filter, 20%
cohort failure rule, minimum HC/motif counts, and matched-permutation veto are
not part of the primary analysis and do not stop it.

Key entry points:

- `process_raw_atac_sample.sh` and `run_raw_atac_cloud.sh`: process all frozen
  PDX and cell-line raw libraries;
- `aggregate_raw_atac_qc.py`: verify processing completeness and emit the
  all-sample analysis manifest while retaining diagnostic metrics;
- `compute_promoter_gate_and_annotations.py`: apply the 158-TF promoter and
  all-three-negative NES logic;
- `prepare_hc_motif_catalog.py`, `run_homer_motif_gate.sh`, and
  `parse_homer_motif_gate.py`: run full-reference HOMER tests and summarize
  sample/system/three-system motif support;
- `run_figure2_post_raw_cloud.sh`: build Figure 2, Supplementary Figure 3,
  Supplementary Figure 4A-C, the descriptive report, and independent
  verification as one run.

Large data and computation remain on the cloud server. Git contains code,
manifests, compact result tables, figures, logs and receipts only. The earlier
`FAIL_DATA` artifacts are retained under `execution/luad-figure2` as historical
evidence of the superseded non-anchor protocol; new runs write the current
original-style result over the active compact output paths.

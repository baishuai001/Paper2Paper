# LUAD Figure 2 chromatin gate

This directory transfers the anchor Figure 2 chain to the 97 externally
replicated LUAD-vs-LUSC TFs frozen by Figure 1.  Figure 2 and Supplementary
Figures 3 and 4A-C are one atomic deliverable: they share the same sample
manifests, thresholds, TF checksum and run receipt.

Current manifest builders:

- `audit_tcga_luad_atac.R`: resolves 44 GDC technical-replicate count columns
  to 22 independent primary LUAD tumors and joins published ATAC QC.
- `audit_gse269746_atac.py`: resolves the 24 public ATAC PDX runs into 13 LUAD,
  4 transformed-SCLC and 7 de-novo-SCLC models; only the 13 LUAD models enter
  the primary gate.
- `fetch_ddbj_cellline_manifest.py`: resolves 28 DRA submission batches to 23
  unique cell lines and emits the one DMSO baseline library per canonical
  Dataset 2 model.
- `classify_ddbj_luad_celllines.py`: freezes 19 strict LUAD cell lines using
  DepMap/Cellosaurus disease and identity evidence.

The full pre-registered analysis, adaptation disclosures and stopping rules are
in `analysis/luad-figure2-gate-protocol.md`.  Large data and computation remain
on the cloud server; Git receives code, manifests, compact tables, figures,
logs and receipts only.

Execution helpers:

- `setup_figure2_references.sh`: freeze genome, annotation, blacklist and tool
  versions on the cloud server.
- `process_raw_atac_sample.sh` and `run_raw_atac_cloud.sh`: process the frozen
  PDX and strict-LUAD cell-line libraries.
- `aggregate_raw_atac_qc.py`: freeze sample-level QC, final artifact paths and
  the PDX/cell-line raw-data gate verdict.  If one fully processed system has
  already made `FAIL_DATA` irreversible, use
  `--early-stop-locked-system <PDX|cell_line>`; unattempted samples in the
  other system are then reported as `NOT_RUN_AFTER_LOCKED_DATA_GATE`, never as
  observed QC failures.  A started sample with an explicit early-stop receipt
  is reported separately as `ABORTED_AFTER_LOCKED_DATA_GATE`.
- `run_ataqv_raw_qc.sh`, `compute_tn5_tss_enrichment.py` and
  `parse_ataqv_metrics.py`: generate and freeze sample-level mapping,
  peak-overlap and TSS QC.  Because ataqv HQAA requires proper pairs, its TSS
  score is retained for paired PDX libraries while a common Tn5 insertion-site
  profile is reported for both layouts.  Neither metric alters the frozen hard
  QC thresholds.
- `run_figure2_post_raw_cloud.sh`: generate Figure 2, Supplementary Figure 3
  and Supplementary Figure 4A-C as one verified atomic result, then stop.  For
  an irreversible raw-data failure, set
  `EARLY_STOP_LOCKED_SYSTEM=cell_line` (or `PDX`); the same driver writes the
  terminal `FAIL_DATA` report and exits before biological-signal steps.

The HOMER parser reads the post-extraction target/background denominators from
`knownResults.txt` and retains GC-weighted decimal motif counts.  A pre-result
smoke test covers both frozen database naming conventions.  If the raw-data
gate fails, promoter/motif analysis and the atomic figure bundle are not run;
the data-failure report is the intended terminal artifact.
`verify_figure2_data_failure.py` independently verifies the full 19-line cell
audit, the early-stop accounting, the terminal verdict, and the required
absence of Figure 2/Supplementary Figure 3/4A-C signal artifacts.

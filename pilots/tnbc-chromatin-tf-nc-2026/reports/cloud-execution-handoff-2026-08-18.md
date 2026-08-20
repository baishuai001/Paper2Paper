# LUAD cloud execution handoff

Last updated: 2026-08-18 (Asia/Shanghai)

## Purpose

This file makes the current LUAD state/motif audit recoverable without the local
Windows computer or an active Codex session. All large downloads, intermediate
data, computation, receipts, and final tables live on the cloud server.

## Authoritative locations

- Repository: `/media/desk16/iy13202/projects/Paper2Paper`
- Figure 1 run: `/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/luad-figure1`
- Figure 2 run: `/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/luad-figure2`
- State audit: `/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit`
- Frozen protocol: `pilots/tnbc-chromatin-tf-nc-2026/analysis/luad-state-and-background-sensitivity-protocol.md`

No password or other credential is stored in this file or in the execution
scripts.

## Immutable result that must not be overwritten

The all-LUAD B0 analysis produced 31 HC-TFs, 20 motif-testable HC-TFs, and three
TFs with motif support in all three systems: FOXA3, NFATC4, and XBP1. B1-B3 and
TRU analyses are sensitivity/diagnostic layers and must not replace B0.

## Completed state work

- TCGA LUAD RNA classifier: 516 cases; TRU 186, PP 169, PI 161.
- RNA-ATAC matched TCGA subset: TRU 6, PP 8, PI 7.
- Frozen ATAC proxy validation: balanced accuracy 0.7667; 10,000-permutation
  one-sided p=0.0308; transfer to PDX was authorized.
- GSE269746 PDX ATAC proxy: TRU 5, PP 3, PI 5. The 13 PDX are not one state.

## Cloud jobs and automatic recovery

The cloud keeper checks every five minutes. If a job disappears before its
receipt is written, it restarts the resumable command without deleting completed
results:

```bash
cd /media/desk16/iy13202/projects/Paper2Paper
nohup pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2/monitor_luad_cloud_jobs.sh "$PWD" \
  > tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit/logs/cloud_keeper.log 2>&1 &
```

It protects:

1. DRA001846 download and Salmon quantification for the frozen 19-cell-line panel.
2. Automatic Wilkerson classifier execution after all 19 quantifications finish.
3. B1-B3 shared-background HOMER runs for JASPAR2024 and CIS-BP2.00, parsing,
   and final sensitivity summary.

A separate detached waiter starts the prespecified TRU diagnostic only after the
cell-state and B1-B3 receipts both exist. It recounts the immutable 31 HC-TFs;
it does not overwrite B0 or claim a new TRU-specific discovery.

The user crontab contains one idempotent `@reboot` entry invoking
`start_luad_cloud_pipeline.sh`, so a server reboot also restores the monitor and
TRU waiter. The launcher checks receipts and running processes before starting
anything.

The quantifier reuses completed `quant.sf` plus sample receipts and resumes SRA
downloads. The HOMER runner skips a sample only when both `.complete` and
`knownResults.txt` exist.

## One-command status check

```bash
cd /media/desk16/iy13202/projects/Paper2Paper
pilots/tnbc-chromatin-tf-nc-2026/code/luad_figure2/cloud_status_luad.sh "$PWD"
```

Primary logs:

- `tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit/logs/cloud_keeper.log`
- `tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit/logs/dra001846_quant_driver.log`
- `tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit/logs/dra001846_state_classifier_driver.log`
- `tmp/tnbc-chromatin-tf-nc-2026/luad-figure2/logs/motif_background_sensitivity/`
- `tmp/tnbc-chromatin-tf-nc-2026/luad-state-audit/logs/tru_state_diagnostic_driver.log`

Completion receipts:

- `luad-state-audit/audit/dra001846_quant/dra001846_quantification_receipt.json`
- `luad-state-audit/audit/state_classifier/state_classifier_receipt.json`
- `luad-figure2/audit/motif_background_sensitivity/final_receipt.json`
- `luad-state-audit/audit/tru_state_diagnostic/final_receipt.json`

## Scientific next step after both receipts exist

Use the frozen labels to report whether the 19 cell lines are state-balanced and
to define the prespecified NKX2-1-high/TRU-concordant subsets. Recount the
existing 31 HC-TFs and their motifs within TRU subsets in the three systems. This
is a diagnostic re-evaluation, not a new TRU-specific discovery analysis. A
formal TRU-specific discovery requires rerunning Figure 1 as TRU versus rest
under a separately frozen protocol.

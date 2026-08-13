# CRC ARACNe3/VIPER regulatory gate

This directory implements the pre-registered `M-vs-rest` entry gate in `analysis/regulatory-gate-m-vs-rest-protocol.md`.

The architecture deliberately separates reusable, phenotype-independent artifacts from a small phenotype manifest:

1. `prepare_tcga_crc.R` downloads GDC TCGA-COAD/READ primary-tumor STAR counts and prepares patient-level TPM.
2. `extract_pango_regulators.py` records the Methods claim of 2,139 PAN-GO genes, audits the 2,138 actual supplement rows and writes their 2,059 unique gene symbols.
3. `run_aracne3.sh` runs one fixed-seed ARACNe3 subnetwork, matching the anchor Code Ocean command's default `-x 1` while making its otherwise time-derived seed reproducible.
4. `build_pseudobulk.py` applies any compatible phenotype manifest and creates patient Cancer-cell pseudobulks.
5. `viper_msviper.R` builds the CRC regulon, computes patient VIPER NES and permutation-based msVIPER results.
6. `statistics.py` performs independent-study effects, REML meta-analysis, leakage-safe leave-one-study-out validation, bootstrap and within-study permutation. For phenotype portability, the selected validation column is normalized to the generic internal field `dataset`; the frozen M-vs-rest manifest maps that field from `study_id`.
7. `decide_gate.py` produces the fail-closed verdict.

All real data and full outputs stay on the cloud server. Git contains source, protocol, compact receipts and reports only. Run the complete cloud workflow with `bash run_cloud.sh`; paths may be overridden through the environment variables documented at the top of that script.

The anchor Code Ocean v1.0 code (commit `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`) explicitly converts `subnet1`, not the consolidated network, and uses `minsize=1` plus 1,000 null permutations. This implementation follows those choices, retains the first real edge that the author script accidentally deletes after `header=TRUE`, and adds no second consensus BH filter.

`viper_msviper.R` rejects fewer than 50 permutations. In `viper` 1.38.0, the empirical-tail helper used by msVIPER can fail to terminate when the null has too few distinct values; the frozen real run uses the anchor's 1,000 permutations.

To switch phenotype after a failed route, add a new JSON manifest with the same schema. The TCGA expression, ARACNe3 network, converted regulon and patient-level VIPER matrix are reusable; only cohort selection, msVIPER contrast, meta-analysis, LODO and decision must be rerun under a newly frozen route-specific protocol.

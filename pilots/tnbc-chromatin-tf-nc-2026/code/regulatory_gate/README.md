# CRC ARACNe3/VIPER regulatory gate

This directory implements the pre-registered `M-vs-rest` entry gate in `analysis/regulatory-gate-m-vs-rest-protocol.md`.

The architecture deliberately separates reusable, phenotype-independent artifacts from a small phenotype manifest:

1. `prepare_tcga_crc.R` downloads GDC TCGA-COAD/READ primary-tumor STAR counts and prepares patient-level TPM.
2. `extract_pango_regulators.py` records the Methods claim of 2,139 PAN-GO genes, audits the 2,138 actual supplement rows and writes their 2,059 unique gene symbols.
3. `run_aracne3.sh` runs the fixed official ARACNe3 commit with 100 subnetworks.
4. `build_pseudobulk.py` applies any compatible phenotype manifest and creates patient Cancer-cell pseudobulks.
5. `viper_msviper.R` builds the CRC regulon, computes patient VIPER NES and permutation-based msVIPER results.
6. `statistics.py` performs dataset effects, REML meta-analysis, leakage-safe LODO, bootstrap and within-dataset permutation.
7. `decide_gate.py` produces the fail-closed verdict.

All real data and full outputs stay on the cloud server. Git contains source, protocol, compact receipts and reports only. Run the complete cloud workflow with `bash run_cloud.sh`; paths may be overridden through the environment variables documented at the top of that script.

To switch phenotype after a failed route, add a new JSON manifest with the same schema. The TCGA expression, ARACNe3 network, converted regulon and patient-level VIPER matrix are reusable; only cohort selection, msVIPER contrast, meta-analysis, LODO and decision must be rerun under a newly frozen route-specific protocol.

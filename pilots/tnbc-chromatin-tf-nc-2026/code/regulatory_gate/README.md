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

Per the user-defined stop rule, the hard scientific verdict is based on the cross-study reproducible TF program (at least 10 reproducible TFs, at least 5 directionally confirmed by msVIPER). Leave-one-study-out classification is always run and reported as a secondary generalization diagnostic, but is not allowed to veto an otherwise present TF program.

All real data and full outputs stay on the cloud server. Git contains source, protocol, compact receipts and reports only. Run the complete cloud workflow with `bash run_cloud.sh`. `REGULATORY_GATE_MANIFEST` and `REGULATORY_GATE_WORK` select a frozen phenotype and its isolated result directory; `REGULATORY_GATE_SHARED_NETWORK_OUTPUTS` can point later phenotypes at the already-built TCGA/PAN-GO/ARACNe3 artifacts.

`run_cloud.sh` takes non-blocking locks for both the route work directory and the shared network directory, and exits with code 73 on contention. This prevents two launchers from corrupting ARACNe3 or route-specific outputs.

The anchor Code Ocean v1.0 code (commit `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`) explicitly converts `subnet1`, not the consolidated network, and uses `minsize=1` plus 1,000 null permutations. This implementation follows those choices, retains the first real edge that the author script accidentally deletes after `header=TRUE`, and adds no second consensus BH filter.

The anchor TCGA helper averages expression rows mapping to a duplicate gene symbol and applies no low-expression filter. The CRC implementation does the same, averages repeat aliquots at participant level, and uses the `viper` default method (`"none"`) because the anchor call omits `method`.

`viper_msviper.R` rejects fewer than 50 permutations. In `viper` 1.38.0, the empirical-tail helper used by msVIPER can fail to terminate when the null has too few distinct values; the frozen real run uses the anchor's 1,000 permutations.

Within-study Welch signatures use an algebraically equivalent vectorized sample-variance formula. This avoids thousands of slow row-wise `apply(var)` calls during the 1,000 frozen label permutations without changing the statistic.

To switch phenotype after a failed route, add a new JSON manifest with the same schema and use a new work directory. The TCGA expression and ARACNe3 network are reused through `REGULATORY_GATE_SHARED_NETWORK_OUTPUTS`; patient selection, VIPER/msVIPER, meta-analysis, LOSO and the manifest-derived decision route are rerun without editing source code.

# Pilot outcome

## Paper-side outcome

The NPJ HCC framework can be imitated by changing cancer type only by replacing
its **evidence roles**, not by renaming HCC in the title. No target paper has
yet been selected and no route has been promoted.

The user-provided CRC-atlas materially changed the evidence balance. It is now
the strongest candidate for the next bounded real-data experiment because this
audit verified:

- a public 30.9 GB processed H5AD with 3,790,266 cells, 588 donors and 1,553
  samples;
- a parsed 650-patient/1,670-sample Table S1 with normal, polyp, primary,
  liver-metastasis, treatment, response, RECIST and genotype fields;
- public Xenium/IMC files with a selective first-pass download plan;
- a substantial author code repository at a fixed commit.

CRC-atlas does not complete the target paper. It lacks a ready malignant
epithelial cNMF module, whole-transcriptome primary-to-CRLM spatial validation,
a single coherent treatment cohort and Geneformer. Its central integrated
neutrophil/spatial claims are already published and cannot be presented as new.
The decision is therefore **continue CRC into a small real-data/code precheck,
not select or promote it yet**.

PDAC and breast cancer are next for file-level audit. HGSC is currently more
valuable as a high-quality spatial/perturbation code donor than as a generic
target-cancer substitution. ccRCC and LUAD remain treatment- and progression-
module candidates whose public/controlled access and bridge cohorts require
further checking.

## Workflow-side outcome

The audit showed four concrete Paper2Paper capability gaps: large public
resources need selective-download and memory plans; targeted spatial panels
need feature-coverage receipts; a large GitHub repository can still lack its
advertised entrypoint and portable paths; heterogeneous treatment labels need
a cohort/regimen-specific estimand contract.

These findings were recorded locally in this Pilot. No Paper2Paper core or
general test suite was enlarged during this audit. Product changes remain
deferred until the real CRC precheck demonstrates exactly which guard prevents
a false feasibility decision.

## Reuse boundary

This Pilot parsed the local CRC-atlas supplementary Table S1 and gene-program
table, checked the CELLxGENE/Zenodo/BioImage resource metadata, downloaded and
checksummed only the 887-byte Xenium sample sheet and 205,363-byte gene panel,
and inspected the author repository at commit
`82c15ecc2ea36baa0c4cffe8818bf58e21dcfc48`.

It did **not** download or execute the 30.9 GB H5AD, rebuild the atlas, parse
GSE225857, run cNMF/spatial/Geneformer, or reproduce a CRC Figure. Therefore it
does not validate the CRC route, unrelated cancer routes or Paper2Paper's full
execution capability.

## Human decision and next boundary

No human promotion decision has been requested or recorded. All three routes
remain candidates at `direction_audited`; `selected_route_id` remains empty.

The next decision boundary is a server-side, bounded CRC experiment: capacity
receipt, backed H5AD read, donor-balanced malignant subset, GSE225857 structure
parse, one real cNMF stability spike, one coherent treatment-cohort contract and
one Geneformer smoke test. Only those results can support selection, refinement,
reduction to a partial imitation, another cancer search, or stopping.

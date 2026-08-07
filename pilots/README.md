# Pilot projects

These workspaces test the generic Paper2Paper workflow against real anchor papers.
They are deliberately isolated from the product core: evidence, decisions, code,
data, and results from one paper must not become defaults for another paper.

A pilot may verify a narrowly defined reusable control, such as checksum
enforcement or a cohort-role conflict check. It does not validate unrelated
modalities, research questions, formulas, expected results, or the workflow as
a whole. Reuse across pilots requires a new compatibility and evidence check.

- `spp1-tam-jitc/`: colorectal cancer liver-metastasis and SPP1+ TAM anchor.
- `gastric-nrrs/`: gastric-cancer nerve-related risk-signature anchor.

Both pilots begin at `anchor_audit`. Previous route recommendations were not
imported. They must be reconstructed from traceable evidence under the current
workflow.

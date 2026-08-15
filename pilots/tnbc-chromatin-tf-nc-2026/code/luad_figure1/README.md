# LUAD Figure 1 transfer

This directory implements the TNBC anchor's Figure 1 regulatory chain for a
pre-frozen **LUAD versus LUSC** comparison:

1. TCGA primary-tumor TPM -> PAN-GO regulators -> ARACNe3.
2. TNBC-compatible `rowTtest` signature, 1,000-permutation `ttestNull`,
   `msviper(minsize=1)`, and per-sample `viper(minsize=1)`.
3. Independent patient replication with GSE81089 and its own ARACNe3 network.
4. Projection of the frozen TCGA regulon to NCI PDMR PDX and DepMap 22Q2 cell
   lines.
5. Figure 1A-E panels, source tables, checksums, receipts, and a hard
   biological verdict.

`render_anchor_style_figure1.R` is a render-only layer for the publication
layout. It recreates Figure 1A and Figure 1C in the anchor paper's visual
grammar and emits separate Supplementary Figure 1A/1B cohort-flow panels. It
reads the frozen activity matrix and manifests, checks all displayed counts,
and does not rewrite upstream numerical results.

All large downloads and computation are designed to run on the project cloud
server. The top-level entry point is `run_luad_figure1_cloud.sh`; it will stop
if any frozen input is missing or any upstream command fails.

For `viper` 1.38.0, `run_viper_anchor.R` contains a narrowly scoped aREA
compatibility patch that restores the expected `1 x n` matrix shape for a
single-target regulon. It preserves `minsize=1` and all numeric definitions;
the triggering failed attempt and the patch are both disclosed in the audit.

`CollecTRI`, ULM, and ATAC are deliberately absent: they are not part of the
anchor Figure 1 chain. ATAC begins at Figure 2.

`render_academic_workflows.R` is the final publication-style rendering pass for
the same three workflow panels. It replaces large enclosing frames and decision
diamonds with an open evidence chain and a CONSORT-like cohort-selection spine.

Design references used for this render-only layer:
- BioRender scientific workflow templates: https://www.biorender.com/templates
- Bioicons open SVG library: https://github.com/duerrsimon/bioicons
- ggconsort cohort-flow grammar: https://github.com/tgerke/ggconsort

No third-party icon file is embedded in the output, so the figures remain fully

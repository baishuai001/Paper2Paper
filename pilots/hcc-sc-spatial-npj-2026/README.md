# HCC single-cell/spatial cancer-switch pilot

This Pilot is the first direct test of whether Paper2Paper can take a long,
multi-scale cancer paper and **change the cancer type without pretending that a
title-level substitution is an executable study**.

The anchor combines six layers: multi-cohort single-cell integration,
malignant-cell meta-programs, immune/stromal ecology, spatial localization,
in-silico perturbation, and clinical/therapeutic validation. A target cancer
must support the same causal and evidential roles, or the target paper must
explicitly state which part of the framework is being omitted.

Current conclusion: a cancer switch is scientifically possible, but no route is
yet qualified as executable. Colorectal cancer is now the leading route for a
**decisive real-data precheck**, because the CRC-atlas audit verified a public
3.79-million-cell processed object, patient/sample metadata, author code and a
selectively downloadable spatial layer. This is not a claim that CRC is already
the best target paper: the atlas cannot directly replace malignant epithelial
cNMF, whole-transcriptome primary-to-liver spatial validation, a homogeneous
treatment cohort or Geneformer. Gastric adenocarcinoma remains a comparator;
ICC remains a valid simple substitution but lacks a complete public role chain.

The CRC conclusion is now resolved panel by panel rather than by whole Figure:
of 92 main panels, 48 have data capable of filling the same evidence role but
remain unrun, 38 are modular only, one must be redesigned, and five are blocked.
The 63 supplementary panels are also mapped: 41 unrun/full-role, 18 modular,
one redesign, and three blocked. The blocked panels depend on a qualifying
treatment-response spatial cohort.

Key deliverables:

- `anchor/audit.md`: anchor framework, public reproduction boundary and flaws;
- `anchor/figure-map.tsv`: one row per main-figure evidence role;
- `anchor/main-figure-panel-map.tsv`: one row per main-text panel, with CRC
  data, replacement level, code source, repair and verification gate;
- `anchor/supplementary-figure-panel-map.tsv`: the same audit for all 63
  supplementary panels;
- `anchor/defect-repairs.tsv`: flaws converted into target-paper repair tests;
- `reports/cancer-switch-assessment.md`: conditional feasibility judgment;
- `reports/crc-atlas-replacement-audit.md`: file-, metadata-, code- and
  availability-level CRC audit;
- `reports/crc-subfigure-replacement-audit.md`: 92-panel data replacement
  judgment;
- `reports/crc-supplementary-panel-audit.md`: 63-panel supplementary evidence
  judgment;
- `reports/crc-code-reconstruction-map.md`: author, official and related-paper
  code reconstruction boundaries;
- `reports/other-cancer-data-screen.md`: preliminary PDAC, breast, HGSC,
  ccRCC and LUAD resource screen;
- `reports/paper2paper-gap-assessment.md`: what this real paper shows the
  workflow can and cannot currently do.

Raw data, PDFs, model weights and third-party code remain outside Git. No
candidate route may be promoted until a representative real-data run closes
the decisive data and code gaps.

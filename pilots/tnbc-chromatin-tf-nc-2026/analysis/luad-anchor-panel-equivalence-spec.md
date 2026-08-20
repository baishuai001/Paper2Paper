# LUAD publication rebuild: anchor-panel equivalence specification

## Purpose

This specification freezes the second-pass publication build. It separates three questions that were blurred in the first pass:

1. **Analytical equivalence** — whether the LUAD panel uses the same type of input, statistic, and selection rule as the TNBC anchor.
2. **Biological adaptation** — whether the breast-specific contrast or covariate has been replaced by a defensible lung-cancer analogue.
3. **Visual equivalence** — whether the panel communicates the same evidence hierarchy without copying the anchor artwork.

The anchor is used as a scientific and editorial reference, not as a source of values or graphical assets. Anchor panels in the comparison atlas are unmodified crops used only for critique.

## Frozen LUAD evidence chain

- Figure 1 starts from 158 TCGA LUAD-specific TFs. Independent patient evidence is reported separately for GSE81089 (95/158), GSE41271 (70/158), and the intersection of both external cohorts (44/158).
- Figure 2 starts from all 158 discovery TFs and applies the anchor promoter-accessibility and mean-activity rules, yielding 31 HC-TFs. Motif evidence is an annotation layer: 20 HC-TFs are motif-testable; FOXA3, NFATC4, and XBP1 meet the strict three-system rule; ETV1 and ZNF75D join them only at q <= 0.05 with support in at least 50% of samples.
- Figures 3–5 use all 31 HC-TFs, matching the anchor's use of all 94 HC-TFs rather than only its three-system motif subset.
- Figure 4 and Figure 5 remain explicitly exploratory where cross-cohort prognosis or in-vivo drug validation is weak. Layout improvements must not imply stronger evidence.

## Canonical supplementary numbering

| Canonical item | Anchor role | LUAD rebuild |
|---|---|---|
| Supplementary Figure 1 | Patient-cohort inclusion | TCGA, GSE81089, and GSE41271 inclusion flows |
| Supplementary Figure 2 | Independent TF discovery, overlap, activity coherence, per-sample distributions | GSE81089 and GSE41271 discovery panels, three-way overlap, TCGA/PDX/cell-line coherence, cross-system distributions |
| Supplementary Figure 3 | ATAC peak counts, saturation, correlations, genomic annotations | Patient/PDX/cell-line ATAC QC with asymptotic saturation and 50/90/95/99% markers |
| Supplementary Figure 4 | HC-TF selection, activity categories, motif availability, MKI67 | 31-HC-TF selection plus four-system MKI67 correlations |
| Supplementary Figure 5 | Cross-system heterogeneity and TME confounding | PDX/cell-line skewness plus TCGA/GSE41271 ESTIMATE correlations |
| Supplementary Figure 6 | Uni-/multivariable survival and permutation nulls | Eight LUAD Cox volcano panels and eight 5,000-permutation null panels |
| Supplementary Figure 7 | Drug coverage, examples, expanded associations | Drug overlap, replicated-pair examples, pathway/TF matrix, representative scatter |
| Extended Supplementary Figure 8 | LUAD-specific method audit | Cell-line mapping, drug eligibility, and classifier robustness |
| Extended Supplementary Figure 9 | In-vivo validation audit | PDX drug coverage and the complete negative association screen |

## Supplementary data packages

Nine supplementary data groups mirror the anchor's data roles while retaining LUAD-specific audit fields:

1. Cohort manifests, histology labels, clinical annotations, and inclusion counts.
2. TF catalogues, ARACNe/VIPER network summaries, differential activity, and sample-correlation matrices.
3. The 158-TF discovery set, 31 HC-TFs, motif-support categories, and signed regulon network.
4. ATAC sample manifest, raw-QC summaries, peak counts, saturation fits, pairwise correlations, and genomic annotations.
5. Promoter accessibility, motif availability, JASPAR-priority/CIS-BP-fallback mapping, strict and sensitivity motif results.
6. Uni- and multivariable Cox models, Kaplan–Meier receipts, endpoint overlaps, and permutation summaries.
7. Pharmacogenomic dataset audit, drug identity mapping, all concordance tests, replicated pairs, and pathway annotations.
8. PDX inventory, classifier scores, drug-response coverage, all PDX association tests, and validated-pair table (including zero-row result).
9. Cell-line and PDX identifier provenance, tissue/histology eligibility, mapping method, and authentication/provenance fields available from the source resources.

## Visual rules

- Use a consistent LUAD palette across all figures: LUAD blue, LUSC rust, PDX olive, cell-line teal, risk red, favorable blue, non-significant gray.
- Preserve panel-specific scales and always show sample size, test, correction method, and evidence status in the panel or legend.
- Heatmaps use symmetric, quantile-clipped activity scales and an explicit neutral midpoint. Large saturated red blocks from non-centered scales are prohibited.
- UpSet plots must show set membership and intersection size; bar charts that merely resemble UpSet plots are not acceptable.
- Network panels must encode node meaning, edge meaning, community, and degree/weight in legends. GO dot plots may accompany but not replace the target-gene community network.
- Negative panels use structured evidence diagrams or complete null distributions; blank boxes or placeholder prose are prohibited.
- Main figures are single-page, landscape, and readable at journal-column reduction. Dense supplementary figures may span multiple pages, but panel lettering and legends remain stable.

## Statistical rules retained from the anchor or frozen LUAD protocols

- Figure 1 differential TF activity: BH-adjusted inference; effect and activity axes shown together.
- Figure 2 promoter accessibility: at least half of samples in each system; mean NES exclusion only when all three systems are negative.
- Figure 2 motif: HOMER adjusted p/q < 1e-5 and at least half of samples in each system for the primary result; JASPAR first, CIS-BP only when JASPAR is unavailable. q <= 0.05 with at least half of samples is a labeled sensitivity analysis.
- Figure 3 skewness: normal-approximation z test with BH correction.
- Supplementary Figure 4D and 5B correlations: Pearson correlation, BH FDR, with the anchor's |r| > 0.4 display threshold reported separately from significance.
- Figure 4: the frozen LUAD primary Cox screen dichotomizes TF activity at the cohort/endpoint median; the displayed Kaplan-Meier cutpoints use maximally selected rank statistics and are descriptive. This is disclosed as an adaptation rather than described as identical to the anchor code, which uses outer tertiles for the METABRIC `NewSelection` analysis and maximally selected cutpoints for TCGA. Both nominal p and BH FDR are reported. Five-thousand-permutation enrichment tests are retained; weak empirical support is not hidden.
- Figure 5: within-dataset BH FDR <= 0.05 and consistent direction in at least two datasets define a replicated drug–TF pair. PDX validation is reported as absent when the required matched data do not exist or no pair passes.

## Official-code audit cautions

The Code Ocean capsule is tagged `v1.0` at commit `edf5314`. Its scientific workflow is used as the anchor, but two verified variable-selection defects are not propagated:

- Figure 1 Supplementary Figure 2F constructs a PDX long table from the TCGA object.
- Figure 1 cell-line PAM50 Fisher output writes the PDX confusion matrix object.

The capsule's run script also relies heavily on precomputed inputs and leaves several upstream steps commented out. The LUAD rebuild therefore records every generated table, input checksum, and executable command rather than claiming raw-to-figure equivalence solely from the capsule.

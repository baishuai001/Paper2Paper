# LUAD Figure 3 protocol: HC-TF regulon architecture and heterogeneity

## Frozen scientific input

Figure 3 starts from all TFs retained by the anchor-equivalent Figure 2
promoter/activity gate.  It does **not** select only TFs with a three-system
motif.  The expected frozen set contains 31 LUAD HC-TFs.  Motif support is an
annotation on this set, not an eligibility filter.

The TCGA ARACNe3 regulon and per-sample VIPER activity are the same objects used
for Figure 1.  Only TCGA primary LUAD samples are used to infer LUAD regulatory
architecture and heterogeneity; LUSC is used solely as a reference in the
representative density panel.

## Anchor-aligned analyses

1. Extract every signed TF-target edge for the 31 HC-TFs.  `tfmode > 0` is an
   activating edge and `tfmode < 0` is a repressing edge.  A target shared by at
   least two HC-TFs is labelled shared.
2. Report per-TF activated/repressed and shared/private target counts.  Annotate
   the three-system motif subset separately.
3. Build a complete TF collaboration graph from pairwise overlaps among
   activated targets.  Louvain communities use all positive-overlap edges and
   the shared-target count as the edge weight.  A high-weight subset is used
   only to keep the network panel legible; no edge is removed from the numeric
   community analysis.
4. For each community, run GO Biological Process over-representation analysis
   on its union of activated targets, using all HC-TF regulon targets as the
   tested universe and BH correction.
5. Calculate Pearson correlations between all HC-TF activities in TCGA-LUAD.
6. Calculate per-TF mean activity and sample skewness in TCGA-LUAD.  As in the
   anchor, use `z = skewness / sqrt(6/n)`, a two-sided normal p-value, BH
   correction, and call strong skew when FDR <= 0.05 and |skewness| > 1.
7. Display representative positively and negatively skewed TFs in LUAD and the
   LUSC reference samples.

## Supplementary heterogeneity analyses

- Correlate each HC-TF activity with MKI67 expression in TCGA-LUAD.
- Recalculate HC-TF skewness in PDMR LUAD PDXs and DepMap LUAD cell lines and
  compare directions with TCGA-LUAD.  These model-system analyses are context
  annotations and do not redefine the 31 HC-TFs.

## Integrity rules

- Exactly 31 Figure 2 HC-TFs must be recovered and all must occur in the TCGA
  regulon and activity matrix.
- Numeric results use all eligible observations; display-only edge thinning is
  labelled and recorded.
- No minimum number of communities, enriched GO terms, skewed TFs, or replicated
  model-system skew directions is imposed as a project-stopping threshold.
- Every table, figure, input hash, software version and independent verification
  result is retained in the Figure 3 execution directory.


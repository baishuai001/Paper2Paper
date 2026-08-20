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
2. Report per-TF activated/repressed and shared/private target counts.  In the
   Figure 3A motif tile, follow the anchor code and mark an HC-TF when its motif
   passes the original enrichment/prevalence rule in at least one of the three
   systems.  Keep the stricter three-system motif subset as a separate Figure 2
   result and supplementary annotation; it does not define Figure 3 entry.
3. Reproduce the anchor Figure 3B target-gene co-regulation network from
   activated shared targets.  An activated shared target is represented in at
   least two HC-TF regulons; an edge is retained only when a pair of target
   genes shares at least three activating HC-TFs.  As in the official
   `Figure3/05-Generate_Regulon_Network.R` script, the plotted graph contains
   only endpoints of qualifying edges.  Louvain communities and node degree
   are calculated on that graph.  Targets without a qualifying edge are
   retained in tables but are not added to the network as artificial isolates.
   GO Biological Process labels are assigned only to connected communities
   containing enough genes for a valid enrichment test; the frozen LUAD graph
   contains three two-gene components and therefore receives no GO label.
4. Reproduce the anchor Figure 3D HC-TF collaboration network from every TF
   pair sharing at least one activated target.  All positive-overlap edges are
   displayed: the 80th-percentile display thinning used in the first LUAD pass
   is removed.  Node size is TF-network degree and edge width is the number of
   shared activated targets, matching the official code and the interpretation
   used in Figure 3E.
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
- Figure 3B uses the frozen minimum of three shared activating HC-TFs per
  target-gene pair and reports the resulting 6-node/3-edge graph without
  lowering the rule for appearance.
- Figure 3D displays all positive TF-TF overlaps (31 nodes and 134 edges); no
  display-only edge thinning or induced isolate is permitted.
- No minimum number of communities, enriched GO terms, skewed TFs, or replicated
  model-system skew directions is imposed as a project-stopping threshold.
- Every table, figure, input hash, software version and independent verification
  result is retained in the Figure 3 execution directory.


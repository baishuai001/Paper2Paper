# LUAD Figure 5 frozen protocol

## Aim

Test whether the 31 chromatin-accessibility-prioritized LUAD HC-TFs predict
pharmacological response.  This transfers the anchor paper's Figure 5 logic,
while replacing the unavailable GRAY breast panel with PRISM 2020.

## Cell-line pharmacogenomics (Figure 5A-F)

- Datasets: GDSC2, CTRPv2 and PRISM, each analysed separately.
- TF activity: the frozen TCGA-LUAD/LUSC ARACNe3 regulon projected to DepMap
  22Q2 RNA-seq.  Only the 76 histologically annotated LUAD lines and the 31
  Figure 2 HC-TFs are eligible.
- Cell identity is joined in the following priority order: DepMap ID,
  Cellosaurus accession, COSMIC ID, then an unambiguous normalized cell-line
  name.  Ambiguous or unmatched records are excluded and reported.
- Sensitivity is AAC.  Recomputed AAC is preferred, followed by published AAC;
  if only AUC exists, AAC is `1-AUC`.  Duplicate cell-drug experiments are
  collapsed by the median.
- A drug enters testing only when at least 10 LUAD lines have paired activity
  and AAC and at least one paired line has AAC > 0.2, matching the anchor's
  requirement that the drug show sensitivity in the tested lineage.
- For every eligible drug and HC-TF, use `survcomp::concordance.index` with
  `surv.time=1-AAC`, all events set to one, and the Noether variance estimator.
  Thus CI > 0.5 denotes higher TF activity associated with greater sensitivity.
- Benjamini-Hochberg correction is applied across all drug-TF pairs within each
  dataset. FDR <= 0.05 is significant.
- Cross-dataset replication requires the same canonical compound, the same TF,
  FDR <= 0.05 in at least two datasets, and the same CI direction relative to
  0.5.  Compound identity uses PubChem CID first, InChIKey second, and a
  normalized drug name only when structural identifiers are unavailable.

## In-vivo PDX analysis (Figure 5G-H)

- Dataset: canonical ORCESTRA PDXE v2 (DOI 10.5281/zenodo.19354439), using
  baseline RNA-seq and model-level Best Average Response.
- The frozen TCGA regulon is projected to PDXE baseline RNA-seq.  Because PDXE
  labels lung models as NSCLC rather than LUAD/LUSC, a frozen LUAD-versus-LUSC
  TF-activity classifier trained only in TCGA is applied; PDXE models are
  reported as *LUAD-like*, not as pathology-confirmed LUAD. The frozen feature
  sets are the 158 TCGA LUAD-specific and 193 TCGA LUSC-specific TFs. Activity
  is row-z standardized without using labels within each external lung cohort,
  and the score is mean LUAD-signature activity minus mean LUSC-signature
  activity. The decision threshold is fitted once in TCGA and then fixed.
- Eligible treatments are single agents, exclude untreated controls, have at
  least 10 LUAD-like PDXs with paired baseline activity and response, and occur
  in at least one cell-line pharmacogenomic dataset. All eligible drug-TF tests
  are reported; BH correction is applied across the complete eligible PDX
  drug-TF screen, which is stricter than within-drug correction.
- Figure 5G-H is generated only when the same canonical drug-TF pair was
  significant with a concordant direction in at least two cell-line datasets
  and then reaches global PDX FDR <= 0.05 with at least 10 LUAD-like PDXs.
  Selection is deterministic by FDR, sample size and name. Otherwise G-H is
  explicitly marked unsupported rather than replaced by a nominal association.

## Recorded classifier amendment before response testing

The initially implemented per-sample percentile-rank classifier showed a
cross-study score-location shift (balanced accuracy 0.545 in PDMR and 0.500 in
DepMap). Before inspecting any PDX treatment outcome, five label-free transform
variants were compared against the independent histology labels in PDMR and
DepMap. The cohort row-z LUAD-minus-LUSC signature was retained because it gave
balanced accuracy 0.940 in TCGA, 0.904 in PDMR and 0.690 in DepMap, while
preserving the pre-defined TF signatures and TCGA-fitted threshold. All variant
results are retained in `SupplementaryFigure8_classifier_variant_audit.tsv`;
this model-selection step makes the PDX state label an exploratory validated
classifier rather than a pristine preregistered classifier.

## Interpretation boundary

Cell-line associations are discovery/replication evidence, not causal drug
targets. PDXE is an independent in-vivo validation layer, but the LUAD-like
label is computational. No panel may be described as clinical treatment
selection without prospective validation.

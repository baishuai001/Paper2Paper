# LUAD Figure 4 protocol: prognostic relevance of the 31 HC-TFs

## Anchor-aligned design

All 31 Figure 2 HC-TFs are tested.  The exposure in each Cox model is high
versus low TF activity at the cohort/endpoint-specific median, matching the
anchor's survival and permutation implementation.  Patients are censored at
10 years and a two-sided Cox p-value of 0.05 defines a nominal prognostic TF.
Raw p-values define the anchor-style Figure 4 counts; BH FDR is also reported
and is never hidden.

Two independent LUAD patient cohorts are analysed separately:

- TCGA-LUAD RNA-seq/VIPER: overall survival (OS) and disease-free survival
  (DFS, the available recurrence endpoint).
- GSE41271 microarray: OS and recurrence-free survival (RFS).  Its independently
  inferred one-subnet ARACNe3 network lacks five of the 31 HC-TFs.  The frozen
  primary external-cohort analysis therefore projects the TCGA regulon onto
  GSE41271, preserving the complete 31-TF feature space; the 26 TFs available
  from the independently inferred GSE41271 regulon are retained as a sensitivity
  analysis.  This adaptation is disclosed and is not presented as identical to
  the anchor's cohort-specific-network survival implementation.

For each TF and endpoint, both univariable and multivariable Cox models are fit.
The LUAD-adapted multivariable model contains age, sex, pathologic stage and
smoking history.  This replaces breast-specific lymph-node/tumour-size/
menopause covariates while preserving the role of established clinical factors.
Complete-case counts and events are frozen for every model.

## Main and supplementary outputs

- Figure 4A-D: multivariable forest plots for GSE41271 OS/RFS and TCGA OS/DFS.
- Figure 4E/G: overlap of favorable and adverse TFs between the two endpoints
  within each cohort.
- Figure 4F/H: 10-year Kaplan-Meier curves for representative favorable and
  adverse TFs; maximally selected rank statistics determine plotted cutpoints.
- Supplementary Figure 7: univariable forests, 5,000-permutation null
  distributions and complete-case/event audit.

## Permutation test

For every cohort, endpoint and model type, survival times alone are shuffled
5,000 times while event status, TF activity and clinical covariates remain fixed.
All 31 median-dichotomized TFs are refit and the number with p <= 0.05 is
recorded.  Empirical p is the fraction of null counts at least as large as the
observed count; fold enrichment is observed divided by the null mean.

## Interpretation and stopping

Figure 4 is computationally executable once both patient cohorts pass endpoint
and covariate audit.  The absence of significant or cross-endpoint prognostic
TFs is a valid negative biological result, not a software failure.  No KM curve
is presented as confirmatory unless its TF and selection rule are disclosed.

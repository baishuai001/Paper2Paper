# Pilot minimal-run specification: GSE62254 NRRS reconstruction

This specification records the bounded real-data test that was actually run.
It does not freeze a future manuscript analysis and does not turn this training
route into a manuscript candidate.

## Cohort and exclusions

- Use all 300 primary gastric-tumor patients in GSE62254/ACRG for whom the GEO
  expression identifier links exactly to the ACRG clinical workbook.
- Require a unique GEO ID, complete overall-survival months and 0/1 death
  status, and expression for every locked score gene after annotation.
- Fail rather than silently remove or impute a patient or required gene.
- Treat GSE62254 as development-exposed because the anchor included it in the
  outcome-guided GEO-meta feature-screening pool. It is a reproduction cohort,
  not untouched external validation.

## Statistical and biological unit

The patient is the unit for scoring and survival inference. Microarray probes
are measurement features, not independent observations. The run contains 300
patients and 152 deaths.

## Exposure, comparison and outcome

- Exposure: continuous reconstructed eight-gene NRRS.
- Prespecified display comparison: deterministic cohort-median high versus low
  NRRS groups, with ties ordered by GEO ID to keep grouping reproducible.
- Primary outcome: overall survival in months with death as the event.
- No additional clinical covariates are introduced in this bounded spike.

## Data resources and identity linkage

Expression comes from the GSE62254 RMA Series Matrix, probe annotation from
GPL570, and survival fields from Supplementary Data 2 of DOI
`10.1038/s41467-018-04179-8`. The three files must match the byte counts and
SHA256 values in `config/gse62254_resources.tsv`. GEO IDs must be unique in
each source, their sets must match exactly, and unmatched counts must be zero
before any score or survival model is calculated.

## Locked score contract

Feature order is `AGT;EPHB3;GNAI1;LPAR2;LRRC4C;NPY1R;NRP1;SEMA6A` and the
published coefficients are `0.040;-0.040;0.018;-0.102;0.014;0.072;0.175;0.006`.
Only GPL570 probes unambiguously assigned to one target gene are retained. The
main value for a gene is the arithmetic mean of all retained probes; each gene
is then standardized within GSE62254 using sample SD (`ddof=1`) before applying
the coefficients in the locked order. A missing target gene stops the run.

The outcome-blind sensitivity analysis replaces the mean-probe value with the
highest-variance retained probe for each gene. It reports score correlation,
group agreement and survival direction; it is not used to choose a rule based
on survival.

## Models, software and reproducibility

The run uses Python 3.12 with the versions recorded in `execution/runs.tsv` and
`code/requirements-gse62254.txt`. Cox models estimate the continuous score and
high-versus-low hazard ratios; Kaplan-Meier and log-rank analyses provide the
grouped display. There is no stochastic feature selection, so no random seed is
required. Resource mismatch, identifier mismatch, missing genes or change in
the 16 reviewed reference metrics fails closed.

## Outputs and acceptance

Required outputs are the per-patient source table, probe-mapping table,
machine-readable summary, noninteractive Kaplan-Meier PNG and run metadata.
Acceptance requires 300 unique patients, 152 events, all eight genes mapped,
finite scores, deterministic grouping, survival estimates derived from the
source table, and the expected-result regression passing.

The patient-level derivative remains local until redistribution review. The
repository retains its checksum receipt and uses the tracked aggregate summary
as the auditable result/Figure source table; this does not turn the aggregate
table into a substitute for patient-level recomputation.

## Claim boundary and falsifier

The strongest permitted claim is transparent retrospective reconstruction of
an adverse survival association under explicit implementation choices. The run
does not establish exact reproduction, independent validation, clinical
utility or single-patient deployment. Failure to compute the locked score or
failure to reproduce the adverse direction would falsify this bounded route.

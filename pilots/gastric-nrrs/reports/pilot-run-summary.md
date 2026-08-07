# GSE62254 NRRS pilot run summary

## Outcome

The bounded reconstruction completed on all 300 GSE62254 patients, including
152 deaths. All 300 expression identifiers matched the clinical workbook
exactly, and all eight score genes were represented by 19 unambiguous GPL570
probes.

Under the locked main rule (mean of all unambiguous probes per gene,
within-cohort sample-SD z-score, published coefficients, deterministic cohort
median split):

- continuous NRRS HR per SD: 1.456 (95% CI 1.239-1.711; p=5.17e-06);
- high versus low NRRS HR: 1.529 (95% CI 1.109-2.107);
- log-rank p=0.00906;
- concordance index=0.596.

The risk direction agrees with the anchor: higher reconstructed NRRS is
associated with poorer overall survival. Discrimination is modest, so the run
supports a retrospective association rather than a strong clinical predictor.

## Mapping sensitivity

Replacing the mean-probe rule with the highest-variance probe per gene was
outcome-blind and retained the adverse survival direction. The two scores had
Spearman rho 0.972 and 95.3% high/low agreement: 14 of 300 patients changed
groups. This shows why an undocumented probe rule cannot be treated as a minor
implementation detail.

## What this run does not prove

- It is not exact figure reproduction because the anchor omits essential score
  computation choices.
- It is not untouched external validation because GSE62254 participated in the
  anchor's outcome-guided feature-screening pool.
- It is not a deployable single-patient model because the reconstructed score
  uses cohort-specific standardization and cutoff.
- It is not a new manuscript direction by itself. It is a worked example for a
  fixed bulk-signature and survival-analysis module, not a benchmark for other
  paper types or for Paper2Paper as a whole.

## Controls exercised and their limits

The pilot exercised field-level resource provenance, cross-source identity
checks, a cohort analysis-role ledger, a computable signature contract,
target-environment path evidence, structured issue recording, checksum pinning,
and separate execution/manuscript decisions. These controls are reusable only
at the level directly covered by core tests. They do not transfer the gastric
data, formula, expected results or scientific conclusions to another pilot.
Details and paper-specific defects are separated in `pilot-findings.md` and
`evidence/issues.tsv`.

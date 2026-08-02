# Workflow and gates

Paper2Paper uses cheap evidence first and expensive execution only after a
route survives triage. Every gate has a concrete exit condition.

## G0 - Scope and anchor decomposition

Tasks:

- define the intended paper type, audience, resources, and stopping rule;
- audit the anchor at claim, figure, dataset, method, and code-module level;
- generate all plausible route modes, including simple substitutions;
- record an explicit adaptation map for every non-reproduction route.

Exit: a reviewable route portfolio without an approved route.

## G1 - Data, code, figure, and overlap triage

Tasks:

- verify dataset identifiers, access, required metadata, cohort overlap, and
  minimum downloadable files;
- locate author code and independent code donors, then assign provisional
  A-D/U grades;
- map every planned main figure to data, code, and an acceptance test;
- search nearest publications and classify overlap.

Exit: every live route has one complete assessment and rule-derived P0-P3.

## G2 - Executable spike and route approval

Tasks:

- download and parse representative data;
- install the environment and run noninteractive smoke tests;
- produce a representative source table or figure fragment;
- resolve critical metadata, license, dependency, or duplicate risks;
- obtain human approval for one P0 route.

Exit: one approved route with verified required datasets, A/B code mappings,
a complete minimum figure map, a non-duplicate overlap record, and an approving
review.

## G3 - Analysis specification

Tasks:

- freeze cohorts, exclusions, biological units, comparisons, outcomes,
  covariates, parameters, random seeds, and signature formulae;
- separate discovery, validation, sensitivity, and exploratory analyses;
- define module input/output contracts and scientific-invariant tests;
- freeze the manuscript claim and figure plan version.

Exit: executable specifications contain no silent analyst choices.

## G4 - Execution

Tasks:

- build adapters and execute approved modules;
- record commit, environment, configuration, data manifest, timestamps, logs,
  and artifacts;
- verify source tables before styling figures;
- register failures and negative results rather than overwriting them.

Exit: required figures and tables can be regenerated from recorded inputs.

## G5 - Result and claim audit

Tasks:

- compare results with prespecified claims and claim ceilings;
- assess robustness, leakage, pseudoreplication, cohort heterogeneity, and
  alternative explanations;
- trigger continue, refine, reroute, or stop through change requests;
- draft Results, Methods, and Limitations from verified artifacts.

Exit: every manuscript claim has adequate evidence or is removed/downgraded.

## G6 - Release

Tasks:

- rerun the submission configuration from a clean environment where practical;
- freeze manifests, source tables, figures, code versions, and checksums;
- update the publication-overlap search;
- complete human review of claims, limitations, and reproducibility package;
- tag a reviewed release.

Exit: a manuscript and reproducibility package are ready for submission.

## Work relevance at every gate

Before adding work, answer:

1. Which claim, figure, section, qualification gap, or review does it serve?
2. What is the minimum deliverable?
3. What is the stopping condition?
4. Would omission materially weaken execution, scientific validity, or review?

If these cannot be answered, the work does not enter the active scope.

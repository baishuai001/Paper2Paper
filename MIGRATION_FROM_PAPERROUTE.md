# Migration from PaperRoute

## Boundary

The Paper2Paper repository imported PaperRoute history through commit
`54050fef58f3938e7e7ba7dfd20dcae2470a28e5` on 2026-08-02. The lightweight
tag `paperroute-import-2026-08-02` marks that boundary.

The import is retained for provenance. It is not a compatibility promise and
does not make PaperRoute policy active in Paper2Paper.

## Retained engineering components

- Git history and review provenance;
- isolated paper workspaces;
- TSV parsing and registry contract validation;
- global ID uniqueness and reference integrity;
- WorkItem, Claim, Module, Run, Result, Review, Decision, ChangeRequest, and
  Dependency concepts;
- result-driven change requests and downstream impact traversal;
- CLI, unit-test, and continuous-integration scaffolding.

## Rebuilt product core

- Python package: `paperroute` -> `paper2paper`;
- CLI: `paperroute` -> `paper2paper`;
- schema: 0.2 -> 0.3;
- direction objects -> execution-qualified routes;
- general resource notes -> typed dataset and code qualification records;
- prose-only figure expectations -> route-linked figure plans;
- broad nearest-neighbor judgment -> explicit publication-overlap dimensions;
- free-form candidate ranking -> rule-derived P0-P3 execution priorities.

## Deliberately not migrated into active state

- former innovation thresholds, fields, scoring, approval conditions, or
  recommendation language;
- old direction assessments and candidate portfolios;
- statements that rejected a route merely because it substituted a marker,
  gene set, cell type, cancer type, scope, or signature;
- old workspace conclusions derived from that policy.

The prior material remains available in Git history and the import tag. Only
neutral provenance facts such as titles, DOIs, accessions, public URLs, file
checksums, access states, and observed code availability may be copied into a
new Paper2Paper workspace. Candidate routes must be generated and assessed
under schema 0.3.

## Breaking-change rule

No automatic migration converts a PaperRoute direction decision into an
approved Paper2Paper route. Every imported paper begins at intake or route
review and must pass data, code, figure, scientific-validity, overlap, and
human-review gates again.

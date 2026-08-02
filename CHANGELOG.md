# Changelog

## [Unreleased]

### Changed

- Rebuilt the product as Paper2Paper 0.3 with a new package, CLI, schema, and
  project discipline.
- Replaced the former direction-audit policy with execution-first route
  qualification.
- Made marker, gene-set, cell-type, cancer-type, pan-cancer, signature, and
  combined substitutions first-class route modes.
- Separated publication duplication from route merit; only substantive
  duplication is a hard publication-overlap stop.
- Added dedicated dataset, dataset-route, code-source, code-module,
  publication-overlap, and figure-plan registries.
- Added transparent P0-P3 execution priorities and hard approval checks for
  scientific validity, verified data, qualified code, figure closure, and
  human review.
- Reset the SPP1 and gastric NRRS workspaces to factual intake records so old
  policy conclusions cannot remain active.

### Removed

- The `paperroute` package and CLI from the active product.
- Legacy novelty fields, assessments, approval gates, initialization policy,
  tests, and active G0 recommendation reports.

## Import boundary

Paper2Paper was copied from PaperRoute at commit `54050fe` on 2026-08-02.
Earlier PaperRoute releases and decisions remain recoverable from Git history;
they are not active Paper2Paper policy.

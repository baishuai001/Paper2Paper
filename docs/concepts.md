# Concepts

## Three levels

### Workflow template

The reusable process, schema, validation rules, route templates, and review
gates in the repository root.

### Paper workspace

One isolated paper-adaptation project containing its anchor, candidate routes,
datasets, code donors, figure plan, decisions, runs, and results.

### Workflow run

One immutable execution snapshot identified by code commit, configuration,
data manifest, environment, timestamps, and output artifacts.

## Source of truth

Paper2Paper separates machine-readable state from generated prose:

1. `PROJECT_DISCIPLINE.md` defines binding product policy.
2. `PROJECT.json` stores project identity, allowed route modes, execution
   policy, current gate, and active route.
3. `registry/*.tsv` stores facts, entities, references, qualification states,
   decisions, and reviews.
4. Code reads approved configuration and registry state.
5. Runs record what actually executed.
6. Figures point to source tables and acceptance tests.
7. Markdown reports are views generated from the sources above.

A decision that exists only in chat is not a project decision. A threshold
that exists only inside a script is not an approved parameter. A figure
without a source table and run record is not a complete result.

## Core entities

- `Paper`: anchor, supporting, method, data, code, or overlap literature.
- `Route`: one candidate manuscript plan and its adaptation mode.
- `RouteAdaptation`: an explicit retain, repair, substitute, extend, or drop
  operation on a named axis.
- `RouteAssessment`: independent execution fields used to derive P0-P3.
- `Dataset`: an accession-level resource with access, metadata, download, and
  provenance status.
- `DatasetRouteMap`: why a route needs a dataset and whether the required
  fields and independence are verified.
- `CodeSource`: a versioned code donor or reconstruction graded A-D/U.
- `CodeModuleMap`: a route module mapped to a concrete code source, entry
  fragment, adaptations, and tests.
- `PublicationOverlap`: a dated comparison with the nearest publication across
  seven overlap dimensions.
- `FigurePlan`: a manuscript output mapped to claims, datasets, code, source
  tables, and acceptance tests.
- `Claim`: a testable statement and its evidence ceiling.
- `WorkItem`: manuscript-linked work with a minimum deliverable and stop rule.
- `Module`: an executable analysis specification linked to a route and claim.
- `Decision` and `Review`: proposed and human-approved judgments.
- `Run` and `Result`: concrete execution and its evidential implication.
- `ChangeRequest`: a versioned adjustment triggered by evidence or results.
- `Dependency`: a directed edge for traceability and impact analysis.

## Lifecycle versus validity

Lifecycle describes progress. Validity describes whether an entity can still
support downstream work. A module can remain `implemented` but become `stale`
when its route, data, code, method, or claim changes. It must be rerun before
its results return to `current`.

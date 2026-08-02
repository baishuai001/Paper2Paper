# Figure planning

## Figure-first execution map

Every required manuscript output is registered before full implementation.
One `figure_plan.tsv` row connects:

```text
anchor figure or module
  -> target manuscript role
  -> claim
  -> dataset(s)
  -> analysis module
  -> qualified code mapping(s)
  -> source table
  -> acceptance test
  -> generated artifact
```

The anchor figure can be left blank when the output is new. The route cannot be
approved if its minimum required figures are missing their data, code, source
table, or acceptance test.

## Status meanings

- `idea`: desired output without a complete execution path;
- `mapped`: data and code donors identified but not yet qualified;
- `ready`: all dependencies qualified and the output contract frozen;
- `generated`: artifact produced from a recorded run;
- `verified`: source table, rendering, labels, statistics, and interpretation
  reviewed;
- `blocked`: an indispensable dependency failed;
- `dropped`: intentionally removed with rationale.

## Acceptance tests

An acceptance test should be falsifiable and local to the output. Examples:

- every point corresponds to one patient-level effect estimate;
- feature order exactly matches the frozen signature formula;
- discovery samples never appear in external validation;
- panel labels and source-table contrasts have the same direction;
- row counts equal the prespecified cohort manifest after exclusions;
- rerunning from the run manifest regenerates the same source table within a
  declared numerical tolerance.

Visual attractiveness is reviewed after data and statistical correctness.

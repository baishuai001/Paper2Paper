# Code qualification

## Code donor grades

| Grade | Meaning | Route use |
| --- | --- | --- |
| A | version and license fixed; environment and entrypoint clear; noninteractive; no hardcoded paths/private inputs; install and full tests pass | approved route |
| B | core module runs noninteractively; environment and entrypoint clear; install and smoke test pass; only bounded adaptation remains | approved route |
| C | useful fragments or processing details, but substantial reconstruction is required | P2 donor only |
| D | blocked, unusable, unlicensed for intended reuse, or dependent on unavailable private inputs | cannot execute route |
| U | not yet examined | triage only |

The grade describes observed executable quality, not journal prestige or code
style. An official package vignette may be A/B for one module while an anchor
repository is C/D for the full paper.

## Qualification procedure

1. Record canonical URI and source type.
2. Pin a commit, tag, release, or package version.
3. record license, language, runtime, system dependencies, and environment.
4. Identify a noninteractive entrypoint and its input/output contract.
5. Check hardcoded paths, embedded credentials, private inputs, manual GUI
   steps, and undocumented global state.
6. Install from a clean or documented environment.
7. Run the smallest real-data smoke test.
8. Run available tests and add adaptation-specific tests.
9. Map each required route module to the exact source path or function.

## Reconstructed code

Missing or poor author code does not automatically kill a good analysis idea.
Paper2Paper may combine:

- scientific logic from the anchor;
- preprocessing from a data-generating paper;
- statistical implementation from official package documentation;
- engineering patterns from a high-quality related repository;
- newly written adapters and tests.

The result must be registered as `reconstructed`; provenance for every donor
module remains explicit. It starts at U/C and becomes A/B only after the same
installation and testing requirements are met.

## Minimum tests

- unit tests for transformations, identifiers, formulas, and edge cases;
- smoke test using a representative real input;
- integration test from input manifest to source table;
- scientific invariants such as patient-level independence, no outcome leakage,
  stable feature order, correct contrast direction, and expected row counts;
- deterministic seeds and parameter capture where stochastic methods are used.

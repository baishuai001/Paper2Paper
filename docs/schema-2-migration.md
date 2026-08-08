# Schema 2.0 migration

Schema 2.0 is a deliberate break from the earlier model in which one workspace
could move from anchor audit directly into manuscript writing and a Pilot issue
could be labeled “verified in core”. That model mixed paper analysis, reusable
code and product maturity, so it could overstate what one real-paper example had
proved.

## Workspace migration

Every workspace now declares `workspace_kind`:

- `pilot` uses `anchor_audit -> route_generation -> verification -> pilot_review
  -> pilot_complete`;
- `manuscript_project` uses `specification -> execution -> interpretation ->
  writing -> complete` and requires source Pilot, route and passed-run IDs.

Do not convert a Pilot by editing its stage. Remove its empty manuscript
template, complete `reports/pilot-outcome.md`, record a user decision, and use
`paper2paper promote` to create a separate project. Promotion copies auditable
direction evidence but resets code checks, Figure status, runs and results so
formal outputs must be regenerated.

## Route field migration

| Before | Schema 2.0 |
| --- | --- |
| `status` | `decision_status` |
| implicit readiness | `evidence_stage` |
| no reusable scope | `capability_ids` |
| `science_status` | `scientific_review` |
| `manuscript_eligible` report | `promotion_evidence_complete` |

`execution_ready` remains a technical/scientific minimum-run conclusion. It is
not a user decision. `promotion_evidence_complete` additionally requires a
manuscript-candidate role and explicit user review of the evidence package; a
separate user `promote` decision is still required.

## Finding migration

Legacy `scope` and states such as `resolved_in_pilot`,
`partially_resolved_in_core` and `verified_in_core` are removed. Each issue now
records:

- where it was observed (`observed_layer`);
- the proposed reuse boundary (`candidate_scope`);
- the issue type;
- the local issue status;
- affected capabilities and an optional central promotion ID.

Local resolution and cross-paper maturity are intentionally independent.
Promotion maturity is stored in `registries/promotions.tsv`; module implementation
maturity is stored in `registries/module_releases.tsv`.

## Evidence migration rules

- Preserve original observations, data facts, failures, results and scientific
  limitations.
- Never invent a user decision while translating an AI recommendation.
- Do not upgrade `minimal_real_run` to `figure_loop_closed` unless every required
  Figure has a verified result from a passed run.
- Runs now declare `run_kind`, `input_provenance` and `exit_code`. Only a
  `real_data` run with `status=passed`, `exit_code=0`, valid timestamps, a real
  manifest, a retained log and every declared artifact present can satisfy
  scientific execution readiness. Unit and synthetic runs remain separate.
- One reference Pilot cannot become cross-paper confirmation; independent
  transfer requires another project, DOI and dataset fingerprint.
- Frozen expected values are exact-input change detectors, not evidence that a
  method transfers to another paper type.

## Repository checks after migration

```bash
paper2paper validate pilots/<pilot>
paper2paper validate-registry .
paper2paper coverage .
python -m unittest discover -s tests -v
```

These commands report structural/contract and recorded regression evidence.
They do not replace a new data download, complete scientific review or the
user's route and manuscript decisions.

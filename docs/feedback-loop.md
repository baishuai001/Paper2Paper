# Result feedback and human review

Paper2Paper assumes that execution can invalidate the plan.

```text
run -> result -> route effect -> change request -> impact analysis
    -> human review -> new version or rejection -> rerun
```

## Route effects

- `none`: no route implication;
- `continue`: compatible with the approved plan;
- `refine`: adjust data, code, methods, parameters, a figure, or a bounded
  claim while answering the same central question;
- `reroute`: change the central question, target disease/object, primary
  outcome, or evidence chain so the approved route no longer answers it;
- `stop`: continuing cannot produce a defensible manuscript.

Every `refine`, `reroute`, or `stop` result requires a ChangeRequest.

## Change control

A change request records the triggering result, work item, affected entities,
proposal, severity, gate to reopen, and approval. High-impact changes cannot
be implemented without a completed approving review.

Use:

```bash
paper2paper impact PROJECT_DIR ENTITY_ID
```

to list active downstream dependencies that may become stale. Impact reporting
does not silently mutate project state; invalidation is a reviewable decision.

## Protection against hindsight bias

1. preserve the original specification and result;
2. record the observation before changing the plan;
3. create a change request;
4. inspect downstream claims, figures, modules, data, and code mappings;
5. obtain review when required;
6. create a new version;
7. label post-result analyses exploratory until independently validated.

Feedback is not permission for unlimited optimization. Every follow-up needs a
manuscript destination, minimum deliverable, and stop condition.

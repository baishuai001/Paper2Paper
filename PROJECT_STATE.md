# Paper2Paper current state

Updated: 2026-08-07

## Product boundary

Paper2Paper is an alpha workflow for a beginner and AI to audit a real anchor
paper, compare bounded imitation routes, verify data and code with real inputs,
and create a separate formal manuscript project only after user approval. Its
success criterion is a scientifically defensible paper package, not continual
platform expansion.

The current schema is 2.0.0 and the package version is 0.5.0-alpha.1. Core,
reusable modules, real-paper Pilots and formal manuscript projects are separate
layers. Local issue resolution, module maturity and cross-paper rule promotion
are separate records.

## Canonical server workspace

The server-local project working directory follows:

`/media/desk16/<server-account>/projects/Paper2Paper`

The concrete account path is kept outside this public repository. Credentials,
hostnames and tokens must never be stored here.

The core workflow currently passes its repository suite on the canonical
server's Python 3.10 runtime. The gastric scientific module was reference-run
under its recorded Python 3.12 environment; its exact lock is not currently
resolvable from the server package index because `lifelines==0.30.3` is absent.
The failed isolated install was not replaced with a nearby version, so the
server must not be described as module-verified until the exact artifacts are
available or a new environment release is rerun and reviewed.

## Real-paper Pilots

- `P2P-GASTRIC-NRRS`: a bounded training/reference reconstruction. GN-R01 has
  a passed 300-patient real-data run and a separate focused-test run. It is not
  a manuscript candidate and does not test SPP1+TAM, single-cell, spatial or the
  whole workflow. It is awaiting a user retain/close decision; none is inferred.
- `P2P-SPP1-TAM-JITC`: an independent SPP1+TAM JITC anchor Pilot. It remains at
  anchor audit and must develop and verify its own route, datasets, code modules
  and cross-scale claims.

## Reusable evidence

- Resource integrity and patient-ID linkage each have one real-paper reference;
  they are not cross-paper confirmed.
- The fixed bulk-signature implementation is reference-verified on GN-R01 and
  has focused unit/adversarial coverage; it is not transfer-verified.
- Frozen GN-R01 values are exact-input change detectors only.
- The dedicated negative identity-linkage case is registered but not run; that
  promotion remains observed rather than being overstated as provisional.
- The gastric patient-level derived table remains local-only pending a
  redistribution review. A tracked receipt records its byte size, SHA256,
  dimensions and columns; repository-visible claims point to tracked aggregate
  source tables instead of depending on that ignored file.

## Formal manuscript projects

None exists. A formal project must be created with `paper2paper promote` after
a manuscript-candidate route has complete evidence review and explicit user
approval. Pilot results are not copied as formal manuscript results.

In this alpha, `run_kind=real_data` remains a declared provenance label rather
than a cryptographic proof of dataset identity. Reviewers must compare each run
manifest with the registered candidate, resource and cohort-usage records.
Likewise, registry case-to-run and assertion records are auditable declarations,
not cryptographic proof that every asserted property was executed by the bound
run. Automated checks enforce identities, statuses, real-data command and
manifest equality, independent anchor/dataset keys and declared unit-test
coverage; a reviewer must still inspect native logs and artifacts before
accepting transfer evidence.

## Immediate boundary

The schema/registry migration is complete. Keep the automated checks green and
use the canonical server directory for continued work. The next scientific task
is the SPP1+TAM Pilot's candidate-route and availability work; no additional
generic feature should be built unless that or another real paper exposes a
concrete blocker.

# Paper workspaces

Each immediate subdirectory is one isolated Paper2Paper instance. Reusable
workflow code and policy live at the repository root; paper-specific facts,
routes, data, code donors, figures, decisions, runs, and reviews live here.

| Workspace | State | Active route | Purpose |
| --- | --- | --- | --- |
| `spp1-tam-jitc/` | G0 intake | none | Regenerate SPP1+TAM adaptation routes under execution-first rules. |
| `bmc-cancer-2025-gastric-nerve-model/` | G0 intake | none | Regenerate gastric NRRS/template adaptation routes under execution-first rules. |

Rules:

1. Every paper receives a separate directory and `project_id`.
2. Entity IDs are workspace-scoped; use `project_id + entity_id` when comparing
   across workspaces.
3. Candidate routes and decisions become active only after validation and
   human approval in their own workspace.
4. Only neutral provenance facts may be imported without reassessment.
5. Source PDFs, credentials, controlled data, and large analysis files remain
   outside Git and are represented by URIs, checksums, and access records.
6. Paper-specific adapters stay inside a workspace until at least two real
   projects demonstrate that they belong in `src/paper2paper/`.
7. A workspace cannot mark a route approved until schema 0.3 validates its
   data, code, figure, scientific, overlap, and review requirements.

# Publication assembly input contract

The assembly consumes completed plotting outputs; it does not execute a
statistical analysis or recreate a panel. Every supported LUAD input must be a
valid PDF inside one of the four declared official-port roots. The assembler
copies its PDF content stream with `pypdf.merge_transformed_page` and records
the source SHA-256. Raster input is allowed only for a TNBC reference crop in
the comparison atlas.

## Frozen roots and roles

| Root key | Expected run | Eligible content |
|---|---|---|
| `fig12` | `official-code-port-v3/run-fig12-final` | Figure 1/2 and Supplementary Figure 2/4 official atoms |
| `fig3` | `official-fig3-port/run-v3` | Figure 3 and Supplementary Figure 4D/5 atoms |
| `fig45` | `official-code-port-v3/run-fig45-v2-p0fix-20260821-r2` | Figure 4/5 and Supplementary Figure 6/7 atoms |
| `supp3` | supplied when complete | Supplementary Figure 3A-L atoms plus one `verification_receipt.json` |

Within `fig12`, publication assembly accepts primary atoms only below
`adapter_final_corrected/results/visuals/`. GSE41271 microarray occupies the
primary external-validation slot. The GSE81089 RNA-seq replay is accepted only
below
`adapter_final_corrected/secondary_validation/GSE81089/results/visuals/` and
is placed on a separately titled secondary-validation page.

Supplementary Figure 2F must carry
`ANCHOR_BUGFIX (official one-line data-object correction)`. No selector points
to `adapter_strict_replay` or `runtime_strict_replay`; the
`CAPSULE_BUG_RETAINED` execution is audit-only. Figure 2C carries
`CAPSULE_PRINT_MISMATCH`. Figure 1A, Figure 2A, Supplementary Figure 1A/B, and
Supplementary Figure 4A remain `MANUAL_REQUIRED`; Figure 5G/H and
Supplementary Figure 7D remain `UNSUPPORTED`.

## Cardinality and failure behavior

`assembly_spec.json` freezes, per atom, whether exactly one PDF or all matching
PDFs are required. Zero matches, an invalid PDF, an out-of-root source, or more
than one match for a `match: one` atom is written to `missing_inputs.tsv`.
Strict QA fails while any required entry remains. It never substitutes a
different cohort, a strict replay, or a newly drawn panel.

Path and filename agreement do not establish official provenance. The assembly
therefore requires root-specific official-port receipts/source locks and binds
each resolved PDF to upstream SHA-256/byte evidence in
`source_identity_manifest.tsv`. Figure 1/2 runtime manifests must report full
commit `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`; Figure 4/5 verification must
report `edf5314`. Figure 3 currently has a complete 19-PDF hash manifest;
Supplementary Figure 3 binds 12 PDFs to receipt byte counts while recording
their observed SHA-256. Figure 1/2 and Figure 4/5 are bound respectively to
their completed `pdf_sha256_manifest.tsv` and `atomic_pdf_sha256.tsv`; the old
unhashed Figure 4/5 package is not accepted merely because filenames match.

Supplementary Figure 3 is additionally release-gated by its upstream receipt.
The receipt must identify exactly 3A-L, report use of official plotting code
with zero constructor rewrites, and bind to the resolved PDF names and byte
counts. Its current `publication_ready=false`, 3A `PARTIAL`, and warnings for
3D/3F/3K/3L force the aggregate state
`STRICT_OFFICIAL_PORT_REVIEW_ONLY`. Technical file validation cannot override
that upstream scientific/visual boundary.

The production runner writes exactly seven PDFs to
`PILOT_ROOT/execution/luad-official-code-port/final/`. Use
`run_assembly.sh --dry-run ...` to materialize the machine-readable
`input_contract.tsv` and check readiness without creating those PDFs.
The default runner is review mode. `--require-publication-ready` is the explicit
non-zero release gate and currently prevents final PDF generation.

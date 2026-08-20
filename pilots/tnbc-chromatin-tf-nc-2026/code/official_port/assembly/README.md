# Strict official-code PDF assembly

This directory assembles the TNBC official-code LUAD atoms into seven stable
PDF deliverables. It performs layout only. It does not run an analysis and it
does not contain a replacement plot constructor.

## Output contract

The run writes exactly:

1. `Figure1.pdf`
2. `Figure2.pdf`
3. `Figure3.pdf`
4. `Figure4.pdf`
5. `Figure5.pdf`
6. `Supplementary_Figures.pdf`
7. `TNBC_vs_LUAD_comparison_atlas.pdf`

The production runner fixes all seven paths under
`PILOT_ROOT/execution/luad-official-code-port/final/`; it does not create a
second competing output directory. This stable path is not a publication-ready
claim. Default runs create seven review artifacts and visibly stamp the
aggregate release state on every page.

LUAD atomic PDFs are imported with `pypdf.merge_transformed_page`. That method
copies the source PDF content stream and resources into the destination page;
it does not rasterize the plot. Assembly may translate and uniformly scale an
atom to fit its publication slot. It never changes the atom's font, palette,
geometry, statistic, clustering, edge set, legend, or internal labels.

TNBC comparison material can be read from the source article PDFs or from the
already exported reference PNGs. A TNBC PNG remains a raster reference crop;
this exception never applies to a LUAD official atom. `embedding_manifest.tsv`
records the method used for every placed object.

## Frozen publication boundaries

- Figure 1A, Figure 2A, Supplementary Figure 1A/B, and Supplementary Figure 4A:
  `OFFICIAL CODE NOT PROVIDED / MANUAL_REQUIRED`.
- Figure 2C: `CAPSULE_PRINT_MISMATCH`. The untouched official `06C` capsule
  output is placed; it is not rearranged to imitate the printed panel.
- Supplementary Figure 2F: `ANCHOR_BUGFIX (official one-line data-object
  correction)`. The final panel is read only from
  `adapter_final_corrected/results/visuals/`. The uncorrected
  `CAPSULE_BUG_RETAINED` strict replay remains an audit artifact and is never
  eligible for publication assembly.
- Figure 5G, Figure 5H, and Supplementary Figure 7D: `UNSUPPORTED` because the
  frozen LUAD PDX eligibility chain has no qualifying input.
- A fig45 constructor with a receipt of `SUPPORTED_ZERO` receives a neutral
  `SUPPORTED_ZERO` placeholder. Scientific absence is not treated as a render
  error and is not replaced by another plot.

Only these neutral boundary placeholders have borders. Real atomic panels are
placed on a white page without dashboard cards, colored banners, or newly
authored statistical graphics.

## Parameterized inputs

The expected cloud roots are:

- Figure 1/2: `official-code-port-v3/run-fig12-final`
- Figure 3: `official-fig3-port/run-v3`
- Figure 4/5: `official-code-port-v3/run-fig45-v2-p0fix-20260821-r2`
- Supplementary Figure 3 A-L: a root containing the atoms and exactly one
  `verification_receipt.json`

Every atom is resolved recursively by the frozen patterns in
`assembly_spec.json`. An unusual filename can be supplied with
`--source-override source_overrides.tsv`. The override schema is shown in
`source_overrides.example.tsv`; normalized crop coordinates use a top-left
origin. Overrides do not permit raster LUAD atoms.
Every resolved source, including an override, must remain inside its declared
official-port root and inside any non-glob directory prefix frozen in the
selector. This prevents a filename override from importing the strict
Figure 1/2 bug replay or an unrelated PDF.

Filename scope is only the first identity check. Before assembly, each root is
also bound to its official-port receipts/source manifests and resolved atomic
PDF evidence. Figure 3 uses `official_port_receipt.json`, its official source
lock, and `official_output_sha256.tsv`; Supplementary Figure 3 uses its
materialization/source locks plus verification-receipt byte counts; Figure 1/2
requires both corrected runtime manifests at official commit
`edf5314ce5b9ee0e2f88b2310e7c2df5619ad888` and an output hash manifest;
Figure 4/5 requires its PASS verification at `edf5314` and its per-PDF hash
manifest. The old Figure 4/5 package is intentionally not grandfathered
without hashes.

Supported output-hash filenames are `official_output_sha256.tsv`,
`atomic_output_sha256.tsv`, `output_sha256.tsv`, `pdf_sha256.tsv`,
`pdf_sha256_manifest.tsv`, or `atomic_pdf_sha256.tsv`. Rows may name
`relative_path`, `artifact`, `path`, `output`, `file`, `pdf`, or `filename`,
with `sha256`/`hash` and optional `bytes`/`size`. A mismatch or ambiguous
basename is a required-input failure.

For the Figure 1 external-validation blocks, GSE41271 microarray occupies the
frozen primary METABRIC-analogue slot. Supplementary Figure 2A/B on the main
Supplementary Figure 2 page therefore resolve only from
`adapter_final_corrected/results/visuals/`. GSE81089 RNA-seq is not silently
substituted for that cohort: its same-constructor official replay resolves only
from `adapter_final_corrected/secondary_validation/GSE81089/results/visuals/`
and is placed on a separately titled secondary-validation page. The combined
supplementary PDF consequently has eight pages while the output-file contract
remains seven PDFs.

The comparison atlas uses audited top-origin crop coordinates in
`tnbc_panel_map.tsv`. The preferred inputs are the main article and supplement
PDFs. Alternatively, pass the local reference directory containing:

- `s41467-026-76385-8_reference_36.png` through `_40.png` for Figure 1-5;
- `41467_2026_76385_MOESM1_ESM_01.png` etc. for supplementary figures.

The explanation block is separate from both images and comes from
`panel_explanations.tsv`. To replace it with Chinese text, pass a replacement
TSV with the same columns and `--font-path` pointing to a Unicode TTF. The run
fails instead of emitting missing-glyph squares when non-ASCII text is used
without such a font.

## Run on the cloud server

```bash
bash run_assembly.sh \
  /media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/official-code-port-v3/run-fig12-final \
  /media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/official-fig3-port/run-v3 \
  /media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/official-code-port-v3/run-fig45-v2-p0fix-20260821-r2 \
  /path/to/official-supplementary-figure3-atoms \
  /path/to/TNBC-main-article.pdf \
  /path/to/TNBC-supplement.pdf \
  - \
  /media/desk16/iy13202/projects/Paper2Paper/pilots/tnbc-chromatin-tf-nc-2026
```

To use the exported TNBC PNGs, pass `-` for both TNBC PDFs and pass their
directory in argument seven.

Python requires `pypdf`, `reportlab`, and `Pillow`. Final QA additionally
requires Poppler `pdftoppm` on `PATH`.

## Upstream release gate

The assembler discovers the Supplementary Figure 3
`verification_receipt.json`, hashes it, checks that it describes exactly 3A-L,
and binds its rows to the resolved PDF names and byte counts. The aggregate
decision is written to `assembly-audit/release_status.json`.

The current receipt is valid but explicitly has `publication_ready=false`:
3A is `PARTIAL`, while 3D, 3F, 3K, and 3L have official-constructor visual
overlaps. Those atoms are not edited. Consequently the current aggregate state
is `STRICT_OFFICIAL_PORT_REVIEW_ONLY`.

Default `run_assembly.sh` operation is review mode: seven PDFs may be generated
and technically verified, but the result is never reported as
`PUBLICATION_READY` or as a fully green strict release. Add
`--require-publication-ready` before the positional arguments to activate the
release gate; with the current receipt it exits non-zero before generating
final PDFs. It can be combined with `--dry-run` for a no-output release
preflight.

Before all atomic outputs are available, run the same command with
`--dry-run` as its first argument. Dry-run opens and validates every resolved
PDF/reference, writes `input_contract.tsv`, `resolved_atoms.tsv`, and
`missing_inputs.tsv` under
`execution/luad-official-code-port/assembly-audit/`, creates no final PDFs, and
exits non-zero if any required atom, reference crop, or atlas explanation is
missing. Thus a pending Figure 1/2 result is reported as a missing contract
item rather than being silently replaced by a provisional plot.

`validation_receipt.json` records the static checks and the current partial
dry-run snapshot. In that snapshot, the shared Figure 3, Figure 4/5, and
Supplementary Figure 3 roots resolve with zero required selector failures;
the non-zero dry-run exit is expected solely because Figure 1/2 atoms are not
yet present.

## Audit and QA

Assembly writes:

- `resolved_atoms.tsv`: exact source path, page count, SHA-256, and boundary;
- `input_contract.tsv`: atom-level root, frozen path scope, filename selector,
  match cardinality, vector-only policy, publication boundary, and consumer;
- `missing_inputs.tsv`: missing, invalid, or ambiguous required inputs;
- `declared_boundaries.tsv`: manual, unsupported, and supported-zero states;
- `embedding_manifest.tsv`: source-to-output placement and embedding method;
- `atlas_reference_crops.tsv`: TNBC source/page/crop provenance;
- `final_output_sha256.tsv` and `assembly_summary.json`;
- `release_status.json`: upstream receipt hash, partial/warning panels,
  blockers, aggregate release status, and publication-ready boolean.
- `source_identity_manifest.tsv` and `source_identity_summary.json`: exact
  receipt/commit/hash evidence for every resolved official atomic PDF.

`verify_assembled_pdfs.py` reopens every PDF, checks stable page counts and
boundary labels, re-hashes every atomic source, rejects non-vector LUAD atom
embeddings, renders every page with Poppler, calculates blank/black-page
metrics, and creates contact sheets under `qa/contact_sheets/`. Its strict mode
fails when any required atom or TNBC reference crop is missing, and also fails
if the output directory is not the frozen
`execution/luad-official-code-port/final/` publication path.
It also requires `release_status.json`, checks that its state is visibly stamped
into all seven PDFs, and records separate technical and release states. Review
mode reports `STRICT_OFFICIAL_PORT_REVIEW_ONLY`, not `PASS`; the explicit
release gate converts a negative upstream receipt into a hard failure.
QA additionally rejects any embedded source whose identity row is absent,
failed, or whose observed SHA changes after assembly.

The contact sheets are QA intermediates, not publication deliverables. Final
visual approval must inspect the latest rendered pages, especially crowded
heatmaps and graph labels, before any later release. The current seven-PDF set
is for review only.

# LUAD Figure 4–5 strict TNBC CodeOcean port

This directory is the only source of truth for the **strict official-code**
rendering pass of LUAD Figure 4, Figure 5, and Supplementary Figures 5–7.
It is deliberately separate from `code/publication_rebuild/`: no custom
renderer is imported here.

## Contract

1. The upstream source is TNBC CodeOcean capsule `7227095`, version 1,
   commit `edf5314`, released under the MIT licence.
2. `source_manifest.tsv` freezes the SHA-256 and source-line ranges used by
   every panel.
3. `materialize_official_blocks.py` verifies every upstream checksum before
   extracting the named R blocks. Verbatim blocks remain byte-for-byte equal
   to the corresponding upstream line range. Label-patched blocks receive a
   unified diff and are rejected if a change is outside the allow-list.
4. `prepare_luad_inputs.R` is a data adapter only. It renames LUAD columns to
   the objects expected by the official code; it does not change thresholds,
   select extra observations, or invent positive evidence.
5. `render_official_fig45.R` evaluates the materialized official blocks and
   supplies only LUAD data objects, output paths, cohort/endpoint labels, the
   four pre-frozen Kaplan–Meier TF choices, and four frozen LUAD drug–TF
   identifiers for the four published Supplementary Figure 7B call slots.
6. Figure 5G, Figure 5H, and Supplementary Figure 7D are not drawn. The LUAD
   frozen PDX screen contains 423 tests, but zero tests enter validation from
   a same-direction, cell-line-replicated drug–TF pair. Consequently neither
   the official PDX concordance panel, response waterfall, nor PDX AAC example
   has eligible input. Each is recorded as `UNSUPPORTED` in the render
   manifest and in `unsupported_panels.tsv`; no replacement graphic is
   permitted.

## Strict output scope

- Figure 4A–H
- Figure 5A–F
- Supplementary Figure 5A–B
- Supplementary Figure 6A–H
- Supplementary Figure 7A–C
- Figure 5G–H and Supplementary Figure 7D: manifest-only `UNSUPPORTED`

The TNBC capsule names only Supplementary Figure 6A–H. The previous custom
LUAD A–P expansion is therefore excluded from the strict port.

## P0 fidelity locks

- Supplementary Figure 5B executes official lines 87–247 verbatim, including
  Helvetica, official significance filtering, dynamic device dimensions, and
  the FDR-star footer. No font compatibility substitution is permitted.
- Supplementary Figure 6D/H retain the frozen `RFS` input object and path, but
  apply an audited visible-label-only `RFS` to `DFS` substitution for TCGA.
- Figure 4A–D retain the official `Multivariate` title qualifier.
- Supplementary Figure 7B executes all four official call/save blocks. The
  published default 7 by 7 inch `cairo_pdf()` device is unchanged.
- `audit/executed_official_blocks.tsv` records every block actually executed;
  `atomic_pdf_sha256.tsv` records every rendered atomic PDF and SHA-256.

## Run

```bash
python3 materialize_official_blocks.py \
  --official-code-root /path/to/TNBC_CodeOcean_7227095_v1/code \
  --output-dir /path/to/work/generated

Rscript prepare_luad_inputs.R \
  /path/to/luad-figure1 \
  /path/to/luad-figure4 \
  /path/to/luad-figure5 \
  /path/to/luad-publication-rebuild \
  /path/to/work/adapted

Rscript render_official_fig45.R \
  /path/to/work/generated \
  /path/to/work/adapted \
  /path/to/work/results

python3 verify_official_port.py \
  --code-dir . \
  --generated-dir /path/to/work/generated \
  --results-dir /path/to/work/results
```

The render script stops on schema or provenance failure. Biological absence
is not a runtime failure: it is recorded explicitly as `UNSUPPORTED`.
`run_official_fig45.sh` refuses an existing work or result directory so an
audited rerun cannot overwrite or silently mix with a superseded run.

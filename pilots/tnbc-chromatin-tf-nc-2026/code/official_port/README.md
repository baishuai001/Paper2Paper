# LUAD strict official-code port (TNBC CodeOcean edf5314)

This directory exists because the TNBC authors **did publish plotting code**.
The publication rebuild must therefore begin from that code, not from a new
house style.

## Non-negotiable boundary

- Files below `upstream_edf5314/code/` are byte-identical copies from
  CodeOcean capsule commit
  `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`.
- `adapt_luad_to_official_schema.R` contains data/schema adaptation only; it
  contains no plotting constructor.
- `materialize_official_runtime.py` changes only paths, LUAD/LUSC or dataset
  labels, frozen sample counts (22/13/19), and explicitly declared adapter
  object names. Every occurrence is checked and written to a machine-readable
  manifest and unified diff.
- The Figure 1 block runner evaluates the exact official line blocks. The
  Figure 2 runner calls the materialised official scripts directly.
- No theme, font size, colour value, geometry, statistical layer, clustering
  option, legend layout, or scale is redesigned in this port.
- A panel without a capsule constructor is `MANUAL_REQUIRED` and is not drawn.

## Publication boundaries discovered in the capsule

- Figure 1A: `MANUAL_REQUIRED`.
- Supplementary Figure 1 (all workflow/inclusion diagrams):
  `MANUAL_REQUIRED`.
- Figure 2A and Supplementary Figure 4A: `MANUAL_REQUIRED`.
- Figure 2C: `CAPSULE_PRINT_MISMATCH`. Official script `Figure2/06C` emits a
  four-block, TF-on-y heatmap (HC status, mean activity, per-sample promoter,
  per-sample motif). The printed paper uses a different three-block horizontal
  composition. This port preserves the capsule 06C output and does **not**
  silently rearrange it. The TCGA portion is explicitly paired multiomic
  `n=21/22`: `TCGA-44-A47F` has ATAC but no RNA/VIPER activity in the frozen
  inputs, so it is excluded only from Figure 2C without filling or imputation;
  Figure 2B and the promoter gate retain all 22 ATAC patients.
- Formal Figure 2C-E motif input is the frozen
  `motif_primary/anchor_consensus_author` run: 20 motif-testable HC-TFs,
  zero JASPAR-priority violations, and exactly three triple-system TFs
  (`FOXA3`, `NFATC4`, `XBP1`). The `candidate_anchor_equivalent` result with
  zero triple-system TFs is retained only as a `VERSION_BOUNDARY` sensitivity
  audit and never occupies the formal figure slot.
- Supplementary Figure 2F: the official source constructs its PDX long table
  from the TCGA object. The strict port retains and flags this capsule defect;
  it is not silently corrected.
- Supplementary Figure 3J-L: the official constructor exists, but the required
  genome-wide 500-bp annotation baseline is currently absent. No proxy is used.
- Supplementary Figure 4C: the LUAD JASPAR-only set is genuinely empty. The
  capsule's scalar database label produces a 0-versus-1 `data.frame` error for
  an empty set. Final-corrected runtime changes the three database labels to
  length-matched `rep` calls (`ANCHOR_EMPTY_SET_BUGFIX`); no dummy TF is added,
  and the eulerr constructor and palette are unchanged.
  Its JASPAR/CIS-BP inputs are isolated complete database-availability files;
  these are not mixed into the author-priority Figure 2C-D motif inputs.
- Figure 1E: ZNF737 is absent only from the 350-TF cross-system heatmap adapter
  because its 35 TCGA-regulon targets have zero overlap with the 19,221 DepMap
  expression genes. Figure 1B remains the full 351-TF discovery result; no
  activity value is filled or imputed.

See `panel_port_manifest.tsv` for every Figure 1/2 and Supplementary 1-4 panel.

## Run on the cloud server

```bash
bash run_official_port.sh \
  /path/to/luad-figure1-v2 \
  /path/to/luad-figure2 \
  /path/to/luad-publication-rebuild \
  /path/to/strict-official-port-run
```

Outputs are written only below the requested work root:

- `adapter/`: LUAD files expressed in the official expected schemas;
- `runtime/code/`: materialised official sources;
- `runtime/results/`: PDFs written by official constructors;
- `runtime/audit/runtime_source_manifest.tsv`;
- `runtime/audit/runtime_substitution_manifest.tsv`;
- `runtime/audit/runtime_manifest.json`;
- `runtime/audit/runtime_unified.diff`;
- panel execution receipts and logs.

## Static verification

```bash
python3 -m py_compile materialize_official_runtime.py
python3 materialize_official_runtime.py /tmp/adapter /tmp/runtime
```

Materialisation fails if any whitelisted source token occurs a different
number of times than recorded. This prevents accidental drift from the frozen
official source.

# Supplementary Figure 3A-L: strict official-code port

This port starts from the TNBC authors' CodeOcean capsule commit `edf5314`.
It does **not** replace the authors' plotting functions with a custom theme or
new chart implementation.

## Exact official mapping

| LUAD panel | Official source | Official plot block |
|---|---|---:|
| Supplementary Figure 3A-C | `Figure2/02J2-Plot_Number_Peaks.R` | 31-60 |
| Supplementary Figure 3D-F | `Figure2/02C-Peak_Saturation_Plot.R` | 44-155 |
| Supplementary Figure 3G-I | `Figure2/03D-Pearson_Correlation_Accessibility_ConsensusPeakSet.R` | 43-96 |
| Supplementary Figure 3J-L | `Figure2/04B-GenomicAnnotation.R` | 58-168 |

`OFFICIAL_SOURCE_LOCK.tsv` freezes every source SHA-256.  Runtime
materialisation permits only the absolute `plotting_aesthetics.R` path and the
official 04B argument-name typo (`consensus_peakset` versus
`consensus_peakset_dir`).  `runtime_unified.diff` and
`runtime_substitution_manifest.tsv` expose every changed byte.

## Data boundary

- Patient panel 3A is `PARTIAL`: TCGA LUAD uses the public fixed open-count
  peak matrix; it is not a newly called per-sample MACS2/IDR object.
- PDX and cell-line panels 3B-C use frozen MACS2 `q=0.01` peak counts. IDR is
  `N/A`, because the frozen cohorts do not contain qualifying true replicate
  sets. The official constructor's hidden fill label is not interpreted as
  evidence that IDR was performed.
- Panels 3D-F use the frozen 1000-permutation saturation results, expressed in
  the exact data objects expected by official 02C.
- Panels 3G-I use the frozen consensus-accessibility Pearson matrices.
- Panels 3J-L use the real 6,062,095-bin hg38 500-bp genomic baseline at
  `luad-figure2/results/supplementary_figure3/genome_annotation_baseline.tsv`.

See `PANEL_STATUS.tsv` for the panel-by-panel scientific boundary.

## Cloud run

```bash
bash run_official_supp3_cloud.sh \
  /path/to/TNBC_CodeOcean/code \
  /path/to/luad-figure2 \
  /path/to/this/supp3 \
  /path/to/new/official-supp3-run
```

The verifier requires all 12 one-page PDFs, 12 Poppler renders, the frozen
22/13/19 cohort sizes, and zero plotting-constructor rewrites.

## Visual QA boundary

The strict port intentionally does not move labels when LUAD values differ
from TNBC. Manual PNG review therefore records four capsule-layout conflicts:
3D and 3F have official saturation-label overlaps, while 3K and 3L have the
official fixed `n` label overlapping the promoter distribution. These outputs
are valid evidence that the official constructors ran, but they are explicitly
`publication_ready: false`. See `VISUAL_QA.tsv`; resolving those conflicts is a
separate, declared post-port layout task and must not be mislabeled as exact
official output.

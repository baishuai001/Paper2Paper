# LUAD Figure 3 strict official-code port

This directory does not implement replacement plotting functions. It executes
the TNBC authors' checksum-locked CodeOcean Figure 3 scripts after converting
the frozen LUAD inputs to the schemas those scripts expect.

Official source:

- CodeOcean capsule `7227095`, v1;
- frozen repository commit `edf5314`;
- MIT license in the capsule root;
- expected code root: a directory containing `Figure3/` and
  `Static_Scripts/plotting_aesthetics.R`.

`OFFICIAL_SOURCE_LOCK.tsv` records the exact official source file, SHA-256,
line count, panel and constructor/data-flow line ranges. At runtime,
`materialize_official_scripts.py` refuses to proceed if any source checksum has
changed. It writes every unified diff to `OUTPUT_ROOT/audit/official_diffs/` and
the source/patch manifests beside them.

## Main panels

| LUAD panel | Official executable source | Preserved behavior |
|---|---|---|
| Figure 3A | `04-Summarize_Regulons_HC_TR.R` | Four signed/shared regulon classes, motif tile, net-effect tile, official palette, geometry and dimensions |
| Figure 3B | `05-Generate_Regulon_Network.R` | Transposed target-by-TF incidence; target pairs require at least 3 shared activating HC-TFs; `graph_from_data_frame(edge_df_targets, directed=FALSE)` receives no explicit vertices; official Louvain hull/network constructor is the main panel |
| Figure 3C | `03-Pearson_Correlation_TRActivity.R` | Pearson correlation, official magma ComplexHeatmap, and official top/bottom 10 bar constructor |
| Figure 3D | `05-Generate_Regulon_Network.R` | All HC-TF pairs sharing at least 1 activated shared target; no percentile edge deletion, no display-only filter, no artificial isolates |
| Figure 3E | `05-Generate_Regulon_Network.R` | Regulon size versus unique partners from the same complete Figure 3D edge table |
| Figure 3F | `06-Activity_Score_Heterogeneity_Across_Patients.R` | Official skewness z approximation, BH FDR, point mappings and dimensions |
| Figure 3G | `06-Activity_Score_Heterogeneity_Across_Patients.R` | Official density/facet/Wilcoxon constructor; TNBC's four hard-coded example names are replaced only by the four highest-degree LUAD TFs from the complete Figure 3D graph |

The frozen LUAD topology is verified independently: Figure 3B must contain 6
endpoint nodes and 3 qualifying edges; Figure 3D must contain 31 nodes, 134
edges, one connected component and zero isolates. The sparse Figure 3B is not
replaced by a compact/custom diagram in the main output.

## Linked supplementary panels

- Supplementary Figure 4D: official `07B` volcano and `07C` heatmap
  constructors across TCGA, GSE41271, PDMR and DepMap.
- Supplementary Figure 5A: official `06` PDX/cell-line skewness constructor.
- Supplementary Figure 5B: byte-identical official `08D` and `08E` ESTIMATE
  volcano/heatmap constructors for TCGA and GSE41271.

The MKI67 and ESTIMATE correlations are frozen upstream statistical results;
the adapter changes their column/cohort names only. The official plotting code
reapplies `|r| > 0.4` and BH `FDR < 0.05` exactly as in the capsule.

## Permitted runtime diffs

The materializer permits only:

1. replacing the CodeOcean absolute `/code` source path with an environment
   path to the byte-identical official aesthetics file;
2. changing TNBC-facing labels to LUAD-facing labels;
3. uncommenting the authors' own table writers needed by the next official
   script and the independent verifier;
4. replacing the fixed denominator `94` with `nrow(HC_TRs)` in an unused
   reporting column;
5. replacing four hard-coded TNBC example TF names with four deterministic
   LUAD names selected from the official Figure 3D edge table;
6. supplying four LUAD cohort file paths to the otherwise unchanged official
   Supplementary Figure 4D heatmap constructor.

No `ggplot`, `ggraph`, `ComplexHeatmap`, scale, theme, layout, geometry, edge
threshold or graph constructor is added or changed by this port.

## Cloud invocation

```bash
bash run_official_fig3_port_cloud.sh \
  OFFICIAL_CODE_ROOT \
  FIGURE1_ROOT \
  FIGURE2_ROOT \
  FIGURE3_ROOT \
  PUBLICATION_ROOT \
  OUTPUT_ROOT \
  R_LIBRARY
```

The run fails closed on a checksum mismatch, missing 31-TF input, altered 3B or
3D topology, missing official panel PDF, invalid PDF header, incomplete diff
manifest or any plot constructor marked as changed.

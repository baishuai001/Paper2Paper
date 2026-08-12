# GSE234129 metadata file check

- Retrieved: 2026-08-08 from `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234129&file=GSE234129_meta.tsv.gz&format=file`.
- Local audit filename: `GSE234129_meta.tsv.gz` (raw file excluded from Git).
- Bytes: 137672.
- SHA256: `172F3041DFE5D36AB090DCE502E6598DD8F7F834E12C5C17B4F929915E70D029`.
- Parsed shape: 19,488 rows × 4 columns.
- Columns: `cell_barcodes`, `patient`, `sample`, `celltype`.
- Unique identities: 6 patients; 17 samples; 62 cell-state labels.
- Key finding: the labels cover immune and stromal compartments but contain no epithelial/malignant state. The resource can support a gastric TME layer, not the anchor's malignant cNMF layer by itself.
- Metastatic sampling: one ovarian sample (`MDA_Pt5-Ov`, 1,477 cells) and one liver sample (`MDA_Pt9-Li`, 1,118 cells); these are inadequate for a donor-level general metastasis claim.

# LUAD publication rebuild: scientific narrative, figure legends, and TNBC panel audit

## Scope and central claim

The LUAD study asks whether a transcription-factor activity state discovered from large patient transcriptomes can be prioritized by chromatin accessibility, retained across patient tumors, PDX models, and cell lines, and then related to intertumoral heterogeneity, outcome, and treatment response.

The frozen evidence chain is:

1. TCGA-LUAD versus TCGA-LUSC identifies 158 LUAD-specific TF activities by ARACNe3/msVIPER.
2. Independent patient cohorts reproduce 95/158 TFs in GSE81089 RNA-seq, 70/158 in GSE41271 microarray data, and 44/158 in both external cohorts.
3. The TNBC-style promoter-accessibility and three-system activity rule prioritizes 31 chromatin-informed TFs (HC-TFs).
4. Of the 20 HC-TFs with a testable motif, FOXA3, NFATC4, and XBP1 satisfy the strict three-system motif rule. ETV1 and ZNF75D join this set only in the declared sensitivity analysis (HOMER q<=0.05 and support in at least half of samples in every system).
5. All 31 HC-TFs—not only the motif-positive subset—enter network, heterogeneity, survival, and pharmacogenomic analyses, matching the anchor paper's use of all 94 HC-TFs.
6. Clinical and treatment results are asymmetric: prognosis is mainly exploratory and only TCGA OS–ZNF75D survives the frozen FDR screen; five cell-line drug–TF associations replicate across at least two pharmacogenomic resources, but none obtains matched public PDX validation.

The defensible manuscript theme is therefore not “a broad motif program conserved in every LUAD model.” It is:

> A chromatin-informed TF-activity framework defines a LUAD regulatory network with a selectively conserved cross-system motif core and reveals testable links between regulatory heterogeneity, patient outcome, and therapeutic response.

The words “patient outcome” and “therapeutic response” must remain qualified as exploratory until independent clinical and in-vivo validation becomes available.

## What was learned from the official TNBC code

The Code Ocean capsule `v1.0` (`edf5314`) was read panel by panel. Its operative visual and statistical grammar was retained:

- Figure 1 separates network discovery, an independently reconstructed patient network, and projection into PDX/cell-line systems.
- Figure 2 defines HC-TFs at the promoter/activity gate; distal motif enrichment is a subsequent orthogonal annotation and is not the entry rule for Figures 3–5.
- Figure 3 uses signed regulons, shared activated targets, target-gene communities, TF collaboration, skewness, MKI67 correlations, and microenvironment-confounding checks.
- Figure 4 in the official code uses dichotomized TF activity with script-specific cutpoints: the METABRIC `NewSelection` analysis compares the outer tertiles and the TCGA analysis uses maximally selected rank cutpoints (with a median fallback). The frozen LUAD primary screen instead uses a cohort/endpoint median split for reproducibility and keeps maximally selected cutpoints only for descriptive Kaplan–Meier displays. This adaptation is explicit; complete screens and 5,000-permutation nulls remain supplementary.
- Figure 5 performs within-resource correction, cross-resource replication, drug/pathway summaries, example scatterplots, and independent PDX response validation.

Two verified coding defects in the capsule were not copied: the PDX long table for one Figure 1 supplement is built from a TCGA object, and the cell-line PAM50 Fisher output writes a PDX object. The LUAD implementation records inputs and outputs explicitly instead.

## Main Figure Legends

### Figure 1 | A transferable LUAD TF-activity framework

**A,** Study design. Separate ARACNe3 networks are reconstructed from TCGA RNA-seq and GSE41271 microarray data; GSE81089 provides a second RNA-seq patient validation, while the frozen TCGA regulon is projected into PDMR PDX and DepMap cell-line transcriptomes. The workflow reports cohort sizes and the 158-TF discovery set, 95- and 70-TF external support, and 44 TFs reproduced in both external patient cohorts. **B,** Relationship between LUAD-versus-LUSC differential gene expression and msVIPER TF activity in TCGA. TF activity is not treated as a surrogate for TF mRNA abundance. **C,** Patient-level activity of the 158 LUAD-specific TFs in TCGA-LUAD/LUSC with lung-cancer clinical and molecular annotations. Activity is centered and symmetrically clipped so the neutral midpoint remains visible. **D,** Projection of the LUAD TF program into independently curated PDMR PDX models. **E,** Projection into DepMap lung-cancer cell lines. Columns are independent biological models; annotations disclose diagnosis/state and available provenance rather than implying equivalence to patient tumors.

### Figure 2 | Chromatin prioritization and selective motif conservation

**A,** ATAC design across 22 primary patient tumors, 13 LUAD PDX models, and 19 LUAD cell lines. **B,** Intersection of the 158 discovery TFs whose promoters are accessible in at least half of samples within each system. **C,** Aligned evidence for the 31 HC-TFs: mean LUAD TF activity, the fraction of samples with an accessible TF promoter, and strict distal motif support in patients, PDX, and cell lines. A TF is excluded by activity only when its mean NES is negative in all three systems. **D,** Cross-system intersections of motif-supported HC-TFs under the primary HOMER rule. **E,** Sample prevalence for the strict three-system motif set (FOXA3, NFATC4, XBP1) and the explicitly labeled q<=0.05/50% sensitivity set (adding ETV1 and ZNF75D). JASPAR motifs are used first; CIS-BP is used only if no JASPAR motif exists for that TF.

### Figure 3 | Network architecture and intertumoral regulatory heterogeneity

**A,** Activated/repressed and shared/private target composition of every 31-HC-TF regulon; strict three-system motif TFs are annotated rather than used as a filter. **B,** Target-gene graph reconstructed according to the anchor logic. Thirty activated targets are regulated by at least three HC-TFs, but only three target pairs share at least three regulators; consequently the graph contains three two-node components and 24 isolated nodes. This sparse topology is the result and does not support a broad convergent target-gene community. GO enrichment is not performed for one- or two-gene components. **C,** Pairwise activity correlation of the 31 HC-TFs across TCGA-LUAD tumors. **D,** TF collaboration network based on shared activated targets. **E,** Relationship between regulon breadth and the number of collaborating TFs. **F,** Skewness of HC-TF activity across TCGA-LUAD, with BH-adjusted skewness inference. **G,** Representative activity distributions illustrating positively skewed, negatively skewed, and approximately symmetric regulatory states.

### Figure 4 | Exploratory clinical associations of the LUAD HC-TF network

**A-D,** Compact multivariable Cox forest plots for GSE41271 overall survival, GSE41271 recurrence-free survival, TCGA overall survival, and TCGA disease-free survival. All 31 HC-TFs were screened after dichotomizing activity at the cohort/endpoint median; only nominally associated TFs (or the prespecified top five when none is nominal) enter the compact main panels, while subtitles disclose the full denominator and FDR count. **E, G,** Overlap of outcome-associated TFs across OS/RFS in GSE41271 and OS/DFS in TCGA, respectively. **F, H,** Representative Kaplan–Meier curves in GSE41271 and TCGA using maximally selected rank cutpoints. These display cutpoints differ from the frozen primary median split and the curves remain descriptive. Complete Cox screens and permutation nulls appear in Supplementary Figure 6.

### Figure 5 | Replicated pharmacogenomic associations and the in-vivo validation gap

**A-C,** Concordance-index screens relating HC-TF activity to drug response in GDSC2, CTRPv2, and PRISM, with BH correction within each resource. **D,** True UpSet representation of significant drug–TF memberships across resources. **E,** Effect-direction heatmap for the five drug–TF pairs that reproduce in at least two resources with a consistent direction. **F,** Drug/pathway/TF bubble plot for the replicated associations. All five indicate higher TF activity associated with resistance: 5-fluorouracil–ZNF254, 5-fluorouracil–ZNF540, PHA-793887–WWC2, vorinostat–ZNF254, and dactolisib–ZNF254. **G,** PDX validation funnel, from state-classified NSCLC PDX models through evaluable drugs/tests to validated pairs. **H,** Complete PDX concordance screen. Zero validated pairs are displayed as a result and a data-coverage limitation, not replaced by an empty placeholder or rephrased as evidence of biological absence.

## Supplementary Figure Legends

### Supplementary Figure 1 | Patient-cohort inclusion

**A-C,** Inclusion and exclusion flows for TCGA, GSE81089, and GSE41271. Each panel distinguishes downloaded profiles, independent tumor samples, accepted histology, matched expression/clinical information, and the final LUAD/LUSC comparison. GSE81089 contains 108 LUAD and 67 LUSC tumors after correcting the two sample-name aliases that were silently lost in the first pass.

### Supplementary Figure 2 | Independent TF discovery and cross-model transfer

Independent GSE81089 and GSE41271 differential-expression/msVIPER panels, overlaps with the 158-TF TCGA discovery set, within-cohort sample-activity coherence for TCGA, PDX, and cell lines, and cross-cohort TF effect distributions. RNA-seq and microarray values are never merged into one expression matrix; networks are reconstructed or projected within the declared platform-specific design.

### Supplementary Figure 3 | ATAC quality, saturation, and genomic distribution

**A-C,** Canonical peak counts and available raw-QC metrics for patient, PDX, and cell-line ATAC data. **D-F,** One-thousand-permutation peak-saturation curves with asymptotic fits and the sample numbers corresponding to 50%, 90%, 95%, and 99% of the fitted asymptote. **G-I,** Pairwise accessibility correlations within each system. **J-L,** Fractions of promoter, exonic, intronic, and distal peaks under one hierarchical annotation definition. Diagnostic QC metrics are reported but are not used as an unauthorized project-stopping rule.

### Supplementary Figure 4 | HC-TF selection, motif availability, and proliferation

**A,** Frozen 158-to-31 filtering flow. **B,** Combinations of promoter accessibility and mean-activity categories across the three systems. **C,** Motif availability and JASPAR-priority/CIS-BP-fallback status for HC-TFs. **D,** Pearson correlations between the activity of all 31 HC-TFs and MKI67 expression in TCGA patients, GSE41271 patients, PDMR PDX, and DepMap cell lines; stars denote |r|>0.4 and within-cohort BH FDR<=0.05. The additional GSE81089 results remain in the table as a sensitivity cohort.

### Supplementary Figure 5 | Regulatory heterogeneity and tumor-microenvironment confounding

**A,** Cross-system HC-TF activity skewness in TCGA, GSE41271, PDMR PDX, and DepMap cell lines; stars denote BH-adjusted skewness tests. **B,** Pearson correlations between HC-TF activity and ESTIMATE stromal score, immune score, and tumor purity in TCGA and GSE41271 LUAD tumors. This panel tests whether patient associations could be explained primarily by nonmalignant-cell content.

### Supplementary Figure 6 | Complete outcome screens and permutation tests

**A-H,** All 31 HC-TFs in univariable and multivariable median-dichotomized Cox screens for OS/RFS in GSE41271 and OS/DFS in TCGA. Purple denotes BH FDR<=0.05, orange denotes nominal p<=0.05, and gray denotes nonsignificant results. **I-P,** Null distributions from 5,000 outcome-label permutations for each cohort/endpoint/model combination; the red line is the observed number of nominally significant TFs, and each panel reports the empirical enrichment p value.

### Supplementary Figure 7 | Expanded pharmacogenomic evidence

**A,** Eligible-drug overlap across GDSC2, CTRPv2, and PRISM. **B,** Scatterplots for every dataset that contributes to each of the five replicated drug–TF pairs; panels report model count, Spearman correlation, and within-resource FDR. **C,** Complete bubble matrix of all within-resource BH-significant drug–TF associations, with drug pathway/mechanism annotations. The figure may span pages to preserve readable labels rather than suppress associations for layout convenience.

### Extended Supplementary Figure 8 | Model and classifier audit

Cell-line identifier mapping, drug eligibility, frozen LUAD/LUSC classifier variants, and classifier scores for patient, PDX, and cell-line systems. This figure documents the assumptions needed to call a model LUAD-like and prevents histology/state relabeling after outcome inspection.

### Extended Supplementary Figure 9 | PDX validation audit

PDX drug coverage, model counts, all evaluable PDX concordance tests, and the zero-row validated-pair table. It distinguishes “no matched/evaluable public data,” “tested but not significant,” and “opposite direction,” which are biologically different failure modes.

## Panel-by-panel interpretation and visual audit

The comparison atlas contains 40 audit pages: 33 pages compare every main Figure 1–5 subpanel, and seven additional pages compare the complete canonical Supplementary Figures 1–7. Each page places the unmodified TNBC crop beside the independently generated LUAD panel. The concise main-panel conclusions are:

| Panel | TNBC biological role | LUAD interpretation | Explanatory/visual difference and required correction |
|---|---|---|---|
| 1A | Defines discovery, independent validation, PDX, and cell-line branches. | Defines the same analytical roles for LUAD. | LUAD has two smaller external patient cohorts rather than METABRIC's 1,980 tumors. Equal-width stacked branches and restrained result badges prevent the workflow from dominating the page. |
| 1B | Shows TF activity contains information beyond TF expression. | Tests the same claim in LUAD versus LUSC. | The logic is equivalent; the LUAD panel must be enlarged and use symmetric axes/category counts. |
| 1C | Displays patient heterogeneity with receptor/PAM50 annotations. | Displays LUAD/LUSC activity with stage, smoking, sex, age, driver/state annotations. | Lung annotations are biologically analogous, not identical. A centered, quantile-clipped color scale is required to avoid the former red saturation block. |
| 1D | Transfers the TNBC program to PDX. | Transfers the LUAD program to PDMR PDX. | LUAD PDX sample size and diagnosis granularity are weaker; donor independence and state labels must be visible. |
| 1E | Transfers the program to breast cell lines. | Transfers it to DepMap lung lines. | Cell-line identity is less clean than the anchor's curated breast panel; model eligibility and lineage state must be exposed. |
| 2A | Introduces the three-system ATAC design. | Introduces 22/13/19 LUAD patient/PDX/cell-line models. | The LUAD PDX arm is smaller; model counts and independent-unit definitions must be prominent. |
| 2B | Intersects promoter-open TFs across systems. | Applies the same half-sample promoter rule to 158 TFs. | A real UpSet is needed; bars that merely resemble an UpSet do not explain membership. |
| 2C | Defines the 94 TNBC HC-TFs. | Defines 31 LUAD HC-TFs. | The first LUAD panel showed all 158 TFs and was unreadable. The rebuild aligns activity, promoter fraction, and motif status only for the 31 selected TFs while preserving the full 158-TF gate in tables. |
| 2D | Intersects motif-supported TFs across systems. | Shows selective LUAD motif conservation. | LUAD has only three strict triple-system TFs, so the figure must emphasize exact sets and denominators rather than graphical area. |
| 2E | Shows per-system motif enrichment for selected TFs. | Shows strict and q<=0.05/50% motif prevalence. | The former LUAD panel mixed many noncentral motifs. The rebuild focuses on the three strict and two sensitivity-only TFs and prints n/N. |
| 3A | Quantifies signed/private/shared regulons. | Quantifies the 31-HC-TF regulons. | The evidence is structurally equivalent but smaller; common color semantics and a motif annotation avoid treating motif as the inclusion gate. |
| 3B | Reveals target-gene communities. | Applies the same >=3-regulator node and >=3-shared-regulator edge rules, yielding 30 nodes but only three edges. | A GO dot plot cannot substitute for topology, but the corrected topology is biologically weak: three two-node components and 24 isolated nodes. The rebuild shows that sparsity rather than lowering the edge rule for appearance; GO is not interpreted for components with fewer than three genes. |
| 3C | Shows HC-TF co-activity structure. | Shows LUAD TF co-activity in TCGA tumors. | The LUAD matrix is smaller and should be square, clustered, and readable at final size. |
| 3D | Shows TF collaboration through shared targets. | Shows the LUAD collaboration network. | Sparse edges are a biological result; node/edge legends and a declared display threshold are required. |
| 3E | Links regulon size and collaboration. | Tests whether broad LUAD regulons collaborate more. | Similar explanatory grammar; labels must be restricted to informative outliers. |
| 3F | Identifies intertumoral activity skewness. | Identifies LUAD HC-TFs active in subsets of tumors. | LUAD must report BH correction and avoid making skewness synonymous with a discrete subtype. |
| 3G | Illustrates representative activity distributions. | Shows ELF3/MAGED2/PHC2 and other representative shapes. | Density panels need matched axes, sample counts, and explicit reasons for representative selection. |
| 4A | Shows adjusted OS associations in the anchor external cohort. | Summarizes the 31-TF median-dichotomized GSE41271 OS screen. | The LUAD rule is a reproducible median split, not the anchor METABRIC code's outer-tertile selection. The compact forest shows nominal hits while the subtitle and supplement retain the 31-TF denominator. |
| 4B | Shows adjusted recurrence associations in the anchor external cohort. | Summarizes the 31-TF median-dichotomized GSE41271 RFS screen. | External recurrence evidence is weak; matching axes and evidence colors with 4A prevent visual overstatement. |
| 4C | Shows adjusted OS associations in the anchor TCGA cohort. | Summarizes the 31-TF median-dichotomized TCGA-LUAD OS screen. | Only ZNF75D survives BH correction, and the LUAD split differs from the anchor code's maxstat cutpoint. The single FDR result is labeled without hiding null TFs in the supplement. |
| 4D | Shows adjusted recurrence associations in the anchor TCGA cohort. | Summarizes the TCGA-LUAD DFS screen. | The screen is largely null. A compact panel and explicit FDR count are more truthful than a visually dominant all-null forest. |
| 4E | Shows outcome-signal overlap in the anchor external cohort. | Shows direction-aware GSE41271 OS/RFS overlap. | Exact membership, expected overlap and permutation support are needed because the sets are small. |
| 4F | Illustrates selected external-cohort survival differences. | Shows GSE41271 KM curves with maximally selected display cutpoints. | The curves are descriptive and use different cutpoints from the primary median-dichotomized Cox screen; risk tables, censoring and selection method are therefore mandatory. |
| 4G | Shows outcome-signal overlap in the anchor TCGA cohort. | Shows direction-aware TCGA OS/DFS overlap. | Exact counts and the 31-TF denominator prevent a small overlap from being read as a broad outcome-stable program. |
| 4H | Illustrates selected TCGA survival differences. | Shows TCGA-LUAD KM curves with maximally selected display cutpoints. | Curves sit below the adjusted Cox evidence, include risk tables, and are not described as independent replication. |
| 5A | Screens drug–TF relationships in the first pharmacogenomic resource. | Screens LUAD associations in GDSC2. | The common CI axis, BH rule and line/drug/test denominators make this panel comparable to 5B/C. |
| 5B | Screens the second pharmacogenomic resource. | Screens LUAD associations in CTRPv2. | Dataset coverage differs, so identical scale and annotation rules are needed before judging signal density. |
| 5C | Screens the third pharmacogenomic resource. | Screens LUAD associations in PRISM. | PRISM has a different response distribution; only same-direction cross-resource replication advances to 5D/E. |
| 5D | Shows cross-resource intersection. | Identifies replicated LUAD drug–TF pairs. | The first pass used a membership bar chart. A true UpSet restores set logic. |
| 5E | Summarizes replicated effects. | Shows five reproducible resistance associations. | Direction and dataset identity must be visible simultaneously; clustering must not imply more pairs than observed. |
| 5F | Places associations in drug/pathway context. | Links LUAD replicated pairs to mechanisms. | A generic integration bar chart was insufficient. The rebuild uses drug/pathway/TF bubbles with effect and FDR encodings. |
| 5G | Validates a selected association in PDX. | Audits whether any public LUAD PDX pair is testable and validated. | LUAD lacks an anchor-equivalent successful experiment. A structured funnel honestly explains where validation is lost. |
| 5H | Shows PDX response ordering/waterfall. | Shows the complete PDX null/negative screen. | A blank placeholder is unacceptable. The rebuilt panel plots every evaluable test and labels zero validation explicitly. |

## Output locations

Final PDFs are written under:

`pilots/tnbc-chromatin-tf-nc-2026/execution/luad-publication-rebuild/results/publication/`

The complete tabular evidence package and workbook are written under:

`pilots/tnbc-chromatin-tf-nc-2026/execution/luad-publication-rebuild/results/tables/`

Raw cloud computation receipts remain under the server work root:

`tmp/tnbc-chromatin-tf-nc-2026/luad-publication-rebuild/`

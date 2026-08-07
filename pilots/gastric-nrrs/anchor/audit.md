# Anchor audit

## Paper identity

- Full citation: Qiu L, Yao S, Yang Z, et al. Development and validation of a novel nerve-related prognostic model for gastric cancer based on bulk and single-cell RNA sequencing data. *BMC Cancer*. 2025;25:1738.
- DOI: `10.1186/s12885-025-15202-9`
- Public article: `https://pmc.ncbi.nlm.nih.gov/articles/PMC12604206/`
- Local source material: PDF, page images and Supplementary Material 1 are stored outside Git.
- Public data statement: TCGA, UCSC Xena and GEO are named generically; the paper gives bulk accessions in Methods/Table S1. The authors' scRNA-seq cohort is deposited as `HRA011727` in GSA-Human. No code-availability statement or code repository is provided.

## Plain-language reconstruction

The paper starts with 441 genes assembled from KEGG nerve-related pathways. It keeps genes associated with overall survival at nominal `P < 0.05` in both TCGA-STAD and a merged GEO cohort, takes the 29-gene intersection, and fits a LASSO Cox model in TCGA-STAD. The published eight-gene score is:

```text
NRRS = 0.040*AGT - 0.040*EPHB3 + 0.018*GNAI1 - 0.102*LPAR2
       + 0.014*LRRC4C + 0.072*NPY1R + 0.175*NRP1 + 0.006*SEMA6A
```

The authors convert expression to a gene-wise z-score within every dataset, divide patients into high and low NRRS groups, test survival and clinical associations, and add pathway, immune, mutation, drug-response and single-cell localization analyses.

## Scientific grammar

- Population: patients with gastric cancer in TCGA-STAD and ten GEO cohorts; two immunotherapy cohorts, one gastric and one urothelial; an eight-patient institutional gastric scRNA-seq cohort plus four published scRNA-seq collections.
- Biological unit: patient for survival, mutation and response claims; donor for tissue-level single-cell comparisons; cell only for localization and descriptive expression.
- Focal object: an eight-gene linear nerve-related risk score derived after outcome-guided screening.
- Main comparison: high versus low NRRS, with an undocumented cutoff; continuous NRRS is used in some correlations and Cox models.
- Primary outcome: overall survival. Recurrence-free survival, immunotherapy response and inferred drug sensitivity are secondary or exploratory.
- Central published claim: NRRS is a stable prognostic model and may identify gastric-cancer patients likely to benefit from chemotherapy or immunotherapy.
- Defensible public-data claim ceiling: a locked eight-gene score may show retrospective association with survival in some public cohorts. The available design cannot establish treatment benefit, mechanism, prospective individual prediction or a clinically transportable cutoff.
- Falsifier: after locking gene identifiers, expression transformation, missing-gene handling and cutoff, the score has no patient-level association or calibration in a cohort untouched by screening, fitting and threshold selection.

## Data inventory and role audit

### Bulk cohorts

| Cohort | Reported n | Platform | Stated role | Audit note |
| --- | ---: | --- | --- | --- |
| TCGA-STAD | 348 | RNA-seq | screening, LASSO fitting, nomogram, downstream biology | training and discovery, not validation |
| GSE62254 | 300 | GPL570 | component of GEO-meta and separate validation | reused after contributing to feature screening |
| GSE15459 | 192 | GPL570 | component of GEO-meta and separate validation | reused after contributing to feature screening |
| GSE57303 | 70 | GPL570 | component of GEO-meta and separate validation | reused after contributing to feature screening |
| GSE29272 | 126 | GPL96 | OS validation | apparently not used for feature selection |
| GSE84437 | 433 | Illumina HT-12 v3 | OS validation | apparently not used for feature selection |
| GSE14210 | 123 | GPL571 | OS validation | apparently not used for feature selection |
| GSE34942 | 56 | GPL570 | OS validation | apparently not used for feature selection |
| GSE13861 | 65 | Illumina WG-6 v3 | OS and RFS validation | endpoint definitions need source verification |
| GSE26942 | 202 | Illumina HT-12 v3 | OS and RFS validation | endpoint definitions need source verification |
| GSE26253 | 432 | Illumina DASL | RFS validation | not an OS cohort in this paper |

The stated total of 2,347 equals the sum of all 11 cohorts, but the paper's use of the term "independent validation" is too broad. `GSE62254`, `GSE15459` and `GSE57303` were merged into GEO-meta for outcome-guided feature screening and then shown individually as validation datasets. They are independent patients from TCGA but not independent of model development.

### Single-cell and treatment-response data

- `HRA011727`: authors' eight-patient matched gastric cohort. Exact public/controlled access status and downloadable files require live verification. Exact reproduction of Fig. 9 and parts of Supplementary Figs. S10-S11 depends on it.
- Public scRNA-seq validation: Sathe cohort plus `GSE167297`, `GSE163558` and `GSE183904`. Table S3 gives no accession for the Sathe cohort and no per-donor mapping or processing manifest.
- `PRJEB25780`: 45 gastric-cancer patients treated with PD-L1 inhibition; response comparison requires expression-to-clinical ID mapping and the original response definitions.
- `IMvigor210`: urothelial, not gastric, cancer. It can test cross-cancer association but cannot directly validate a gastric-cancer treatment-benefit claim.
- TIDE, SubMap, IPS and GeneMANIA are dynamic or externally hosted resources. Query date, version, parameters and returned files are not reported.

## Computation contract gaps

The published formula is not by itself computable without undocumented choices:

1. The high/low cutoff is not stated. Figure group sizes are not consistent with a simple median split.
2. "Expression converted to a z-score" does not specify orientation, denominator, handling of missing values or whether scaling occurred before or after cohort merging and batch correction.
3. Microarray probe-to-gene mapping and multi-probe aggregation are not specified.
4. TCGA count-to-TPM implementation, gene annotation version and duplicate-gene handling are not specified.
5. LASSO folds, seed, standardization, censor coding and exact package/version are absent.
6. Time-dependent ROC, C-index, calibration, decision-curve and nomogram implementations and parameters are absent.
7. DEG thresholds, GSEA ranking statistic, gene-set versions, multiplicity control and online-service snapshots are incomplete.
8. Single-cell sample-to-patient mapping, doublet handling, normalization, variable features, clustering resolution, marker list and donor-aware testing are absent.

Because every cohort is standardized separately, a patient's NRRS depends on the other patients in that cohort. This supports retrospective cohort stratification but not a deployable individual-patient score unless a reference distribution is frozen.

## Figure-to-evidence map

| Figure | Manuscript role | Required inputs | Main reproducibility limit |
| --- | --- | --- | --- |
| Fig. 1 | feature screening and eight-gene formula | TCGA, GEO-meta, 441-gene list | outcome leakage roles, missing LASSO seed/folds and probe mapping |
| Fig. 2 | survival and multivariable association | five displayed cohorts plus clinical covariates | cutoff and covariate coding undocumented; three cohorts participated in screening |
| Fig. 3 | score distribution, PCA and time ROC | TCGA, GEO-meta, GSE62254 | cutoff and ROC implementation undocumented; PCA separation is partly constructed by grouping on the score |
| Fig. 4 | nomogram, calibration and DCA | TCGA clinical data | package, resampling and validation details absent |
| Fig. 5 | DEG, network and enrichment | TCGA score groups, gene sets, GeneMANIA | thresholds, versions and mutable service output absent |
| Fig. 6 | immune infiltration | TCGA and GEO-meta expression | signatures, package versions, multiplicity and patient-level source tables absent |
| Fig. 7 | immunomodulators, mutation, TMB and IPS | TCGA expression/MAF plus IPS | exact sample matching and IPS source/version absent |
| Fig. 8 | inferred and observed immunotherapy response | TIDE/SubMap, PRJEB25780, IMvigor210 | dynamic services, cross-cancer validation and response mapping |
| Fig. 9 | scRNA localization | HRA011727 eight-patient cohort | access restriction and undocumented single-cell pipeline choices |
| Figs. S1-S12 | batch, extra cohorts, drugs and scRNA extensions | multiple public/private inputs | many are conceptually reproducible, but exact outputs require unavailable choices/code |

## Module disposition

- **Retain:** published eight-gene coefficients; Supplementary Table S4; exact public bulk accessions; patient-level survival as the statistical unit.
- **Repair:** dataset-role leakage labels; probe/gene mapping; preprocessing contract; fixed cutoff; continuous-score evaluation; resampling; multiplicity; cohort-independent external validation; donor-aware single-cell inference.
- **Substitute:** use maintained, versioned local implementations for deprecated or dynamic services where possible; use publicly downloadable scRNA-seq data when private HRA access is unavailable.
- **Extend:** save cohort manifests, raw-to-analysis ID maps, locked score specification, source tables, negative results and sensitivity analyses.
- **Drop from a minimum paper:** decorative PCA/alluvial panels, network graphics and drug-sensitivity inference when they do not support the central claim.
- **Blocked for exact reproduction:** private/controlled scRNA input, undocumented author choices, original code, exact dynamic-service outputs and any wet-lab evidence.

## Reproduction boundary

The public materials are sufficient to attempt an analysis reproduction of the fixed NRRS in selected bulk cohorts and to test whether the reported survival direction can be recovered. They are insufficient for pixel-identical figures, exact model reconstruction, an uncontaminated validation of the development procedure, or exact reproduction of the institutional single-cell and online-service modules. The example should therefore begin with one public-cohort score-computation spike and treat all stronger claims as conditional.

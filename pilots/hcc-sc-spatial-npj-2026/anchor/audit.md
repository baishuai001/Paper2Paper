# Anchor audit

## Paper identity

- Title: *Large-scale single-cell analysis and in silico perturbation reveal dynamic evolution of HCC: from initiation to therapeutic targeting*.
- DOI: `10.1038/s41698-026-01307-2`; published 30 January 2026 in *npj Precision Oncology*.
- Local sources reviewed: the 18-page article PDF, the 19-page supplementary PDF, and Supplementary Data 1–6 workbooks in `WorkSpace/HCC_sc_spatial_npj_2026/PDF与图片/`.
- Public data statement: six GEO single-cell accessions, four older HCC spatial samples at a Lifeome link, GSE238264 spatial treatment samples, GSE109211 and TCGA-LIHC. Supplementary workbooks expose sample descriptions, program genes, subtype DEGs and Geneformer-ranked results.
- Public code statement: custom analysis code and other data are available only from the corresponding author on reasonable request. No runnable author repository or frozen environment was found in the supplied materials or article page.

## Scientific grammar

- Problem: connect recurrent malignant-cell states, tumor-microenvironment states and spatial niches to HCC progression and computational therapeutic vulnerabilities.
- System: a cross-sectional integration of public normal liver, primary HCC, portal-vein tumour thrombus (PVTT) and metastatic lymph-node (MLN) samples, plus independent bulk and spatial cohorts.
- Main comparisons: Normal versus Tumor versus PVTT versus MLN; malignant meta-programs; primary versus advanced tissue; treatment responder versus non-responder; source-state versus desired-state virtual perturbations.
- Primary outcomes: recurrent malignant meta-program identity/usage, cell-state abundance, spatial co-localization and communication, survival/treatment association, and ranked virtual knockout candidates.
- Central claim: recurrent tumor-intrinsic programs and stromal/immune niches track HCC progression, and convergent in-silico perturbation identifies HSP90B1 as a candidate vulnerability.
- Honest claim ceiling: cross-sectional state association and computational target prioritization. The public evidence cannot prove temporal evolution, tumor initiation, causal cell-cell interaction or therapeutic efficacy.
- Falsifiers: programs disappear under donor/dataset holdout; stage differences vanish with donor-level statistics; spatial associations fail across samples or predefined regions; virtual-perturbation rankings are unstable to split/seed/state definition or do not outperform simple baselines; candidates lack tumor-selective perturbational support.

## Anchor framework assets

The paper's valuable asset is the **ordered evidence chain**, not any particular HCC marker:

1. Assemble multiple single-cell cohorts representing clinically meaningful states and make an explicit sample manifest.
2. Identify malignant cells, derive recurrent within-tumor expression programs, consolidate them into interpretable meta-programs and test robustness across donors/cohorts.
3. Place those programs in a cellular ecosystem by resolving immune, endothelial and fibroblast states and their donor-level changes.
4. Project programs and cell states into spatial data, then test whether specific tumor–stroma neighborhoods recur across specimens.
5. Link the programs to independent prognosis and treatment cohorts without reusing discovery outcomes.
6. Define biologically explicit source and target states, run virtual perturbation, intersect robust candidates across relevant compartments, and validate prioritization with orthogonal dependency/expression/treatment evidence.
7. End with a claim hierarchy: atlas/state description → spatial association → clinical association → candidate dependency. Each layer has a different evidential strength.

This framework is transferable to another cancer only if each retained layer has data of the required role. A route may deliberately retain only a subset, but it must then be described as a partial framework imitation rather than the full anchor paper.

## Data and cohort anatomy

- Six scRNA-seq studies contributed 183 described raw samples. The supplementary sample sheet lists 19 GSE125449, 5 GSE134355, 21 GSE149614, 46 GSE151530, 58 GSE156625 and 34 GSE189903 samples.
- The same sheet contains 160 Tumor, 19 Normal, 2 PVTT, 1 MLN and 1 adjacent non-tumor row. Pathology labels identify at least 39 clearly non-HCC/ICC-related rows. The paper reports 115 retained pure-HCC samples after pathology and QC, but does not publish the exact retained 115-sample list or an exclusion reason for every discarded sample.
- Supplementary Data 2–5 provide top cNMF genes and endothelial/fibroblast DEG tables. Supplementary Data 6 provides three Geneformer result tables, not tokenized inputs, trained models, code, seeds or split manifests.
- Most processed single-cell GEO files and the bulk cohorts are publicly retrievable. Some raw data are controlled or absent. The old Lifeome spatial link was not reachable during this audit, so the complete spatial layer is not currently reproducible from the stated public link alone.

## Figure-to-evidence map

`anchor/figure-map.tsv` is only an eight-Figure index. It must not be used to assign one feasibility label to a multi-panel Figure. The CRC cancer-switch audit is resolved at the 92-main-panel level in `anchor/main-figure-panel-map.tsv` and the 63-supplementary-panel level in `anchor/supplementary-figure-panel-map.tsv`, with human-readable reviews in `reports/crc-subfigure-replacement-audit.md` and `reports/crc-supplementary-panel-audit.md`, and the code reconstruction plan in `reports/crc-code-reconstruction-map.md`. Main-figure logic:

- Figure 1: integrated atlas and tissue composition.
- Figure 2: malignant cNMF programs, meta-programs, pseudotime and cross-reference comparisons.
- Figure 3: regulators, trajectory, survival, sorafenib association and predicted drug sensitivity.
- Figure 4: myeloid/macrophage ecosystem.
- Figure 5: endothelial/fibroblast states and prognostic associations.
- Figure 6: cell-cell communication changes across primary and advanced sites.
- Figure 7: spatial programs, niches and neighborhood modeling.
- Figure 8: Geneformer virtual knockout, 11-gene intersection and HSP90B1 validation.

## Analysis-step disposition

- **Retain:** multi-role cohort architecture; recurrent malignant programs; explicit tumor–TME–space–clinical–perturbation evidence order; orthogonal target prioritization.
- **Repair:** immutable retained-sample manifest; donor-level composition and differential testing; cohort/stage confounding audit; within-dataset and leave-one-dataset-out program robustness; spatial replication across all samples; predeclared source/target states; independent perturbation controls; multivariable clinical validation.
- **Substitute:** use official cNMF/Geneformer/cell2location/COMMOT/MISTy implementations or reliable related-paper code where author code is unavailable; choose each target-cancer cohort for a stated paper task rather than by accession count.
- **Extend:** require negative controls and simple baselines for foundation-model perturbation; add tumor-versus-normal selectivity for candidate dependencies; use external or prospective treatment response only when outcomes are explicitly defined.
- **Drop unless evidence appears:** literal “initiation” and “dynamic evolution” wording from cross-sectional normal/primary/metastatic samples; causal language for inferred communication/spatial correlation; therapeutic-target language without experimental or strong perturbational validation.
- **Blocked for exact reproduction:** the authors' complete preprocessing choices, exact 115-sample selection, original scripts/environment, trained Geneformer artifacts and the currently unavailable old spatial download.

## Defect-to-repair contracts

The complete contracts are in `anchor/defect-repairs.tsv`. The most important are:

1. **Cross-sectional data are not temporal evolution or initiation.** A target paper must either obtain true premalignant/longitudinal/lineage evidence or use “state continuum/association” language.
2. **Stage is nearly confounded with dataset and tiny advanced groups.** In the supplied manifest PVTT has two samples and MLN one. The target must report per-group donors, require bridge cohorts and run dataset/donor holdouts before progression claims.
3. **Cell counts are not independent patients.** Cell fractions and marker differences must be tested at donor/sample level with compositional or mixed models.
4. **cNMF input provenance is underspecified.** Clipping negative Harmony residuals is not a sufficient reproducible contract. The target will use a documented nonnegative count-derived matrix, donor-aware program discovery and stability analysis.
5. **Spatial niche selection risks circularity.** All specimens or prespecified/held-out regions must be analyzed; co-enrichment cannot both select and validate the same ROI.
6. **Geneformer output is a hypothesis, not a knockout experiment.** Freeze model revision, tokenization, splits, states, seeds, nulls and baseline comparisons; require orthogonal selectivity/perturbation evidence before using “dependency” or “target”.
7. **Clinical associations are under-specified.** Use explicit arm/outcome manifests, multivariable models, proportional-hazards checks, multiple-testing control and external validation where available.

## Reproduction boundary

Exact end-to-end reproduction is not currently possible from public materials. The public record permits independent reconstruction of many modules, but not verification that the reconstruction matches the authors' implementation. Missing elements include the exact retained-sample manifest, author code, lockfiles, random seeds, trained model/checkpoints, train/validation/test splits, complete perturbation background/null construction, several spatial-processing choices and experimental validation of HSP90B1. These limits do not invalidate the framework; they define which modules must be rebuilt and which claims must be lowered in the target paper.

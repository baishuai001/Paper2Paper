# Data qualification

## Why a separate registry exists

A public accession can still be unusable because the required samples,
biological units, endpoints, batch fields, patient IDs, or processed matrices
are missing. `datasets.tsv` records dataset facts; `dataset_route_map.tsv`
records whether those facts satisfy one route's actual needs.

## Qualification sequence

1. Resolve the accession and canonical download URI.
2. Record access class: public, controlled, request-only, local, unavailable,
   or unknown.
3. Identify platform, disease, cohort size, biological unit, data level, and
   indispensable metadata.
4. Download metadata and inspect a representative data file.
5. Record file size and checksum when a stable file is obtained.
6. Build patient/sample/source overlap groups before assigning discovery and
   validation roles.
7. Map each route to required datasets, fields, and independence assumptions.

## Download states

- `not_checked`: only a citation or link is known;
- `metadata_verified`: accession and metadata were inspected;
- `sample_verified`: a representative file was downloaded and parsed;
- `downloaded`: the required route files were obtained;
- `checksum_verified`: required files were obtained and integrity recorded;
- `blocked`: access or file state prevents use.

`metadata_verified` alone is insufficient for route approval. Every required
dataset needs at least `sample_verified`, an accepted dataset state, and a
`verified` route-map verdict.

## Independence rules

Never infer independence from different accession names. Check paper, author,
institution, recruitment period, patient identifiers, sample descriptions,
and deposited object provenance. Record `independent`, `overlap`, `unknown`, or
`not_applicable` per route use.

Cells, spots, image patches, and repeated samples are not independent patients.
The biological unit in the analysis specification must match the unit that
supports the manuscript claim.

## Stop conditions

Stop or replace a dataset when an indispensable field is absent, controlled
access cannot be obtained within the project window, cohort identity cannot be
resolved, or the only validation data overlap the discovery cohort. Lowering a
claim is acceptable; fabricating metadata or treating cells as replicates is not.

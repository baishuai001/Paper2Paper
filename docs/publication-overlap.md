# Publication-overlap audit

## Purpose

Paper2Paper does not rank routes by perceived originality. It does require
evidence that the planned manuscript is not a substantive duplicate of an
existing paper. This is a narrower, auditable question.

## Required dimensions

Compare the nearest paper with the planned route across:

1. disease or population;
2. marker, gene set, cell type, pathway, or signature;
3. primary outcome;
4. discovery and validation datasets;
5. analysis skeleton;
6. main-figure sequence and narrative;
7. central manuscript claim.

Record the search query, source URI, date, dimension-level judgments, explicit
distinguishability, and proceed/stop decision in `publication_overlap.tsv`.

## Classification

- `clear`: no material match found after a documented search;
- `adjacent`: some dimensions overlap but the manuscript's central claim and
  evidence remain clearly different;
- `high_overlap_distinguishable`: many dimensions overlap; continuation is
  allowed only with a concrete distinction visible in title, aim, figures, or
  decisive evidence;
- `duplicate`: central claim and principal evidence are substantively the same.

Only `duplicate` is an automatic stop. “Same cancer” or “same marker” alone is
not duplicate. Conversely, renaming a score while reusing the same disease,
datasets, pipeline, figure story, and claim can still be duplicate.

## Timing

Run the audit twice:

- during G1 before investing in full implementation;
- immediately before manuscript release because the literature may have changed.

The second audit may force a bounded distinction, reroute, or stop even after
the analysis has run.

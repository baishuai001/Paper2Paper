# Route selection standard

## Purpose

Route selection answers a practical question: which adaptation of this anchor
paper is most likely to become a correct, reproducible, reviewable manuscript
with the data, code, time, compute, and skills actually available?

It does not ask candidates to prove a novelty threshold. A simple substitution
is legitimate when the scientific design remains valid, its execution path is
real, and the final central claim is not a substantive publication duplicate.

## Step 1: decompose the anchor

Describe the anchor at figure and module level:

- scientific question, population, biological unit, comparison, and outcome;
- each main and supplementary figure and its manuscript role;
- dataset, metadata, code, method, parameters, and output behind each figure;
- which modules can be retained, repaired, substituted, extended, or dropped;
- claims that require unavailable private data, wet-lab intervention, or
  undocumented author decisions.

The output is an anchor-module map, not a verdict on what the new paper must be.

## Step 2: generate a route portfolio

Generate every locally plausible route mode before ranking:

| Mode | Typical change | Default value |
| --- | --- | --- |
| faithful reproduction | same question and data | training, baseline, code qualification |
| marker substitution | replace one marker or target gene | fast template-preserving test |
| gene-set substitution | replace the central program or gene collection | common signature or pathway route |
| cell-type substitution | keep cancer and workflow, change focal cell | common single-cell adaptation |
| cancer-type substitution | keep object and workflow, change disease | common cross-disease adaptation |
| pan-cancer extension | expand one-cancer design to several cancers | broader public-data route |
| signature construction | construct and validate a transcriptomic score | common bulk-to-single-cell route |
| combined substitution | change multiple explicit axes | allowed when internally coherent |
| custom | bounded route not covered above | requires plain-language explanation |

Do not reject a route because it is “only” one of these transformations. Do
reject or repair it if the comparison, biological unit, outcome, validation,
or claim ceiling is invalid.

## Step 3: perform execution triage

For each route, fill independent fields rather than one total score:

- scientific validity;
- data readiness and missing fields;
- code readiness and donor grade;
- figure coverage;
- implementation and beginner burden;
- anchor reuse;
- publication overlap and distinguishability.

The machine derives priority:

- `P0`: pass + verified data + qualified code + complete figure map +
  non-duplicate overlap;
- `P1`: verified data and adaptable code, with a bounded spike remaining;
- `P2`: usable data but critical code must be rebuilt;
- `P3`: a hard dependency remains unverified, blocked, invalid, or duplicate.

P0-P3 is not a measure of scientific prestige. It predicts execution readiness.

## Step 4: run the minimum spike

Before approval:

1. download and parse a representative sample from every indispensable data
   source;
2. verify the fields that implement the central comparison and outcome;
3. install every required code donor in a pinned environment;
4. run at least one noninteractive end-to-end smoke test;
5. generate a small representative output and its source table;
6. update overlap search and document the distinguishable central claim.

Textual confidence cannot replace this spike.

## Step 5: human selection

Only P0 routes can be approved. When several P0 routes exist, the project owner
selects one after reviewing the explicit trade-offs. P1-P3 routes remain
backup, reconstruction, training, or stopped routes. AI recommendations stay
`proposed` until recorded in `reviews.tsv`.

## Route-change boundary

Changes compared before approval are candidate adaptations. After approval,
changing the central question, target disease/object, primary outcome, or
evidence chain is a `reroute`; changing a dataset, parameter, adapter, or
bounded analysis while answering the same question is normally `refine`.

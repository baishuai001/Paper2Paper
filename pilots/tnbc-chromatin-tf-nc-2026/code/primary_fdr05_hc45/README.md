# LUAD FDR05 -> HC45 primary branch

This directory implements the user-frozen replacement analysis:

- Figure 1 LUAD entry: BH FDR <= 0.05, NES > 0 and expression logFC > 0.
- Figure 2: the unchanged promoter/activity and motif principles applied to the
  resulting 187 candidates, yielding 45 HC-TFs before motif annotation.
- Figure 3B: an edge requires at least two shared activating HC-TFs.
- Figure 4: the primary display uses an unadjusted log-rank P <= 0.05 after a
  maximally selected cutpoint (10-year administrative truncation and at least
  20% of samples per group). Adjusted values remain audit columns.

The completed 158/31 run remains immutable. Every new run is written to a new
execution root and carries a machine-readable receipt.

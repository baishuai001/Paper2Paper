# Bulk fixed-signature contract

This directory defines the smallest reusable contract learned from the
`gastric-nrrs` GN-R01 reference exercise. The current release points to the
pilot implementation because no second independent paper has yet demonstrated
transfer. It is therefore `reference_verified`, not `transfer_verified`.

The contract covers deterministic calculation of a published fixed signature,
patient-level source tables, explicit preprocessing and mapping choices, and
fail-closed behavior. It does not claim that a signature is biologically valid,
independently validated, deployable for one new patient, or portable to another
assay or disease.

The frozen GN-R01 numerical results are a change detector for one exact resource
manifest and implementation. They are deliberately registered separately with
`oracle_source=reviewed_freeze` and cannot advance module maturity on their own.

The release reaches `reference_verified` only by combining the accepted
GN-R01 real-data run GN-RUN05 with the separately logged GN-RUN06 focused unit
and adversarial tests. GN-RUN05 is not treated as proof that those tests ran. The
identity-linkage negative test remains an explicit `not_run` registry gap and
therefore its proposed core promotion stays `observed`.

A future release may be marked `transfer_verified` only after a transfer case
uses a different pilot, anchor paper, and dataset while satisfying this contract.

The release registry pins LF-canonical SHA256 values for the implementation,
dependency lock, this contract, and the unit-test file. Semantic text drift
without a new release record fails strict registry validation, while CRLF/LF
checkout differences do not create a false stale-release failure. Individual
run metadata still retains the raw hash of the bytes actually executed.

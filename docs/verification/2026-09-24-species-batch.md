# Frozen composition validation in the scalar group solve

The scalar inner iteration assembles each independent species equation against
one unchanged composition tuple. Previously every assembly scanned every
independent species in the trial, accepted and previous layers. This repeats a
whole-tuple scan once per equation, giving quadratic work in the species count.

`SpeciesCompositionBatch` validates that tuple once after the iteration's halo
and physical boundary update. The row and retained-diagonal assemblers can use
this private, scoped validation result. They retain all material, source, flux,
alias, boundary and equation checks. The validation loops and summation order are
unchanged; no physical field, stopping criterion, source policy or model default
changes. No case configuration is needed.

The batch binds the plan, box and complete field-view identities, revisions and
layouts of all three history layers. A different binding falls back to full
validation. Outputs must be disjoint from every borrowed composition layer. A
failed preparation invalidates the previous result. The scalar owner invalidates
the batch before updating the iterate; a new batch is built on every iteration,
solve and retry. Statistical assembly and ordinary public assembly continue to
validate independently. This is not a global or cross-step cache.

## Validation

- The existing numerical species suite passes, including common-face transport,
  IBM species matrices, conservation and public-assembler atomic rejection.
- New production-assembler checks compare cached and uncached coefficients,
  check all three history revisions, storage identity, box extent, output aliases,
  and invalidation after a failed prepare. Focused launches pass with 1 and 2
  ranks; the analytical fixture itself uses `MPI_COMM_SELF`.
- Production reaction tests pass with 1 and 2 MPI ranks, including interval PaSR,
  missing interval-provider rejection, one integration per active cell/attempt,
  and injected partial-source rollback.
- A private compressible reacting case with 6,099,456 allocated cells, 128 MPI
  ranks and 47,652 cells per rank uses the same initial state and `dt=3e-7 s`.
  Both programs accept the first attempt after 23 outer iterations and 93
  pressure iterations. All 330,382,720 checkpoint values compare numerically
  equal: current/previous fields, accepted/previous rates, current/previous face
  fluxes and pressure/mass references. Metadata matches. Mass, energy and element
  budgets pass without changing their thresholds.

The runtime build uses the changed translation units over the V1.1.0 archive.
Tests link the same current common-face kernel overlay as the runtime. An initial
numerical-test link omitted that overlay and reproduced the archived kernel's
plateau-coefficient failures; adding the runtime overlay resolves those failures.
This is focused changed-module validation, not a fresh full-repository build.

## Performance scope

The measured host is an AMD EPYC 7763 system. Runs use 128 processes, core binding,
socket mapping and one OpenMP/BLAS thread per process. Builds, tests and simulations
are serialized by one work lock. Baseline source is `9e90d30`.

Two sequential fixed-state comparisons give:

| Pair | Baseline step (s) | Batch step (s) | Reduction |
| --- | ---: | ---: | ---: |
| Initial | 72.4357 | 69.9669 | 3.41% |
| Adjacent repeat | 71.8434 | 70.4013 | 2.01% |

These are two single-step observations, not long-run throughput statistics.
This is a small scalar-path improvement;
it does not reduce the 23 outer iterations or resolve overall compressible
coupling performance. Detailed timing and native comparison evidence remain with the private case.
There is no same-case executable COAST speed-parity claim and no long-test
acceptance claim. Private geometry, configuration and result files are not part
of this change.

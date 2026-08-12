# Task 11 science-closure mutation-sensitive RED design

Status: authorized design; no product repair has been authorized or performed
by this record. The current Task 11 candidate remains rejected.

Accepted base:
`0db56e463470dd1a605709ba05d8bd6a900f496b`

Takeover HEAD:
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8`

Controlling public records:

- `docs/superpowers/specs/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale-design.md`
- `docs/superpowers/specs/2026-08-05-hundun-flow-task11-force-consistency-authority-addendum.md`
- `docs/superpowers/plans/2026-08-05-hundun-flow-stage3-science-first-execution-amendment.md`

This design preserves the current tracked and untracked A22/A3/R1 chain. It
does not treat those candidate records, the current dirty tree, or any old
24/48/96 log as accepted capability.

## 1. Question to close

Task 11 cannot proceed by changing another numerical formula until executable
REDs distinguish these questions:

1. What sign does one product immersed wall-row pressure/viscous contribution
   have relative to the force exerted by the fluid on the solid?
2. Is `operator_reaction_N` a physical body-force report or the opposite-sign
   balance reaction used to close the summed fluid residual?
3. Does the consistency report compare two same-orientation physical forces or
   two opposite-orientation balance terms?
4. Do pressure and viscous consumers use one shared boundary-row authority, or
   can a per-link reconstruction silently diverge from the row actually used
   by the operator?
5. Is the observed first-order pressure/total consistency produced by exact
   spatial measure/reconstruction or by solved-state pressure/PISO error?

The REDs must answer the questions independently. A screen result alone is not
a sign oracle, and a sign oracle alone is not an order proof.

## 2. Mathematical definitions

Let

```text
n_s = unit normal from solid into fluid
n_f = outward unit normal of the active fluid domain at the immersed wall
n_f = -n_s
sigma = -pi I + tau
```

The force exerted by fluid on solid at one true-surface quadrature patch is

```text
F_surface = sigma n_s A_surface
F_surface_pressure = -pi n_s A_surface
F_surface_viscous  = tau n_s A_surface
```

For the product momentum residual convention

```text
R = transient + convective - div(sigma),
```

the wall term predicted by the divergence theorem for a matching physical
patch is

```text
R_wall = -sigma n_f A = sigma n_s A.
```

The opposite term used to close a global residual sum is

```text
B_wall = -R_wall.
```

For a Ghost-Cell/full-background-cell operator, one link's discrete operator
measure is not assumed equal to one true-triangle quadrature weight. Let
`F_operator` denote the independently evaluated shared-row functional using
its own frozen projected/background measure. The sign prediction is

```text
R_wall = F_operator                    (fluid-on-solid orientation)
B_wall = -F_operator                   (fluid balance orientation)
```

while `F_operator - F_surface` is required to converge at least second order
only after the approved global accumulation/refinement. This RED must not add
a cut-cell-like pointwise measure equality that Stage 3 does not implement.

Current code places `-evaluated.residual` in
`operator_reaction_N` and then adds that report to `surface_traction`. That is
algebraically coherent only if the report is `B_wall`, but it conflicts with
the frozen text that also treats the report as the same physical body force as
the surface result.

The derivation is not allowed to repair the product by assertion. RED-S1 must
observe the actual product row. If the observed raw product residual uses a
different convention, stop and revise this derivation before any product edit.

The intended separation to test is:

```text
operator physical force = independently accumulated R_wall/F_operator
surface physical force  = independently integrated sigma n_s A_surface
physical consistency    = operator physical force - surface physical force
budget reaction         = -R_wall
budget closure          = summed residual + budget reaction
```

No operator physical force may be copied from surface quadrature. No budget
oracle may redefine the public physical-force orientation.

## 2.1 RED-S0: signed-force non-degeneracy preflight

Before building RED-S1 or launching another formal sequence, evaluate only
candidate-independent analytic data for the existing asymmetric MMS on the
approved body/mapping sequences:

```text
analytic pressure, viscous and total force vectors
integral of componentwise absolute traction and traction L1
frozen pressure/viscous/total normalization scales
48/96 independent analytic quadrature agreement
sign-reversal and omitted-component mutation separation
```

Use the existing fixed `max_abs(vector)` force-error norm; do not select a
component after seeing results. For viscous force define

```text
T_abs = integral_surface ||tau n_s||_1 dA
D_net = max_abs(F_viscous_reference) / T_abs
```

The analytic reference records the three independent positive component
integrals `A_abs = integral_surface abs(tau n_s) dA` and defines
`T_abs = A_abs.x + A_abs.y + A_abs.z`. The scalar must equal that fixed sum;
an L2 norm, a signed component integral or a single-coordinate substitute is
not the same evidence.

The componentwise absolute value has kinks where a traction component changes
sign, so raw 48/96 Gauss agreement of `T_abs` is not used as a smooth-reference
certificate. Instead use the pointwise norm inequality and Cauchy--Schwarz:

```text
T_abs
  <= sqrt(surface_area * integral_surface ||tau n_s||_1^2 dA)
  <= sqrt(3 * surface_area * integral_surface ||tau n_s||_2^2 dA)
  = T_cert

T_cert = sqrt(3) * viscous_traction_rms_force
D_cert = max_abs(F_viscous_reference) / T_cert <= D_net
```

`viscous_traction_rms_force` is the existing smooth positive-weight analytic
quadrature quantity. Both its 48/96 agreement and the net force-vector 48/96
agreement are checked independently. This certificate avoids accepting a
poorly resolved absolute-value integral and avoids rejecting a valid fixture
merely because an absolute-value kink is not spectrally converged. `T_abs` and
the RMS input use the same positive-weight rule at each level, for which the
discrete Cauchy--Schwarz inequality also holds. No floating
`T_abs <= T_cert*(1+roundoff)` assertion with an unspecified tolerance is
introduced.

and require, before any solved candidate value is read,

```text
for both the 48- and 96-point references:
  each component of F_viscous_reference is finite
  each A_abs component >= 4096*eps*viscous_force_scale
  T_abs is the exact fixed component sum and is finite
  max_abs(F_viscous_reference) >= 4096*eps*viscous_force_scale
  max_abs(F_viscous_reference)
    >= 1024*sqrt(eps) * T_cert
48/96 independent net-force and T_cert reference differences
  <= 1e-13*max(1, viscous_force_scale)
```

The direct force-component finiteness check is mandatory. It must not be
inferred from `max_abs` or `max_abs_difference`: a non-finite component can
make a subtraction non-finite and a comparison helper based on ordered
`std::max` can otherwise hide that component. The component floors together
with the exact finite `T_abs` sum imply finite positive `A_abs` components, so
there is no separate redundant absolute-component finiteness predicate.

Because `D_net >= D_cert`, the certified `1024*sqrt(eps)` lower bound still
bounds the cancellation condition number before the formal error is measured;
it is a fixture-discrimination precondition, not a relaxed order threshold.
The implementation records `D_net` and `D_cert` only after their primitive
inputs have passed. The acceptance comparison is the division-free
`F_net >= 1024*sqrt(eps)*T_cert`, so no redundant ratio-finiteness predicate or
overflowing division participates in classification. Pressure and total
references retain their existing frozen scale/floor checks. All quantities and
mutations are evaluated without a product candidate.

The initial choice is the existing approved asymmetric MMS, not a tuned new
case. If it fails this analytic preflight, stop the formal matrix and create a
separate fixture addendum containing the exact new field coefficients, body,
component/vector norm, analytic references and mutation results. That addendum
must be reviewed before any product run. No candidate error or observed order
may select the replacement fixture.

The preflight extends the test-side `verified_force_reference()` contract.
Its signed-viscous classifier must be the same helper used by both the formal
reference path and a focused behavioral test; the test may not mirror the
classifier. The analytic accumulator and that focused test must likewise share
one pure componentwise-absolute increment primitive. The hand-derived input
`traction = {-2, 3, -4}`, `positive_weight = 0.25` must produce exactly
`{0.5, 0.75, 1.0}`; a signed or omitted-component implementation must fail.

The TDD transition is fixed. First perform a behavior-preserving extraction
under the existing green oracle, with default-zero new evidence fields, the
old classifier behavior and a temporary zero-returning increment primitive.
Then add the focused selector and observe an executable assertion RED from the
zero primitive/real zero L1 evidence and the old classifier accepting the
controlled negative references. A compile or link failure is not the required
RED. The selector must compute and print every named observation before one
final aggregate assertion; it must not lose later observations to the test
harness's first-failure exception. Full-wrapper exceptions are caught and
recorded as failed observations. Only after the complete RED vector is present
may the minimal absolute-component arithmetic, fixed L1 finalization and new
classifier predicates be implemented.

For each real full-wrapper reference the independently named structural
observation is

```text
real_l1_valid = finite(T_abs) && T_abs > 0
                && T_abs == (A_abs.x + A_abs.y) + A_abs.z.
```

The positivity term is required in the initial RED: default-zero scaffold
fields satisfy zero-equals-zero composition but are not valid L1 evidence. If
a wrapper throws, its wrapper, L1, omitted-error and reversed-error
observations default to false; the selector still prints all four.

Before GREEN, the focused test supplies candidate-independent copies of a
valid analytic reference and requires rejection after each of these independent
mutations:

```text
non-finite viscous net-force component at both levels
one zero absolute component with the other components and L1 positive
finite absolute components whose fixed L1 sum overflows
positive L1 scalar formed from L2 or one component instead of the fixed sum
net force below the roundoff floor while D_cert is above its floor
net force above the roundoff floor but below the D_cert floor
48/96 net-force mismatch
48/96 T_cert mismatch
```

It also supplies a valid reference whose fixed `max_abs` component is not `x`,
so a single-coordinate norm mutation fails. The low-`D_cert` fixture is chosen
analytically relative to the frozen threshold, not from product results, and
must distinguish omission of the `sqrt(3)` certificate factor. Each synthetic
negative keeps its absolute components in fixed positive proportions summing
to `T_abs` and above the component floors unless that specific component or
sum predicate is the mutation under test. Until the behavioral RED is
observed, the helper is minimally implemented, the negative fixtures are
rejected, the real full-preflight wrapper has exercised every unique approved
body, and actual one-at-a-time classifier/arithmetic mutants fail, no
signed-force screen or formal acceptance may start.

## 3. RED-S1: signed one-link authority

### Fixture

Add one focused MPI-capable test executable, provisionally named
`test_task11_signed_force_authority`. It must use product geometry,
classification, GhostStencilPlan, boundary-row assembly and wall quadrature.
The fixture must contain at least one selected active row with exactly one
selected wall link and an oblique `n_s`; a coordinate-aligned or closed-body
net-zero assertion is not sufficient.

The test inspects the selected link/row contribution before global closed-body
cancellation. It runs three exact cases:

1. pressure-only: nonzero pressure, zero viscosity/velocity gradient;
2. viscous-only: zero pressure, divergence-controlled linear velocity with
   nonzero normal and tangential traction;
3. combined: both contributions nonzero and not parallel.

The surface expected values are computed directly from
`sigma n_s A_surface`, without calling WallForceIntegrator. The operator
expected values are computed from independently queried topology, geometry,
exact stress and the frozen LFP signed-measure formula, without reading the
product's final shared-row descriptor, calling the product row evaluator or
copying its final snapshot. RED-S1 is intentionally limited to downstream
direction/report semantics; RED-S2 owns row-authority construction. The
fixture must prove that operator and surface measures are both nonzero; it
should keep them non-identical when practical so a surface-copy mutation
cannot pass.

### Required assertions

For pressure and viscous parts separately:

```text
raw product wall-row contribution == analytic F_operator
surface quadrature contribution    == analytic F_surface
reported physical operator force   == analytic F_operator
budget reaction                    == -analytic F_operator
physical consistency               == operator - surface
summed residual + budget reaction  == 0 for the isolated balance
```

The test does **not** assert `F_operator == F_surface` for one link. Their
global difference is an order oracle in RED-S3/S4, not a pointwise identity.

The first RED may use current `operator_reaction_N` as the declared physical
operator report. On the current candidate it is expected to fail the physical
orientation assertion; that failure is useful evidence that one report is
serving two incompatible meanings. A later semantic repair may separate or
rename the budget quantity, but only after the RED result is recorded and the
public science wording is amended.

### Mutations that must die

| ID | Actual product mutation site | Independent failing oracle |
|---|---|---|
| S1-M1 | `ImmersedOperatorAdapter::accumulate_momentum`, reverse the report/raw-row sign relation | topology/geometry analytic operator direction and raw row snapshot |
| S1-M2 | `WallForceIntegrator`, reverse `-pi n_s A_surface` | direct Cauchy pressure traction |
| S1-M3 | `WallForceIntegrator`, reverse `tau n_s A_surface` | direct analytic viscous stress |
| S1-M4 | `collect_final_force`, exchange addition/subtraction | expected physical consistency vector, not only its norm |
| S1-M5 | `collect_final_force`, copy surface result into operator result | independent operator-only and surface-only perturbations |
| S1-M6 | report assembly, swap/merge pressure and viscous components | pressure-only and viscous-only cases |

The surface-copy mutation is detected by two follow-up perturbations: one
changes an operator-row coefficient while keeping the surface reconstruction
fixed; the other changes a surface quadrature input while keeping the raw row
snapshot fixed. Exactly one independent result may change in each probe.

Analytic-oracle self-mutations such as reversing `n_s` are useful test-quality
checks but are not claimed as product mutants. Product mutants above are built
one at a time in a disposable exact-candidate worktree; the focused RED binary
must fail at the listed assertion. The main worktree remains unchanged and no
mutation hook enters a non-test build.

## 4. RED-S2: shared boundary-row authority

Primary host:
`tests/mpi/test_immersed_operator.cpp`. Product rows and
`last_boundary_row_evaluations` are observed values only. The expected donor
union is rebuilt independently from accepted GhostStencilPlan link donor IDs,
topology shared-face endpoints and geometry, then passed through the accepted
QuadraticReconstruction primitive. The expected oracle must not call
`ImmersedBoundaryAuthorityAccess::row_reconstruction`; that access is used only
to compare the product authority with the independent construction.

The fixture must contain a multi-link active row for which:

```text
union(row donors) != donors(any single link)
shared-row reconstruction result != legacy per-link reconstruction result
```

The test first proves this non-degeneracy. It then independently evaluates
pressure and viscous exact polynomial data from the independently constructed
row and requires the product authority/snapshot/report to match it.

Current pressure tests already distinguish the legacy per-link path. The
viscous oracle around `independent_viscous_hybrid_row` is insufficient because
it calls `ghost_plan.reconstruction(link.id)`, the same per-link authority as
the suspected product path. Replace or complement it with a truly independent
shared-row viscous oracle before claiming RED.

Mutations that must die:

| ID | Actual product mutation site | Independent failing oracle |
|---|---|---|
| S2-M1 | `build_boundary_row_reconstruction`, keep only first-link donors | independently sorted donor union and fingerprint |
| S2-M2 | pressure/viscous affine plan, use `ghost_plan.reconstruction(link.id)` | independent shared-row polynomial value/gradient |
| S2-M3 | complete A22 row, omit diagonal or neighbour defect | exact polynomial residual plus independent source-term inventory |
| S2-M4 | final boundary row, omit or repeat background removal | before/background/removed/wall/after decomposition identity |
| S2-M5 | viscous wall derivative, reuse stale/different authority | exact shared-row derivative with non-polynomial donor probe |
| S2-M6 | row execution, evaluate affine row twice | exact evaluation/write counters and final vector |
| S2-M7 | link ordering or rank decomposition affects assembly | canonical donor IDs/fingerprint and 1/2/4-rank bitwise result |

Each S2 mutant is compiled one at a time in the same disposable-worktree
protocol as S1. A copied product descriptor is never accepted as the expected
construction.

Downstream regression remains the existing shared-row predictor/history
coverage in `tests/mpi/test_immersed_piso.cpp`; it does not replace the direct
operator RED.

Proposed focused registrations:

```text
test_task11_shared_boundary_row_1_rank
test_task11_shared_boundary_row_2_rank
test_task11_shared_boundary_row_4_rank
labels: mpi;numerical;stage3;task11;fast
```

## 5. RED-S3: exact-state force-consistency decomposition

Primary hosts:

- `tests/support/stage3_mms.hpp`
- `tests/support/stage3_mms.cpp`
- `tests/numerical/test_laminar_ibm_order.cpp`

For each force part `q in {pressure, viscous, total}`, define at level `h`:

```text
O_h(U) = independently accumulated operator physical force
S_h(U) = true-surface force
C_h(U) = O_h(U) - S_h(U)
U*_h   = exact manufactured cell-average state on the product mesh
U_h    = solved product state

C_exact_h = C_h(U*_h)
C_solved_h = C_h(U_h)
C_state_h = [O_h(U_h) - O_h(U*_h)] - [S_h(U_h) - S_h(U*_h)]
```

The test first verifies, to a frozen roundoff-scaled tolerance,

```text
C_solved_h = C_exact_h + C_state_h.
```

This equality is a bookkeeping closure, not an independent RED: it follows
algebraically from the definitions. Root classification requires the separate
grid orders of `C_exact_h` and `C_solved_h` plus the frozen S1/S2/S3 product
mutations. With those observations, use this decision tree:

- `C_exact_h` first order: geometry/reconstruction/operator/surface authority
  is the root; do not patch PISO.
- `C_exact_h` second order but `C_solved_h` first order: pressure/PISO/state
  coupling is the root; do not tune surface quadrature.
- pressure passes but viscous fails: inspect shared viscous row/derivative;
  do not change pressure A22/A3.
- both exact and solved consistency pass: proceed to the complete formal gate.

Add selector `exact_force_consistency_fast` to the existing product numerical
executable. Proposed scope:

```text
12³: sphere-uniform pressure/viscous/total sign and decomposition
12³: warped-prism shared-row viscous discrimination
12³/24³: direction-only order screen using the same product path
```

Register the stable fast selector as
`test_laminar_ibm_exact_force_consistency_fast` with labels
`numerical;stage3;task11;fast`.

The exact-state path may inject approved exact cell averages and analytic wall
data through existing private test access. It must still use product mesh,
classification, reconstruction coefficients, immersed operator and wall-force
integration. It may not implement a replacement operator in test support.

Mutations that must die include operator sign, surface sign, wall-link
association, pressure measure, projected defect, shared-row-to-per-link
viscous reconstruction and removal of one pressure/viscous component.

## 6. RED-S4: solved fast and screen

Only after RED-S0 passes, RED-S1/S2/S3 fail for the intended reason and a
minimal repair makes them green:

```text
fast:
  12³/24³ sphere uniform
  12³/24³ warped prism
  direct two-corrector count and conservative-flux regressions

screen:
  24³/48³ sphere uniform
  24³/48³ warped prism
  24³/48³ translated sphere
```

All 15 formal categories remain hard. Their concrete signed-force fixture and
component are enabled only after RED-S0; screen then uses the frozen `>=1.8`
criterion to decide whether to launch 96³ and is not final acceptance.

If one consistency row fails, return to the exact-state decomposition. Do not
change threshold, case coefficients, filter, damping or corrector count.

## 7. RED-M: MPI, collective failure and rollback

These tests are independent of the force-sign decision and may be implemented
by a bounded worker after their exact fixtures are frozen:

1. Extract one internal production helper with inputs
   `(MpiContext, uint64_t local_count, uint64_t maximum, context)` that first
   classifies the count collectively and only then permits a checked cast.
   Its 2-rank test passes scalar counts `INT_MAX + 1` and `1`; it never
   allocates an `INT_MAX`-sized array. The initial test-only RED may be a
   compile failure for the missing helper; after minimal implementation it
   must return the same classified failure on every rank.
2. Build an underprovisioned Halo fixture whose execution reaches the product
   missing-donor guard rather than an earlier fixture validation.
3. Build an in-domain near-periodic donor row and require the product periodic
   authority to select it identically on 1/2/4 ranks.
4. Use one authoritative trial-state equality helper for committed/history,
   pressure, velocity, face flux, inactive canonical zero, PISO report and
   force report.
5. Inject failure in corrector 1, corrector 2 and final-force collection; each
   failed trial must roll back bitwise on 1/2/4 ranks.

Any local `size > INT_MAX` check that can strand peers is correctness work and
cannot be deferred as extreme-scale optimization. After GREEN, a disposable
mutant that restores a rank-local throw must fail under a bounded MPI timeout.
Both ActivePressureOperator and ActiveMomentumOperator must call the helper
before cast/allocation; codegraphf/caller review and a tests-off build must
prove that no synthetic-count override or test hook enters production.

## 8. RED-E: selector and evidence harness

Before a new long run:

- register every affected QR/Ghost/wall quadrature/reconstruction/operator/
  pressure/PISO/force/CG test with `stage3;task11` labels;
- gate exact-A22 and row snapshot diagnostics behind explicit selectors;
- prove selector-off runs do not retain large diagnostic vectors;
- install a detached runner that preserves the child exit code and records the
  manifest required by the 2026-08-05 execution amendment;
- test the runner with success, nonzero exit and signal termination using a
  small non-product command before entrusting it with a numerical job.

The currently running candidate2b sphere-uniform job is not runner acceptance
because its shell expanded the exit-code variable before execution.

## 9. Initial T11-S-only test-file proposal

Freeze the exact list again immediately before RED implementation. The maximum
T11-S-only proposal is:

```text
CMakeLists.txt
tests/mpi/test_task11_signed_force_authority.cpp        (new)
tests/mpi/test_immersed_operator.cpp
tests/mpi/test_immersed_piso.cpp                        (registration/reuse)
tests/numerical/test_laminar_ibm_order.cpp
tests/support/stage3_mms.hpp
tests/support/stage3_mms.cpp
finite_volume/src/immersed_operator_test_access.hpp     (only if unavoidable)
flow/src/stage3_flow_test_access.hpp                    (only if unavoidable)
```

No implementation file belongs in the RED commit/diff. If a test cannot be
made mutation-sensitive without a new private test seam, the main agent must
review that seam separately and prove it is excluded from non-test builds.
RED-M and RED-E are separate acceptance clusters and require separate frozen
allowed-file packets; this list does not authorize their transaction, Halo,
periodic, state-equality, CMake-label or runner hosts.

## 10. RED acceptance packet

Before product repair, record:

```text
exact candidate identity and test-only diff SHA-256
RED command, ranks, binary SHA-256 and log SHA-256
the exact failed assertion and why it is the intended failure
each independent mutation and the oracle that killed it
codegraphf changed-symbol and affected-test inventory
proof that no product implementation changed
main-agent requirements and test-quality review
```

Only then freeze the minimal product allowed-file list. RED evidence does not
authorize a broad cleanup or reopen A22/A3/R1 beyond the proven root.

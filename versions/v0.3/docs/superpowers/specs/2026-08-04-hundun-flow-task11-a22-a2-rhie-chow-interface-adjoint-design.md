# HUNDUN-FLOW Task 11 A22-A2 Rhie--Chow Interface-Adjoint Design

Status: approved by the user on 2026-08-04. This is the additive scientific
authority for the A2 repair only. It does not reopen the rejected full
normal-equation Route A and does not authorize G1/G2/G3, Task 12 or Stage 4.

Accepted Task 10 base: `0db56e463470dd1a605709ba05d8bd6a900f496b`

Candidate HEAD at A2 freeze:
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8`

Rejected predecessor:
`docs/superpowers/specs/2026-08-04-hundun-flow-task11-a22-adjoint-pressure-coupling-design.md`

Stopping evidence:
`.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-route-a-adjoint-normal-equation-stopping-report.md`

## 1. Root cause and selected correction

The rejected Route A formed the pressure operator from the complete centered
mechanical-pressure derivative. On a periodic collocated grid that derivative
annihilates both the constant mode and an odd/even checkerboard mode. The
result is a structural nullspace, not a tolerance, solver, MPI or coefficient
bug. The first integrated pressure solve consequently reached the frozen
500-iteration limit before corrector one.

A2 retains the accepted compact, time-consistent Rhie--Chow operator `K_RC`
as the volume-wide checkerboard authority. It adds only a positive
semidefinite interface correction derived from the difference between the
accepted A22 mechanical boundary row and its already-declared background row.
It changes no threshold, case parameter, corrector count, momentum equation,
wall flux, public API or Stage 1/2 contract.

The alternatives are closed as follows:

- full `C_A22^T M C_A22`: rejected by the recorded spectral proof;
- `K_RC` plus the complete A22 row: rejected because it repeats ordinary
  shared-face content and has no unique no-double-count interpretation;
- an augmented coupled projection: explicitly deferred unless A2 fails its
  algebraic or numerical gates and the user approves a larger amendment.

## 2. Unique background and interface-difference authority

For an owned active interface cell and one momentum component, let:

- `C_shared` contain the ordinary active--active shared-face pressure terms
  and any unchanged body-fitted physical-boundary term;
- `C_wall_bg` be the donor derivative in
  `BoundaryRowPlan::background_pressure_constrained`;
- `C_wall_A22` be the donor derivative in
  `BoundaryRowPlan::pressure_constrained`;
- `B_bg` and `B_A22` be their mechanical wall-normal-gradient derivatives.

The two complete mechanical rows are

```text
C_background = C_shared + C_wall_bg
C_A22        = C_shared + C_wall_A22
```

and the only permitted A2 difference is

```text
E = C_A22 - C_background = C_wall_A22 - C_wall_bg
B = B_A22 - B_bg.
```

The current declared background is constructed from the unconstrained wall
value functional, so `B_bg` is exactly zero and `B=B_A22`. The implementation
nevertheless computes the difference from both immutable plans so a future
change cannot silently add the background affine term twice.

`pressure_unconstrained` is not `C_background`. It is the older unconstrained
LFP replacement alternative. Using it in `E` would compare two replacement
schemes rather than compare A22 with the declared physical background and is
forbidden.

Only rows with one or more immersed links may contribute to `E`. Ordinary
active rows export no A2 term. Shared active faces and unchanged physical
patch terms are never inserted into the difference plan. Tests must prove
these exclusions by identity and by mutation, not by source-text matching.

For every component, both complete rows reproduce constant pressure, hence
`E 1=0` within the frozen scaled FP64 construction bound. The immutable plan
may canonicalize only one stable anchor coefficient when the raw closure
defect is at most

```text
64*epsilon*max(1, sum(abs(E_row))).
```

The anchor is the momentum cell when present, otherwise the lowest donor
global ID. A larger defect is a collective, non-retryable construction error.
The raw difference, pre-closure defect, anchor and canonical value are part of
the tests-only evidence. The mechanical A22 row itself is not changed.

## 3. Hybrid pressure algebra

Let:

- `p` be pressure correction in Pa;
- `q=S p`, with `S=diag(sqrt(V_i))`;
- `M=diag(1/a_P)` use the exact current assembled momentum diagonal for each
  owned active cell and component;
- `D` map one canonical owner-oriented active-face mass flux to owned-cell net
  outward mass flux;
- `H_E` satisfy the interface-only discrete adjoint identity
  `E=-H_E^T D^T`;
- `delta_g` be the corrected-minus-current mechanical wall-normal-gradient
  increment for the same pressure corrector.

Define

```text
E_s = E S^-1
L_E = rho E_s^T M E_s
s_g = rho S^-1 E^T M B delta_g
L_A2 = K_RC + L_E.
```

The normalized pressure equation is

```text
L_A2 q = -S^-1 D mdot_star - s_g.
```

`K_RC` is not claimed to equal `C_background^T M C_background`. It remains
the accepted time-consistent collocated projection operator. `L_E` is an
interface-consistency energy term only; it contains no cross term and no
ordinary shared-face coefficient. This distinction is part of the capability
claim and prevents overclaiming the hybrid as a complete mechanical normal
equation.

Because `K_RC` is symmetric positive semidefinite, `L_E` is symmetric
positive semidefinite and `E 1=0`, their sum preserves the constant mode. The
existing `K_RC` checkerboard energy remains present everywhere. In a closed
connected domain the tests must prove that only the constant mode survives;
with an active pressure outlet the existing `K_RC` reference remains the sole
ground and A2 adds no duplicate reference.

The exact Jacobi diagonal is

```text
diag(L_A2)_j = diag(K_RC)_j
  + rho/V_j * sum_d M_d E[d,j]^2.
```

The compact and interface terms use the same dependency revision. A stale
base operator, stale interface mobility, nonpositive diagonal or mismatched
layout is rejected collectively before a solve.

## 4. Conservative face correction and affine sign

After solving, form

```text
w = E_s q + B delta_g
delta_m_E = -rho H_E M w.
```

Then

```text
S^-1 D delta_m_E = L_E q + s_g.
```

The final correction is the sum of:

1. the existing reciprocal time-consistent Rhie--Chow face correction, whose
   divergence is `K_RC q`; and
2. the A2 canonical interface correction `delta_m_E`.

Therefore the corrected shared flux closes exactly the equation solved by
`L_A2`. The A2 contribution updates only canonical active--active faces. Each
stored face mass-flux increment also updates the face-velocity normal
component by

```text
delta_u_face = delta_m_E * S_face / (rho*|S_face|^2),
```

so `rho*dot(u_face,S_face)` and the shared mass-flux field remain identical.
Tangential predictor velocity is untouched. Immersed wall velocity and mass
flux remain exact positive zero. Periodic pairs receive reciprocal values
through one stable canonical pair record; no physical patch rule is bypassed.

The full cell correction remains the accepted equation

```text
A_u delta_u + C_A22 delta_p + B_A22 delta_g = 0
```

solved by the three existing momentum-correction solves. A2 does not replace
that equation with diagonal mobility and does not overwrite its result.

## 5. Immutable private plan and execution

One private factory consumes the already-built A22 and background affine
plans. It does not reconstruct a quadratic stencil. It produces stable,
sorted records for:

```text
E:   (momentum global cell, component, pressure donor global cell, coefficient)
B:   (momentum global cell, component, immersed link, coefficient)
H_E: (canonical active face/pair, momentum global cell, component, coefficient)
```

The plan stores raw/canonical closure evidence, row fingerprints, the active
layout, exact momentum mobility, reusable forward/transpose/face workspaces
and all MPI counts/offsets required by the existing active layout. It is
private to Stage 3. No public header gains a callable API; the existing
private friend declaration remains the only permitted public-header change.

The minimal implementation may reuse the already-existing active global-ID
exchange layout and preallocated buffers. Matrix apply and face update may not
allocate, rebuild a row, rebuild a route or create `FieldStorage`. Every
additional collective is explicit and covered by 1/2/4-rank deterministic
tests. A later peer-only optimization may replace the private communication
mechanism without changing the algebra, but it is not allowed to create a
second coefficient authority.

For each interface momentum DOF, `H_E` is constructed on the canonical active
face graph. Donor values route to the stable anchor over the shortest active
path; ties are broken by canonical global face ID. Missing paths, solid-face
routing, inconsistent periodic pairs, duplicate ownership or a non-closable
row are collective, non-retryable construction failures.

## 6. Transaction and failure behavior

The following are collective non-retryable input/construction failures:

- duplicate or nonfinite `E/B/H_E` identity;
- raw closure defect above the frozen bound;
- missing donor, wall link, active path, periodic reciprocal or owner;
- nonpositive/nonfinite `rho`, `a_P`, volume or hybrid diagonal;
- stale structure/numeric revision;
- MPI operation failure.

Linear non-convergence, numerical breakdown, nonfinite trial state and final
physical residual failure keep the frozen recoverable classification. A
failed attempt changes no committed/history/trial value, wall datum, face
field, pressure revision, controller state or fingerprint. Tests-only probes
are read-only and change no solver counter.

## 7. TDD and proof gates

Product integration is forbidden until independent tests-only oracles prove:

1. the exported raw rows equal `C_wall_A22-C_wall_bg` and
   `B_A22-B_bg`, including ordinary-row emptiness;
2. using `pressure_unconstrained`, repeating `C_shared`, omitting the
   background term or changing the stable anchor is detected;
3. `E 1=0`, `E=-H_E^T D^T`, symmetry, nonnegative energy and the exact
   additive diagonal;
4. `K_RC` retains nonzero checkerboard energy and the hybrid has the frozen
   constant/reference behavior;
5. homogeneous and affine operator/RHS/face-divergence identities agree;
6. wall flux remains positive zero and face velocity/mass flux agree;
7. the full mechanical cell correction identity remains unchanged;
8. repeated probes are deterministic and state-neutral on 1/2/4 ranks.

The direct pressure/PISO test must first fail for the missing A2 interface
term and then pass after the minimal implementation. Only after all algebraic
and direct Debug gates are green may the coordinator run 12/24 localization,
24/48 screening and the one final 24/48/96 Release matrix. A failed A2 proof
or screening stops this hypothesis; it does not authorize damping, filtering,
threshold relaxation, case tuning or the augmented projection.

## 8. Scope and capability statement

A2 claims only: a compact time-consistent Rhie--Chow pressure projection plus
one conservative A22 interface-difference correction, with the complete A22
row retained for cell momentum correction. It does not claim that the full
mechanical gradient is itself a checkerboard-safe pressure operator.

No public schema, Restart, diagnostic stable ID, plugin ABI, output format or
Stage 1/2 behavior changes. No new dependency, Python runtime, GPU
implementation, WALE, moving body, cut-cell method, thermal wall, chemistry,
private legacy source access, publication or push is authorized.

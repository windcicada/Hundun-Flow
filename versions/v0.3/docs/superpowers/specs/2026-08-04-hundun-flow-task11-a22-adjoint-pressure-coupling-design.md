# HUNDUN-FLOW Task 11 A22 Adjoint Pressure Coupling Design

Status: approved Route A hypothesis, rejected at its pre-grid spectral proof
gate on 2026-08-04. It is retained as an exact record of the tested
hypothesis, not as an implementation authority. The stopping evidence is in
`task-11-a22-route-a-adjoint-normal-equation-stopping-report.md`.

Accepted Task 10 base: `0db56e463470dd1a605709ba05d8bd6a900f496b`

Candidate HEAD at design freeze:
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8`

Superseded Task 11 pressure component: the compact active-cell pressure
operator and its independently interpolated face correction. The accepted A22
complete six-face quadratic mechanical-pressure row remains authoritative.

## 1. Problem and fixed boundaries

The A22 mechanical row already satisfies its polynomial, affine and multi-link
contracts, but the pressure solve still applies an unrelated compact
face-difference operator. The resulting pressure-operator/flux mismatch is
observable in the direct R4 RED and the formal pressure Linf sequence.

Route A replaces only that split pressure-correction coupling. It keeps:

- one `hundun` executable and all Stage 1/2 public contracts;
- the existing C++17/MPI-3 backend-neutral linear interfaces;
- the existing A22 row-level LFP reconstruction and its mechanical wall datum;
- the existing SPD momentum predictor operator;
- exactly two PISO correctors;
- the existing full cell correction equation
  `A_u delta_u + C delta_p + B delta_g = 0`;
- positive-zero immersed wall mass flux, rollback, gauge and Restart behavior;
- all frozen Task 11 numerical thresholds and formal sequences.

It adds no public method, public data type, dependency, filter, damping, case
parameter, third corrector, post-solve overwrite or alternative Stage 3
capability. One private friend-access forward declaration may be added to the
PIMPL class; it changes neither object layout, exported symbol nor callable
public API and must remain unavailable to public consumers.

## 2. Algebra and signs

Let:

- `p` be the owned active-cell mechanical pressure in Pa;
- `q = S p`, where `S = diag(sqrt(V_i))`;
- `C` map active pressure values to the three integrated cell momentum
  pressure residuals in N;
- `B` map immersed mechanical wall-normal gradients to the same residuals;
- `M = diag(1/a_P)` be the positive pressure-correction mobility, using the
  exact assembled momentum diagonal for the current step and component;
- `D` be the signed owner/neighbour incidence that maps one canonical shared
  face mass flux to owned-cell net outward mass flux;
- `H` map cell velocity changes to canonical shared-face *volume* flux and
  satisfy the exact discrete adjoint factorization

```text
C = -H^T D^T .
```

For a wall-gradient increment `delta_g`, define

```text
C_s = C S^-1
L   = rho C_s^T M C_s
s_g = rho S^-1 C^T M B delta_g.
```

The normalized pressure equation is

```text
L q = -S^-1 D mdot_star - s_g.
```

After the solve,

```text
w         = C_s q + B delta_g
delta_u_R = -M w
delta_m   = rho H delta_u_R
```

and therefore

```text
S^-1 D delta_m = L q + s_g.
```

The corrected face flux closes the same equation used by the pressure solve.
The separately retained cell correction solves

```text
A_u delta_u = -(C_s q + B delta_g)
```

with the accepted momentum solver. Thus `M` is the standard positive
Rhie--Chow/PISO mobility; it is not substituted for `A_u` in the cell momentum
identity. The rejected A22 experiment used a non-adjoint face interpolation;
it did not construct `H` from the accepted `C` coefficients and is not this
design.

For a closed connected active domain, `C_s(S 1)=0`, so `L` has exactly the
constant-pressure nullspace. With a prescribed-pressure outlet, the boundary
ground route removes that nullspace and no duplicate mean constraint is added.

## 3. One immutable coefficient authority

One private factory builds a `PressureCouplingPlan` from the same immutable
row construction used by `ImmersedOperatorAdapter`. It contains sorted stable
records for:

```text
C: (momentum global cell, component, pressure donor global cell, coefficient)
B: (momentum global cell, component, immersed link, coefficient)
H: (canonical global face, momentum global cell, component, coefficient)
```

The factory must not independently reimplement the quadratic reconstruction.
The complete six-face constrained row coefficients are produced once by the
existing A22 builder and consumed by both mechanical residual evaluation and
pressure coupling. Ordinary active rows add the existing shared-face pressure
terms; interface rows add the accepted constrained replacement terms.
Existing physical zero-normal-gradient pressure faces add their
owner-pressure area coefficient. Fixed physical pressure values are affine
constants and do not enter `C`.

The plan also owns:

- active layout and owned/ghost coefficient offsets;
- canonical face ownership and owner-oriented signs;
- precomputed forward and transpose MPI peer schedules;
- a structure fingerprint independent of rank decomposition;
- a numeric revision dependent on density, exact assembled `a_P`, geometry,
  active layout and wall-condition fingerprint;
- reusable `C`, `C^T`, `H` and face-update workspaces.

No matrix apply may rebuild a `FieldStorage`, row reconstruction, MPI schedule
or coefficient map, and no apply may allocate.

## 4. Conservative construction of H

For each momentum degree of freedom `d=(cell,component)`, the corresponding
row of `C` supplies the target pressure-cell divergence

```text
t_j = -C[d,j].
```

The factory constructs a deterministic local graph flow `h_d` such that

```text
D h_d = t
```

and stores `h_d` as the column of `H`.

The graph contains active-cell adjacency through canonical shared faces. For
a closed row, exact constant reproduction requires `sum_j t_j=0`. Coefficients
are accumulated in stable global-ID order; a bounded roundoff closure is
applied to the authority's stable anchor coefficient only when the pre-closure
defect is at most `64*epsilon*max(1,sum(abs(C[d,j])))`. A larger defect rejects
the plan. Both mechanical evaluation and pressure coupling consume the closed
coefficient, so this is algebraic canonicalization, not a pressure-only patch.

The stable anchor is the momentum cell when present, otherwise the lowest
donor global ID. Each other donor routes to the anchor over a shortest active
path; ties are broken by global face ID. A row with nonzero constant response
is legal only when it has an explicit prescribed-pressure boundary route; its
net flow terminates at that boundary ground face. Missing paths, solid-face
routing, duplicated canonical faces or an unexplained nonzero sum are
collective non-retryable construction errors.

Every internal face has one owner-oriented stored value. Its neighbour sees
the bitwise negated incidence contribution. Immersed wall faces never appear
in `H`, so their corrected flux remains positive zero. The transpose peer
schedule returns remote donor contributions to the owning active pressure DOF
without a global replicated sparse matrix.

## 5. Pressure operator and face update

The replacement remains a private implementation of `linear::LinearOperator`.
For every apply it performs:

```text
owned q
-> exchange required pressure donors
-> z = C S^-1 q
-> z = M z
-> y = rho S^-1 C^T z
-> collective finite check
```

Its exact Jacobi diagonal is derived from the same sparse coefficients:

```text
diag(L)_j = rho / V_j * sum_d M_d C[d,j]^2.
```

Remote contributions use the precomputed transpose schedule. The diagonal is
positive for every constrained owned row, finite, and revision-bound. It is
never borrowed from the retired compact operator.

The affine wall contribution is evaluated once per corrector as
`s_g = rho S^-1 C^T M B delta_g` and subtracted from the predictor-divergence
RHS. After the pressure solve, exactly the same `C`, `B`, `M` and `H` update the
canonical shared face mass flux. The face velocity receives the unique normal
increment

```text
delta_u_face = delta_m * S_face / (rho * |S_face|^2),
```

so `rho*dot(u_face,S_face)` and the stored shared face mass flux remain
identical. Tangential predictor velocity is not overwritten. Physical inlet,
symmetry and pressure-outlet rules remain owned by the existing boundary
registry. Immersed wall face velocity and flux remain exact positive zero.

## 6. Failure and transaction contracts

The following are collective, deterministic and non-retryable before a trial
is committed:

- invalid/duplicate coefficient identity;
- unexplained constant-mode defect;
- disconnected conservative route;
- missing or duplicated peer ownership;
- nonpositive/nonfinite `M` or operator diagonal;
- stale structure/numeric revision;
- MPI operation failure.

Linear non-convergence and nonfinite trial results retain the frozen
recoverable classifications. A failed attempt changes no committed/history
pressure, wall datum, face flux, controller state, revision or fingerprint.
Read-only test probes and diagnostics change no state or counter.

## 7. Proof gates before numerical grids

Product integration is forbidden until a tests-only sparse oracle proves:

1. `C=-H^T D^T` for every stored coefficient and canonical face;
2. symmetry `x^T L y = y^T L x` within the scaled FP64 bound;
3. nonnegative `x^T L x`, with a mutation that makes the oracle fail;
4. closed-domain constant nullspace and outlet reference behavior;
5. `G(p+dp,g+dg)-G(p,g)=C dp+B dg`;
6. exact operator/RHS/face-divergence identity for homogeneous and affine wall
   data;
7. exact diagonal reconstruction from `C` and `M`;
8. stable fingerprints and 1/2/4-rank decomposition invariance;
9. positive-zero wall flux and unchanged state under repeated probes.

Any algebraic failure rejects Route A before a 12/24 localization run. Route
B (an explicitly augmented coupled projection) requires a new user-approved
scientific amendment; it may not be entered silently.

## 8. Numerical closure order

After all algebraic gates are green:

1. Debug direct operator/PISO/force/transaction tests on 1/2/4 ranks;
2. Release tests-only 12/24 localization;
3. Release 24/48 screening for every affected frozen sequence;
4. one final Release 24/48/96 nine-sequence matrix, both adjacent orders
   `>=1.8`;
5. approved decomposition and long engineering gates;
6. Task 11 full requirements, code-quality and main-agent exact-HEAD closure.

No 96-grid task may run before steps 1--3 are green, and only one large
numerical task may run at a time with at most 96 logical CPUs.

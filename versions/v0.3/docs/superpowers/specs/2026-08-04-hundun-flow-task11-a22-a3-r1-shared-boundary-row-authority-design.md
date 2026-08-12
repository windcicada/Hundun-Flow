# HUNDUN-FLOW Task 11 A22-A3/R1 Shared Boundary Row Authority

Status: coordinator-approved minimal scientific repair under the user's
delegated Task 11 authority. It is additive to A22 and A3; Task 11 is not
accepted by this design.

## 1. Problem statement

For one active cell `P`, all immersed links `L(P)` must consume one immutable,
row-owned quadratic cell-average reconstruction. The present candidate uses a
joint reconstruction only for the A22 mechanical pressure row, while the wall
predictor, pressure-history reconstruction and wall-force path evaluate
separate link-owned fits. Near-coincident parallel links therefore supply
incompatible exact Neumann constraints and create an unbounded curvature
response.

## 2. Frozen mathematical authority

For the sorted row donor union `D(P)`, build one scaled quadratic basis and
one deterministic factorization

```text
q_P(x) = b(x)^T a,
a_0 = argmin_a ||A_D a - q_D||_2.
```

For wall-normal data `g_l`, form the joint constrained solution

```text
a_g = argmin_a ||A_D a - q_D||_2
      subject to c_l^T a = g_l for every independent link l in L(P).
```

The constraint rows are not dropped or regularized by a geometry threshold.
Exact rank deficiency is a layout failure. Close but independent constraints
remain legal because their data are generated from the same `q_P`, making the
divided difference consistent and bounded for smooth fields.

The same authority supplies:

- pressure-free predictor values at every row wall intercept;
- BE/BDF2 history wall values;
- unconstrained and constrained mechanical-pressure wall gradients;
- the A22 complete cell pressure functional;
- correction-pressure affine increments;
- pressure and velocity wall-point values/gradients used by force integration.

There is one donor order, one basis frame, one factorization identity and one
row/link association fingerprint. Per-link geometry remains distinct; shared
authority does not mean averaging wall points, normals or fluxes.

## 3. Conservation and boundary contracts

Each wall face retains its own predictor mass flux and the frozen identity

```text
delta_g_l = mdot_l_star / D_l,
mdot_l_star - D_l delta_g_l = +0.
```

The joint polynomial consumes those link-local data simultaneously. Active
shared-face flux remains reciprocal and is evaluated once. A3 continues to
apply the exact maps

```text
A delta_u = -(C delta_p + B delta_g),
delta_m = R delta_u + Q delta_p + q_g delta_g.
```

No compact substitute, filter, damping, extra corrector or solved-field
overwrite is introduced.

## 4. Minimal private architecture

Factor the existing A22 row donor-union, frame, unconstrained QR, constraint
nullspace and particular map into one private immutable
`BoundaryRowAuthority` (final name may vary). It is constructed once from
topology, geometry, domain and Ghost plan and keyed by stable active-cell ID.

`ImmersedOperatorAdapter`, the Stage 3 predictor and wall-force adapter consume
the same object through private detail interfaces. Public headers, schema,
Restart, diagnostic IDs and Stage 1/2 behavior do not change. Hot kernels do
not gain virtual calls, shared ownership or runtime donor selection.

## 5. Direct RED/GREEN contracts

1. A two-link row with close parallel normals and a non-polynomial donor vector
   must show that independent link fits produce a divided-difference response
   different from the shared row authority.
2. One shared row authority exactly reproduces constant, linear and all ten
   quadratic cell-average modes at both wall points and in the complete A22
   pressure functional.
3. Smooth non-polynomial data under geometrically similar refinement retain a
   bounded dimensionless wall-functional response; halving link separation
   cannot double an independent-fit mismatch.
4. Link permutation, donor permutation and 1/2/4-rank ownership produce the
   same stable IDs and approved FP64 values.
5. Missing/duplicate link constraints, wrong normal, wrong row association and
   mixing one old link reconstruction fail mutation-sensitive oracles.
6. A3 momentum, routed face divergence, bitwise positive-zero wall flux,
   rollback and exactly two correctors remain unchanged.
7. Release sphere `12^3/24^3` must pass both velocity L2 and Linf order `>=1.8`
   before any 24/48 screen.

## 6. Stop conditions and scope

Stop before product integration if the shared donor union cannot reproduce
quadratics within the approved bound, if the joint constraints are rank
deficient within the approved Halo reach, or if shared-face conservation is
lost. Do not compensate with a threshold, tuned grouping, high-pass,
regularization, diagonal damping, a third corrector or post-processing.

Only Task 11 private immersed/finite-volume/flow force adapters, tests-on
access, focused tests and Task 11 evidence records may change. G1/G2/G3,
Task 12, Stage 4, private legacy code, publication and push remain forbidden.

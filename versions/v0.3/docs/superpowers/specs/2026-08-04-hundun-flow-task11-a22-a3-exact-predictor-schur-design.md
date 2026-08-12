# HUNDUN-FLOW Task 11 A22-A3 Exact-Predictor Schur Design

Status: coordinator-approved within the user's delegated Task 11 repair
authority on 2026-08-04, after the A2-P direct rejection. This additive
scientific amendment does not accept Task 11 or authorize later tasks.

Accepted Task 10 base: `0db56e463470dd1a605709ba05d8bd6a900f496b`

Candidate repository HEAD at freeze:
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8`

Root evidence:
`.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-a2p-predictor-consistency-result.md`

## 1. Selected mathematical authority

A3 replaces the rejected A2 product experiment with the exact matrix-free
Schur complement of one discrete A22 momentum/predictor linearization. It does
not add a pressure filter or a second pressure authority.

For one corrector, write the already assembled momentum equation and face
predictor as

```text
A u = h - C p - B g
m = R u + Q p + q_g g + m_history.
```

Here:

- `A` is each actual implicit momentum operator, not its diagonal;
- `C/B` are the complete A22 mechanical pressure and wall-gradient rows;
- `R` is the exact active-face interpolation used by
  `assemble_face_predictor()`;
- `Q` is the exact explicit Rhie--Chow pressure-defect map, including the A22
  reconstructed cell gradients used by that same function;
- immersed wall mass flux remains exact positive zero;
- all maps use the current coefficients, geometry, active layout and time
  history revision.

For a pressure increment `delta p` and known wall increment `delta g`, solve

```text
A delta_u = -(C delta_p + B delta_g)
delta_m = R delta_u + Q delta_p + q_g delta_g.
```

With `q=S delta_p`, define

```text
L_A3 q = S^-1 D delta_m(q, delta_g=0)
s_g    = S^-1 D delta_m(delta_p=0, delta_g)
L_A3 q = -S^-1 D m_star - s_g.
```

The solved `delta_u` and the exact same `delta_m` are the only accepted cell
and face corrections. No separate three-solve post-processing and no A2 face
route may overwrite them.

## 2. Why this is the minimal general repair

A2-P proved that neither the compact operator nor `K_RC+L_E` is the exact
Jacobian of the reassembled predictor. A3 removes that approximation boundary
without changing the A22 reconstruction, momentum discretization, boundary
condition, corrector count or numerical threshold. It is algebraically the
velocity-eliminated form of the coupled momentum/continuity block system; a
future monolithic block solver may implement the same map without changing
the contract.

Nested momentum solves are permitted in the CPU-reference implementation only
as the initial correctness backend. They must use preallocated workspaces,
the existing backend-neutral linear interfaces and a frozen inner residual
contract. Inner failure is a collective recoverable linear failure. It is not
legal to hide failure by loosening the outer residual or returning an
approximate success.

## 3. Single predictor-Jacobian authority

Before product integration, one private linearization helper must be factored
from `assemble_face_predictor()` so that product predictor assembly, A3 apply,
A3 affine source and A3 face update use the same formulas and stable face/link
identity. Copying the formula into a second implementation is forbidden.

The helper consumes cell velocity increment, pressure increment, A22 pressure
gradient increment, wall-gradient increment, current face mobility and the
same owner-oriented geometry. It produces one canonical face-velocity and
mass-flux increment. Physical boundary behavior remains the frozen Stage 2/3
rule; immersed wall values are always bitwise positive zero.

The helper is private. No public header, schema, Restart field, diagnostic ID
or plugin ABI is added.

## 4. Direct proof before integration

The first A3 RED is a tests-only Jacobian identity. For deterministic pressure
directions on uniform and warped fixtures it compares:

1. the derivative obtained by applying complete A22 pressure residual,
   solving the three actual momentum correction equations and passing the
   increments through the factored predictor map; and
2. the desired `L_A3` action and routed face divergence.

It proves elementwise action, conservative face divergence, cell momentum
closure, wall positive zero, face velocity/mass consistency, superposition,
sign and state neutrality. The current A2 action must fail this oracle by a
non-negligible amount. A test that merely compares two copies of the new
helper is invalid.

Constant pressure remains the closed null mode. The explicit Rhie--Chow defect
must give nonzero odd/even energy. Bilinear symmetry is required within the
scaled FP64/inner-solve bound on the accepted adjoint discretization. If the
literal Jacobian is measurably nonsymmetric, A3 must use BiCGStab and record
the reason; it may not symmetrize the product map to retain CG.

## 5. Inner and outer solve contracts

The inner momentum control is derived from the existing Task 11 momentum
control and is never weaker. Each inner solve must independently recompute its
FP64 residual. The A3 apply reports/accumulates inner iterations, matvecs,
preconditioner calls and reductions without modifying committed flow
diagnostics until the trial succeeds.

Outer success still uses the frozen independent FP64 residual contract. The
final continuity equation is recomputed from the stored corrected face flux,
not from the Krylov recurrence. The full cell identity is recomputed from
`A delta_u+C delta_p+B delta_g`.

## 6. Transaction, MPI and performance

All A3 workspaces and immutable routes are created before iteration. Apply and
face update allocate no `FieldStorage`, route or coefficient plan. Any-rank
inner/outer failure produces one collective classification and lowest failing
rank. Failed attempts restore state, wall data, pressure revision, face fields
and metadata bitwise.

The initial correctness backend may be slower than A2. Task 11 acceptance
still requires the approved performance evidence; optimization may reuse the
same exact Jacobian through a block or better Schur preconditioner, but cannot
replace it with a different numerical map.

## 7. TDD gates

1. tests-only exact-predictor Jacobian oracle RED against current A2;
2. factored private predictor-increment helper GREEN and mutation closure;
3. A3 homogeneous/affine operator and face-divergence identity;
4. constant/reference, checkerboard energy, symmetry/nonsymmetry decision;
5. one-rank frozen checkerboard with exactly two correctors and parity
   `<=1e-8`;
6. 1/2/4-rank direct algebra, rollback and checkerboard;
7. 12/24 localization, 24/48 screening, then the single final Release
   24/48/96 matrix only after all direct gates are green.

Failure at a gate stops A3. No threshold relaxation, damping, filtering,
per-case tuning or extra corrector is authorized.

## 8. Scope

Only Task 11 pressure/predictor private implementation, tests-on access,
focused tests and Task 11 records may change. G1/G2/G3, Task 12, Stage 4,
private legacy source, publication and push remain forbidden.

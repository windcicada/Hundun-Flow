# HUNDUN-FLOW Task 11 A22-A2-P Predictor-Consistency Design

Status: approved diagnostic continuation on 2026-08-04. This document is an
additive Task 11 authority only. It does not accept the rejected A2 candidate,
authorize an augmented projection implementation, or enter G1/G2/G3, Task 12
or Stage 4.

Accepted Task 10 base: `0db56e463470dd1a605709ba05d8bd6a900f496b`

Unchanged candidate HEAD at freeze:
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8`

Predecessor evidence:
`.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-a2-checkerboard-stopping-report.md`

## 1. Question and bounded hypothesis

The A2 operator, affine source, conservative face correction and complete A22
cell-momentum correction close their direct identities, yet the next PISO
predictor leaves a pressure parity amplitude of `4.6065766642193698e-5` after
the frozen two correctors. The only A2-P question is:

> Does `assemble_face_predictor()` reproduce the same current-pressure Schur
> action that `correct_active_pressure()` solves, or does reassembly omit,
> duplicate or sign-reverse the A2 interface action?

A2-P is tests-only until this question has a mutation-sensitive answer. It
does not authorize a current-pressure filter, damping, case tuning, threshold
change, third corrector or selectable alternate operator.

## 2. Exact diagnostic algebra

Let `p` be the pressure at the start of one corrector, `q=S p`, and let

```text
b0 = -S^-1 D mdot_star
r  = b0 - s_g
Kp = K_RC q
Ep = L_E q
```

where `r` is the exact right-hand side passed to the pressure solver after the
single affine-source subtraction and any existing closed-domain projection.
For a pure current-pressure mode to be removed by `delta p=-p`, the reassembled
predictor and solved operator must satisfy

```text
F_hybrid  = r + Kp + Ep = 0.
F_compact = r + Kp.
```

All vectors are also reported in physical per-cell-volume scaling by dividing
the normalized vector by `sqrt(V)`. The diagnostic stores both vector values
and parity/L2/Linf summaries; summaries never replace the vector oracle.

The missing-interface hypothesis is proved only if all of the following hold
on the frozen checkerboard case:

1. `F_compact` closes at the scaled FP64 accumulation bound;
2. `F_hybrid` equals `Ep` at the same bound and is non-negligible;
3. changing the interface sign, omitting it or duplicating it is detected;
4. the observed next-corrector pressure reintroduction has the same signed
   parity authority, rather than merely a similar norm.

If these conditions do not hold, A2-P rejects the small repair hypothesis and
the next design must be an explicitly approved augmented velocity-pressure
projection.

## 3. Smooth-field and no-double-count evidence

The diagnostic also evaluates constant, linear and quadratic cell-average
pressure fields with analytic wall-normal gradients. It records separately:

- the compact `K_RC` response;
- the A2 interface `L_E` response;
- the affine `s_g` response;
- the complete current predictor-divergence response;
- the residual of the independently derived hybrid identity.

No A2-P product term may be accepted merely because it removes a parity mode.
A proposed term must be derived as the exact missing part of the predictor
identity and must not duplicate the full A22 cell-gradient contribution
already present in the momentum predictor. Constant response must be zero at
the scaled FP64 bound. Linear and quadratic probes use their derived discrete
identity and grid-convergence behavior; they are not required to vanish unless
the derivation says so. A mutation that adds the same smooth contribution
twice must fail.

## 4. Tests-only record and state neutrality

`Stage3CellPressureCorrectionRecord` may gain tests-only vectors for:

```text
compact_pressure_action_per_volume
interface_pressure_action_per_volume
hybrid_pressure_action_per_volume
affine_wall_source_per_volume
compact_predictor_defect_per_volume
hybrid_predictor_defect_per_volume
```

The values are captured inside the same corrector after coefficient
preparation and before pressure mutation, so corrector-one evidence cannot be
silently evaluated with corrector-two coefficients. The record remains absent
from tests-off products and creates no public API or object-layout change.

The probe must not change committed/history/trial fields, pressure revision,
allocation identity, solver reports, controller state or numerical counters.
Repeated queries must be bitwise identical. Nested-size, ordinary-value,
aggregate, sign, omission and duplication mutations must all be detected by
one authoritative helper.

## 5. Disposition rules

- **P1 proved:** freeze a minimal product addendum that inserts exactly the
  missing conservative current-interface predictor action, followed by a new
  RED, direct one-rank checkerboard gate, 1/2/4-rank gate and smooth-field
  order screen.
- **P1 disproved or ambiguous:** make no product change; freeze an augmented
  coupled-projection proposal for separate approval.
- **Any algebra/state-neutrality failure:** repair the diagnostic before
  drawing a numerical conclusion.

No 12/24/48/96 numerical matrix is run during diagnosis. Large numerical work
remains forbidden until a focused direct product gate is green.

## 6. Scope

Allowed files are the tests-only portions of `flow/src/stage3_flow.cpp`,
`flow/src/stage3_flow_test_access.hpp`, `tests/mpi/test_immersed_piso.cpp` and
Task 11 evidence/design/plan records. Product behavior requires a separate
post-diagnostic addendum. No public schema, Restart, diagnostic stable ID,
plugin ABI, dependency, private legacy source, publication or push is changed.

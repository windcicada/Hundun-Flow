# HUNDUN-FLOW Task 11 signed-force semantic amendment

Status: frozen for RED-S1 GREEN implementation on 2026-08-06
(Asia/Shanghai). This amendment repairs one contradictory Stage 3 reporting
contract. It does not accept RED-S1, Task 11, or Stage 3.

## Authority and scope

This document is narrower than, and otherwise preserves, the Stage 3 design
and its science-first execution amendment. It is authorized by the accepted
RED-S1 executable evidence and supersedes only the signed-force vocabulary and
comparison clauses enumerated below.

Controlling inputs are:

- `docs/superpowers/specs/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale-design.md`,
  SHA-256 `9220bc03da8df13d55f9b8a0a311879af1191e052f8aae2dbd443eb4be3f6de3`;
- `docs/superpowers/plans/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale.md`,
  SHA-256 `cafce7917c27a77835a8649320b4279a38a0ea7931972a3da4287148c8828b56`;
- `docs/superpowers/specs/2026-08-05-hundun-flow-task11-force-consistency-authority-addendum.md`,
  SHA-256 `94d4d5655f774cedf1191f69caf846013fcbb7fdb0c77f0cc74c01e935a3d6aa`;
- `docs/superpowers/plans/2026-08-05-hundun-flow-stage3-science-first-execution-amendment.md`,
  SHA-256 `0cb469e3f8aa1310d175b7c87dc6403c4ccd8a4b8cb0c14eb7de755700262c60`;
- `docs/superpowers/plans/2026-08-05-hundun-flow-task11-red-s1-signed-force-authority.md`,
  SHA-256 `b2514f3697d134499d84f9d3e55871c4f2acd98784513701610c5127ca43b767`;
  and
- `.superpowers/task-11-red-s1-coordinator-acceptance-2026-08-05.md`,
  SHA-256 `6ecfa95d47987b26c4ff6f8d09d7194be1531f35b8387551827fb4dc797c9378`.

Historical handoffs, reviews, failed candidates and diagnostic records remain
evidence. This amendment does not rewrite or retroactively accept them.

## Proven sign identity

Let `n_s` point from solid into fluid and let the fluid stress be

```text
sigma = -p I + tau.
```

For each pressure, viscous and total component, define:

```text
F_operator = sum of the raw immersed wall-row residual contribution
             = physical discrete force exerted by fluid on solid

B_budget   = sum of the actual momentum-budget reaction accumulator
             = -F_operator

F_surface  = integral_surface (sigma n_s) dA
             = physical true-surface force exerted by fluid on solid

C          = F_operator - F_surface.
```

The accepted one-link RED proves independently that the product wall row has
the `F_operator` orientation, that every actual accumulator delta has the
`B_budget` orientation, and that the real wall quadrature has the `F_surface`
orientation. The aggregate probe proves that the reduced adapter report is
`B_budget = -sum(raw rows)`. Therefore the current candidate's assignment of
that adapter report to a field presented as the physical operator result is a
reporting error; it is not evidence that the row or surface formula is wrong.

Pressure, viscous and total are separate contracts:

```text
F_operator.total = F_operator.pressure + F_operator.viscous
B_budget.total   = B_budget.pressure   + B_budget.viscous
F_surface.total  = F_surface.pressure  + F_surface.viscous
C.total          = C.pressure          + C.viscous.
```

## Public attempt-report contract

`ForceAttemptReport` is attempt-local derived evidence and shall expose four
unambiguous fields:

```cpp
struct ForceAttemptReport final {
  immersed::ForceComponents operator_force;
  immersed::ForceComponents budget_reaction;
  immersed::ForceComponents surface_traction;
  immersed::ForceComponents consistency;
};
```

- `operator_force` is `F_operator`, in the same physical fluid-on-solid
  orientation as `surface_traction`.
- `budget_reaction` is `B_budget`, the opposite-sign quantity used to close
  the discrete fluid momentum budget.
- `surface_traction` is `F_surface`.
- `consistency` is `operator_force - surface_traction` component by component.

There is no `operator_reaction` compatibility member and no alias that
presents `budget_reaction` as a physical force. Stage 3 is unaccepted, the
report is not serialized, and no accepted checkpoint or restart ABI contains
this type, so the explicit source-contract repair is required instead of a
deprecated ambiguous name.

The public lower-level adapter report is renamed without changing its layout
or accumulation:

```cpp
struct ImmersedOperatorReport final {
  // unchanged counters
  ImmersedResidualParts budget_reaction_N;
};
```

`budget_reaction_N` is the authoritative rank-local actual `B_budget`
accumulator. `ForceAttemptReport::budget_reaction` is its collective MPI sum.
The old public member name `operator_reaction_N` is removed; keeping it would
retain exactly the operator/budget ambiguity this amendment closes. No
compatibility member or alias is permitted. This source-contract rename does
not alter the adapter row, sign, reduction, storage layout or numerical
behavior.

## Single conversion authority

The shared final report assembler used by both `collect_final_force()` and the
guarded test seam is the only conversion point:

```text
report.budget_reaction = reduced adapter budget
report.operator_force  = -report.budget_reaction
report.surface_traction = reduced physical surface traction
report.consistency      = report.operator_force - report.surface_traction.
```

No caller may negate again, reconstruct a second raw-row sum, substitute the
surface result for the operator result, or maintain a second report algorithm.
All finiteness checks cover all four fields.

## Precisely superseded Stage 3 wording

Only the following statements in the 2026-07-27 Stage 3 design are replaced:

1. Section 9.2's name `operator_reaction_force` for the negative sum of the
   immersed row residual becomes `budget_reaction`. The positive raw-row sum
   is `operator_force`. Its convergence comparison against the surface
   quantity is `operator_force - surface_traction_force`.
2. Section 12's required force report now contains `operator_force`,
   `budget_reaction`, and `surface_traction_force`; its three consistency
   residuals are `operator_force - surface_traction_force`.
3. Section 18.2's formal pressure, viscous and total consistency rows use
   `operator_force - surface_traction_force`. Their two-segment `>= 1.8`
   thresholds and non-degeneracy rules are unchanged.
4. Section 18.5's transaction/decomposition determinism requirement applies
   independently to all four attempt-report fields. A rejected attempt still
   publishes no force report and persists none of them.
5. The Stage 3 implementation plan's frozen `ForceAttemptReport` interface at
   lines 665--669 is replaced by the four-field interface above, and its
   `ImmersedOperatorReport::operator_reaction_N` interface at lines 2060--2065
   is renamed to `budget_reaction_N`.
6. Section 3 of the 2026-08-05 force-consistency authority addendum is now
   resolved by the accepted executable RED and this vocabulary. Its hard
   gates, independent-authority rules and closure order remain unchanged.

Every other Stage 3 clause remains controlling, including the unique residual
chain, conservative pressure/flux authority, two PISO correctors, rollback,
collective failure, MPI decomposition and all formal convergence thresholds.

## Forbidden interpretations and changes

This amendment does not authorize:

- changing the sign or coefficients of an immersed pressure or viscous row;
- changing the adapter's budget accumulation;
- reversing a wall normal or changing wall quadrature;
- changing Ghost-Cell/LFP reconstruction, pressure closure, Rhie--Chow,
  Poisson, PISO or trial/rollback behavior;
- deleting or weakening a force observation or convergence row;
- changing a scientific threshold, case parameter, selector algorithm,
  damping or corrector count; or
- publishing, pushing, or accessing private reference source or research
  data.

The amendment becomes accepted capability only when RED-S1 GREEN, its formal
mutations, independent reviews and coordinator acceptance all pass on one
frozen candidate. Until then it is an implementation contract, not a Task 11
verdict.

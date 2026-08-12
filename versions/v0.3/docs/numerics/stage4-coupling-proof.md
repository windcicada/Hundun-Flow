# Stage 4 partitioned C-T-C coupling proof

## Frozen schedule

For one accepted step from `t_n` to `t_n + dt`, Stage 4 uses:

1. chemistry on `[t_n, t_n + dt/2]`;
2. one conservative scalar transport over `dt`;
3. PISO corrector 1 on the midpoint predictor;
4. chemistry on `[t_n + dt/2, t_n + dt]`;
5. PISO corrector 2 consuming the post-chemistry-2 thermodynamic state.

There are exactly two chemistry calls and exactly two PISO corrections. No
damping, filter, third corrector, case-specific threshold, or endpoint-rate
source is part of this schedule.

## Non-commuting analytic proof

The focused proof uses

`C=[[0,1],[0,0]]` and `T=[[0,0],[1,0]]`.

Both subflows are exact because the matrices are nilpotent, while the combined
solution is `exp((C+T)t)`, whose entries are `cosh(t)` and `sinh(t)`.
The symmetric composition
`exp(C dt/2) exp(T dt) exp(C dt/2)` has measured local slope above 2.9
and repeated-step global slope above 1.95. The tests explicitly reject
`C-T`, `T-C`, full-duration chemistry stages, and incorrect stage time.

## Fixed-two-PISO flux decision

A midpoint predictor flux alone is insufficient after chemistry stage 2:
composition and temperature change density and therefore the final divergence
constraint. Reusing the predictor would make the transported all-species sum
consume a flux with a different thermodynamic epoch.

Stage 4 therefore selects a conservative predictor-to-final delta flux inside
PISO 2. The delta is applied to the existing accepted face-mass-flux authority;
it does not create a reacting-only face flux. PISO 2 must report the same epoch
as post-chemistry-2 state and publish the final corrected provenance. This is
the minimal fixed-two-corrector choice that closes total mass, species sum, and
final divergence simultaneously.

The rejected candidate is “midpoint predictor flux with no final delta.” It is
not retained as a runtime option.

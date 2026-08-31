# HUNDUN-FLOW v0.4 Adaptive IBM Reconstruction Plan

Date: 2026-08-31

Status: implementation candidate complete; focused verification complete;
coarse-grid product admission must be recorded by an immutable Evidence V6
run and is not a v0.4 release gate.

## Scope and root cause

The original IBM compiler admitted only a full three-dimensional quadratic
stencil. That is a useful accuracy oracle, but it makes the whole case fail when
one local interface group cannot uniquely determine the ten polynomial
coefficients. Coarser grids and combustion-chamber geometry make that local
event more likely through narrow gaps, holes, thin walls, corners, adjacent
surfaces, strong curvature, special grid alignment, or an MPI partition cutting
the available donor neighborhood.

The production correction is a local, explicit, auditable order hierarchy. It
does not relax donor-side geometry or material semantics and does not turn an
invalid stencil into a nearest-cell copy.

## Accepted policy contract

Two runtime policies are accepted:

- `strict_quadratic` is the verification and high-accuracy mode. It preserves
  the original 14--32 donor, three-normal-band, tangential-coverage, full-rank,
  condition-number and reach gates. Any quadratic failure aborts compilation.
- `adaptive_order` is the engineering mode. For each group it tries the same
  quadratic at the standard search reach, expands the admissible donor search
  through the existing maximum reach, and then locally retries a full-rank
  three-dimensional linear basis `[1,n,t1,t2]`. The linear row requires at
  least six donors, two normal bands, reachable tangential coverage, positive
  wall-normal fluid-side donors, and the same `1e8` condition ceiling.

Both policies fail closed if their last permitted reconstruction cannot be
compiled. Constant and nearest-point-copy reconstruction are forbidden.
Reduced-dimensional quadratic reconstruction is not selected from numerical
rank deficiency. It remains a separate future feature requiring an explicit
two-dimensional or axisymmetric geometry certificate.

The selected order is shared by the compiled affine row family. Momentum,
enthalpy, species and passive-scalar boundary equations therefore consume the
same geometry/order decision while retaining their declared value, gradient,
thermal and composition boundary functionals.

## Immutable plan and audit

The input policy participates in case serialization, semantic identity and the
case fingerprint. Each stencil group records reconstruction order, search
reach, condition estimate, functional L1 norm and, for a fallback, the first
quadratic failure category: donor count, coverage, rank, or condition.

The product plan aggregates, with MPI ownership semantics:

- total, quadratic and linear group counts;
- groups that required expanded search;
- fallback counts by failure category;
- maximum accepted condition estimate and functional L1 norm.

Zero-group ranks publish a valid zero-count local audit so a distributed case
does not fail merely because one partition contains no immersed surface. The
thin-domain runner writes the global policy and group counts to `RUN.meta` and
prints a concise warning-style summary when IBM is active. These records are
diagnostic and do not weaken any runtime physical gate.

## Falsifiable verification

The implementation is accepted only with all of the following:

1. A deliberately rank-deficient quadratic donor set must fail atomically in
   `strict_quadratic`, then compile a rank-four linear stencil in
   `adaptive_order` and exactly reproduce a linear field.
2. A stencil recoverable by expanding search reach must remain quadratic; it
   may not prematurely change order.
3. Invalid `nearest_copy` input must be rejected. Explicit strict and adaptive
   cases must have different immutable fingerprints, while the legacy omitted
   policy must remain strict.
4. The product must compile the same immersed plan at 1, 2 and 4 MPI ranks, and
   a rank with no local IBM groups must retain a valid zero-count audit.
5. A previously rejected coarse Re=3900 cylinder case must pass a single-job
   128-rank compile/preflight with all ranks bound one-to-one to physical cores.
   The run must publish the actual quadratic/linear/expanded counts and pass
   positivity, EOS, continuity, energy, closed-mass, gauge, dual-CFL and AFC
   health gates.

The coarse-grid preflight is an engineering admission for the selected mesh,
not proof that every combustor geometry is reconstructible. A geometry that
cannot form even the certified linear stencil still fails closed and must be
re-meshed or supplied with a future explicit reduced-dimensional certificate.

## Benchmark execution guard

Only one CFD case may run on the host at a time. Every MPI launch must acquire
the shared benchmark lock, reject an already-running `mpirun`, set threaded
math libraries to one thread, and use at most 128 ranks with one rank bound to
one physical core. HUNDUN and COAST comparisons run sequentially from matched
physical conditions; no concurrent throughput measurement is admissible.

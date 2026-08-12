# HUNDUN-FLOW Task 11 R1 Authority Equivalence Addendum

Status: coordinator-frozen clarification of the approved R1 mathematical
contract. This addendum narrows the meaning of "same authority" without
changing A22, A3, Task 11 thresholds, or any public interface.

## 1. Authoritative object

For an active immersed row `P`, the authority is the unique least-squares
quadratic polynomial defined by:

- the sorted donor union `D(P)`;
- the cell-average moment matrix on that donor order;
- the row-owned scale and polynomial space; and
- the stable row/link association.

Two private reconstruction objects represent the same authority when they use
that same ordered donor projection and differ only by an invertible change of
quadratic basis (translation, rotation, or scale). They need not share a C++
object address, cached QR storage, or a basis frame. Within the approved FP64
conditioning bound they must reproduce the same value and gradient at a common
physical point.

This is a mathematical equivalence contract, not permission to select a new
donor set, perform an independent per-link fit, or average link geometry.

## 2. Consumer responsibilities

- Predictor and history paths evaluate every compatible wall value or gradient
  from the one unconstrained row polynomial `q_P`.
- A22 consumes all link-normal data `g_l` jointly because the complete cell
  pressure row is the constrained functional of the whole row. It may own the
  constraint nullspace/particular factorization required for that functional.
- Wall-force integration remains an independent true-surface quadrature. It
  re-expresses the same donor projection at each actual quadrature point and
  normal, then consumes the associated authoritative pressure gradient. It
  must not reuse cell-face area, merge wall points, or replace force by an
  operator reaction.
- A3 consumes the resulting A22 maps and retains its exact predictor/Schur
  identities.

## 3. Required evidence

Acceptance requires mutation-sensitive tests proving all of the following:

1. every consumer uses the exact sorted row donor union;
2. basis-re-expressed objects agree at a common point on non-polynomial donor
   data within a condition-scaled FP64 bound;
3. the old independent per-link fits disagree on a discriminating multi-link
   row, so the shared-authority oracle cannot pass vacuously;
4. the complete A22 reaction is independently reproduced from its full
   constrained row coefficients and authoritative wall gradients;
5. wall force remains a distinct true-surface integral; and
6. stable fingerprints, 1/2/4-rank behavior, conservation, rollback, and the
   approved convergence thresholds remain unchanged.

No relaxed threshold, tuned grouping, damping, filter, third corrector, or
post-processing is authorized by this clarification.

# HUNDUN-FLOW WallQuadrature feasible-owner addendum

Date: 2026-08-06

Status: proposed Task 11 blocker repair; supersedes the owner-selection
sentence in Stage 3 plan Task 5 Step 6 and supplements the Halo-provisioning
contract after independent review.

## Evidence and problem

The accepted Stage 3 design requires every triangle to be integrated exactly
once, with decomposition-independent surface quadrature and a collective
failure when the required fixed Halo is unavailable.  The implementation plan
made the execution owner the lowest rank owning any of the triangle's three
associated active rows, tied by active global cell ID.

Task 11's rotated, partition-conforming fixture exposed that rank priority is
not an executability proof.  Its first 2-rank failure has three association
owners `0,0,1`.  The frozen tuple rule selects rank 0, but the complete point
and pressure-authority donor union contains logical cell `(10,1,6)`, outside
rank 0's expanded half-open x interval `[-4,10)`.  The independently
recomputed feasible-owner mask is `0x2`: rank 1 can execute the complete
triangle with the existing reach four.  The disposable diagnostic was removed
and the product source and binary were restored byte-for-byte.

Therefore the failure is neither an unavailable-Halo fixture nor a validation
bug.  It is a product capability gap caused by choosing owner rank before
checking the already-required donor coverage.

## Revised ownership contract

For triangle `t`, construct all three existing point reconstructions and their
existing pressure-authority reconstructions without changing association,
donor selection, weights, or reach.  Define the complete donor set

```text
D_t = union over q in {0,1,2} of
      (D_point(t,q) union D_pressure_authority(t,q)).
```

For decomposition rank `r`, let `B_r` be its owned logical-cell box and let
`expand(B_r,4)` use the existing componentwise half-open predicate:

```text
B_r.begin - 4 <= donor < B_r.end + 4.
```

The feasible set is

```text
F_t = { r | every valid donor in D_t lies in expand(B_r,4) }.
```

The triangle execution owner is the numerically lowest rank in `F_t`.  If
`F_t` is empty, construction still fails collectively with the existing
owner-Halo classification.  Selection is made from the complete global rank
set, not only the association-owner set: execution ownership means that the
rank has all required field data, while link pressure authority remains owned
by the link's unchanged fluid-row owner.

All three quadrature points of a triangle retain the same execution owner.
Every triangle is still integrated exactly once.  Per-point ownership is not
introduced.

## Wall-plan Halo contract

`expand(B_r,4)` is an executability proof only when every wall-force input
field on rank `r` actually has four populated logical Halo layers.  The
link-local `GhostStencilPlan::maximum_halo_reach()` is an exact measured reach
and may be smaller than four; it is not sufficient for an execution owner
chosen relative to a different owned box.

Stage 3 therefore gives `WallQuadraturePlan` its own fixed required reach:

```text
WallQuadraturePlan::maximum_halo_reach() == 4
```

Whenever a wall plan is present, field descriptors and `ExchangePlan` use

```text
max(GhostStencilPlan::maximum_halo_reach(),
    WallQuadraturePlan::maximum_halo_reach()).
```

An existing caller that provisions the literal width four is already exactly
conforming, not an exception: the GhostStencilPlan contract is at most four
and the wall-plan contract is exactly four, so its allocation equals this
maximum.  Such fixed-four callers remain regression-only in Task 1D; they do
not need a mechanical rewrite to spell the same value dynamically.  The
dynamic ghost-only callers exposed by the transaction and Task 11 fixtures
spell the explicit maximum for forward correctness even though the frozen
current fixtures both measure four.  Runtime construction and integration
checks remain authoritative if a future real GhostStencilPlan measures less.

The current fixed wall reach deliberately favors correctness and a usable
global owner rule over minimizing Halo storage.  Measuring a smaller safe
wall reach is deferred performance work and may not change ownership or
numerical results.

`FixedStepStage3Flow` rejects an under-provisioned supplied Halo during
construction.  `WallForceIntegrator` independently checks the actual input
`FieldView` ghost widths collectively before any donor dereference, so a
standalone or future caller cannot silently bypass the contract.

## Preserved contracts

This addendum does not change:

- triangle geometry, degree-2 three-point quadrature, weights, normals, stable
  IDs, association link selection, or link tie-breaking;
- point or pressure-authority donor construction, QR, LFP, pressure boundary,
  force traction, or any numerical coefficient;
- the fixed maximum reach four, owner boxes, half-open containment predicate,
  or no-fallback collective failure; the new wall-plan reach contract only
  ensures that the already-authorized four layers are actually provisioned;
- link boundary-authority owner/catalog semantics;
- the rule that one rank owns all three points and integrates each triangle
  once; or
- the decomposition-independent wall-plan fingerprint.  Raw execution rank is
  not added to that fingerprint.  The separate future Checkpoint/static
  `triangle_ownership` identity remains authoritative for same-process-grid
  owner-map compatibility and is not removed or weakened.

Checkpoint v2 payload, manifest, fingerprint, and restore behavior remain
unchanged: they do not currently consume a wall plan.  Execution owner is
recomputed during plan construction and is not serialized as cached plan
state.  This addendum does not claim immersed Restart coverage; the frozen
future Checkpoint v3 work must rebuild and compare its separate
`triangle_ownership` identity.

The implementation must construct the three candidate point reconstructions
before selecting the execution owner.  It must not reject the lowest rank
early and must not discard donors to make a rank feasible.

## Required RED and acceptance evidence

The frozen Task 1C snapshot is the integration RED.  A Task 11 test-owned
oracle must also
prove on the 2-rank fixture that triangle 0 has:

```text
point authority owners = [0,0,1]
complete donor logical x range = [4,10]
rank 0 reach interval = [-4,10)   -> infeasible
rank 1 reach interval = [2,16)    -> feasible
triangle execution owner = 1
```

The oracle must include point and nested pressure-authority donors, retain the
existing exact three-point/one-owner coverage checks, and fail under the old
lowest-associated-owner mutation.

The unique production owner-selector helper must additionally be called by
focused synthetic tests that are independent of the Task 1C geometry:

1. the only feasible rank is outside a test-annotated association-owner set;
2. at least two ranks are feasible and the numerically lowest is selected;
3. no rank is feasible and the helper reports no selection; and
4. point and pressure-authority donor sets deliberately differ, so ignoring
   either dependency changes the result.

Product integration must convert case 3 into the same collective owner-Halo
classification and lowest-failing-rank contract.  A wall-force input with
only three ghost layers must be rejected collectively before donor access.
The mutation packet must include the old lowest-associated rule,
association-only scanning, highest-feasible selection, defaulting when the
feasible set is empty, and omission of the divergent pressure-authority donor.

The Halo evidence must be mutation-sensitive without inventing a numerical
fixture.  A read-only survey measured
`GhostStencilPlan::maximum_halo_reach() == 4` for every existing clean
integrated candidate, including the exact two-rank A1/A2 Task 11 fixture.  A
natural `ghost_reach <= 3` comparison is therefore unavailable in the current
product test space and must not be manufactured by coordinate, threshold, or
selection tuning.

Use two replacement product-path REDs:

1. On one rank, all donors of the existing wall-force plan lie in the owned
   global box.  The pre-change integrator therefore succeeds with three-layer
   fields.  After the contract change, the identical plan and field values
   must be rejected collectively before donor access solely because
   `3 < WallQuadraturePlan::maximum_halo_reach() == 4`; width four must pass.
   Removing the pre-dereference validator must restore the old success.
2. In the existing immersed transaction, an explicitly three-layer Halo
   currently reaches the older generic ghost-plan-insufficient failure.  The
   new wall-specific collective construction check must run first and expose
   its exact distinct classification; width
   `max(ghost_reach, wall_reach) == 4` must construct successfully.  Removing
   or moving the new check after ImmersedReconstruction must restore the old
   classification and fail the oracle.

The Task 11 path separately freezes its actual ghost reach at four.  Dynamic
callers still spell `max(ghost,wall)`, but acceptance does not falsely claim a
currently natural `ghost<wall` case.  The public wall reach and the two runtime
validators make future lower-reach callers fail closed.  The first registered
real GhostStencilPlan that measures below four must add the natural combined
Halo regression before it can be accepted.

The real mixed-owner Task 11 force consumer must exercise the supplied
pressure-authority-gradient path, not only the independent point-gradient
fallback.  For triangle 0, rank 1 executes points whose link 77 authority
remains on rank 0.  Supplying the complete link-gradient catalog from the
unchanged authority ranks must succeed and preserve force, moment, and global
point count; moving link 77's datum to execution rank 1 must fail with the
existing wrong-provider collective contract.

Finally, an approved-decomposition comparison must first prove that at least
one triangle's execution owner/local distribution differs and then require the
wall-plan fingerprints to remain bitwise equal.  Hashing raw owner rank or a
decomposition-specific feasible mask must be killed by this test.

After GREEN, require serial 1/2/4 Task 11 preflights, the focused existing
wall-quadrature and wall-force MPI controls, then Task 11 behavioral RED on
1/2/4 ranks.  No force result may influence owner selection.  Tests-off proof
and an explicit old-rule restoration mutation remain mandatory before RED-S1
acceptance.

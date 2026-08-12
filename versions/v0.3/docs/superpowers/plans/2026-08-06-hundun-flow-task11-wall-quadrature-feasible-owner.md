# Task 11 WallQuadrature feasible-owner implementation plan

Date: 2026-08-06

Authority: `docs/superpowers/specs/2026-08-06-hundun-flow-wall-quadrature-feasible-owner-addendum.md`

## Scope

Repair only the WallQuadraturePlan execution-owner decision exposed by Task
11.  Keep the A1 fixture and A2 loader repair frozen.  Do not change the
maximum allowed reach four, donor selection, association selection, QR, force
formulas, LFP, PISO, thresholds, or process grids.  This task does add the
missing wall-plan required-Halo contract and caller validation.

Implementation ownership is limited to:

```text
immersed/include/hundun/immersed/ghost_stencil_plan.hpp
immersed/src/ghost_stencil_plan.cpp
immersed/src/ghost_stencil_plan_internal.hpp
immersed/src/wall_force.cpp
flow/src/stage3_flow.cpp
tests/unit/test_ghost_stencil_header_contract.cpp
tests/mpi/test_wall_quadrature_plan.cpp
tests/mpi/test_wall_force.cpp
tests/mpi/test_immersed_transaction.cpp
tests/mpi/test_task11_signed_force_authority.cpp
```

Coordinator evidence, reviews, and plan records may be added under the
existing Task 11 RED-S1 SDD directory and `.superpowers/`.

The other known integrated test/support callers provision a literal four:

```text
tests/mpi/test_immersed_closed_transient.cpp
tests/numerical/test_laminar_ibm_engineering.cpp
tests/numerical/test_laminar_ibm_order.cpp
tests/support/stage3_mms.cpp
```

They are deliberately outside the Task 1D edit boundary.  Because the Ghost
plan is bounded by four and the new Wall plan reach is exactly four, their
literal allocation is already equal to `max(ghost,wall)`.  They remain normal
regression/acceptance consumers and must not be mechanically rewritten in
this focused repair.  Only the dynamic ghost-only transaction and Task 11
callers require source changes.

## Task 1D.1: Freeze the RED and add the independent oracle

1. Record HEAD, addendum/plan hashes, current product/test hashes, Task 1C
   report/log hashes, focused binary hash, and empty index.
2. Extend `PlanFixtureObservations` in the Task 11 test.  On two ranks,
   independently collect triangle 0's three point-index authority owners and
   both the point and nested pressure-authority donor logical bounds.
3. Require authority owners `[0,0,1]`, donor x bounds `[4,10]`, rank 0
   infeasible under the unchanged half-open reach-four box, rank 1 feasible,
   and one execution owner equal to rank 1.  Preserve the existing three-point
   mask and one-owner oracles.
4. Build the focused target with `-j2` and reproduce the existing 2-rank
   owner-Halo RED before editing product code.  The new post-plan observation
   is unreachable under the old rule; the exact existing collective failure
   remains the RED.
5. Through the unique production owner-selector helper, add synthetic focused
   cases for: a unique feasible rank outside a test-annotated association set;
   multiple feasible ranks selecting the lowest; no feasible rank; and
   deliberately divergent point/pressure donor sets whose union changes the
   result.  These tests call product logic and do not copy it.
   The no-feasible case must traverse the same helper throw and collective
   classification used by `WallQuadraturePlan::create`, and must lock the
   exact owner-Halo message plus lowest failing rank 0.
6. Add an approved-decomposition oracle to the existing wall-quadrature test:
   collect the complete triangle owner maps for the default and alternate
   4-rank grids, prove at least one owner/local distribution differs, and only
   then require bitwise-equal wall-plan fingerprints.
7. Freeze the Halo-discriminator adjudication.  The read-only fixture survey
   and exact two-rank A1/A2 GDB observation all report ghost reach four; record
   the binaries, source hashes, commands, exits, and absence of temporary
   residue.  Do not tune a fixture to force a lower value.
8. Add the replacement REDs before product changes.  In the one-rank
   wall-force case, use identical values in three-layer fields and prove the
   old integrator succeeds because every donor is owned; the new expected
   pre-dereference rejection must therefore be RED.  In the immersed
   transaction, supply an explicit three-layer Halo and require the new exact
   wall-specific collective error; the current path's distinct generic
   ghost-plan-insufficient error is RED.  Retain width-four success controls.
   Freeze Task 11's actual ghost reach equal to four.  These replace only the
   unavailable natural `ghost<=3<wall` comparison.

## Task 1D.2: Implement the minimal feasible-owner rule

Implement the owner rule in one internal production helper and call it from
`WallQuadraturePlan::create`:

1. leave the three point/link associations and normal validation unchanged;
2. build all three point-authority candidates from the already-built
   pressure-authority donor sets before assigning any execution owner;
3. evaluate every rank in ascending numeric order using the existing
   `owner_boxes`, `kMaximumReach`, donor global IDs, cell logical coordinates,
   and `within_expanded_box` predicate;
4. select the first rank covering every point and pressure-authority donor;
5. if no rank is feasible, throw the existing owner-Halo error collectively;
6. publish all three points with that one owner while preserving each point's
   independent link authority owner; and
7. do not hash execution rank or alter point/reconstruction ordering.

No helper may maintain a second owner or donor algorithm outside this product
path.

Add `WallQuadraturePlan::maximum_halo_reach() noexcept`, fixed to four for
Stage 3.  Do not change the exact measured
`GhostStencilPlan::maximum_halo_reach()`.  Any construction path that supplies
both plans must size its cell fields and `ExchangePlan` to the maximum of the
two reaches.  `FixedStepStage3Flow::create` must reject a supplied Halo below
the wall reach, and `WallForceIntegrator` must reject actual input
`FieldView`s below that reach collectively before reading any donor.

Raw execution rank remains absent from the decomposition-independent
wall-plan fingerprint.  Do not remove or weaken the separate future static
triangle-ownership identity.

Checkpoint v2 files, schemas, fingerprints, and restore code are a frozen
non-change surface.  Do not serialize WallQuadraturePlan or execution owner in
this task.

## Task 1D.3: Focused GREEN and regression order

Build only the affected targets with `-j2`:

```text
test_ghost_stencil_header_contract
test_task11_signed_force_authority
test_wall_quadrature_plan
test_wall_force
test_immersed_transaction
```

Run serially:

1. `test_ghost_stencil_header_contract`;
2. `fixture_preflight_1_rank`;
3. `fixture_preflight_2_rank`;
4. `fixture_preflight_4_rank`;
5. registered `test_wall_quadrature_plan` success controls on 1/2/4 ranks;
6. registered `test_wall_force` focused ownership/force controls on 1/2/4
   ranks, including the registered 2/4-rank alternate process-grid controls,
   selected by exact CTest names;
7. the under-provisioned three-layer standalone wall-force failure control;
   and
8. the affected immersed-transaction construction/control test proving its
   fields and ExchangePlan use `max(ghost,wall)`.  This control must freeze
   `ghost_plan.maximum_halo_reach() == 4` and
   `wall_plan.maximum_halo_reach() == 4`, require width four to construct,
   and require an explicit width-three Halo to fail with the new exact
   wall-specific collective classification before the older generic
   ghost-plan failure.  It does not claim a natural lower-reach fixture.

On the Task 11 2-rank mixed-owner fixture, also run the real supplied
pressure-authority-gradient consumer.  Build the complete zero-gradient
catalog from each link's unchanged authority rank, require success, and match
the constant-pressure reference force, moment, and exactly-once global point
count.  Then move link 77's datum from authority rank 0 to execution rank 1
and require the existing wrong-provider collective failure.  The ordinary
fallback path remains a separate control.

Stop on the first failure.  Do not run sanitizers, broad suites, or 24/48/96
matrices in Task 1D.

## Task 1D.4: Resume RED-S1

Only if every preflight and focused regression is GREEN:

1. run Task 11 `behavioral_red` on 1, 2, and 4 ranks serially;
2. require every fixture/raw/surface/trace/budget/closure/rank observation to
   be true and only the plan-declared physical-report/consistency observations
   to be false;
3. run the tests-off library build proving both test seams disappear; and
4. generate exact prestate-to-current deltas, hashes, named vectors, commands,
   environments, and exit statuses.

## Task 1D.5: Mutation and review

Run disposable patches restored byte-for-byte afterward for:

1. the old lowest-associated-owner rule, which must restore the frozen
   2-rank owner-Halo failure;
2. restricting candidates to association owners;
3. choosing the highest rather than lowest feasible rank;
4. defaulting to a rank when the feasible set is empty;
5. ignoring either the divergent point donor set or nested
   pressure-authority donor set;
6. hashing execution owner or a decomposition-specific feasible mask into the
   wall-plan fingerprint; the alternate-grid owner-map oracle must kill it;
   and
7. removing the standalone wall-force pre-dereference width validator; the
   one-rank three-layer case must revert to its frozen old success; and
8. removing or moving the Stage3 wall-specific Halo check after the existing
   GhostStencilPlan check; the transaction must revert to the distinct old
   generic error classification.

Mutations 2--5 must be killed by the focused product-helper cases; mutation 6
must be killed by the owner-map/fingerprint oracle; mutations 7--8 must be
killed by the replacement product-path Halo controls.  Restore source and
binary hashes after the packet.  Do not claim a mutation-sensitive natural
`ghost<wall` caller until the product supplies such a real plan; the runtime
checks make that future case fail closed in the meantime.

Then perform, in order:

1. independent requirement/science review;
2. independent code-quality/evidence review;
3. resolution of every Critical or Important finding; and
4. main-agent full diff, log, hash, and fresh 1/2/4 behavioral RED review.

Acceptance of Task 1D means the owner capability blocker is repaired and the
RED-S1 scientific contradiction is exposed.  It does not accept the later
signed-force semantic GREEN.

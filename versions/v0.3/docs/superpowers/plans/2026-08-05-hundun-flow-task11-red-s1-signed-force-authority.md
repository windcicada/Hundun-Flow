# Task 11 RED-S1 Signed Force Authority Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish a mutation-sensitive, MPI-capable one-link executable that independently fixes the pressure/viscous wall-row direction, surface-traction direction, budget-reaction direction, and final force-report comparison before any Task 11 semantic repair.

**Architecture:** A new 12^3 test constructs one closed, rotated Cartesian box through the normal product geometry/classification/GhostStencilPlan/operator/wall-quadrature path. Existing operator test access supplies observed row values and is extended only with the actual per-row accumulator delta; independent expectations are rebuilt from test-owned STL triangles, `ImmersedDomain::links()`, topology, geometry, and analytic stress. Three narrowly guarded observation surfaces trace the real row accumulator and WallForceIntegrator point path and call the same pure report assembler used by `collect_final_force`; none changes a non-test numerical path or supplies an alternative algorithm.

**Tech Stack:** C++20, CMake/CTest, MPI 1/2/4 ranks, project-owned mesh/immersed/finite-volume/flow libraries, Debug and Release Clang/libc++ builds.

## Global Constraints

- Work only in `/home/wyf/code_dev/.worktrees/hundun-flow-stage3`; preserve every inherited tracked and untracked change and keep the index free of unrelated files.
- Accepted Task 10 base remains `0db56e463470dd1a605709ba05d8bd6a900f496b`; RED-S0 is accepted at current HEAD `77eb15cd31a28fab055fc2af18e8f923492e2bd6`.
- The initial RED-S1 phase may add test code, CMake registration, test-only observations, and a behavior-preserving extraction of the current report assembly. It must not change a pressure, viscous, LFP, PISO, force, consistency, tolerance, corrector, or solver formula.
- Product geometry, classification, `GhostStencilPlan`, `ImmersedOperatorAdapter`, and `WallForceIntegrator` must execute. No mock, copied row evaluator, cut-cell replacement algorithm, or alternative force integrator is allowed.
- Expected operator values must not read the final operator row descriptor or use `last_boundary_row_evaluations()` as an oracle. That snapshot is the observed product result only.
- Expected surface values must not call `WallForceIntegrator`, its test trace, or any existing wall-force oracle. Triangle normals, areas, barycentric points, and weights are rebuilt from the test-owned transformed STL triangles; product quadrature geometry is an observed value compared against that fixture, not the oracle.
- The pressure-only, viscous-only, and combined cases share one selected one-link row with an oblique solid-to-fluid normal. Both operator and surface measures must be nonzero and non-identical.
- The test prints every named observation before one final aggregate assertion. A compile/link failure or an early assertion is not the required behavioral RED.
- No threshold relaxation, case-result tuning, filtering, damping, added PISO corrector, private-source access, publication, or push is authorized.
- Long or mutation runs must record exact HEAD, test-only diff SHA-256, binary SHA-256, command, environment, log SHA-256, exit status, and immutable candidate seed.
- Use one implementation worker at a time. The main agent owns the mathematical judgment, full diff review, semantic amendment, mutation adjudication, and final acceptance.
- All eight RED files overlap the inherited, unaccepted dirty Task 11 candidate. Implementation workers must not stage or commit them: doing so would capture unrelated pre-existing hunks. The main agent freezes read-only pre-task copies/hashes, reviews an exact before/after delta, and keeps the index empty. A later coordinator-approved consolidation may commit only after the complete inherited candidate scope is reviewed.

## Frozen baseline

At plan creation the index is empty and the whole tracked patch SHA-256 is:

```text
9fdab54d2fb44919f3bb9cd42c1b657b52157a395feb677349b8a3f9655e740c
```

The RED-S1 input hashes are:

```text
9e2e7590c0993de98b80f97a59e81a1ec69c37b1a2fc01e3c89916b9db3e5eca  CMakeLists.txt
6d362b9114fe4ffbb688851d112fe5b258f48568f385418b690d1c9a1f5b34d9  immersed/src/wall_force.cpp
ac93542abc8444bf95aa31062ce8af237f97cdd986a37098131477fd336ceaf4  immersed/src/wall_force_detail.hpp
5a0eba74f361a1db23cb129eee20fdbbe7d2cc0cc9526e56902b2f34a474acb1  finite_volume/src/immersed_operator_test_access.hpp
21dd801fd69ee01452e8946593ef7ce862a65c8117cad8bb5837527a5adc91d7  finite_volume/src/immersed_operator.cpp
74601d8ba93f981a3ba1af4a2030d7a55a925c0b25ade3f1613ed46ce8574e24  flow/src/stage3_flow.cpp
27b90a43f8a1b9842b58db13a99149ecf56475ddb29f4544b751016d966b2d98  flow/src/stage3_flow_test_access.hpp
94176948bce9b970ccf57a8eb5ee3f1f6ecf7fc6410ad13c0fcda20cb1895f0f  tests/mpi/test_immersed_operator.cpp
1f9676daa3d7da99652be0f5b16c636ba4c56f819a113b85baa002c5b93e42f4  tests/mpi/test_wall_force.cpp
```

`codegraphf` has four pending modified files and is intentionally not synced. Its indexed result for `accumulate_momentum` (finite-volume adapter, two Stage3 flow callers, and the operator MPI test) is advisory; exact `rg` results over the dirty tree control this plan. `collect_final_force` and anonymous wall-force helpers are absent from the stale graph and are controlled by exact source inspection.

The older science-closure design and signed-force audit brief record takeover
HEAD `28f5dd541a0e3ce9ecf852e53d83981add3a5be8`. That identity is historical.
RED-S0 coordinator acceptance at HEAD
`77eb15cd31a28fab055fc2af18e8f923492e2bd6`, SHA-256
`93fc89f23a119898a2c0561b274f0ed66759376731d4180fb671f81d65988bf3`,
supersedes it for RED-S1. All evidence hashes are refreshed from this plan's
actual start HEAD.

## Frozen file boundary

### RED test and observation scaffold

Only these implementation files may change before the behavioral RED is accepted:

```text
CMakeLists.txt
tests/mpi/test_task11_signed_force_authority.cpp             (new)
immersed/src/wall_force.cpp
immersed/src/wall_force_detail.hpp
finite_volume/src/immersed_operator.cpp
finite_volume/src/immersed_operator_test_access.hpp
flow/src/stage3_flow.cpp
flow/src/stage3_flow_test_access.hpp
```

The finite-volume changes are test-only observation: extend
`BoundaryRowEvaluationSnapshot` with the actual per-row delta of the real
`local_reaction` accumulator. Capture the accumulator before and after each
row, then store their difference; do not recompute `-wall_contribution` into a
second budget algorithm. Existing `test_immersed_operator.cpp`,
`test_wall_force.cpp`, and `test_immersed_piso.cpp` are regression targets, not
RED-S1 edit hosts.

Coordinator-only evidence and planning records may be added under `.superpowers/task-11-red-s1-*` and this plan's SDD workspace. No other source file is authorized until the behavioral RED is reviewed and a semantic amendment freezes the GREEN list.

## Fixed mathematical fixture

The surface starts from `tests/support/stage3_stl_fixture.hpp::outward_cube()`. Refine each triangle twice, rotate each centered vertex with the right-handed orthonormal columns

```text
q1 = { 1/3,  2/3,  2/3}
q2 = { 2/3,  1/3, -2/3}
q3 = {-2/3,  2/3, -1/3}
```

then scale by `0.36` and translate to `{0.5, 0.5, 0.5}` metres. Its analytic axis-aligned bounding box is `[0.2,0.8]^3`: each projection half-width is `(0.36/2)*(5/3)=0.30`, leaving `0.20 > 2/12` separation on every side. Load the resulting closed STL with scale `1.0`. Use a uniform `12 x 12 x 12` mesh on `[0,1]^3`, `ImmersedFluidSide::outside`, a non-periodic decomposition, and six `no_slip_wall` outer patches. The ghost width is the product GhostStencilPlan requirement. Process grids are `{1,1,1}`, `{2,1,1}`, and `{2,2,1}` for 1/2/4 ranks.

### Coordinator fixture amendment A1 — partition-conforming tessellation

The first implementation of the paragraph above serialized the 192 rotated
triangles directly. Its force-free 1-rank preflight passed, but construction
of the real `WallQuadraturePlan` on `{2,1,1}` failed collectively with
`wall quadrature donor exceeds owner Halo reach`. This happened before link
selection, field fill, force integration, or report assembly. The preserved
failure is a fixture-plan RED, not a force result.

Stage 3 fixes one owner for all three quadrature points of a triangle: the
lowest rank owning an associated active row, tied by active global cell ID.
Every donor used by any of those three points must be visible through that
single owner's fixed reach-four Halo. The rotated 192-triangle tessellation
contains triangles that strictly cross the approved `x=0.5` process plane;
its points may therefore associate rows on both ranks even though each
link-local reconstruction is valid. The fixed Stage 3 owner rule and Halo
rejection are not changed by RED-S1.

Before STL serialization, split the already rotated/scaled/translated
intermediate triangles along `x=0.5`, then `y=0.5`. This is a tessellation-only
operation: it must preserve the exact closed surface, winding, file normal,
area vector, total area, and `[0.2,0.8]^3` bounding box.

For each plane, retain a triangle once if its vertex coordinates do not
strictly straddle the plane. Otherwise:

1. clip its ordered polygon into the low and high closed half-spaces;
2. compute every edge intersection from the strictly lower-coordinate
   endpoint toward the strictly higher-coordinate endpoint, then set the
   split-axis coordinate exactly to `0.5`, so a shared edge produces the same
   floating-point vertex in either orientation;
3. remove only adjacent bitwise-identical vertices, including an identical
   cyclic first/last pair;
4. fan-triangulate each resulting polygon without changing winding and copy
   the source triangle's analytic file normal; and
5. reject a non-finite or non-positive emitted triangle before writing the
   STL.

No tolerance, vertex perturbation, body translation, grid change, Halo
increase, owner-rule change, or force-dependent choice is permitted. The
test-owned oracle must verify before product construction that no emitted
triangle strictly crosses either process plane; every child geometric normal
has positive dot product with its source file normal; each source triangle's
children preserve its scalar area and oriented area vector; and the complete
tessellated surface preserves total scalar area, oriented area vector, and
AABB within the already frozen componentwise arithmetic bound. A bitwise
undirected-edge catalog must also prove that every emitted edge occurs exactly
twice with opposite directions, including split-plane seams. After plan
construction it must additionally verify that the loaded product triangle
count equals the emitted test-owned triangle count, every local point has
`point.owner_rank == mpi.rank()`, and globally every triangle has exactly
three quadrature memberships with point-index mask `{0,1,2}` whose concrete
owner-rank identity is unique across all three points.

The 1/2/4 force-free preflights are then rerun serially. If any still fails,
or if no common eligible link is selected, stop again and return the complete
fixture evidence to the coordinator. No second fixture adjustment is
authorized from a force observation.

### Coordinator prerequisite amendment A2 — legal shared-vertex contact

A1 emitted 244 triangles and passed every test-owned manifold, orientation,
per-source/global area, AABB, split-plane, and bitwise opposite-edge oracle.
The 1-rank product load nevertheless rejected the surface before query/domain
construction as self-intersecting. Read-only branch instrumentation, restored
bit-for-bit immediately afterward, identified emitted triangles 4 and 88 and
the legacy `coplanar_contact_is_forbidden` path for triangle 88 edge 0.

The independent exact fixture geometry is:

```text
P = (0.47, 0.26, 0.62)          shared bitwise vertex
C = (0.50, 0.32, 0.68)          triangle 4 intersection-line endpoint
D = (0.44, 0.20, 0.56)          triangle 88 intersection-line endpoint
C - P = P - D = (0.03, 0.06, 0.06)
```

The triangles lie on distinct rotated-cube faces. Their plane-intersection
segments are `[P,C]` and `[D,P]`; therefore their closed intersection is
exactly the allowed welded topology vertex `{P}` and their interiors are
disjoint. The legacy implementation trims an allowed shared endpoint by only
`4*coincidence/segment_length`, then tests the still tolerance-expanded point
as though it were an interior overlap. For this legal child geometry the trim
remains inside the expansion and produces a false positive.

A2 supersedes the earlier RED-S1 source-file prohibition only for this loader
prerequisite. Do not change welding, minimum area, coincidence, broad-phase
boxes, full coplanar-triangle overlap, non-coplanar segment hits, shared-edge
handling, or any immersed/force numerical path. Replace the heuristic endpoint
trim inside `coplanar_contact_is_forbidden` with a conservative shared-vertex
tangent-cone predicate:

1. require the allowed endpoint to be the exact welded shared vertex ID in
   both indexed triangles, not merely a coordinate within `coincidence`;
2. project the target triangle, shared vertex, and direction away from that
   vertex along the triangle's existing dominant projection axis;
3. for the two half-spaces incident on the shared vertex, evaluate the
   directional derivative of the consistently oriented half-space function
   in higher precision (`long double` is higher precision where the platform
   provides it, not an exact predicate);
4. take the orientation sign from the nonzero dominant component of the
   already validated geometric normal with the existing projection parity:
   `sign(n.x)` when dropping x, `-sign(n.y)` when dropping y, and `sign(n.z)`
   when dropping z. For projected edge `e` and direction `d`, bound the
   determinant derivative error by
   `64*double_epsilon*(abs(e.x*d.y)+abs(e.y*d.x))`; classify the segment as
   immediately leaving the closed triangle only if at least one signed active
   derivative is strictly below the negative bound. An ambiguous,
   tangent-inward, or interior direction remains forbidden; and
5. by convexity, a segment that provably leaves at its shared endpoint cannot
   re-enter the triangle. The existing exact shared-edge exclusion remains
   authoritative before this helper is called.

No positive-width intersection interval, however small, is ignored. If an
optional half-space clip is used as an additional cross-check, `lo < hi` is
always forbidden; only an exact singleton (`lo == hi`) at the canonical shared
endpoint may be allowed. The old `256*epsilon` parameter window may not be
used to turn positive width into a singleton.

The original A1 1-rank loader failure is the RED. The existing true
self-intersection and shared-vertex-intersection negative fixtures must remain
red/green controls, while the A1 surface must load. Restoring the heuristic
trim after GREEN must make the A1 preflight fail again. No fixture coordinate,
triangulation, threshold, owner, Halo, or force observation changes under A2.

Before any force field is filled, run selector `fixture_preflight` on 1, 2, and
4 ranks. It uses no WallForceIntegrator result. The selected link ID and all
geometry below must agree within the frozen componentwise arithmetic bound.
Failure to find the same valid link is a fixture-plan failure returned to the
main agent; an implementation worker may not change the body, grid, point
count, or selector after observing a force result.

Select the globally smallest link ID satisfying all of the following before reading a force result:

```text
the owned operator row contains exactly this one link
min(abs(n_s.x), abs(n_s.y), abs(n_s.z)) >= 0.2
at least three WallQuadraturePoints carry this boundary-authority link
norm(sum(n_triangle * A_triangle/3)) / sum(A_triangle/3)
  is within 4096*epsilon of one
dot(normalize(sum(n_triangle*A_triangle/3)), n_fixture) >= 1-4096*epsilon
signed_background_measure = -dot(A_background, n_s) > 0
analytic operator pressure and viscous vectors are nonparallel
analytic surface pressure and viscous vectors are nonparallel
```

`WallQuadraturePoint` does not publicly carry an authority link. The test may
use the already existing internal read-only
`detail::boundary_authority_link(point.reconstruction)` only to decide which
public point belongs to the selected link. The expected point position,
normal, and weight are recomputed from that point's public `triangle` and
`point_index` into the test-owned transformed STL triangle. The product
position, normal, and weight are separately required to match these analytic
fixture values; neither trace values nor product quadrature geometry feed the
expected force.

The selected `ImmersedLink` is read on its owner from
`ImmersedDomain::links()` and broadcast. Its `triangle` selects the independent
fixture normal `n_fixture`; with outside fluid, the product
`solid_to_fluid_normal` must match `n_fixture`, and the wall intercept must lie
on that triangle plane. The expected background face is independently found
as the unique topology face joining the link's `fluid_cell` and `solid_cell`;
orient `A_background` outward from the fluid cell using
`MeshGeometry::face_area_vector_m2`. The row snapshot may identify the
selected row and supply the observed raw contribution, but its copied area,
normal, signed measure, pressure quadrature, coefficients, and descriptors are
not used to compute the expectation.

All ranks execute the same collective sequence. They first all-reduce the
maximum link ID seen in local domain links, rows, and quadrature membership,
then allocate the same dense per-link table. They all-gather row eligibility
and all-reduce point count, independent weight, independent weighted normal,
and product-versus-fixture geometry checks by link ID. Every rank chooses the
same minimum eligible ID from that global table. Exactly one domain-link owner
and one row owner are required; owner-side link/row/face metadata are then
broadcast. Every rank calls each trace and report-path collective
unconditionally and in the same order. No candidate-owner-only call is
allowed. Each selected `(triangle, point_index)` membership must occur exactly
once globally; duplicate or missing membership is a preflight failure.

Let `n_q` be the independently rebuilt common fixture normal,
`A_surface = sum(A_triangle/3)` over selected memberships, `n_s=n_fixture`,
and `x_wall` be exactly the selected `ImmersedLink::wall_intercept_m`:

```text
p0 = 2.75
mu = 0.35
t = normalize(cross(e_least_aligned, n_s))
a = 0.70*t + 0.40*n_s
G = a outer n_s
u(x) = G * (x - x_wall)
tau = mu * (G + transpose(G) - (2/3)*trace(G)*I)
sigma = -p0*I + tau
```

`e_least_aligned` is the coordinate axis with the smallest absolute component
of `n_s`, with fixed tie order `x`, then `y`, then `z`.

The selected planar wall satisfies `u=0`. `tau*n_s` has independently nonzero normal and tangential parts. The fixed expectations are

```text
F_operator_pressure =  p0 * A_background
F_operator_viscous  = -tau * A_background
F_surface_pressure  = -p0 * n_q * A_surface
F_surface_viscous   =  tau * n_q * A_surface
```

and totals are componentwise sums. This is the Ghost-Cell full-background-cell residual measure, not a cut-cell pointwise equality. The signed projected measure

```text
A_signed = -dot(A_background, n_s) > 0
A_background = (A_background + A_signed*n_s) - A_signed*n_s
```

fixes the background-face orientation relative to the true wall normal; it is
not silently replaced by `-n_s*A_surface`. The test requires both force
measures nonzero and `F_operator != F_surface` by a pre-result,
epsilon-scaled separation check.

The constant/linear LFP cancellation used by this RED is explicit. Constant
pressure gives zero affine remainder at the face, wall, diagonal, and all
neighbour samples, so every transformed-minus-background pressure defect is
zero. For `u(x)=G(x-x_wall)`, origin-constrained quadratic reconstruction is
linear-exact and gives the same constant `G` at the wall, diagonal, and every
neighbour sample. The frozen LFP constant-preservation invariant gives
`source_transformed=source_background=0` and
`sum(diagonal + six neighbours + source)_transformed =
 sum(diagonal + six neighbours + source)_background`; therefore all
constant-stress viscous defects cancel and the wall row is exactly
`-tau*A_background`. The fixture preflight separately records the complete
affine-row sum, the zero source, linear-gradient reconstruction exactness, and
constant-stress defect cancellation before reading force results. If any
exactness/cancellation observation fails, that is a non-semantic RED and
product repair stops for RED-S2 adjudication.

All algebraic comparisons use one frozen componentwise bound

```text
2^18 * epsilon * max(1, max_abs(actual), max_abs(expected))
```

chosen before any candidate result. It is an exact-path arithmetic bound, not a convergence threshold.

## Required observation seams

### Wall force trace

Under `HUNDUN_IMMERSED_ENABLE_TEST_ACCESS`, add to `wall_force_detail.hpp`:

```cpp
struct WallForcePointSnapshot final {
  ImmersedLinkId link{};
  TriangleId triangle{};
  std::uint32_t point_index{};
  runtime::Real3 position_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double weight_m2{};
  runtime::Real3 pressure_force_N{};
  runtime::Real3 viscous_force_N{};
  runtime::Real3 total_force_N{};
};

struct WallForceTrace final {
  WallForceSample reduced;
  std::vector<WallForcePointSnapshot> local_points;
};

WallForceTrace trace_wall_force_for_test(
    const WallForceIntegrator &,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell);
```

The implementation sets a scoped, test-only local snapshot sink and calls `WallForceIntegrator::integrate`. The existing `integrate_local` loop appends the already-computed point geometry and force values to that sink immediately after the real pressure/viscous/total calculation. It must not recompute traction, bypass canonical ordering, filter the plan, or alter the reduction. The trace's sum over **all** local snapshots, reduced with MPI in the test, must equal the full-surface `trace.reduced.surface_traction`. The selected-link product surface contribution is instead the MPI sum of traced point forces filtered by `snapshot.link`; that observed subset is compared with the independently rebuilt selected-triangle analytic sum. The full-surface reduced value is never compared with the one-link analytic value.

### Per-row budget trace

Add this field to the existing `BoundaryRowEvaluationSnapshot`:

```cpp
BoundaryResidualPartsSnapshot budget_reaction_delta;
```

In `accumulate_momentum`, under the existing finite-volume test-access guard,
snapshot `local_reaction` immediately before a wall row and subtract it from
the accumulator immediately after that row. Store that actual pressure and
viscous delta. The selected-row budget assertions and the input to the real
final report assembler use this observed delta. Summing all per-row deltas
must reproduce `ImmersedOperatorAdapter::report()`; this is independent of the
raw-row snapshots and kills a report/raw sign reversal.

### Final report assembler

Extract the current `collect_final_force` assignment block into one anonymous-namespace pure function in `stage3_flow.cpp`:

```cpp
ForceAttemptReport assemble_candidate_force_report_from_budget(
    const immersed::ForceComponents &budget_reaction,
    const immersed::ForceComponents &surface_traction);
```

Before RED acceptance this helper must preserve the candidate exactly: copy the budget input to `operator_reaction`, copy the surface input, and add them into consistency. `collect_final_force` must call this helper after the existing reductions and wall integration. Expose it through a deliberately different member name:

```cpp
static ForceAttemptReport assemble_force_attempt_report_from_budget(
    const immersed::ForceComponents &budget_reaction,
    const immersed::ForceComponents &surface_traction);
```

on `Stage3FlowTestAccess`. This seam is the actual report path; it is not a second report algorithm.

## Required named observations

The rank-zero output prints every boolean and relevant vector before one aggregate `HUNDUN_CHECK`. At minimum it prints:

```text
fixture.selected_one_link_row
fixture.oblique_normal
fixture.surface_points_present
fixture.planar_link_patch
fixture.domain_normal_equals_triangle_normal
fixture.wall_intercept_on_triangle_plane
fixture.quadrature_geometry_equals_fixture
fixture.lfp_full_affine_row_sum_conserves_constant
fixture.lfp_source_is_zero
fixture.linear_gradient_reconstruction_exact
fixture.constant_stress_defects_cancel
fixture.operator_measure_nonzero
fixture.surface_measure_nonzero
fixture.operator_surface_measures_distinct
trace.reduced_equals_point_sum

pressure.raw_equals_analytic_operator
pressure.row_update_uses_positive_raw_wall_term
pressure.surface_equals_analytic_surface
pressure.budget_equals_negative_raw
pressure.physical_report_equals_analytic_operator
pressure.consistency_equals_operator_minus_surface
pressure.isolated_budget_closure

viscous.normal_traction_nonzero
viscous.tangential_traction_nonzero
viscous.raw_equals_analytic_operator
viscous.row_update_uses_positive_raw_wall_term
viscous.surface_equals_analytic_surface
viscous.budget_equals_negative_raw
viscous.physical_report_equals_analytic_operator
viscous.consistency_equals_operator_minus_surface
viscous.isolated_budget_closure

combined.operator_pressure_viscous_nonparallel
combined.surface_pressure_viscous_nonparallel
combined.raw_equals_analytic_operator
combined.surface_equals_analytic_surface
combined.physical_report_equals_analytic_operator
combined.consistency_equals_operator_minus_surface

aggregate.raw_wall_sum_nonzero
aggregate.adapter_budget_equals_negative_raw_sum
perturb.operator_only_isolated
perturb.surface_only_isolated
perturb.assembly_depends_on_operator
perturb.assembly_depends_on_surface
components.pressure_viscous_separate
```

The current candidate is expected to make only the physical-report and physical-consistency observations false. If any geometry, raw operator, surface, trace, budget, traction, separation, or component observation is false, stop: do not repair the product and do not reinterpret the failure as the expected RED.

For every selected case, `budget_equals_negative_raw` requires the actual
per-row accumulator delta to satisfy
`budget_reaction_delta == -F_operator == -raw_wall_contribution`; it is not
constructed by negating the raw snapshot in the test. The row-update
observation independently checks

```text
residual_after = (residual_before + background) - removed_background
                 + raw_pressure + raw_viscous
```

so changing the sign with which a raw wall term enters the fluid row cannot be
hidden by an unchanged snapshot.

The aggregate nonzero probe uses the fixed affine pressure
`p(x)=2.75+{0.31,-0.17,0.23} dot (x-x_wall)` in an additional call, while the
three signed cases retain the fixture above. It requires both the sum of
actual per-row budget deltas and the sum of raw row snapshots to be nonzero,
then compares the former with the adapter report and the latter with the
opposite-sign budget. It exists solely to make the aggregate
`accumulate_momentum` sign relation mutation-sensitive; it is not a
surface-force oracle.

The operator-only and surface-only probes call the real report assembler with independent copies of the already verified analytic vectors. The first changes only the budget/operator input by `{0.125, -0.25, 0.375}` times the case force scale; the second changes only surface input by `{-0.30, 0.20, 0.10}` times that scale. Each probe requires only the corresponding output dependency to change.

---

### Task 1: Build the test first and add only the three observation seams

**Files:**
- Modify: `CMakeLists.txt`
- Create: `tests/mpi/test_task11_signed_force_authority.cpp`
- Modify: `immersed/src/wall_force_detail.hpp`
- Modify: `immersed/src/wall_force.cpp`
- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Modify: `finite_volume/src/immersed_operator.cpp`
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Modify: `flow/src/stage3_flow.cpp`

**Interfaces:**
- Consumes: current product `ImmersedDomain::links()`, topology/geometry, operator snapshots, WallForceIntegrator, and ForceAttemptReport.
- Produces: the guarded `WallForceTrace` seam, actual per-row budget delta, guarded report-assembly seam, `fixture_preflight_{1,2,4}_rank`, and `test_task11_signed_force_authority_{1,2,4}_rank` CTests.

- [ ] **Step 1: Write the test against the wished-for seams**

Implement the fixed fixture, analytic stress, deterministic link selection, independent geometry reconstruction, all named observations, and one final aggregate assertion. Register one executable linked with `hundun_flow`, `stage3_stl_fixture`, `hundun_options`, `hundun_warnings`, and `MPI::MPI_CXX`; add the root include directory and `HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1` to that target.

- [ ] **Step 2: Prove the seam is initially missing**

Configure and build only the new target. Expected at this scaffold checkpoint: compile or link failure naming `trace_wall_force_for_test`, `budget_reaction_delta`, and/or `assemble_force_attempt_report_from_budget`. Capture it as TDD scaffold evidence, not as the scientific RED.

- [ ] **Step 3: Add the WallForceIntegrator observation without a second algorithm**

Append snapshots inside the existing point loop and implement the scoped trace wrapper by calling the real integrator. Capture the finite-volume per-row accumulator delta before/after the actual updates. Do not move or duplicate the pressure/viscous formulas and do not change their signs.

- [ ] **Step 4: Extract and expose the current report assembly exactly**

Move only the current assignment/addition block into the distinctly named pure helper, call it from `collect_final_force`, and forward it through `Stage3FlowTestAccess` without an unqualified same-name recursive call. Run the existing 1-rank wall-force, immersed-operator, and immersed-PISO regressions before the new scientific test.

- [ ] **Step 5: Run the force-free fixture preflight on 1/2/4 ranks**

All ranks execute the frozen gather/reduction/broadcast protocol and print the
same selected link/row/triangle geometry. They do not call a force integrator
or report assembler. If any rank count lacks an eligible link or the selected
geometry differs beyond the frozen arithmetic bound, stop and return a
fixture-plan failure to the main agent.

- [ ] **Step 6: Build successfully and run the 1-rank behavioral RED**

Run the new 1-rank CTest verbosely. Required: the executable reaches the final aggregate assertion, prints the full vector, all non-semantic observations are true, and only the declared physical-report/consistency observations are false.

- [ ] **Step 7: Reproduce the same semantic RED on 2 and 4 ranks**

Run only the new target on 2 and 4 ranks, one at a time. The selected global link, true/false observation vector, and failing assertion category must match exactly; analytic and observed floating-point vectors must match within the already frozen componentwise bound, not bitwise, because reduction grouping may change rounding.

- [ ] **Step 8: Prove the test guards disappear from a tests-off build**

Configure `HUNDUN_BUILD_TESTS=OFF` in a separate focused tree and build the
affected `hundun_immersed`, `hundun_finite_volume`, and `hundun_flow` libraries
with `-j2`. No test-access declaration or symbol may enter that build.

- [ ] **Step 9: Record the exact tests-only RED package**

Save source hashes, full allowed-file diff, build/log/binary hashes, commands, environment, exit statuses, selected geometry, analytic vectors, observed vectors, and the unchanged whole-tree/index invariants under `.superpowers/task-11-red-s1-*`.

- [ ] **Step 10: Self-review and freeze the exact unstaged delta**

Review the eight-file boundary, run `git diff --check`, verify no numerical
formula changed, and generate one exact patch against the main agent's frozen
pre-task file copies. Record its SHA-256. Do not run `git add` or `git commit`;
verify the index remains empty. The intended future consolidation subject is:

```text
test: expose Task 11 signed force authority RED
```

### Task 1A: Apply coordinator fixture amendment A1 and resume Task 1

**Files:**
- Modify: `tests/mpi/test_task11_signed_force_authority.cpp`
- Update evidence only under the existing RED-S1 SDD workspace

**Interfaces:**
- Consumes: the preserved unsplit 2-rank fixture failure, amendment A1, and
  the existing Task 1 test/seams.
- Produces: a partition-conforming version of the same analytic surface,
  explicit tessellation/triangle-owner observations, refreshed exact delta,
  and resumed Task 1 Steps 5--10.

- [ ] **Step 1: Preserve the unsplit failure and freeze the repair prestate**

Record the current test SHA-256, Task 1 exact delta SHA-256, binary SHA-256,
1-rank PASS, and 2-rank reach-four collective failure. Do not overwrite the
original Task 1 report.

- [ ] **Step 2: Implement only the deterministic plane splitter and oracles**

Keep the original twice-refined rotated cube as the intermediate geometry.
Add the exact `x=0.5`, then `y=0.5` split and the pre-product surface
invariants above. Add global post-plan checks for three points and one owner
per emitted triangle. Do not change any selector, analytic stress, expected
force, tolerance, process grid, product source, or observation seam.

- [ ] **Step 3: Run the force-free 1/2/4 preflights serially**

Require all three to pass and print the same selected link/row geometry within
the frozen bound. Stop on the first failure; do not run a behavioral mode
after a failed preflight.

- [ ] **Step 4: Resume Task 1 Steps 6--10**

Run the 1/2/4 behavioral REDs serially, then the tests-off guard build,
`git diff --check`, exact allowed-file delta refresh, and evidence report.
Only the declared physical-report/consistency observations may be false.

### Task 1B: Repair the loader's legal shared-vertex false positive

**Files:**
- Modify: `CMakeLists.txt` only to enable existing immersed test access for
  `test_immersed_surface`
- Modify: `immersed/src/immersed_surface.cpp`
- Modify: `immersed/src/immersed_test_access.hpp`
- Modify: `tests/unit/test_immersed_surface.cpp`
- Update evidence only under the existing RED-S1 SDD workspace

**Interfaces:**
- Consumes: the blocked A1 source/delta, exact 1-rank loader RED, restored
  product prestate, pair/branch diagnostics, and amendment A2.
- Produces: the minimal conservative tangent-cone predicate, loader regressions,
  resumed A1 1/2/4 preflights, and a refreshed Task 1 evidence package.

- [ ] **Step 1: Freeze the RED and clean diagnostic restoration**

Record the A1 report/delta/log/binary hashes, triangle-pair proof, diagnostic
log hashes, and the restored product/test SHA-256 values. Rebuild once from
the restored sources and reproduce the generic 1-rank loader RED before any
accepted product edit.

- [ ] **Step 2: Write focused tangent-cone RED probes**

Through the existing tests-only `ImmersedTestAccess`, call the real predicate
with exact shared-vertex geometry. Require: the A1 outward direction is
allowed; an inward direction is forbidden even when the corresponding convex
intersection parameter width is below the rejected `256*epsilon` window; a
second finite-overlap case above that window is forbidden; tangent-inward is
forbidden; and tangent-outward is allowed. Include a dropped-y case that would
fail if projection parity were omitted, plus a coordinate-near endpoint with a
different welded vertex ID that must remain forbidden. The seam contains no
duplicate geometry algorithm and disappears from tests-off builds.

- [ ] **Step 3: Implement the bounded tangent-cone proof**

Change only the legacy coplanar segment/shared-vertex subpath specified by A2.
No multiplier tuning, physical tolerance change, owner/Halo change, or
fixture edit is allowed.

- [ ] **Step 4: Run focused loader GREEN controls**

Build `test_immersed_surface` and the Task 11 executable with `-j2`. Run the
complete existing immersed-surface unit executable, requiring both true
self-intersection controls to remain rejected and all positive surfaces to
pass. Then run only the affected immersed-surface MPI tests serially on their
registered 1/2/4 ranks.

- [ ] **Step 5: Resume the A1 hard gates**

Run `fixture_preflight_1_rank`, then 2, then 4, serially. Stop on the first
failure. If all pass, resume Task 1 behavioral RED 1/2/4 and the tests-off
proof without any additional fixture or product change.

- [ ] **Step 6: Freeze and review the prerequisite delta**

Generate exact prestate-to-GREEN product and updated Task 1 deltas, hashes,
commands, environments, statuses, named observations, and a restoration
mutation record. The index remains empty until coordinator review.

### Coordinator blocker diagnosis A3 and authority amendment A4

Task 1B passed its focused loader and immersed-surface controls.  The resumed
Task 11 one-rank preflight passed, but the two-rank preflight still failed with
`wall quadrature donor exceeds owner Halo reach`.  A disposable, fully
restored Task 1C diagnostic proved that triangle 0 has association owners
`[0,0,1]`: the old tuple rule chooses rank 0, while the complete point and
pressure-authority donor union reaches logical x cell 10.  Rank 0's reach-four
half-open interval ends at 10 and is infeasible; rank 1's interval contains
the complete union.  The feasible-owner mask is `0x2`.

This is a product capability blocker, not a fixture defect and not authority
to discard donors, enlarge reach, or change reconstruction.  The following
documents supersede this plan's earlier lowest-associated-owner sentence and
supplement its Halo-provisioning contract for Task 1D:

```text
docs/superpowers/specs/2026-08-06-hundun-flow-wall-quadrature-feasible-owner-addendum.md
docs/superpowers/plans/2026-08-06-hundun-flow-task11-wall-quadrature-feasible-owner.md
.superpowers/task-11-red-s1-feasible-owner-a4-authority-2026-08-06.md
.superpowers/task-11-red-s1-halo-discriminator-adjudication-a5-2026-08-06.md
```

Task 1D must finish before Task 1 resumes.  It chooses the numerically lowest
globally feasible execution owner from the complete three-point point/pressure
donor union, gives `WallQuadraturePlan` its own fixed required reach four,
keeps link pressure authority independent, and provisions integrated callers
with `max(ghost_reach, wall_reach)`.  Raw execution owner remains outside the
decomposition-independent wall-plan fingerprint and Checkpoint v2 remains a
non-change surface.  No force result may influence ownership.

For Task 1D only, the file boundary in its implementation plan replaces this
plan's earlier eight-file behavioral-RED boundary.  This narrowly authorizes
the owner/Halo product and focused-test files listed there, including the
previously regression-only wall-force, wall-quadrature, and transaction tests.
After Task 1D acceptance the original RED-S1 boundary resumes; no other source
or numerical semantic surface is reopened.

The original Task 1D pre-product hard stop requested a real
`0 < ghost_reach <= 3` fixture.  A read-only survey of every existing clean
integrated candidate, followed by an exact 2-rank GDB observation of the A1/A2
Task 11 fixture, found `ghost_reach == 4` throughout.  No lower-reach fixture
can be obtained without tuning geometry or selection, so that natural
discriminator is unavailable and may not be fabricated.

The coordinator therefore replaces only that unavailable test precondition,
not the wall reach or fail-closed contract.  Before product code changes, the
replacement REDs must prove that a one-rank standalone wall-force case whose
donors are all owned currently accepts three-layer fields, and that a
three-layer transaction reaches the old generic ghost-plan error rather than
the new wall-specific error.  After implementation, the same cases must
reject width three at the wall-contract boundary, accept width four, and emit
the new collective wall-specific construction classification before the old
ghost check.  Task 11 itself is frozen at actual ghost reach four.  This is a
mutation-sensitive product-path discriminator and does not claim that a
natural lower-reach GhostStencilPlan currently exists.

### Task 2: Independent RED requirements and code-quality review

**Files:**
- Create: `.superpowers/task-11-red-s1-requirements-review-2026-08-05.md`
- Create: `.superpowers/task-11-red-s1-code-quality-review-2026-08-05.md`
- Create: `.superpowers/task-11-red-s1-coordinator-acceptance-2026-08-05.md`

**Interfaces:**
- Consumes: Task 1 brief/report, exact diff package, RED logs, hashes, and current worktree invariants.
- Produces: an explicit ACCEPT/REJECT verdict for the RED only; no product semantic verdict is implied.

- [ ] **Step 1: Dispatch an independent requirements/science reviewer**

Require it to check independent operator geometry, analytic stress, surface oracle independence, non-degeneracy, complete observation vector, rank invariance, and expected failure classification.

- [ ] **Step 2: Dispatch an independent code-quality/evidence reviewer**

Require it to check compile guards, absence of duplicate algorithms, real-path forwarding, mutation observability, test cleanup, MPI ownership/reductions, file boundary, and evidence integrity.

- [ ] **Step 3: Resolve every Critical or Important finding**

Return findings to the original implementer for at most three fix rounds, each followed by a scoped re-review. Any fix that changes a product numerical sign/formula is rejected and moved to Task 3 after RED acceptance.

- [ ] **Step 4: Main-agent full diff and evidence acceptance**

The main agent recomputes all hashes, reads the complete before/after task
delta plus the surrounding inherited file contents, reruns fresh 1/2/4 RED
tests, and writes the RED acceptance record. Acceptance means “the executable
correctly exposes the current contradiction,” not “Task 11 is green.”

### Task 3: Freeze the semantic amendment and GREEN implementation packet

**Files:**
- Create: `docs/superpowers/specs/2026-08-05-hundun-flow-task11-signed-force-semantic-amendment.md`
- Create: `.superpowers/task-11-red-s1-green-implementation-packet-2026-08-05.md`

**Interfaces:**
- Consumes: accepted RED observations.
- Produces: one exact semantic contract and a newly hashed GREEN allowed-file list.

- [ ] **Step 1: Apply the stop/branch rule**

If the selected raw pressure/viscous row does not equal the independently
derived one-link `F_operator`, if its actual per-row budget delta does not
equal `-F_operator`, or if the selected traced surface subset does not equal
the one-link `F_surface`, stop and revise the derivation before any product
edit. Continue only if, separately, the nonzero aggregate probe proves that
the sum of every actual per-row budget delta equals the full adapter report
and is the negative of the full raw-row sum, and the shared final report
assembler copies the supplied selected-case budget input into the field it
currently presents as the physical operator result.

- [ ] **Step 2: Freeze one unambiguous vocabulary**

The amendment must distinguish: raw wall residual/physical discrete operator force; opposite-sign budget reaction; true-surface physical force; and physical consistency. It must state that two same-orientation physical forces are subtracted and that the budget reaction is not exposed under a physical-force name or compatibility alias.

- [ ] **Step 3: Re-run `codegraphf impact` plus exact `rg`**

Inventory every caller/test affected by the required rename or separation. The GREEN packet freezes exact source files and hashes only after this inventory; no provisional source host silently becomes authorized.

- [ ] **Step 4: Freeze the minimal GREEN behavior**

The implementation must preserve the operator residual and budget closure, convert the final physical operator report exactly once, compute consistency as physical operator minus physical surface, and keep pressure/viscous/total components separate. No LFP reconstruction or PISO formula is changed in RED-S1 GREEN.

### Task 4: Implement the semantic GREEN and verify focused regressions

**Files:**
- Modify only the exact GREEN list frozen by Task 3.

**Interfaces:**
- Consumes: semantic amendment and GREEN packet.
- Produces: a green signed-force executable without changing row/surface numerical authorities.

- [ ] **Step 1: Run the accepted RED at exact pre-fix HEAD**

Confirm the full expected failing vector once more and record binary/source hashes.

- [ ] **Step 2: Make the smallest naming/orientation repair**

Preserve adapter raw rows and budget closure. Separate/rename the budget quantity as mandated, form the physical operator result from the raw orientation once, and subtract the surface physical result.

- [ ] **Step 3: Run the new 1-rank test GREEN**

Every named observation must be true; no observation may be removed or weakened.

- [ ] **Step 4: Run focused Debug regressions**

Run the new 1/2/4 tests plus existing wall-force, immersed-operator, immersed-PISO, and immersed-transaction tests on 1/2/4 ranks, scheduled so only one MPI test runs at a time.

- [ ] **Step 5: Run focused Release regressions**

Build only affected targets with `-j2`; run the same signed-force 1/2/4 tests and the task-focused existing subset. No 96^3 or full numerical matrix belongs here.

- [ ] **Step 6: Review and freeze the GREEN delta**

After independent requirements and code-quality reviews, freeze the exact
GREEN before/after delta and its hashes without staging inherited dirty files.
The index remains empty until the main agent has reviewed the complete Task 11
candidate and explicitly chooses a consolidation boundary.

### Task 5: Kill the normal/raw/report mutations and S1-M1 through S1-M6

**Files:**
- Create only disposable candidate worktrees/directories and `.superpowers/task-11-red-s1-mutation-*` evidence.

**Interfaces:**
- Consumes: exact immutable GREEN source snapshot and manifest.
- Produces: eight one-at-a-time mutant records and final RED-S1 acceptance.

- [ ] **Step 1: Create one immutable exact GREEN seed**

Record HEAD, every allowed-file hash, the exact RED/GREEN deltas, whole tracked
patch hash, untracked/source manifests, compiler/MPI environment, and focused
baseline binary hash. Copy the complete required candidate files into a
read-only seed directory; the seed need not be a Git commit. Do not
mutate/revert the main worktree.

- [ ] **Step 2: Run S1-M0 normal reversal**

Inside the actual WallForceIntegrator point calculation only, negate the local
normal supplied to both pressure and viscous traction without changing the
WallQuadraturePlan or fixture geometry. Required designated failures:
`pressure.surface_equals_analytic_surface`,
and `viscous.surface_equals_analytic_surface`. The WallQuadraturePlan geometry
checks remain true, proving the failure comes from the downstream product
normal use. Because the oracle normal/area comes from the test-owned STL
triangles, it cannot follow this product mutation.

- [ ] **Step 3: Run S1-M1a raw-row insertion reversal**

Reverse the sign with which `evaluated.residual.pressure` and
`evaluated.residual.viscous` enter the actual boundary-row candidate while
leaving the observed evaluated functional unchanged. Required designated
failures: `pressure.row_update_uses_positive_raw_wall_term` and
`viscous.row_update_uses_positive_raw_wall_term`.

- [ ] **Step 4: Run S1-M1b budget/report reversal**

Reverse both pressure and viscous budget accumulation signs in
`ImmersedOperatorAdapter::accumulate_momentum`. Required designated failures:
the selected per-row `budget_equals_negative_raw` observations and
`aggregate.adapter_budget_equals_negative_raw_sum`.

- [ ] **Step 5: Run S1-M2**

Reverse `-pi*n*A` in the actual WallForceIntegrator point calculation. Required designated failure: `pressure.surface_equals_analytic_surface`.

- [ ] **Step 6: Run S1-M3**

Reverse `tau*n*A` in the actual WallForceIntegrator point calculation. Required designated failure: `viscous.surface_equals_analytic_surface`.

- [ ] **Step 7: Run S1-M4**

Replace physical consistency subtraction with addition in the shared final report assembler. Required designated failures: the pressure, viscous, and combined consistency observations.

- [ ] **Step 8: Run S1-M5**

Inside the shared pure assembler used by both `collect_final_force` and the
test seam, copy surface traction into the physical operator result. A mutation
only in the caller is invalid because it would bypass the seam. Required
designated failures: analytic physical-operator comparisons and
`perturb.operator_only_isolated`/assembly operator sensitivity.

- [ ] **Step 9: Run S1-M6**

Swap pressure and viscous report components before totals. Required designated failures: pressure-only, viscous-only, and component-separation observations.

- [ ] **Step 10: Re-run the unmutated seed and final reviews**

Require all eight mutants to fail through their designated observation,
remove only disposable mutant worktrees/candidates, rerun the exact unmutated
Debug/Release 1/2/4 tests, and obtain independent evidence/science plus
code-quality PASS verdicts.

- [ ] **Step 11: Main-agent RED-S1 acceptance**

Write one aggregate record containing all hashes, observations, mutation vectors, review verdicts, and remaining Task 11 blockers. RED-S1 acceptance authorizes RED-S2; it does not accept Task 11 or the formal force matrix.

## Build and scheduling contract

Use fresh focused build trees and at most two compiler jobs:

```bash
CC=/home/wyf/.local/opt/hundun-toolchain/clang/bin/clang
CXX=/home/wyf/.local/opt/hundun-toolchain/clang/bin/clang++
CXXFLAGS='-stdlib=libc++ -DOMPI_SKIP_MPICXX'
LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu
cmake --build build/task11-red-s1-debug \
  --target test_task11_signed_force_authority -j2
```

Use `build/task11-red-s1-release` for Release and
`build/task11-red-s1-tests-off` for `HUNDUN_BUILD_TESTS=OFF`; do not reuse a
mutation candidate build directory.

Non-MPI configure/build and read-only reviews may overlap. The new 1/2/4 MPI cases, existing MPI regressions, and every mutation executable run are serialized. No large MMS matrix runs until RED-S1 through RED-S3 are green.

## Completion boundary

RED-S1 is complete only when the accepted RED was observed before repair, the semantic amendment was frozen, the focused GREEN passes Debug/Release on 1/2/4 ranks, S1-M0, both S1-M1 submutants, and S1-M2 through S1-M6 die through named observations, and the main agent accepts the full diff/evidence. The next work is RED-S2 shared boundary-row authority; formal 24/48 and 24/48/96 convergence remains prohibited until S2 and S3 close.

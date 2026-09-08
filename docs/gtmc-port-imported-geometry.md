# GTMC imported-marker Cartesian IBM port

## Purpose and root cause

The GTMC cold case cannot preserve its frozen fluid domain by recompiling the
available STL with the current HUNDUN surface path.  This is not an open-mesh
repair problem:

- `gtmc13_si.STL` has 22,848 non-degenerate triangles, 11,410 exact vertices,
  34,272 undirected edges with incidence two, no incidence-one or incidence
  greater-than-two edge, and three edge-connected closed components.
- The three component signed volumes are
  `9.01265160476906e-4`, `3.797018838792935e-5`, and
  `4.693052039073086e-6 m^3`; all have the same orientation.
- HUNDUN's float32-effective y-ray scan has 299 odd rays.  At the shared plane
  `y=-0.0135`, two positive-normal events from different closed components are
  collapsed into one event, changing a raw six-event ray into five events.
- A component-local, nonzero-winding union is deterministic and agrees across
  x/y/z scans, but still differs from the frozen marker in 46,607 cells
  (0.8167437% of the grid).  The best of all 27 one-cell shifts is the zero
  shift, so this is not a coordinate-offset error.
- The frozen marker has 218,175 directed fluid-to-solid links.  Its 250 fuel
  owner links are all present.  Reconstructed STL union geometry has 207,018
  links and preserves only 104 of those 250 fuel links.

Therefore the frozen `global_marker.bin`, not a guessed STL repair or a new
scan, is the geometry authority for this compatibility mode.  Current COAST
source was consulted only for semantic evidence (normal mode 1 forms a
six-neighbour marker direction and the default wall distance is the minimum
half-centre distance).  It is not treated as proof that the current source can
reproduce the historical frozen marker.

## Imported authority contract

`ImportedIbmCompiler::compile` accepts an immutable, globally replicated
marker with these fail-closed requirements:

- x-fastest Cartesian ordering, with total size exactly equal to
  `global_cells.x * global_cells.y * global_cells.z`;
- every byte is exactly `RegionFlag::solid` (`0`) or `RegionFlag::fluid`
  (`1`);
- `marker_source` is a nonzero fingerprint supplied by the case/file layer;
- geometry, process-grid, marker identity, stencil limits, and resource limits
  agree collectively on every rank;
- publication is atomic: topology, boundary stencil, and quadrature outputs
  are moved into caller storage only after every collective and local compiler
  stage succeeds.

File I/O, JSON parsing, and case binding deliberately remain outside this
compiler.  That keeps the mesh module independent of GTMC filenames and case
schema.

## HUNDUN program changes

### `versions/v0.4/include/hundun/v04_ibm.hpp`

Symbols changed:

- added the public `ImportedIbmCompiler` declaration and its single collective
  `compile` entry point;
- added narrowly scoped `friend class ImportedIbmCompiler` access to
  `QuadraticStencilPlan`, `IbmInterfaceMetricPlan`, `EBTopology`,
  `BoundaryStencilPlan`, and `SurfaceQuadraturePlan`.

Reason: the compiler must publish the same sealed immutable plan types used by
the existing solver without adding mutable setters or fabricating an
`ImmersedSurfacePlan`.  No STL compiler API or existing public call signature
was changed.

### `versions/v0.4/src/mesh_ibm_import.cpp`

Symbol added: `ImportedIbmCompiler::compile`.

The implementation:

1. validates the replicated marker and collective compile contract;
2. copies the exact owned marker bytes and fills the width-four region halo
   directly from the same global authority;
3. enumerates every owned fluid cell's six neighbours and creates one stable
   link key `6 * fluid_global_cell + direction` for each in-domain solid
   neighbour;
4. places the wall point on the exact shared Cartesian face, uses the exact
   transverse face area, and assigns the solid-to-fluid unit normal `-delta`;
5. builds `IbmInterfaceMetricPlan` directly from that planar face.  The area,
   area vector, first moment, normal first moment, and normal second moment are
   exact Cartesian planar moments;
6. selects only positive-normal fluid donors from the immutable marker.  Before
   the 32-donor cap is filled by distance, it retains the available tangential
   quadrants and up to three nearest distinct coordinate bands in each of the
   normal and two tangent directions;
7. delegates ordinary fitting, rank, coverage, conditioning, and adaptive-order
   decisions to the existing `QuadraticStencilCompiler`, then builds the same
   four boundary rows (Dirichlet ghost value, zero-normal ghost value, wall
   value, and wall directional derivative);
8. creates one centre quadrature point per Cartesian link.  For a planar face
   with constant axis normal this point exactly represents every moment stored
   by the current metric plan.  Its reconstruction reuses the already
   certified wall rows; no nearest-cell copy path exists.

Imported links use `global_link` as their synthetic surface-element identity,
so `ImmersedLink::triangle`, `IbmInterfaceLinkMetric::source_triangle`, and
`SurfaceQuadraturePoint::triangle` remain bijective without claiming an STL
triangle.  `fluid_side` is stored as the canonical legacy value `outside`, but
normal authority is the explicit marker transition, not STL winding.

Status details reserved by this translation unit are 13801 through 13807 for
input, collective, region, link, donor, metric, and allocation failures.

#### Explicit one-normal-band closure

The frozen marker contains 285 links for which every reachable positive-normal
fluid donor lies in one normal band.  Their direction counts in
`[x-, x+, y-, y+, z-, z+]` order are `[62, 68, 0, 26, 67, 62]`.  Widening the
configured reach cannot create a second band because the next layer is solid.
The regular adaptive compiler therefore rejects these links with exactly
`invalid_plan/1304` (`QuadraticCoverage`).

Only when all three conditions hold -- imported-marker geometry,
`adaptive_order`, and a one-normal-band request independently rejected as
detail 1304 -- the importer publishes an explicit one-dimensional linear
two-point closure.  For fluid-centre distance `d_f` and solid-centre distance
`d_s` from the shared Cartesian face, its four rows are:

```text
u_ghost,Dirichlet = -(d_s/d_f) u_fluid + (1+d_s/d_f) u_wall
u_ghost,zero-n     = u_fluid - (d_f+d_s) (du/dn)_wall
u_wall,reconstruct = u_fluid
(du/dn)_wall       = (u_fluid-u_wall)/d_f
```

Every such group records `order=linear`, `donor_count=1`,
`normal_band_count=1`, and `fallback_reason=quadratic_coverage`.  Strict policy
still returns detail 1304, and every other failure remains fail closed.  This
is a face-normal closure, not a nearest-cell-copy geometry path: the link's
axis normal, wall point, area, and metric moments remain unchanged.

#### 1305 donor-selection root cause

After the 285 genuine one-band links were handled, the real case exposed 698
detail-1305 links, with direction counts `[273, 216, 0, 0, 0, 209]`.  All had
three selected normal bands and 32 selected donors.  The marker itself was not
lower-dimensional.  For example, link `(12,238,13) -> (11,238,13)` has 217
reachable positive-normal fluid candidates with x/y/z cardinalities `5/5/9`
and y range `234..238`; the old distance truncation retained `5/1/8`, all at
y=238.  Thus one tangent column of the linear basis was identically zero and
the apparent rank failure was manufactured by selection.  The symmetric
x-positive example had 176 available candidates (`4/5/9`) but retained
`4/1/8`.

The fix is the coordinate-band retention in step 6 above.  It exposes existing
valid marker cells to the unchanged reconstruction compiler; it does not add a
rank fallback, change the 32-donor limit, widen reach, or alter the marker.

### `versions/v0.4/CMakeLists.txt`

Change: added `src/mesh_ibm_import.cpp` to `hundun_v04_core_sources`.

Reason: compile and link the new explicit geometry mode into both the product
core and the test core.  Source ordering and all existing targets remain
unchanged otherwise.

### `versions/v0.4/tests/mpi/mesh_ibm_import_mpi_test.cpp`

New focused MPI fixture covers:

- a 16-cubed half-domain plane marker with exactly 256 links;
- a 16-cubed fluid-octant/solid-complement corner marker with exactly 192
  links, including three directional links owned by cell `(8,8,8)`;
- an internal one-cell-thick fluid sheet with exactly 512 links, proving strict
  detail 1304, adaptive two-point row coefficients, and explicit linear /
  coverage-fallback audit counts;
- a stretched-grid half-space edge with exactly 272 links.  Before the donor
  selection fix it deterministically failed detail 1305; the passing test
  requires the target reconstruction to retain more than one tangential y
  band;
- bit-identical owned region publication;
- topology/metric/boundary/quadrature one-to-one identities;
- axis normals, Cartesian areas, and global area conservation;
- reuse of the boundary wall value and wall-normal-gradient rows by
  quadrature;
- collective, atomic rejection when rank zero supplies a marker byte `2`.

### `versions/v0.4/tests/CMakeLists.txt`

Change: added `v04_mesh_ibm_import_mpi_test` and registered 1-rank and 2-rank
CTest cases.

Reason: exercise both serial ownership and a decomposed marker without adding
case-specific data files.

## Verification receipts

Core build:

```text
cmake --build /tmp/hundun-gtmc.Itnz3O/build-ibm-source-face \
  --target hundun_v04_core -j2
[100%] Built target hundun_v04_core
```

Focused imported-marker tests:

```text
ctest --output-on-failure -R '^v04_mesh_ibm_import_mpi_[12]$'
v04_mesh_ibm_import_mpi_1  Passed  0.83 sec
v04_mesh_ibm_import_mpi_2  Passed  0.62 sec
100% tests passed, 0 tests failed out of 2
```

Read-only GTMC asset checks independently established the marker dimensions
`153 x 247 x 151`, the 20 j=0 air labels, the paired fuel labels at j=80/j=79,
all 250 fuel owners, and total inlet area `50.625 mm^2`.

Real-case marker-only dry-plan after the two diagnosed reconstruction fixes:

```text
mpiexec -n 16 build/gtmc-tests/versions/v0.4/hundun \
  validate /tmp/hundun-gtmc.Itnz3O/case-exact-marker-probe --dry-plan
VALID case=9228097579569678706 product=14320855744224918633 \
  global=153x247x151 local=77x62x76 fields=57 arena_doubles=68334704 \
  stages=15 correctors=2 pressure_rtol=1e-13 pressure_atol=1e-13 \
  pressure_maxit=400 pressure_restart=12 sealed=1
```

A one-shot read-only `CaseValidationReport::summary` probe of the same 16-rank
compile reported the global boundary reconstruction audit:

```text
groups=218175 quadratic=176204 linear=41971 expanded_search=85329
rank_fallback=41236 condition_fallback=0 coverage_fallback=735 donor_fallback=0
```

The 735 coverage fallbacks include the 285 explicit one-normal-band groups;
the other 450 and all 41,236 rank fallbacks are ordinary decisions made by the
existing adaptive compiler.  `quadratic + linear == groups` exactly.  The
quadrature reconstruction is cloned from this sealed boundary plan and its
fingerprint is checked before publication.

The same 16-rank marker-only summary reported, per rank, 68,334,704 arena
doubles (546,677,632 bytes), a 34,831,104-byte maximum workspace, and
54,554,528 service-staging bytes for the rank-zero local patch
`77 x 62 x 76`.  Their deliberately conservative arithmetic sum is
636,063,264 bytes (606.5972 MiB); it is a sizing bound, not a claim that the
three categories peak simultaneously.

## Compatibility

- The existing STL scan, `ImmersedSurfaceCompiler`, `EBTopologyCompiler`,
  `BoundaryStencilCompiler`, and `SurfaceQuadratureCompiler` paths are
  unchanged and remain the default unless the case layer explicitly selects
  imported-marker geometry.
- Marker cells and link ownership are decomposition-independent.  Each rank
  currently receives the replicated marker by API contract; a future
  distributed reader can be added outside this compiler without changing the
  published link semantics.
- Imported donor search has no implicit periodic images.  At physical-domain
  intersections it records the exact reachable tangential quadrant subset,
  while `QuadraticStencilCompiler` remains the rank/coverage/condition
  authority.
- The one-dimensional closure is unreachable for STL geometry, strict policy,
  multi-normal-band requests, or any failure other than exact detail 1304.
- No diagnostic logging or one-shot audit probe is retained in the product
  source tree.

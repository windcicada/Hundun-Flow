# HUNDUN-FLOW Stage 3 Numerical Contracts

This document is the compact implementation oracle for the approved Stage 3
specification and plan. It does not replace either document. Where this
summary and an approved source differ, the approved specification and then
the approved plan control.

The 2026-08-08 no-96 execution amendment supersedes any earlier `96^3`
execution wording: active Task 11/Stage 3 flow matrices use `12^3/24^3/48^3`.
Historical 96-cubed logs are evidence only and must not be rerun.

## Scope and capability boundary

The public capability names are `profile-1` through `profile-9`; each is
`implemented-and-accepted` in the capability ledger. The source tree is a
`0.2.0 candidate` until S3-V1 accepts the frozen formal matrix.

Stage 3 adds one static, closed, connected STL surface; a sharp-interface
Local-Flow-Pattern Ghost-Cell immersed-boundary path; real-surface pressure,
viscous and total force integration; and one WALE LES model. These additions
compose with the Stage 2 constant, material-density and ideal-gas paths on the
CPU-reference backend.

Stage 3 does not add moving or multiple immersed bodies, cut-cell
conservation, wall functions, configurable thermal walls, general mesh
import, production accelerators, external solver backends, rank-changing
restart, or any Stage 4 physical model.

The nine profiles are constant/material/ideal-gas crossed with body-fitted or
static IBM, with WALE present in the three WALE rows. Checkpoint v3 presence
1--9 is the serialization identity for that same truth table.

## Schema and composition

Schema version 3 preserves `simulation.type=variable_density_flow` and all
schema-v2 common-flow values. Its internal `common_flow.schema_version` is
always 2. The only legal module combinations are:

```text
immersed boundary absent + WALE absent  rejected
immersed boundary absent + WALE         accepted
LFP-GCIBM                + WALE absent  accepted
LFP-GCIBM                + WALE         accepted
```

The immersed geometry format is exactly `stl`; its length scale is positive
and finite; `fluid_side` is explicit; and its case-root path is normalized,
relative and non-escaping. Wall velocity is mathematically zero and is
canonicalized to FP64 positive zero. Enthalpy and every generic scalar use
`zero_normal_diffusive_flux`.

The closed input-safety ranges are:

```text
1.0e-6 <= Cw   <= 1.0
0.1    <= Pr_t <= 10.0
0.1    <= Sc_t <= 10.0
```

Approved numerical cases use `Cw=0.50`, `Pr_t=0.90`, and `Sc_t=0.70`.

## Geometry and surface

Define:

```text
L_ref = max(surface bounding-box diagonal, domain diagonal)
h_i   = cbrt(cell volume i)
h_max = collective max over active-cell h_i
eps   = numeric_limits<double>::epsilon()
```

The derived constants are:

```text
weld tolerance                    = 128*eps*L_ref
minimum triangle area             = 1024*eps*L_ref^2
intersection coincidence bound    = 512*eps*L_ref
minimum surface/nonperiodic-domain separation = 2*h_max
coverage point/link distance      <= 2*h_local
fluid/solid witness distance      <= 2*sqrt(3)*h_local
```

The reader accepts standard little-endian binary STL and strict C-locale
ASCII STL. It rejects non-finite input, size/count overflow, extra records or
tokens, degenerate triangles, open or non-manifold topology, inconsistent
orientation, self-intersection, multiple components, zero enclosed volume,
ambiguous classification, and intersection with an open domain boundary.
File normals are diagnostic only; product normals come from validated vertex
ordering.

For a periodic axis, a closed parent surface may extend beyond both planes of
the fundamental periodic domain instead of satisfying the `2*h_max`
separation. Merely touching the planes is not sufficient, because a cap on a
periodic plane would enter the half-open active window. The surface mesh must
be split at both periodic planes: no
triangle may straddle either plane. Classification still uses the complete
closed parent surface, while coverage, ghost surface measure, wall
quadrature, and force integration consume only triangles whose centroids lie
in the half-open fundamental interval `[origin, origin + length)`. A surface
that crosses only one plane is rejected because it does not define a complete
periodic image. Nonperiodic axes retain the original separation requirement.
For any reconstruction or local-flow-pattern sample crossing a periodic
plane, the canonical global cell ID remains the field identity. Reconstruction
moments, Halo-reach checks and execution ownership use the stencil-anchor-
nearest logical image and its domain-length position shift. Actual field
access uses the equivalent periodic image nearest the owning rank's box; this
can differ from the stencil image while naming the same canonical cell. LFP
sample geometry uses the wall-row-adjacent image. Keeping these three roles
separate preserves weights, ownership and local field indexing without
increasing the four-cell Halo reach.

Active pressure and momentum operators may gather a deterministic global
layout to freeze remote value offsets. The gathered ID vectors and lookup
tables are construction-only and must be released after binding; only the
runtime value buffers and frozen offsets remain. Releasing these tables does
not alter exchange order, operator coefficients, or reductions.

Triangle quadrature is the degree-2 three-point rule:

```text
(2/3,1/6,1/6), (1/6,2/3,1/6), (1/6,1/6,2/3)
weight = triangle area / 3
```

Stable IDs derive from canonical input order and global background IDs, never
from rank-local traversal. Classification, wall intersections, normals and
fingerprints must be invariant across approved decompositions.

## Active domain and immersed links

Cell centres are physically either fluid or solid. “Ghost” is a link-local
algebraic role of a solid neighbour for one active-fluid row; it is not a
third physical cell state. Every fluid--solid background link has one wall
intersection and one stable link record. Solid rows do not enter the linear
unknown vector, residual norm, conservation sum or source balance.

`SurfaceCoverage` requires nonzero global fluid cells, solid cells and
interface links; complete quadrature ownership; fluid- and solid-side
resolution witnesses; and a valid wall intercept for every operator link.
An unresolved body is rejected collectively.

`ActiveBoundaryLayout` retains only background faces owned by active cells.
An enabled inlet or outlet must retain an active face. Closed active domains
use constant-pressure nullspace projection and an active-volume zero-mean
gauge. Periodic active faces remain reciprocal.

## Reconstruction and Ghost constraint

The normalized local basis order is:

```text
1, n, t1, t2, n^2, n*t1, n*t2, t1^2, t1*t2, t2^2
```

The plan constants are:

```text
minimum donors       = 14
maximum donors       = 32
maximum halo reach   = 4 logical cells
required matrix rank = 10
condition limit      = 1.0e8
normal layers        >= 3 distinct positive-normal bands
tangential coverage  >= 4 deterministic quadrants
```

Donors are physical fluid cells ordered by physical distance and global cell
ID. Cell-average moments, not centre-point values, form the design matrix.
Uniform and warped cells use the same deterministic 12-tetra polyhedral
representation as `MeshGeometry`. A column-pivoted Householder QR chooses the
lower basis index when pivot magnitudes tie within
`64*eps*max(1,max_column_norm)`. The accepted condition estimate is
`max(abs(R_ii))/min(abs(R_ii))`.

Reproduction must satisfy:

```text
abs(error) <=
  512*eps*max(1, analytic_scale, weight_l1*data_scale)
```

Rank deficiency, conditioning failure, directional-coverage failure,
cross-wall donors or unavailable Halo reach is a collective initialization
failure. There is no lower-order fallback, implicit regularization, skipped
link or post-solve overwrite.

Each link-local Ghost value has one affine form:

```text
q_G = sum_j(w_j*q_Dj) + c_D*q_wall + c_N*g_wall
```

Velocity uses zero Dirichlet wall data. Enthalpy and generic scalars use zero
normal diffusive flux. Material density uses unconstrained one-sided
fluid-side extrapolation; ideal-gas density remains derived from enthalpy,
constant material properties and uniform thermodynamic pressure.

## Local Flow Pattern and unique residual

The standalone coefficient oracle uses:

```text
[A_N,A_S,A_W,A_E,A_L,A_R,A_P,S] = [2,3,5,7,11,13,41,17]
k0=0.37, k1=0.40, k2=0.65, k3=0.25
```

Its independent matrix oracle tolerance is:

```text
abs(product-oracle) <= 256*eps*max(1,row_l1_norm)
```

Product replacement additionally uses a direct quadratic evaluator. Link
groups are an exact partition and are permutation invariant. A row with one
through six immersed neighbours starts from one immutable background-row
snapshot. Shared active--active face contributions remain reciprocal; each
immersed occurrence is removed once, transformed once and replaced once.
Wall quadrature reads the accepted reconstruction but never writes solver
rows.

## Pressure, PISO and transaction

The Stage 2 two-corrector PISO transaction remains authoritative. Pressure
Ghost data is refreshed from the current attempt density, momentum diagonal
and wall-normal predictor flux before each corrector. Exactly two correctors
run; a failed final residual causes rollback and retry, never a third
corrector.

All attempt-local IBM, pressure-Ghost and WALE data is discarded on failure.
Committed/history state changes only after final residual, conservation, wall
penetration, outlet and force contracts pass collectively. The lowest
failing rank and stable failure classification are identical on every rank.

## Wall traction and force

Pressure, deviatoric viscous and total traction are evaluated on the real STL
surface with the same quadratic wall reconstruction used by the operator.
Quadrature ownership is unique by stable triangle ID. Pressure, viscous and
total operator/surface consistency errors are assessed separately.

Manufactured force scales are:

```text
A_ref                = analytic surface area
pressure force scale = rho_ref*U_ref^2*A_ref
viscous force scale  = mu*U_ref*A_ref/L_ref
total force scale    = max(pressure scale, viscous scale)
```

An independent oracle must show every reference force norm exceeds
`1.0e-6` of its scale before the product solve begins.

## WALE

For velocity gradient tensor `g`, let `S` be its symmetric part and `Sd` the
traceless symmetric part of `g*g`. WALE uses the Nicoud--Ducros homogeneous
ratio with no fixed additive epsilon:

```text
Delta  = cbrt(active background cell volume)
nu_t   = (Cw*Delta)^2
         * (Sd:Sd)^(3/2)
         / ((S:S)^(5/2) + (Sd:Sd)^(5/4))
mu_sgs = rho*nu_t
mu_eff = mu + mu_sgs
Gamma_h,sgs   = mu_sgs/Pr_t
Gamma_phi,sgs = mu_sgs/Sc_t
```

The exact zero-gradient result is bitwise positive zero. Near a canonical
wall, for `y=2^-k`, `k=4..12`:

```text
2.9 <= log(nu_t(y)/nu_t(y/2))/log(2) <= 3.1
```

At each attempt WALE is evaluated exactly once from the approved lagged
velocity and trial density. Startup uses `u^n`; BDF2 history uses
`u^n + r*(u^n-u^(n-1))`. The resulting coefficients remain frozen through
the predictor, two correctors, final transport, residual and force
evaluation. Retry recomputes them from unchanged committed/history state and
the new attempt time stencil.

## Formal accuracy and decomposition

Every formal second-order row uses three grids and two independently checked
orders:

```text
p1 = log(E_h/E_h2)/log(2)
p2 = log(E_h2/E_h4)/log(2)
p1 >= 1.8 and p2 >= 1.8
```

Errors must be finite, strictly positive, strictly decreasing and above the
approved roundoff-discrimination bound. The manufactured matrix uses
`12^3/24^3/48^3`, exact BDF2 history, and
`dt=0.05*h_max^2/(U_ref*L_ref)` for sphere, finite cylinder and oblique
rectangular prism on uniform and `[0.02,-0.015,0.01]` warped mappings, plus
translated-sphere and inside-cavity sequences.

Velocity and gauge-normalized pressure volume `L2/Linf`, near-wall errors,
wall penetration, pressure/viscous/total forces and all three
operator/surface consistency rows each meet both order gates. If every
penetration error is below `8192*eps*max(1,U_ref)`, the exact-enforcement
branch passes without fabricating a `0/0` order.

The formal near-wall velocity and pressure norms use one fixed physical
support across the matrix. For the approved 12/24/48 sequence the thickness
is frozen before execution at the coarse two-cell width,
`2*(L_ref/12)=L_ref/6`. A resolution-dependent `2*h_max` band is diagnostic
only and cannot satisfy or fail the formal row.

Approved decomposition grids are:

```text
1 rank: (1,1,1)
2 ranks: (2,1,1) and (1,2,1)
4 ranks: (4,1,1) and (2,2,1)
max field difference <= 5e-12*max(1,global reference Linf)
```

## Equality, inactive storage and failures

Rollback, failed-trial, Checkpoint continuation and bitwise metadata compare
FP64 bit patterns. Nested fields compare outer size, inner size and every
element. Numerical discretization uses only its approved tolerance.
Inactive owned cell and face slots remain bitwise `+0.0` through
initialization, trial, Halo, rollback, commit and restore.

Configuration, geometry, layout, capability, MPI-operation and file-integrity
failures are not retried. Trial non-finite/non-positive state, a derived
non-finite WALE coefficient, linear breakdown/non-convergence, outlet
backflow and final residual failure are recoverable under the frozen Stage 2
retry controller.

## Checkpoint, diagnostics and performance

Checkpoint v3 preserves the Stage 2 transaction and byte conventions and
adds canonical IBM/LES presence tags plus deterministic static-plan
fingerprints. It supports only identical rank count, process grid, owned
boxes and active layout. Derived query caches, raw addresses, force reports
and instantaneous WALE fields are not authoritative stored state.

Diagnostic module values append without renumbering Stage 2:

```text
immersed_surface=18, ghost_stencil=19, local_flow_pattern=20,
wall_force=21, les=22
```

Absent modules register no provider, instance or business counter. Local
diagnostics perform no implicit collective; collective diagnostics are
explicit, deterministic and read-only. Portable CI hard-checks exact
counters, bytes, messages and collective counts; timing, memory and
throughput remain informative compatible-baseline data.

## Public scientific sources

- Wang et al. (2024), DOI `10.1063/5.0195598`: local-flow-pattern geometry,
  local coordinates and coefficient reconstruction.
- Möller and Trumbore (1997), DOI
  `10.1080/10867651.1997.10487468`: ray/triangle intersection.
- Tseng and Ferziger (2003), DOI `10.1016/j.jcp.2003.07.024`: second-order
  Ghost-Cell boundary conditions on a non-staggered grid.
- Pan and Shen (2009), DOI `10.1002/fld.1942`: velocity/pressure Ghost-Cell
  convergence evidence.
- Seo and Mittal (2011), DOI `10.1016/j.jcp.2011.06.003`: conservation and
  pressure behavior of sharp-interface immersed boundaries.
- Nicoud and Ducros (1999), DOI `10.1023/A:1009995426001`: WALE.

These sources provide public mathematics and validation context. Product
code is independently derived from the approved contracts and tests.

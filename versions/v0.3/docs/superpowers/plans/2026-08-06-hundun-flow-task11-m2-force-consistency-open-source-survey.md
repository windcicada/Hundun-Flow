# Task 11 M2 force-consistency open-source survey and proposed repair

Status: frozen for review on 2026-08-06. This is the solution proposal the
user requested (literature + open-source survey, independent of the blocked
derivation). Product changes still wait for the RED-S3 protocol, but this
provides the concrete repair direction.

## 1. Literature consensus

The sharp-interface immersed-boundary literature consistently computes the
physical body force as a single surface quadrature of RECONSTRUCTED wall
values, using the same reconstruction authority as the momentum/ghost rows:

- Tseng & Ferziger (2003), *J. Comput. Phys.*: ghost-cell IBM; ghost values are
  defined by the wall condition; forces follow from the reconstructed wall
  state.
- Mittal & Iaccarino (2005), *Annu. Rev. Fluid Mech.*: survey; force is the
  integral of the pressure/stress reconstructed at the immersed surface.
- Yang & Balaras (2006), *J. Comput. Phys.*: sharp-interface IBM; force by
  integrating the reconstructed surface stress.
- Wang et al. (2024), DOI 10.1063/5.0195598: the Ghost-Cell/LFP construction
  that HUNDUN-FLOW follows; wall values and momentum rows share one
  reconstruction authority.

The key principle: there is ONE force path (boundary-pressure quadrature), so
an operator-vs-surface consistency problem does not arise.

## 2. Open-source survey (GitHub, read locally under /tmp/hf-refs)

| Project | Activity | Force architecture | Reuse / avoid |
|---|---|---|---|
| Basilisk `embed.h` (Popinet) | very active | `embed_force`: `Fp=-∫ p n dΓ`, `Fmu=-∫2μD·n dΓ`, integrated per cut cell at the fragment barycentre; pressure value from the cut-cell reconstruction | **single force path**; read-only GPL reference |
| **ghost-cell-IBM** (`SpencerSchwart/ghost-cell-IBM`, commit 2026-08-03) | **active** | `ibm_force`: `Fp = Σ area*p_boundary*n` with `p_boundary = extrapolate_scalar()` (inverse-distance × normal-projection weighted extrapolation of the fluid pressure to the boundary point); ghost rows replaced by the wall BC in the WLSQ polynomial | **the canonical single-path implementation**; shows exactly how the boundary pressure is defined for the force |
| AMReX EB (LBNL, daily) | very active | EB face values/forces; EBFluxRegister conservation | EB infrastructure, BSD reference for conservation; not a ghost-cell force |
| NASSLARD2D (Toronto CFD group, 2021) | inactive | `COMP_SURFACE_FORCE`/`COMP_DRAGF` on the ghost-cell IBM | second-order staggered ghost-cell; Fortran, simple |
| 3D-NS-FSI-WLSQ (2023) | low | WLSQ reconstruction; FSI | WLSQ reference; not force-focused |
| gslib `findpts` (Nek5000) | maintained | global donor lookup | donor semantics only |

## 3. Root cause restated from the survey

HUNDUN-FLOW has TWO force paths: the A22 row-sum operator force and the
WallForceIntegrator surface quadrature. The RED-S3 evidence (M2-A..M2-P)
established: every local quantity (wall value 2.46 order, linear row 1.96,
MMS row 1.92, surface 2.0-2.1) is second order, but the INTEGRATED `C_exact`
is first order (1.23). The open-source consensus shows the fix: the operator
row's wall-face pressure term must use the SAME wall-anchored authority
reconstruction (origin at the wall, link-normal frame, the exact object the
surface quadrature integrates) that defines the boundary pressure -- not the
row reconstruction (cell-center origin, identity frame).

## 4. Proposed minimal repair (R8 coupling, literature-backed)

In `build_complete_pressure_boundary_row_plan`
(`finite_volume/src/immersed_operator.cpp`), replace the wall-face pressure
value contribution -- currently implied by the row reconstruction's constraint
solution -- with the per-link wall-anchored authority reconstruction's value at
the wall intercept:

```text
p_wall(link) = boundary_coefficient(g) + sum(donor_weight_i * p_avg_i)
              from the authority reconstruction (origin = wall, frame = link normal)
F_operator_wall = sum over links of p_wall(link) * A_wall * n_link
```

This is exactly the ghost-cell-IBM `ibm_force` structure: the force integrates
the boundary pressure reconstructed by the wall authority. The operator force
then shares the surface force's wall-value authority and its 2nd-order
convergence; the momentum-budget closure (`budget_reaction_N`, the row-sum)
is unchanged.

### File whitelist (M2)

```text
finite_volume/src/immersed_operator.cpp          (A22 row wall-face value)
immersed/src/ghost_stencil_plan.cpp              (expose per-link authority reconstruction to the row build)
immersed/src/quadratic_reconstruction.cpp        (only if a new functional is needed)
tests/support/stage3_mms.{hpp,cpp}               (exact-state oracle stays as RED gate)
tests/numerical/test_laminar_ibm_order.cpp       (decision-tree selector)
tests/mpi/test_immersed_operator.cpp             (row-level RED)
```

### Mutation list

- M2-M1: revert the row wall-face value to the row reconstruction -> `C_exact`
  must return to first order.
- M2-M2: evaluate the authority value at the cell center instead of the wall
  intercept -> first order.
- M2-M3: use the unconstrained (non-wall-gradient) authority value -> the wall
  gradient datum must be lost.
- M2-M4: integrate against the row area/normal instead of the link
  area/normal -> first order.

### Verification gates

1. Exact-state wall-value oracle: the row wall value must equal the authority
   wall value (bitwise/near) and stay second order.
2. `C_exact` (decision tree): >= 1.8 on the formal rows.
3. `functional_selection_fast` pressure consistency >= 1.8.
4. Screen (24/48) and then the formal 24/48/96 matrix.
5. Full focused regression (Release 26 tests + Debug decomposition 5/5).

### Rollback

Disposable candidate trees; frozen M2 pre-change hashes; no threshold, PISO
corrector, damping or per-case tuning. The momentum budget closure and the
RED-S1 four-field semantics are untouched.

## 5. Expected outcome

The operator force and the surface force share one wall-value authority, so
their integrated difference converges at second order, closing the `C_exact`
gate and unblocking RED-S4, the formal matrix, and Task 11.

# Task 11 A22-A2 Rhie--Chow Interface-Adjoint Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to execute task-by-task. The approved Stage 3 protocol assigns the coupled mathematics, cross-module review and final acceptance to the main agent; no worker receives this complete plan.

**Goal:** Add a conservative A22 interface-difference correction to the existing time-consistent Rhie--Chow pressure operator while preserving checkerboard control, the complete mechanical cell correction and every frozen Task 11 contract.

**Architecture:** Keep `ActivePressureOperator` as `K_RC`. Build one private immutable interface plan from `pressure_constrained-background_pressure_constrained`, wrap `K_RC` in an additive `A2PressureOperator`, route the affine source through the same coefficients, and add a canonical conservative interface face correction. The complete A22 momentum residual and three full momentum-correction solves remain unchanged.

**Tech Stack:** C++17, MPI-3, existing HUNDUN mesh/runtime/execution/linear/field interfaces, CMake/CTest, codegraphf, `rg`, Debug/Release/ASan/UBSan. No new runtime dependency and no Python public build or execution path.

## Global Constraints

- Scientific authority: `docs/superpowers/specs/2026-08-04-hundun-flow-task11-a22-a2-rhie-chow-interface-adjoint-design.md`.
- Accepted Task 10 base: `0db56e463470dd1a605709ba05d8bd6a900f496b`.
- Candidate HEAD at freeze: `28f5dd541a0e3ce9ecf852e53d83981add3a5be8`.
- Preserve all pre-existing dirty Task 11 work and inspect the complete diff from the accepted Task 10 base before acceptance.
- Use the main agent for mathematical derivation, cross-module changes, requirements review, code-quality review and exact-HEAD verification.
- Use `codegraphf` for symbol/caller/impact navigation and `rg` for exact and repeated-pattern searches.
- Use `apply_patch` for manual edits.
- Add no callable public API, schema key, Restart field, diagnostic stable ID, plugin ABI, dependency, filter, damping, case tuning, third corrector or post-solve overwrite.
- The sole permitted public-header change remains the private friend forward declaration already present in `immersed_operator.hpp`; it changes no object layout or exported callable API.
- Do not enter G1/G2/G3, Task 12, Stage 4, private legacy source, publication or push.
- Do not run 96-grid work until the algebraic/direct gates and 24/48 screening are green. Run one large numerical job at a time with at most 96 logical CPUs.
- Keep one final Task 11 candidate commit with exact DCO; do not create intermediate commits from the existing dirty worktree.

---

### Task 1: Freeze A2 evidence and establish the missing-interface RED

**Files:**

- Create: `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-a2-evidence-matrix.md`
- Modify: `tests/mpi/test_immersed_operator.cpp`
- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Test: `tests/mpi/test_immersed_operator.cpp`

**Interfaces:**

- Consumes: immutable `BoundaryRowPlan` snapshots already produced by `ImmersedOperatorAdapter`.
- Produces: tests-only raw A22, background and difference row evidence; no product behavior.

- [ ] **Step 1: Record the evidence matrix**

  Map each A2 requirement to its product location, positive test, mutation-sensitive failure, FP64 equality rule, 1/2/4-rank coverage, presets and prohibited adjacent scope. Record design/plan SHA-256 and the rejected Route A logs.

- [ ] **Step 2: Write a compile-time RED for the desired tests-only snapshot**

  Add the wished-for test API before adding its implementation:

  ```cpp
  struct InterfacePressureRowSnapshot final {
    mesh::GlobalCellId momentum_cell{};
    std::uint64_t authority_fingerprint{};
    std::vector<PressureCouplingDonorTermSnapshot> background_donor_terms;
    std::vector<PressureCouplingDonorTermSnapshot> a22_donor_terms;
    std::vector<PressureCouplingDonorTermSnapshot> difference_donor_terms;
    std::vector<PressureCouplingWallTermSnapshot> background_wall_terms;
    std::vector<PressureCouplingWallTermSnapshot> a22_wall_terms;
    std::vector<PressureCouplingWallTermSnapshot> difference_wall_terms;
  };

  static std::vector<InterfacePressureRowSnapshot>
  interface_pressure_rows(const ImmersedOperatorAdapter &);
  ```

  In the test, call `interface_pressure_rows(adapter)` and require at least one nonempty interface difference row. The production change that makes this pass is an immutable export of the declared A22/background difference; the old full-row export cannot satisfy it.

- [ ] **Step 3: Run the RED and bind its evidence**

  Run:

  ```bash
  cmake --build --preset debug -j2 --target test_immersed_operator
  ```

  Expected: compile failure naming the missing `interface_pressure_rows` implementation, not an unrelated header, environment or MPI error. Save the log and SHA-256.

- [ ] **Step 4: Preserve the integrated Route A failure as a separate RED**

  Run the existing one-rank PISO test on the current full-normal-equation experiment and retain the recorded 500-iteration failure. Do not use this failure as the A2 GREEN criterion; it proves only that the rejected product experiment remains unaccepted.

### Task 2: Prove the unique raw interface difference in tests only

**Files:**

- Modify: `finite_volume/src/immersed_operator.cpp`
- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Modify: `tests/mpi/test_immersed_operator.cpp`

**Interfaces:**

- Consumes: `BoundaryRowPlan::background_pressure_constrained`, `BoundaryRowPlan::pressure_constrained`.
- Produces: tests-only raw maps and independent A2 algebra oracle.

- [ ] **Step 1: Add one internal aggregation helper**

  Implement one private helper that aggregates an existing `PressureAffinePlan` without rebuilding reconstruction weights:

  ```cpp
  using PressureDonorMap =
      std::map<std::pair<mesh::GlobalCellId, std::uint32_t>, double>;
  using PressureWallMap =
      std::map<std::pair<immersed::ImmersedLinkId, std::uint32_t>, double>;

  struct InterfacePressureMaps final {
    PressureDonorMap background;
    PressureDonorMap a22;
    PressureDonorMap difference;
    PressureWallMap background_wall;
    PressureWallMap a22_wall;
    PressureWallMap difference_wall;
  };
  ```

  Build `difference=a22-background` in stable key order. Use this helper for tests-only snapshots and later product export; do not keep the current duplicated full-row snapshot implementation.

- [ ] **Step 2: Make the snapshot RED GREEN**

  Export only rows with immersed links. Include exact raw maps before closure canonicalization. Verify background wall coefficients are currently exact zero and the difference wall map equals A22, while computing both sides explicitly.

- [ ] **Step 3: Add independent no-double-count assertions**

  In `test_immersed_operator.cpp`, derive literal map differences and prove:

  ```text
  difference == a22 - background
  ordinary active row count in A2 export == 0
  no shared-face-only donor identity appears twice
  pressure_unconstrained substitution changes the oracle and fails
  omitted background term changes the oracle and fails
  duplicated A22 term changes the oracle and fails
  background/A22 nested-size mutation fails
  ```

  Use scaled FP64 comparisons for coefficients and exact comparisons for row, donor, component and link identities.

- [ ] **Step 4: Prove row closure and mutation sensitivity**

  For each component require raw `sum(E_row)` within
  `64*epsilon*max(1,sum(abs(E_row)))`. Prove a constant vector passes, then mutate one donor coefficient by a scaled nonzero amount and prove closure fails.

- [ ] **Step 5: Run uniform and warped 1/2/4-rank tests**

  Run the registered `test_immersed_operator` uniform/warped cases for ranks 1, 2 and 4. Require identical sorted global row/donor/link identities and coefficients within the frozen decomposition threshold.

### Task 3: Build the private immutable E/B/H_E authority

**Files:**

- Modify: `finite_volume/src/pressure_coupling_detail.hpp`
- Modify: `finite_volume/src/immersed_operator.cpp`
- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Modify: `tests/mpi/test_immersed_operator.cpp`

**Interfaces:**

- Consumes: the Task 2 aggregation helper.
- Produces:

  ```cpp
  struct InterfacePressureDonorTerm final {
    mesh::GlobalCellId pressure_cell{};
    std::uint32_t output_component{};
    double coefficient{};
  };

  struct InterfacePressureWallTerm final {
    immersed::ImmersedLinkId link{};
    std::uint32_t output_component{};
    double coefficient{};
  };

  struct InterfacePressureRow final {
    mesh::GlobalCellId momentum_cell{};
    std::uint64_t authority_fingerprint{};
    std::vector<InterfacePressureDonorTerm> donor_terms;
    std::vector<InterfacePressureWallTerm> wall_terms;
  };

  class PressureCouplingAccess final {
  public:
    static std::vector<InterfacePressureRow>
    interface_rows(const ImmersedOperatorAdapter &);
  };
  ```

- [ ] **Step 1: Write the product-export RED**

  Change the test oracle to call `PressureCouplingAccess::interface_rows` and require equality with the independently captured raw difference. Expected first result: compile failure because `interface_rows` is absent.

- [ ] **Step 2: Implement bounded canonical closure**

  Copy the raw difference map, choose the momentum cell or lowest donor as stable anchor, and replace only that coefficient with the negative sum of the other coefficients when the pre-closure defect is within the frozen bound. Reject larger defects. Preserve raw closure evidence in tests-only snapshots; do not alter the mechanical row.

- [ ] **Step 3: Replace the rejected full-row export**

  Remove the old shared-face/full-`C_A22` accumulation from `PressureCouplingAccess::rows`. Export only canonical interface `E` and difference `B`. Remove the duplicated tests-only full-row implementation and make both product/test snapshots consume the single aggregation helper.

- [ ] **Step 4: Extend the independent H_E oracle**

  Construct the canonical active-face graph from gathered global face IDs and active endpoints. For each `(momentum cell,component)` route each non-anchor donor to the anchor over a shortest active path with global-face-ID tie breaking. Prove coefficient-by-coefficient:

  ```text
  E = -H_E^T D^T
  ```

  Mutate a face orientation, path face, periodic reciprocal and duplicate route; each mutation must fail.

- [ ] **Step 5: Prove additive spectrum and diagonal in the oracle**

  Build literal small-fixture `K_RC`, `E`, positive `M` and
  `K_RC+rho(E S^-1)^T M(E S^-1)`. Prove symmetry, nonnegative energy,
  constant/reference behavior, exact additive diagonal and nonzero 4x4x4
  checkerboard energy. This replaces the old proof that only demonstrated the
  full centered-gradient checkerboard nullspace.

### Task 4: Replace the rejected experiment with the additive A2 operator

**Files:**

- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Modify: `tests/mpi/test_immersed_pressure_operator.cpp`
- Modify: `tests/mpi/test_immersed_piso.cpp`

**Interfaces:**

- Consumes: `ActivePressureOperator` as `K_RC`, `PressureCouplingAccess::interface_rows` as `E/B`.
- Produces a private `A2PressureOperator : linear::LinearOperator` with:

  ```cpp
  void replace(double rho_ref,
               const std::array<std::vector<double>, 3>& momentum_diagonal,
               std::uint64_t dependency_revision);

  void affine_wall_source(
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>&
          wall_gradient_increments,
      execution::VectorView<double> output) const;

  void prepare_face_correction(
      execution::VectorView<const double> scaled_pressure_correction,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>&
          wall_gradient_increments) const;

  void add_interface_face_correction(
      runtime::FaceFieldView<double> face_mass_flux,
      runtime::FaceFieldView<double> face_velocity) const;
  ```

- [ ] **Step 1: Write additive operator RED tests**

  Require product apply values to equal an independent literal
  `K_RC+E^TME` oracle, diagonal to equal `diag(K_RC)+sum(E^2 M)`, constant
  mode to remain zero, parity energy to remain nonzero, repeated apply to be
  bitwise deterministic, and the apply schedule to retain `1234`. Mutate the
  additive sign and omit one E row in test copies; each oracle must fail.

- [ ] **Step 2: Delete the rejected full-normal-equation product path**

  Rename/replace `AdjointPressureOperator`; do not retain a selectable full
  `C_A22^T M C_A22` fallback. Keep `ActivePressureOperator` unchanged as the
  compact face-mobility and checkerboard authority.

- [ ] **Step 3: Implement additive apply and exact diagonal**

  `A2PressureOperator::apply` first invokes `K_RC.apply(x,y)`, then reuses the
  compact operator's completed global scaled-pressure exchange to evaluate
  local `E S^-1 x`, applies the exact component mobility, performs the
  preallocated transpose accumulation, and adds
  `rho S^-1 E^T M E S^-1 x` to `y`. It performs no allocation or plan rebuild.

- [ ] **Step 4: Bind one revision and reference authority**

  Require the base operator revision, A2 revision, density, momentum diagonal,
  active layout and interface authority fingerprint to match. Delegate
  pressure-reference ownership solely to `K_RC`. Reject stale/mismatched input
  collectively.

- [ ] **Step 5: Wire solver and tests-on probes to A2**

  Update preconditioner, solve, independent residual, reference normalization,
  operator-value probe and apply-schedule probe to use `A2PressureOperator`.
  Continue using `ActivePressureOperator` for time-consistent face mobility.

- [ ] **Step 6: Run focused operator tests**

  Run header/contract tests, `test_immersed_pressure_operator` and the
  operator-only subset of `test_immersed_piso` on 1/2/4 ranks. Do not run
  numerical grids while the homogeneous/affine face identity remains RED.

### Task 5: Route the affine wall source and conservative face increment

**Files:**

- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_piso_detail.hpp` only if a private report field is required
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Modify: `tests/mpi/test_immersed_piso.cpp`
- Modify: `tests/mpi/test_immersed_transaction.cpp`

**Interfaces:**

- Consumes: corrected-minus-current wall-gradient increments, scaled pressure correction `q`, one shared face-mass-flux field.
- Produces: exact `s_g` and `delta_m_E` with no committed-state mutation.

- [ ] **Step 1: Write the affine RHS RED**

  Construct one nonzero wall-gradient increment and require the pressure RHS
  difference to equal the independent literal
  `-rho S^-1 E^T M B delta_g`. Check zero increment is exact zero. Wrong sign,
  duplicated `B`, missing link and cross-assigned link mutations must fail.

- [ ] **Step 2: Reorder coefficient preparation before RHS projection**

  Compute the dependency revision and call compact/A2 `replace` after the
  predictor-divergence RHS is assembled but before nullspace projection.
  Build `delta_g=corrected-current`, evaluate `s_g`, subtract it once, then
  apply the existing closed-domain projection. Do not change solve controls.

- [ ] **Step 3: Implement immutable H_E route records**

  Gather the canonical active-face graph only at construction. Store stable
  local records for internal and periodic-pair faces plus preallocated global
  momentum-workspace counts/offsets. Validate owner, reciprocal, active path
  and route identity collectively.

- [ ] **Step 4: Compute and add the interface face correction**

  In `prepare_face_correction`, reuse the compact pressure exchange, form
  `w=E_s q+B delta_g`, multiply by exact `M`, gather the preallocated values
  required by local route records, and compute `delta_m_E=-rho H_E M w`.
  After the existing compact face loop, add each canonical increment once and
  update only the normal face-velocity component with
  `delta_m_E*S/(rho*|S|^2)`. Set periodic partners reciprocally. Leave
  immersed wall face velocity/flux exact positive zero.

- [ ] **Step 5: Make direct homogeneous and affine identities GREEN**

  The authoritative helper must compare outer size, inner size and every
  element. Require on 1/2/4 ranks:

  ```text
  operator(q) + s_g
    == normalized divergence(compact face correction + delta_m_E)
  rho*dot(face_velocity,S) == stored face_mass_flux
  wall mass flux bits == bits(+0.0)
  corrector count == 2
  ```

  Use scaled FP64 tolerance for algebraic accumulation and bitwise semantics
  for positive zero and transactional state.

- [ ] **Step 6: Preserve the full mechanical cell correction**

  Rerun the existing mutation-sensitive record proving
  `A_u delta_u+C_A22 delta_p+B_A22 delta_g=0`. Search all pressure-correction
  call sites with codegraphf/`rg` and prove no diagonal-mobility overwrite,
  third corrector or alternate wall datum was introduced.

### Task 6: Close MPI, failure, lifetime and tests-off evidence

**Files:**

- Modify only the focused test/helper files needed by a demonstrated evidence gap.
- Update: `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-a2-evidence-matrix.md`

**Interfaces:**

- Consumes: final A2 product candidate.
- Produces: algebraic/direct task-candidate evidence, not yet formal numerical acceptance.

- [ ] **Step 1: Run failure and rollback tests**

  Cover invalid closure, donor, link, route, periodic reciprocal, density,
  mobility, diagonal, revision and one-rank MPI failure. Require deterministic
  classification and lowest failing rank. Cover first/second-corrector and
  final rejection rollback with the authoritative bitwise nested-state helper.

- [ ] **Step 2: Prove allocation/rebuild neutrality**

  Add tests-on counters for interface plan builds, route builds and apply/face
  workspace use. Repeated apply and repeated face probes must allocate/rebuild
  zero product plans, return bitwise-identical output and leave state/counters
  unchanged except the explicitly observed apply schedule.

- [ ] **Step 3: Run the focused rank/preset matrix**

  Run Debug unit/MPI/operator/PISO/force/transaction tests in header, unit,
  MPI-1, MPI-2, MPI-4 order. Run focused Release. Run ASan/UBSan lifetime,
  repeated apply, collective failure and rollback tests.

- [ ] **Step 4: Verify private-boundary discipline**

  Build tests-off fully, compile standalone frozen public-header consumers,
  and verify test-access symbols/private headers do not enter installed public
  interfaces. Run source-policy, provenance, `nm`, `ldd` and
  `git diff --check`.

- [ ] **Step 5: Perform main-agent requirements and code-quality reviews**

  Review the complete diff from accepted Task 10 base. Requirements review
  closes every A2 evidence row and prohibited-scope boundary. Code-quality
  review checks coefficient ownership, units/signs, deterministic ordering,
  integer/MPI counts, exception safety, revision/lifetime, allocation-free
  apply, duplicated formulas and all codegraphf/`rg` callers. Repair and repeat
  the affected full review if any finding remains.

### Task 7: Numerical screening and Task 11 closure

**Files:**

- Modify `tests/support/stage3_mms.*` or numerical oracle files only for a proved oracle defect; do not add a case-specific product branch.
- Update Task 11 requirements, code-quality and acceptance records under `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/`.

**Interfaces:**

- Consumes: algebraically accepted A2 candidate.
- Produces: final Task 11 verdict; no Task 12 capability.

- [ ] **Step 1: Run 12/24 localization**

  Run Release tests-only sphere, prism and translated-sphere localization on
  uniform/warped mappings. Inspect maximum-error row type, A22/background/E
  fingerprints, pressure/velocity/force authority and conservation. A failed
  screen stops A2; do not tune a case or add stabilization.

- [ ] **Step 2: Run the complete 24/48 screen**

  Run every affected frozen sequence and require the first adjacent order
  `>=1.8` with all original thresholds. Classify any failure as program defect,
  oracle defect or A2 scientific failure before changing code.

- [ ] **Step 3: Run one formal Release 24/48/96 matrix**

  On one stable exact candidate, run the nine frozen shape/mapping/translation/
  cavity sequences once. Require both adjacent orders `>=1.8`, conservation,
  penetration, force, consistency and nullspace contracts. Bind HEAD, preset,
  toolchain, command, rank/parameters, binary SHA-256, log SHA-256, result and
  elapsed time.

- [ ] **Step 4: Run decomposition and engineering acceptance**

  Run the approved 1/2/4-rank process-grid matrix and long cylinder/sphere/
  transient engineering gates one large task at a time. Require frozen field,
  force, energy, residual, conservation and lowest-rank contracts.

- [ ] **Step 5: Run exact-HEAD final Task 11 verification**

  Run complete Debug on the exact candidate, focused Release/ASan/UBSan,
  tests-off, headers, policy, provenance, linkage and all frozen Task 11 exit
  commands. Any later product/public/build/MMS/oracle edit invalidates and
  reruns the affected evidence.

- [ ] **Step 6: Accept or stop Task 11**

  Confirm requirements compliant, code quality approved, complete diff
  inspected, no skipped tests, no residual worker/test process, no private
  access and no publication/push. Verify the one final commit parent, subject
  and exact DCO, record accepted HEAD/report/log SHA-256, then stop. Only after
  formal Task 11 acceptance may the already-approved G1 -> G2 -> G3 sequence
  begin; this plan itself never enters it.

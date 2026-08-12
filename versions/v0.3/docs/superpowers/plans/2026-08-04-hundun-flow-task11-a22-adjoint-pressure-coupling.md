# Task 11 A22 Adjoint Pressure Coupling Implementation Plan

Status: stopped at the Task 2/Task 5 algebraic-spectral gate. The full
`C^T M C` normal-equation hypothesis is rejected and Tasks 6--11 in this plan
are not authorized without a replacement scientific amendment.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. The current user instruction instead assigns the coupled mathematical design, review and final acceptance to the main agent; a worker may receive only a bounded independent coding task.

**Goal:** Replace the split compact pressure/face correction with one private,
conservative, discrete-adjoint pressure coupling derived from the accepted A22
mechanical row, while preserving all frozen Task 11 and Stage 1/2 contracts.

**Architecture:** Build one immutable private sparse plan containing the A22
pressure derivative `C`, wall derivative `B`, and canonical conservative
velocity-to-face route `H` with `C=-H^T D^T`. Implement
`L=rho (C S^-1)^T M (C S^-1)`, derive its diagonal from the same coefficients,
route the affine wall source through `B`, and update the one shared face-flux
field through `H`. Keep the full SPD momentum solve for cell correction.

**Tech Stack:** C++17, MPI-3, existing HUNDUN execution/linear/mesh/field/Halo
interfaces, CMake/CTest, codegraphf, `rg`, Debug/Release/ASan/UBSan. No new
runtime dependency and no Python public build or execution path.

## Global Constraints

- Accepted base: `0db56e463470dd1a605709ba05d8bd6a900f496b`.
- Candidate HEAD at freeze: `28f5dd541a0e3ce9ecf852e53d83981add3a5be8`.
- Review the complete diff from accepted base, not only A22 edits.
- Preserve all current dirty Task 11 work and unrelated files.
- No public method/type, schema, Restart, diagnostic stable-ID or ABI change.
  The sole permitted public-header edit is a private friend forward
  declaration that exposes no callable API and changes no object layout.
- No G1/G2/G3, Task 12, Stage 4, private legacy source, publication or push.
- No threshold change, case tuning, filter, damping, third corrector or
  post-solve overwrite.
- No 96-grid run until every algebraic/direct and 24/48 screening gate is
  green; one large task at a time, at most 96 logical CPUs.
- Main agent performs the mathematical proof, complete-diff review,
  requirements review, code-quality review and exact-HEAD acceptance.
- Use `apply_patch` for edits and codegraphf before broad manual navigation.
- Keep one final Task 11 commit with exact DCO; do not create intermediate
  candidate commits from the existing dirty Task 11 worktree.

---

### Task 1: Freeze evidence matrix and verify the existing R4 RED

**Files:**

- Read: `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-*.md`
- Modify: `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a22-adjoint-coupling-evidence-matrix.md`
- Test: `tests/mpi/test_immersed_piso.cpp`

- [ ] Record each design requirement, implementation location, positive and
  mutation-sensitive negative evidence, equality semantics, rank coverage,
  preset coverage and forbidden adjacent scope.
- [ ] Use `codegraphf context` for `ActivePressureOperator`,
  `correct_active_pressure`, `evaluate_boundary_row_once`, and all test-access
  callers; use `rg` to enumerate every pressure/face update assertion.
- [ ] Run the existing direct pressure-flux identity test on one rank and save
  the expected mismatch log and SHA-256. Confirm that the failure is the
  operator/routed-divergence assertion, not setup, MPI or environment failure.
- [ ] Run `git diff --check` and record the unchanged accepted base and current
  candidate HEAD.

### Task 2: Add a tests-only sparse algebra oracle

**Files:**

- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Modify: `finite_volume/src/immersed_operator.cpp` only for tests-on immutable
  snapshots; no product behavior yet
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Modify: `flow/src/stage3_flow.cpp` only for tests-on probes; no product
  pressure path yet
- Modify: `tests/mpi/test_immersed_operator.cpp`
- Modify: `tests/mpi/test_immersed_piso.cpp`

- [ ] Add immutable tests-on snapshots for the complete `C` and `B` rows using
  the already-built A22 authority; do not reconstruct coefficients in the
  test-access layer.
- [ ] In the test oracle, construct signed incidence `D`, deterministic route
  `H`, scaled `C_s`, positive diagonal `M`, and dense small-fixture `L` by
  literal loops independent of the future product apply.
- [ ] Add mutation-sensitive helpers: exact copy passes; donor coefficient,
  normal sign, wall coefficient, route face, duplicate route and nested vector
  size mutations fail.
- [ ] Prove constant/linear/quadratic row reproduction, affine `C/B` response,
  `C=-H^T D^T`, symmetry, nonnegative quadratic form, constant nullspace,
  pressure-reference behavior and exact diagonal on one rank.
- [ ] Add 1/2/4-rank fingerprint/decomposition checks with the frozen max-field
  equality semantics.
- [ ] Run the focused tests. If any mathematical gate fails for the accepted
  row coefficients, stop Route A before product integration and write a
  stopping report; do not continue to Task 3.

### Task 3: Extract one private immutable pressure-coupling factory

**Files:**

- Create: `finite_volume/src/pressure_coupling_detail.hpp`
- Modify: `finite_volume/include/hundun/finite_volume/immersed_operator.hpp`
  only for the private friend forward declaration
- Modify: `finite_volume/src/immersed_operator.cpp`
- Modify: `finite_volume/src/immersed_boundary_authority_detail.hpp`
- Modify: `finite_volume/src/immersed_operator_test_access.hpp`
- Test: `tests/mpi/test_immersed_operator.cpp`

- [ ] Move the accepted complete-row coefficient result into one private
  authority object consumed by mechanical residual evaluation and coupling
  export; do not duplicate the quadratic reconstruction or coefficient map.
- [ ] Assemble ordinary shared-face `C` terms, physical zero-normal-gradient
  owner terms, interface replacement `C` terms and wall `B` terms in stable
  `(row,component,input)` order; prescribed pressure values remain affine
  constants.
- [ ] Apply only the bounded constant-mode roundoff closure defined in the
  design, to the shared authority. Reject a defect above the bound.
- [ ] Construct canonical active paths and aggregate `H` route coefficients in
  stable `(global_face,row,component)` order. Reject solid-wall routing,
  disconnected donors, duplicate face ownership and unexplained nonzero sum.
- [ ] Freeze structure fingerprint and verify that residual evaluation and
  coupling export report the same row/plan identities.
- [ ] Run the direct row algebra and mutation tests on 1/2/4 ranks.

### Task 4: Implement reusable distributed C/C-transpose/H schedules

**Files:**

- Modify: `finite_volume/src/pressure_coupling_detail.hpp`
- Modify: `finite_volume/src/immersed_operator.cpp`
- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Test: `tests/mpi/test_immersed_piso.cpp`

- [ ] Build immutable owner/peer schedules once from global IDs and canonical
  face ownership. Validate every send and receive identity collectively.
- [ ] Allocate all product workspaces at flow construction/revision update;
  `apply_C`, `apply_C_transpose` and `apply_H` allocate nothing and rebuild no
  field or reconstruction.
- [ ] Implement interior work before wait and peer-boundary work after wait
  where the schedule permits; preserve deterministic accumulation order.
- [ ] Add tests-on counters proving no per-apply plan build/allocation and
  repeated-call bitwise determinism.
- [ ] Run homogeneous/affine apply and decomposition tests on 1/2/4 ranks,
  plus collective missing/duplicate/nonfinite input rejection.

### Task 5: Replace the compact pressure operator

**Files:**

- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_piso_detail.hpp`
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Test: `tests/mpi/test_immersed_pressure_operator.cpp`
- Test: `tests/mpi/test_immersed_piso.cpp`

- [ ] Change private `ActivePressureOperator` to apply
  `rho C_s^T M C_s` with the exact current assembled momentum diagonals.
- [ ] Derive and revision-bind the Jacobi diagonal from `C` and `M`; remove the
  compact operator's independent diagonal and face-mobility authority.
- [ ] Preserve `LinearOperator` layout/context/revision/event contracts,
  constant-nullspace normalization and pressure-reference behavior.
- [ ] Add direct symmetry/PSD/diagonal/zero-RHS/stale-revision/nonfinite tests.
- [ ] Watch the former R4 operator/route test remain RED until Task 6 wires the
  affine RHS and face update; all new operator-only tests must be GREEN.

### Task 6: Route the affine wall source and canonical face correction

**Files:**

- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_piso_detail.hpp`
- Modify: `flow/src/stage3_flow_test_access.hpp`
- Test: `tests/mpi/test_immersed_piso.cpp`
- Test: `tests/mpi/test_immersed_transaction.cpp`

- [ ] Evaluate `s_g=rho S^-1 C^T M B delta_g` once per corrector and subtract
  it from the normalized predictor-divergence RHS.
- [ ] After the solve, evaluate `w=C_s q+B delta_g` once and update the one
  canonical shared face mass-flux field with `-rho H M w`.
- [ ] Update only the face-velocity normal component needed to make
  `rho*dot(u_f,S_f)` identical to the stored face mass flux; keep wall faces
  exact positive zero and preserve physical-boundary rules.
- [ ] Retain the full three-component momentum correction solve with the same
  `C/B` residual difference and verify
  `A_u delta_u+C delta_p+B delta_g=0`.
- [ ] Make homogeneous and affine pressure-operator/face-divergence identities
  GREEN on 1/2/4 ranks; verify corrector count exactly two.
- [ ] Verify failed first/second corrector, stale revision and final rejection
  restore nested state, mechanical wall data, histories, face fields and
  metadata bitwise.

### Task 7: Align force and pressure-authority fingerprints

**Files:**

- Modify: `immersed/src/wall_force.cpp` only if the existing adapter lacks the
  shared authority fingerprint
- Modify: `immersed/src/wall_force_detail.hpp`
- Modify: `flow/src/stage3_flow.cpp`
- Test: `tests/mpi/test_wall_force.cpp`
- Test: `tests/mpi/test_immersed_piso.cpp`

- [ ] Make the true-surface pressure reconstruction report the same stable
  mechanical polynomial/link/donor identity as the A22 row authority while
  retaining independent surface integration.
- [ ] Prove a wall-gradient mutation changes row, predictor and traction once;
  missing, duplicated and cross-assigned link mutations fail.
- [ ] Verify no surface-force value is backfilled into the cell operator and no
  post-solve overwrite is introduced.

### Task 8: Focused Debug, sanitizer and tests-off closure

**Files:**

- Modify only direct test/oracle helpers when a proven evidence gap remains.

- [ ] Run focused Debug unit/MPI/operator/PISO/force/transaction tests in the
  order header/contract, unit, MPI 1, MPI 2, MPI 4.
- [ ] Run Release-focused direct algebra and PISO tests.
- [ ] Run ASan/UBSan-focused plan lifetime, peer schedule, repeated apply,
  collective failure and rollback tests.
- [ ] Build tests-off fully and run standalone frozen public-header consumers;
  verify the private test seams are absent.
- [ ] Run source-policy/provenance checks, `nm`, `ldd` and `git diff --check`.
- [ ] If two review cycles expose the same category, perform a main-agent
  closure sweep with codegraphf and `rg` before one bounded repair.

### Task 9: Numerical screening without formal 96

**Files:**

- Modify tests/support or numerical oracle files only for a demonstrated
  oracle defect; no case-specific product branch.

- [ ] Run tests-only Release 12/24 pressure localization for sphere, prism and
  translated sphere on uniform/warped mappings.
- [ ] Require the accepted threshold/margin and inspect maximum-error row type,
  coefficient fingerprint, pressure/velocity/force authority and conservation.
- [ ] Run the full affected Release 24/48 screening sequences and require the
  first adjacent order `>=1.8` with no threshold adjustment.
- [ ] If screening fails, classify program defect versus algorithm failure. An
  algorithm failure stops Route A and records evidence; do not launch 96 or add
  an empirical stabilization.

### Task 10: Formal Task 11 numerical closure

**Files:**

- No product edit after binding formal evidence unless the invalidated matrix
  is rerun.

- [ ] On one stable exact candidate, run exactly one Release 24/48/96 formal
  nine-sequence matrix (shape, mapping, translated sphere and cavity) and
  require both adjacent orders `>=1.8` plus all frozen conservation/force
  thresholds.
- [ ] Bind evidence to HEAD, preset, compiler/libstdc++/MPI versions, command,
  parameters, test binary SHA-256, log SHA-256, result and elapsed time.
- [ ] Run the approved 1/2/4-rank decomposition matrix and long engineering
  validation, one large task at a time and at most 96 logical CPUs.
- [ ] Any product, public/build, MMS/oracle or acceptance-helper edit after the
  evidence invalidates and reruns the affected formal matrix.

### Task 11: Main-agent review and final Task 11 acceptance

**Files:**

- Create/update Task 11 evidence, requirements, code-quality and final
  acceptance records under
  `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/`.
- Modify coordinator ledger only after exact-HEAD verification.

- [ ] Inspect the complete diff from accepted Task 10 base to candidate HEAD,
  including every A22 and pre-A22 Task 11 file.
- [ ] Requirements review: close every evidence-matrix row, public-header and
  tests-off boundary, MPI/failure classification, numerical threshold and
  prohibited-scope check.
- [ ] Code-quality review: inspect ownership/lifetime/revision, deterministic
  ordering, overflow, units/signs, MPI request/error paths, allocation-free
  apply, duplicated formulas and all callers found with codegraphf/`rg`.
- [ ] After any product/build/test-registration repair, rerun exact-HEAD full
  Debug and the affected focused Release/ASan/UBSan/formal evidence.
- [ ] Run the frozen Task 11 exact-HEAD final matrix, confirm no skipped tests,
  no residual worker/test process, no private access, no publish/push, and no
  unrelated worktree change.
- [ ] Verify final commit parent, subject and exact DCO, then record accepted
  HEAD and report/log SHA-256. Stop before G1 until Task 11 is formally
  accepted; after acceptance follow the already-approved G1 -> G2 -> G3 order.

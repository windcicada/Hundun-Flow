# v0.4 Tasks 14--19 spec/code completion audit

Date: 2026-08-22

Audit snapshot: clean governance commit
`1e6279fb386b28ef14cfa92550053c5c834fb62c`, tree
`7ed38beb4a4aeaa15390d25d9f76a886dc981c54`.  The subsequently delegated
FGMRES restart experiment is excluded from this baseline audit and must pass
its own frozen matrix before it can replace the audited Krylov path.

This is a main-agent source/spec audit, not a restatement of a short test run.
It compares every Task 14--19 production responsibility in
`docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md`
against public types, production consumers, failure paths and mutation/oracle
coverage.  It does not accept Task 20, a performance candidate or the release.

## Completion result

| Task | source/spec result | inspected production authority |
|---:|---|---|
| 14 unified equations | code-complete at the snapshot | `solver_equations.cpp` owns one coverage/alias-checked `AssemblyEpoch`; scalar Schmidt values are positive, cold-hashed inputs; pressure work distinguishes absolute pressure; IBM equation replacement is part of the same certified assembly path. No persistent `rhoU`/`rhoh` product field or PISO-specific duplicate equation implementation was found. |
| 15 exactly two PISO correctors | code-complete at the snapshot; any FGMRES replacement reopens its numerical integration row | case parsing rejects any pressure-corrector count other than two; `PisoPlan` publishes fixed FGMRES pressure authority; `PisoPressureSolveEpoch` admits exactly correctors one then two and finalizes only after two solves. Corrector one owns trial state, corrector two owns pending final pressure-equation flux, and terminal EOS/continuity/closed-mass/gauge audit precedes collective publication. |
| 16 WALE/Vreman wall treatment | code-complete at the snapshot | case-level model combinations are statically compiled; WALE/Vreman share the gradient/effective-viscosity authority. The IBM equation path evaluates positive bounded wall viscosity, retains resolved normal traction and substitutes wall-law tangential traction from relative velocity. Donor construction retains 14--32 donors, three positive-normal bands and four tangential quadrants rather than a relaxed fallback. |
| 17 IBM pressure/final force | code-complete at the snapshot | the pressure operator removes impermeable fluid-solid mass-flux links while keeping solid identity rows and the exact outer operator. Final force requires final velocity/pressure/gradient/effective-viscosity/final-flux revisions, uses the positive donor-envelope rule, performs status consensus before force reduction, and has independent oracle plus stale/sign/normal/area/provisional-state mutations. No final-force read of `HbyA` or corrector scratch was found. |
| 18 one production seal | code-complete at the snapshot | `ProductCompiler::compile` performs the complete cold registration/allocation/bind/seal transaction and rejects a second seal. The frozen graph contains exactly pressure stages 40/50 and model-contribution seams 12/45. The real hot-step authority exercises 100 repetitions with stable storage and allocation observation. Service snapshot capacity is part of the sealed product rather than added by Task 19. |
| 19 driver/restart/output | code-complete at the snapshot | CLI/init/run use run-directory authority; output and restart consume committed snapshots. Restart uses synced pending files/directories, collective validation, atomic generation/current switch, bounded pruning and root metadata broadcast before rank-remapped publication. VTI/VTR/Visit and structured evidence remain bounded services; startup/retry/restart-recovery rows are ineligible for statistics. |

## Cross-module findings

1. The `FrozenExecutionGraph`/`ContributionPlan` production seam remains
   `C1 -> transport -> PISO1 -> C2 -> PISO2`: graph stages 12 and 45 are the
   two model-contribution seams and stages 40 and 50 are the only pressure
   stages.  This audit neither imports nor changes Stage 5.
2. Native MG `numeric_refreshes` is the transaction that refills the current
   fine coefficients and rebuilds/validates all current coarse coefficients
   before publication.  Therefore Task 19's evidence mapping of both exact
   and coarse numeric refills to this counter is a lifecycle invariant, not a
   measured equality inferred from Re=3900.  `hierarchy_rebuilds` is the
   policy-authorized preconditioner setup count; numeric refreshes minus
   rebuilds is the setup-reuse count.  A future backend that can refresh fine
   and coarse numerics independently must add distinct counters before it can
   use the same evidence schema.
3. The recent MG changes preserve exact/coarse refill identities, operator
   order, transfer, outer FGMRES true-residual authority and hot allocation
   contract.  The rejected strict single-reduction FGMRES experiment did not
   enter this snapshot.  Its restart-on-breakdown successor is a separately
   frozen Task-15 numerical experiment and cannot inherit this audit's
   acceptance without fresh PISO/focused/sanitizer/isolation evidence.
4. Positive-property handling is atomic: nonpositive/nonfinite donors reject;
   a finite quadratic overshoot is projected only to the strictly positive
   donor envelope.  Force evaluation reaches an all-rank status consensus
   before any force reduction, so a rank-local reconstruction failure cannot
   diverge collective order.

## Gate consequence

No missing Task 14--19 production responsibility was found that requires a
baseline source fix.  Their code-completion status at the audited snapshot is
confirmed, but release completion is intentionally not claimed: Task 20
remains active, the current performance candidate is rejected, the
Parnaudeau profile receipt is incomplete, and no new immutable candidate or
formal COAST pairing is authorized.  Any accepted pressure-path modification
must rerun the full downstream evidence matrix and supersede this exact
snapshot explicitly.

## 2026-08-23 superseding addendum: collective fusion and IBM diagonal

This addendum supersedes the source-completion conclusion above for the
current route.  The original audit missed two cross-module defects later
exposed by the Re=3900 full-grid diagnostic: a rank-local failure could leave
different ranks entering different existing collective epochs, and the IBM
momentum replacement retained the exact quadratic traction only in the
residual without restoring a positive owning-row coefficient after removing
the Cartesian fluid-solid face.  The reviewed base is HEAD
`e854a575ddd74a799da46e7e4af86a8d805383a8`, tree
`abe242324395de5352fe9ebb3da1730f6a339a48`, with exact seven-file production
and test diff SHA-256
`6752b0578e96669b6b02cfc4db30859ac52e4a045b38fa2d9b7d634c82ce39a2`.

| Task | amended source/spec result | reviewed current authority |
|---:|---|---|
| 14 unified equations | complete for numerical admission; formal candidate gates remain pending | The exact quadratic IBM traction remains in the residual.  Cold compile stores wall distance and the L1 majorant `1/d + sum_j abs(w_j - delta(j,owner)/d)`.  The hot component solve restores a positive geometric-plus-correction diagonal and covers the absolute component row of the traction tensor; the equilibrium-wall path also covers the complete tangential-projector row.  Recomputing `rhs = diagonal*U - residual` preserves the exact equation at the trial state.  No hot scan, geometry query, allocation, field, message or collective is added. |
| 15 exactly two PISO correctors | complete for numerical admission | Rank-local pressure/terminal failures are passed into the already scheduled common terminal consensus, and rank-local final-gradient/turbulence failures enter the already scheduled stage-60 halo preflight.  The successful path adds no collective, correction, solve, halo or retry.  Corrector count, true-residual authority, pressure cap and terminal tolerances are unchanged. |
| 16 WALE/Vreman wall treatment | complete for the current scope | WALE/Vreman still share one gradient and effective-viscosity authority.  The equilibrium wall function retains exact resolved-normal traction and opposing tangential drag; only the deferred scalar row majorant changes.  Donor, band, quadrant, rank and conditioning contracts are unchanged. |
| 17 IBM pressure/final force | complete and unchanged downstream | The pressure operator and final-force authority are unchanged.  Final force still consumes only terminal-certified final state, reaches status consensus before reduction and never reads `HbyA` or corrector scratch. |
| 18 one product seal | complete | The new per-link coefficients allocate only during cold `IbmEquationInterfacePlan::compile` and are immutable after seal.  Graph stages 12/45 remain the contribution seams and 40/50 the only pressure stages. |
| 19 driver/restart/output | complete | Committed-attempt evidence, restart publication, output bounds and statistics eligibility are unchanged.  No output/restart schema or service capacity is added. |

The rejected predecessor used only the own-component tensor factor.  Its
frozen replay accepted 19 steps, then rejected attempted step 20 before either
pressure solve because two near-wall cells had predicted enthalpies below the
existing 200 K inversion range.  The thermodynamic limit was a detector, not
an authorized clipping or range change.  Changing only the deferred tensor
factor to its absolute component-row majorant produced exactly 20 accepted
rows and 40 converged pressure solves with no retry, restart recovery or hot
allocation.  Step 20 required seven PISO2 terminal audits and accepted at
`9.030541169188891e-7`, so the route is numerically admitted with an explicit
thin-margin risk; this dirty HUNDUN-only run is not a Task 20 gate or timing
sample.

Release focused passed 108/108.  A matching Clang C/C++ UBSan build passed
the frozen 28-test union once under delegated execution and again under root
execution, with no diagnostic.  The registered ASan 13-test manifest likewise
passed one delegated and one root execution against the same binaries with no
diagnostic.  Invalid earlier sanitizer attempts that mixed GCC C with Clang
C++, or omitted the frozen libc++ runtime path, ran no valid test bodies and
are retained only as rejected provenance.

The stable production/Stage-5 compatibility seam remains
`FrozenExecutionGraph/ContributionPlan`,
`C1 -> transport -> PISO1 -> C2 -> PISO2`.  Stage 5 is neither imported nor
modified by this addendum.  The dirty-tree prerequisites permit
candidate-identity closure; they do not accept `focused`, create an immutable
candidate, seal COAST equivalence, authorize pairing or advance any later
Task 20 gate.

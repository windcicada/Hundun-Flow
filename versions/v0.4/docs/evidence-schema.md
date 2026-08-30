<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 runtime evidence schema

`evidence/runtime.jsonl` uses `HUNDUN_V04_EVIDENCE_V5`. Each line is one
collectively agreed committed-step record. Writers reject rank-disagreeing
records and never observe trial state.

V5 retains the V4 pressure--enthalpy refinement lifecycle and makes two
accepted-step policies explicit. `momentum_predictor_limiter` identifies
`scheme:"common_face_afc_v2"`, records the L1 fraction of retained
anti-diffusive correction as `retained_correction_l1_ratio`, and counts global
unique limited faces as `limited_faces`. An inactive limiter must report ratio
one and zero faces; a limited row requires a strictly positive ratio below one
and at least one limited face. The legacy V4 names `theta` and `activations`
are not V5 fields.

The limiter embeds the momentum-advection certificate
`advective_cfl:{present,plan,time_revision,density_revision,face_flux_revision,activity_collective,dt,out_max,abs_max,limit}`.
`present` must be true; the plan and all three revisions must be nonzero; `dt`
and `limit` must be positive. Both CFL maxima must be finite and nonnegative,
and `out_max` must not exceed `limit` apart from a 64-binary64-epsilon
comparison allowance. Its `limit` must equal the terminal committed CFL
`limit`; both are the same configured acceptance ceiling. This certificate observes the face flux used by the
momentum predictor, so its `face_flux_revision` must differ from the committed
`terminal_physical_audit.final_flux_revision` produced after pressure
correction. Missing, invalid, or same-revision certificates reject the row.
The rank-local runtime certificate additionally binds the exact density view
(base, shape, stride, component, replica, field, revision, storage and
revision domain) and all three provisional-flux views (base, extents, stride,
axis, storage and revision domain). These pointers are checked by the product
driver before FGMRES and are deliberately not serialized into V5.

One V5 JSONL file is one immutable run. Its schema and
`build,binary,case,stl,product,cpu_plan` identity are frozen by the first row;
steps must then be contiguous and committed time must increase. For every row
after the first, the advective-CFL `dt` must equal the adjacent committed-time
increment within 128 binary64 epsilons. The momentum plan, IBM activity
collective and configured CFL limit are also static across the file. A zero
STL requires a zero activity collective, while an IBM case requires a nonzero
one. Mixed schemas, identity mutation, skipped steps, time reversal, changed
static authority or time--dt disagreement reject the file.

The terminal audit additionally contains
`committed_convective_cfl:{out_max,abs_max,limit}`. These are the rank-global
maximum outward and absolute convective CFL measures computed from the
committed mass flux, plus the positive configured acceptance limit. All three
values must be finite; both maxima must be nonnegative, and `out_max` must not
exceed `limit` apart from a 64-binary64-epsilon comparison allowance. Because
the values live inside `terminal_physical_audit`, they certify the same
committed final-flux revision as EOS, continuity, energy, mass and gauge. The
advective and committed CFL objects therefore form an explicit dual-revision
record rather than treating predictor and final mass flux as interchangeable.

V4 retained the V3 pressure lifecycle and terminal physical acceptance, and
adds an explicit same-target nonlinear-refinement record. The ordinary
`pressure_solve_calls` remains exactly two: refinements are not extra PISO
correctors. `pressure_energy_refinement_solve_calls` is a prefix count from
zero through six, `pressure_energy_refinement_termination` must be
`component_residuals_converged` for an accepted row, and
`pressure_energy_refinement` contains exactly that many successful solves in
ordinal order. Every entry binds one nonzero target generation and a unique,
rank-invariant collective state/flux lineage; all entries in a prefix share
the same target. Rank-local pressure-state and linear-system identities remain
inside the solver report and are deliberately not serialized. A V3 record
carrying any V4 refinement field is rejected rather than silently interpreted
under the older contract.

V3 made the pressure lifecycle and terminal physical acceptance explicit.
`pressure_solve_contract` is one rank-invariant contract shared by both
correctors: `pressure_continuity` or `continuity_energy_coupled`. A
pressure-continuity C2 row retains the V2 supplemental linear-audit
obligation. A coupled row requires both correctors to report zero supplemental
linear audits, matching the coupled Schur solve's null `selected_audit`; its
acceptance authority is the terminal physical audit below. An invalid or
missing contract, a pressure-continuity C2 without its audit, or a coupled
solve with fabricated audit counts rejects the row.

`terminal_physical_audit` is mandatory and binds the nonzero committed
`final_flux_revision`. It records residual and tolerance pairs for EOS,
continuity, energy, closed mass, and gauge closure. Every value must be finite,
every residual must be nonnegative, and every enabled residual must be no
larger than its positive paired tolerance. The historical pressure-continuity
contract may encode its disabled energy gate as `energy_tolerance:0`; the
finite nonnegative energy residual is then evidence only and is not an
acceptance limit. A continuity-energy-coupled row must use a positive energy
tolerance. Missing proof, zero final-flux identity, NaN/Inf, an over-limit
enabled metric, or a disabled coupled-energy gate rejects the complete row. In
the current v0.4 coupled product the normalized conservative-energy gate
intentionally uses the configured continuity tolerance.

V2 introduced mandatory per-corrector pressure evidence. `pressure_solve_calls` is
exactly two and `pressure` contains correctors one and two in order. Every
entry records status, termination, iterations, initial/final FP64 true
residual, recursive residual and supplemental convergence audit counters. A
cap, failed termination, contract-incompatible audit, or metric above its
limit rejects the row. New writers emit V5; immutable historical V1/V2/V3/V4
evidence artifacts retain their original schema tag and are never rewritten.
Consumers dispatch by the explicit schema tag; a newly constructed V5 runtime
record defaults to an invalid pressure contract until its producer supplies
the contract, terminal audit, advective and committed CFL certificates, and
refinement termination. V4 rows keep their original limiter and terminal-audit
shape; V5-only CFL or limiter fields do not retrofit a V4 row.

Required identity fields are `build`, `binary`, `case`, `stl`, `product`, and
`cpu_plan`. A zero `stl` value means that the case has no immersed surface.
Step identity is `step`, `time`, and the actually executed `bdf_order`.

Timing and memory fields are `launcher_ns`, `max_rank_step_ns`,
`max_rank_rss_bytes`, and `max_node_rss_bytes`. Communication fields are
`structured_messages`, `structured_bytes`, `ibm_messages`, `ibm_bytes`,
`blocking_collectives`, `nonblocking_collectives`, and `reduction_ns`.
Traffic counts are whole-job sums. Logical reduction calls and reduction time
are maximum-rank values; this both exposes the slow rank and makes the
committed record rank invariant.
Solver lifecycle fields are `linear_iterations`, `exact_numeric_refills`,
`coarse_numeric_refills`, `preconditioner_setups`,
`preconditioner_reuses`, and `heap_allocations`. `stages` contains stage id and
min/mean/max-rank nanoseconds without an evidence-only timer barrier.
The frozen stage-id set is exactly `10,12,15,20,30,40,45,50,60,70`; missing,
duplicate or extra stages reject a candidate runtime record.

`startup`, `retry`, and `restart_recovery` classify nonstationary steps. If any
is true, `statistics_eligible` must be false. Runtime evidence records facts;
Task 20 candidate manifests bind compiler, launcher, MPI, CPU/NUMA, affinity,
environment, command, input hashes, and immutable gate receipts.

The machine gate authority is `HUNDUN_V04_CANDIDATE_V1`. Its only legal order
is `focused -> full2 -> frozen -> full20 -> literature -> final`; every receipt
binds the unchanged canonical candidate SHA-256 and new, non-overwritable
evidence paths. `tools/v04_evidence_validate.py` rejects reordered gates,
candidate mutation, reused paths and artifact hash mismatch before replacing
the manifest.

Candidate identity stores absolute paths and SHA-256 values for the tests-off
product, tests-on product and test executable plus every flat case input. The
validator re-hashes those files on every manifest operation. It also validates
Git HEAD/tree object ids, compiler/linker/flags, positive product and CPU-plan
fingerprints, MPI/thread level, 64-rank process grid, CPU/NUMA/affinity facts,
an explicit environment allowlist, argv, output/checkpoint state and the
embedded `HUNDUN_V04_PERFORMANCE_POLICY_V1`.

Every accepted `HUNDUN_V04_GATE_RECEIPT_V1` carries start/end timestamps, an
unchanged candidate hash, new immutable evidence paths and gate-specific
checks. In particular, full2 requires a complete physics/work equivalence
receipt, full20 requires an `ACCEPT` paired-statistics result, and literature
requires a complete primary-data receipt plus physical-accuracy acceptance.
An accepted receipt containing any false check is invalid.

COAST/HUNDUN scientific-work equivalence uses
`HUNDUN_V04_COAST_EQUIVALENCE_V2`. The rule can be validated while the receipt
is `UNSEALED`; formal pairing additionally requires
`v04_evidence_validate.py equivalence --require-sealed`, which rejects a
changed rule, a non-`SEALED` status, or any unresolved capability and identity
field. Native SIMPLE outer iterations are not forced to equal HUNDUN's two
PISO corrections. Both products instead bind the same discrete problem, FP64
true-residual termination, common terminal audit and max-rank hot-step timer.

`HUNDUN_V04_PERFORMANCE_POLICY_V1` freezes `480x480x48/64 ranks`, full2 steps
1--2, full20 warmup steps 1--5, measured steps 6--20, five through at most nine
alternating pairs starting `HC`, maximum-rank hot-step-sum timing, disabled
serialized output, and the full2 directional ceiling. Paired input uses
`HUNDUN_V04_PAIRED_INPUT_V1`; the tool verifies the policy hash, alternating
order, exact measured-step count, positive samples and equality of each hot
total to its step-sample sum before bootstrapping pair ratios.

Literature data use `HUNDUN_V04_LITERATURE_RECEIPT_V1`. Its `complete` flag is
derived from the bound reference JSON rather than asserted by a caller.
`receipt-validate --require-complete` rejects any reference with
`pending_profiles`; this prevents a partial scalar receipt from authorizing a
long statistics run.

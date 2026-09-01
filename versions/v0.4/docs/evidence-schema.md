<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 runtime evidence schema

New writers use `HUNDUN_V04_EVIDENCE_V7`. Each line is one collectively agreed
committed-step record. Writers reject rank-disagreeing records and never
observe trial state. Historical V1--V6 artifacts and validators remain frozen;
V7 fields do not retrofit an older schema. V7 retains the V6 candidate,
time-anchor, CFL and owner-face AFC semantics, but raises the typed
same-target pressure--enthalpy refinement prefix bound from six to twelve.
This bound change is explicit in the schema identity, so a V6 validator cannot
mistake a legal seven-through-twelve-entry V7 prefix for V6 evidence.

V6 and V7 identify the exact running candidate with
`candidate_identity:{schema,head,tree,build_manifest_sha256,executable_sha256,identity_sha256}`.
V7 requires candidate-identity V2, whose canonical hash payload also contains
`evidence_schema=HUNDUN_V04_EVIDENCE_V7`; V6 retains candidate-identity V1.
Consequently, changing only a frozen V6 row's top-level schema string cannot
turn it into V7 evidence.
`head` and `tree` are exact Git object ids embedded at configure time. Git
HEAD/index/ref authorities and all v0.4 executable inputs are configure
dependencies, so an incremental build cannot silently retain identity from a
prior checkout or dirty source state. The build-manifest digest covers the
clean/dirty source state, compiler, flags,
build type and V0.4 feature switches without paths or timestamps. The
executable digest is computed from `/proc/self/exe`; every rank must produce
the same bytes. `identity_sha256` hashes a versioned canonical payload over
the four identities. The legacy numeric `build` and `binary` fields are the
first 64 bits of the build-manifest and executable digests and must agree with
the full identity. Mixed executable/build identities in one JSONL file are
rejected.

Every V6/V7 row carries a common
`run_start:{kind,previous_step,previous_time,restart_manifest_sha256}` and
records `previous_committed_time`. A fresh anchor is the canonical zero
step/time with a null manifest. A restart anchor binds the predecessor
step/time to the SHA-256 of the exact integrity-checked `manifest.bin` loaded
by `RestartReader`. The writer and validator both
require `time - previous_committed_time == advective_cfl.dt` within 128
binary64 epsilons, including the first row. A fresh first row must be step one
with previous time zero. A continuation first row must be the immediate
successor of the anchored snapshot and be marked `restart_recovery`;
subsequent rows preserve the same run-start anchor and bind the preceding
committed time. This closes the former first-record gap without assuming every
run starts at time zero or trusting two mutable time fields in the same row.
Accordingly, external validation of a restart JSONL requires
`runtime <evidence> --run-start-manifest <manifest.bin>`; the validator checks
the manifest format and FNV integrity, recomputes its SHA-256, and decodes the
authoritative predecessor step/time. Restart evidence without that frozen
manifest is rejected. Fresh evidence needs no sidecar and rejects one if
supplied.

The limiter policy is `scheme:"common_face_afc_v3_owner"`. A partition or
periodic face has one negative-side alpha authority; the positive-side copy
consumes the published value. Owner-only aggregation records
`active_correction_faces`, `limited_faces`, `minimum_face_alpha`,
`limited_face_fraction`, and the L1 `retained_correction_l1_ratio`. When no
nonzero correction face exists, `correction_metrics_applicability` is
`"not_applicable"`, both counts are zero, and the three ratio/alpha values are
JSON `null`; default numeric zero is not evidence. Applicable metrics require
positive active count, exact `limited/active` fraction, alpha in `[0,1]`, and
ratio in `(0,1]`. A limited row requires a positive retained ratio below one.

The provisional certificate is
`advective_cfl:{present,plan,time_revision_collective,density_view_collective,face_flux_view_collective,activity_collective,dt,out_max,abs_max,limit}`.
The three collective fingerprints bind every rank-local logical revision/view
identity while keeping the committed JSON row byte-identical across ranks.
The in-memory certificate still binds exact base, shape, stride, replica,
field, storage and revision-domain capabilities. The fingerprints are nonzero,
and the provisional face-view fingerprint must differ from the terminal
final-flux view fingerprint.

The terminal audit contains the full compact committed certificate:
`committed_convective_cfl:{valid,density_revision,final_flux_revision,density_view,face_flux_view,activity_collective,dt,out_max,abs_max,limit,out_winner,abs_winner}`.
The exact density/final-flux revisions must match the terminal audit. The view
objects carry communicator-wide fingerprints of all exact rank-local logical
views. Each deterministic winner binds global cell, rank, `rho*V`, outgoing
mass-flow sum, raw six-face absolute mass-flow sum, `dt`, and both CFL values.
The validator recomputes `Co_out=dt*outgoing/(rho*V)` and
`Co_abs=dt*0.5*absolute/(rho*V)`, verifies the advertised maximum, and rejects
an outward maximum above the shared configured limit (apart from the frozen
64-epsilon comparison slack). On an actual committed-CFL rejection, the same
winner data are also available to the CLI failure path.

One V6/V7 JSONL file freezes schema, full candidate identity,
`case,stl,product,cpu_plan`, momentum plan, IBM activity and CFL limit from its
first row. Steps must be contiguous and the time anchors must form an exact
committed chain. A zero STL requires zero activity; an IBM case requires
nonzero activity. Mixed schemas/candidates, skipped steps, time reversal,
changed static authority, aliased provisional/final views, invalid winners or
time--dt disagreement reject the file.

V5 remains supported only as an immutable historical schema. It retains
`common_face_afc_v2`, raw provisional revision numbers, the adjacent-time
check starting at the second row, and the scalar committed-CFL object. Frozen
V5 oracle data are not regenerated from V6 or V7. V6 likewise remains an
immutable historical schema with a maximum six-entry refinement prefix.

V4 retained the V3 pressure lifecycle and terminal physical acceptance, and
adds an explicit same-target nonlinear-refinement record. The ordinary
`pressure_solve_calls` remains exactly two: refinements are not extra PISO
correctors. `pressure_energy_refinement_solve_calls` is a prefix count from
zero through twelve in V7 and zero through six in frozen V4--V6 records.
`pressure_energy_refinement_termination` must be
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
limit rejects the row. New writers emit V7; immutable historical V1--V6
evidence artifacts retain their original schema tag and are never rewritten.
Consumers dispatch by the explicit schema tag; a newly constructed V7 runtime
record defaults to an invalid pressure contract until its producer supplies
the contract, terminal audit, advective and committed CFL certificates, and
refinement termination. V4 rows keep their original limiter and terminal-audit
shape; V5/V6/V7-only CFL or limiter fields do not retrofit a V4 row.

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

The separate machine-gate candidate manifest stores absolute paths and SHA-256 values for the tests-off
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

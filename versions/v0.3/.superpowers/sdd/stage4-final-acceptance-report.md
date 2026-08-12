# Stage 4 final acceptance report

```text
schema=hundun.stage4.final_acceptance.v1
decision=STAGE4_ACCEPT
version=0.3.0
accepted_code_head=6407cd7c591ce088db7f1dd7e296d77acd18da1c
accepted_code_parent=1674e60b192887a03f37ee03da3bad09052ad550
accepted_code_tree=2791a1cee7ac8114f1696670d30c8951212d6024
accepted_code_diff_sha256=fb60034a28a91d96236103442b463a8c35c8bfedb5bb4db63f98b5150c666c2c
acceptance_evidence_commit=60f40701be1200950e7d70e5f71f5e5ba09c4515
governance_seal=commit-containing-this-report
stage5_status=not-started-user-decision-required
```

## Decision basis

The tested code candidate C4 is clean, signed, and separated from later
governance-only commits. The compact exact-C4 matrix passed:

- Debug Stage 4: 45/45;
- focused Release: 15/15;
- focused ASan: 8/8;
- focused UBSan: 8/8;
- frozen Jammy/GCC11 Cantera backend, thermo, transport, chemistry interval,
  0D and PSR: 6/6;
- schema-v4 CLI validate/print-resolved, tests-off configure, source policy,
  package/RPATH contract and 1/2/4-rank compact MPI selectors: PASS.

The complete log and binary hashes are frozen in
`.superpowers/sdd/stage4-4A-3-acceptance-manifest.md`. No optional detached
screen was launched and no long-test wait occurred.

## Exact-head and governance audit

C4 is commit `6407cd7c591ce088db7f1dd7e296d77acd18da1c`, parent
`1674e60b192887a03f37ee03da3bad09052ad550`, tree
`2791a1cee7ac8114f1696670d30c8951212d6024`. The only commit between C4 and
this seal before the seal itself is `60f40701be1200950e7d70e5f71f5e5ba09c4515`,
which changes only the 4A-3 manifest and Stage 4 capability ledger. The seal
also changes documentation and `AGENTS.md` only. C4 remains the tested product
and test head.

Every one of the 26 commits after the execution handoff baseline and through
C4 contains exactly one `Signed-off-by:` trailer. The 4A-3 evidence commit is
also signed; the seal is committed with DCO. Cantera v3.2.0 provenance remains
tag `v3.2.0`, commit `4a8358eb80cfeb50474386b5f9ec0b3a83519889`,
zero local patches, artifact SHA-256
`093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760`,
and the tracked/transitive license inventory recorded by `UPSTREAM.json` and
`PREBUILT-LINUX-X86_64.json`. No mechanism is redistributed.

## Runtime and package audit

The host Debug/Release executable resolves all dependencies. Its backend-neutral
development RUNPATH contains the pinned local libc++ and configured OpenMPI
locations. In the frozen Ubuntu 22.04 rootfs the C4 Cantera test resolves
`libcantera_shared.so.3` and only the declared runtime closure; the artifact has
SONAME `libcantera_shared.so.3` and no embedded RPATH/RUNPATH. The formal
relocatable-package contract uses relative `$ORIGIN/../lib`; the moved package
consumer evidence was established in 4P-4 and its exact source contract is part
of the 45-test C4 matrix. The tests-off build registers zero tests.

At seal preparation the worktree was clean, the evidence hash audit matched,
and no process executable resolved to an exact-path `hundun` binary in either
the Stage 4 worktree or product repository.

## Accepted capabilities and limits

Accepted Stage 4 capabilities are:

- all-species conservative `rhoY`, total thermochemical enthalpy and typed
  schema-v4 composition/mechanism identity;
- backend-neutral services plus the audited Cantera 3.2.0 runtime, thermo,
  mixture-averaged transport and stiff cell-chemistry interval;
- fixed symmetric C-T-C with exactly two PISO correctors, shared final face
  mass flux, correction-velocity species diffusion and atomic rollback;
- open, closed and partially closed `p0`, body-fitted boundaries, IBM and WALE
  composition, schema-v4 root dispatch;
- Checkpoint v4 validate-then-publish sections and read-only reacting
  diagnostics;
- no-Python normal configure/build/install/test/runtime and relocatable
  Cantera package boundary for the audited Ubuntu 22.04/GCC11 profile.

This is not a generic Linux claim. It does not accept a real-fuel mechanism,
flame/ignition validation, 48-cubed or 96-cubed evidence, TPDF/TCR, ESF, spray,
AMR, moving IBM, rank-changing restart, NativeChemistryBackend, GPU production,
Python integration, or private-code equivalence. Stage 5 remains stopped at
the user decision gate.


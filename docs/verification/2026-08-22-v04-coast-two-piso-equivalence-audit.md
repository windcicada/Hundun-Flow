<!-- SPDX-License-Identifier: Apache-2.0 -->

# COAST SIMPLE / HUNDUN two-PISO equivalence audit

Date: 2026-08-22

Decision: `RULE_FROZEN / RECEIPT_UNSEALED / FORMAL_PAIRING_FORBIDDEN`

This is the root-agent scientific-work decision for the pressure-continuity
architecture at HUNDUN product commit
`393e02287d24020635d4c0892265bc890ee9a683`, product tree
`7a28c3dc62dcb137f08e0c4ee572feb8d3e46726`. It does not freeze a performance
candidate or accept Task 20.

## Comparison unit

The unit remains one complete accepted physical step of the same discrete
problem. HUNDUN executes exactly two PISO pressure corrections. COAST executes
exactly two native SIMPLE outer iterations. Native SIMPLE/PISO, linear-solver
or residual-audit counts are not forced to be equal. Every registered native
operation is charged to its product's max-rank hot-step timer.

HUNDUN's PISO2 supplemental continuity audit is scientifically admissible
under that rule because it:

- runs only after the canonical FP64 `||b-Ax||_2` pressure criterion passes;
- can only continue the same FGMRES solve within the frozen 500-iteration cap;
- cannot add a PISO correction, relax or replace the common L2 pressure rule,
  publish a rejected trial, or make a failed step statistics-eligible;
- uses the active masks and the byte-identical
  `hf_coast_common_terminal_cell_v1` continuity cell formula;
- reports audit count, rejection count, final metric and limit in
  `HUNDUN_V04_EVIDENCE_V2`; and
- leaves every extra Krylov, MG, halo, collective and audit operation inside
  HUNDUN's measured interval.

COAST need not reproduce this native audit cadence. It must complete its
sealed two-SIMPLE schedule and both products must still pass the independent
common post-step terminal audit before accepted-state commit. The in-solver
HUNDUN evaluation does not replace that common terminal authority.

## Receipt disposition

The receipt formerly sealed to rejected candidate
`h53b1075-c3c22e0f-r5` remains recoverable from Git history and from that
immutable candidate's artifacts, but cannot authorize any new HUNDUN binary.
The current machine-readable template is therefore deliberately reset to
`UNSEALED`; all new-candidate HUNDUN identity and direction-screen fields are
null and `formal_pairing_authorized=false`.

- Rule document SHA-256:
  `7d7a12ebee60b908fdf7d6e16c6db98a3267db1dbf6ca645f0d12bffb1022b10`
- Unsealed template file SHA-256:
  `b42c146a4948ecadeb1cf50c22ad686bc56e72e3205e2441b7cfbc1c1a571cff`
- Validator canonical digest:
  `81d75f98b54f080400e5a686ddec25b49ec06ffea5804622d2d61b5a6e33e773`

The ordinary equivalence validator accepts the unsealed rule. The same
validator with `--require-sealed` rejects it with `COAST equivalence receipt
is not SEALED`, which is the required fail-closed state before a new immutable
candidate exists.

## COAST and Stage 5 identities

The nonreacting COAST equivalence source remains the clean isolated worktree
at commit `3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3`, tree
`ba449790e4918e7a3fd5c21e71c4e9f980a4691f`. No COAST file was edited or run
for this audit.

Stage 5 remains outside this lane. Read-only governance inspection found:

- historical Stage 5 framework acceptance: tested code
  `41b2aac97d28da3949a6a4bd629c079d7aa6a8b7`, tree
  `7456632c41ae749b66f2fc5b9325587782e6781d`;
- sealed framework branch head
  `02b57cce45311bfd6c1507b6e6b32eb0617a59cd`, tree
  `32e72d14083f91576c0d8f504529a68f9baae6a9`; and
- later field-validation branch head
  `8ffdf2b5673374fb14639fc4dce09a1b586ee5db`, tree
  `e567af1e7682edc09aaf3f6b543dfdb5242de117`.

Those Stage 5 receipts distinguish framework acceptance from still-pending
field-science/performance claims. This lane neither reinterprets nor waits on
those claims. It imports no combustion code and records only the stable
Cartesian compatibility seam:
`FrozenExecutionGraph/ContributionPlan`,
`C1 -> transport -> PISO1 -> C2 -> PISO2`.

## Next authorization boundary

This decision closes the rule question but not the candidate receipt. Pressure
performance work may now proceed under a frozen diagnostic matrix. Only after
the product path passes focused and bounded full-grid screening may the root
agent create an immutable candidate, populate the null identity fields,
revalidate COAST capability bytes, and reseal the receipt. Formal pairing is
forbidden until that reseal succeeds.

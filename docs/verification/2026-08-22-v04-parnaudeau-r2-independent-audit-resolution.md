<!-- SPDX-License-Identifier: Apache-2.0 -->

# Parnaudeau R2 independent-audit resolution

The R2 development extraction from commit `21996c5` is rejected and must not
be attached to a literature receipt. This decision was made without reading
HUNDUN long-statistics output.

Independent audit A passed marker identity, sequence grouping, mapped values,
render presence and reproducibility. Its report is
`/tmp/parnaudeau-h21996c5-audit-a.md`, SHA-256
`1d5951ee3281029faa78d4719fd8390860b12a0b4d64d14e84212e8cffda89a1`.
Independent audit B repeated those checks but returned the controlling `FAIL`.
Its report is `/tmp/parnaudeau-h21996c5-audit-b.md`, SHA-256
`3f2b51a75e286e0b3baa363141903617f0769c459c7703961ebfbc5ea675ff07`.

Audit B found two method defects:

1. R2 fitted panel-box-derived evenly spaced coordinates rather than actual
   printed SVG tick centres, so its near-zero residual was synthetic.
2. R2 compared two deterministic conversions but had no independent frozen
   clean-source cardinality/hash. A common-mode marker insertion, deletion,
   style change or sequence-boundary shift could therefore pass both copies.

R3 resolves both defects before a replacement authoritative extraction:

- every panel parses and binds its actual SVG axis/frame path and all major
  tick centres in both passes;
- the two tick inventories must be exactly equal;
- clean-source per-sequence counts and canonical marker-inventory SHA-256 are
  mandatory in both passes;
- clean extraction plus duplicate, style, deletion, boundary, ordering and
  axis-tick-shift mutations form the required fail-closed self-test matrix.

Root review of the development implementation produced nonzero tick-fit
residuals from `1.9304729274638444e-05` through
`9.700656854216838e-05`. All six mutations were rejected, both tick and
marker inventories were equal, and all 15 profiles were ordered over a common
interval. Those `/tmp` runs are development evidence only. The authoritative
trace must be generated from a later immutable R3 commit and stored at a new,
non-overwriting external evidence path.

That authoritative R3 extraction was subsequently generated from exact commit
`1adbc9979ae47d76dba4e5b894ff4ae8c04e54e7`, tree
`f06af27cedbf7e7ae5415b836a3280dd98ba7faf`, and committed extractor SHA-256
`c73fb4eb663c450f2b288edbe23d7c4a4233c94f13570ebb790834a986c3dc25`.
Its immutable external directory is
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/parnaudeau-b8a775e5-h1adbc99-r3`.
The bound source, calibration, raw trace and embedded-profile payload hashes
are respectively
`b8a775e5a5078e19fc47d9c5f47e95b81b4a68d4449e4444459088ed9befcdd4`,
`a3cb287d34033808f3cd7f77594b2b9f77cf8a27289ad4a45e4a6eaec6c475ee`,
`7c017aeff21ca05a11dd3f102bd1d5a132353af872af2bc20e1800fe387ce048`
and
`a96b47d9f1b1823950db65af8de760ddd2182a2c3f8d4b66db0e5b45266413b9`.

The main agent reran the complete extraction, all six mutations, JSON/hash
checks and render review. An independent Luna audit reparsed both SVG passes
without importing the extractor, inspected all ten renders and returned
`PASS`; its report is
`/tmp/parnaudeau-r3-immutable-independent-audit.md`, SHA-256
`aa08336f9fd21640e36aadc21ebf92b77f98ac35761d6f0e2d4944f530a2f243`.
The 15 profiles and their controlled-digitization bounds are therefore
accepted as the Parnaudeau profile authority and embedded in the literature
reference. This acceptance is limited to that profile matrix.

This resolution does not complete the literature receipt, authorize full20 or
long statistics, choose a physical-accuracy tolerance, or decide the final
gate. Direct total-drag and finite-span whole-force lift-RMS authority remain
separate fail-closed requirements.

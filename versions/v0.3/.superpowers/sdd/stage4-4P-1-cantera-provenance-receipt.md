# Stage 4 4P-1 Cantera provenance receipt

- Date: 2026-08-11 (Asia/Shanghai)
- Result: `4P-1_PASS`
- Scope: source, dependency, license and binary-input identity only
- Product capability accepted: false

The official Cantera release asset was re-hashed locally as
`a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b`.
Its `v3.2.0` tag resolves to commit
`4a8358eb80cfeb50474386b5f9ec0b3a83519889`. The candidate shared object was
re-hashed as
`093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760`.

The verified fmt, yaml-cpp, SUNDIALS and Eigen archives also matched the hashes
in `third_party/cantera/UPSTREAM.json`. Boost and every builder-only tool remain
excluded from the installed/runtime closure. No local source patch exists.

Legal review retained the Cantera non-endorsement condition and keeps all
dependency notices distinct. Mechanisms and datasets are explicitly outside
the Cantera code-license claim; none is bundled by 4P-1.

Mutation evidence:

- changed Cantera archive hash: rejected;
- omitted dependency license: rejected;
- unapproved Cantera-named product source: rejected.

This receipt does not accept the CMake consumer, relocation, backend behavior,
or reacting-flow science. Those remain owned by 4P-2 through 4A-4.

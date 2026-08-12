# Stage 4 4A-3 low-cost acceptance manifest

```text
schema=hundun.stage4.acceptance_manifest.v1
stage=4
result=PASS
accepted_code_head=6407cd7c591ce088db7f1dd7e296d77acd18da1c
accepted_code_parent=1674e60b192887a03f37ee03da3bad09052ad550
accepted_code_tree=2791a1cee7ac8114f1696670d30c8951212d6024
candidate_diff_sha256=fb60034a28a91d96236103442b463a8c35c8bfedb5bb4db63f98b5150c666c2c
version=0.3.0
worktree=clean
detached_stage4_run=not_started
forbidden_large_tests=not_registered_or_run
```

The tested code candidate is `C4 = 6407cd7c591ce088db7f1dd7e296d77acd18da1c`.
All product and test changes stop at C4. This manifest and the accompanying
ledger update are governance-only and do not relabel their commit as the tested
code head.

## Exact low-cost matrix

The terminal matrix ran from
`/home/wyf/code_dev/.hundun-flow-stage4-acceptance/6407cd7-run2`. Its
`results.tsv` SHA-256 is
`c809e6653ca84efd5b7537ce09ca1af2c679ea452df6aec4c4c4f4f47a3c2740`.

| Selector | Result | Log SHA-256 |
| --- | --- | --- |
| Debug `ctest -L stage4`, 45 tests | 45/45 PASS | `adb6a54a320d06d6393f0f379612160f6c5d21512b336a4c0ffce798a0110d1e` |
| focused Release, 15 tests | 15/15 PASS | `d82438b07b669aa63236adf43a740e4ac20a15df53456de72748ea432ea7d62b` |
| focused ASan, 8 tests | 8/8 PASS | `f253ddedf79d6faf2ba05ce56bc70200e4702186eff176faab508f22c18da663` |
| focused UBSan, 8 tests | 8/8 PASS | `259c1855967e48e0c5a945e86f2f534b4bb67160bd8224bd26d4abc077fa25c2` |
| schema-v4 `--validate` | PASS | `de545cc7e7ff8eaacefae8fcb4a954555a1224f24fe7ce8ff62b9874c54f3fe2` |
| schema-v4 `--print-resolved` | PASS | `5e5f335ca86df5b4cf11228eb88b3b0ad6a5fdeea775a6436e89b93a0cec48cc` |
| Cantera runtime/thermo/transport/interval/0D/PSR | 6/6 PASS | empty-success logs, SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

The Debug selector includes source-policy and mutation fixtures, schema-v4
typed broadcast, source transactions and rollback, C-T-C/MMS, body-fitted and
IBM+WALE paths, closed pressure at 1/2/4 ranks, Checkpoint v4, diagnostics,
Cantera builder and package/RPATH contracts. No 48-cubed or 96-cubed case,
Flame D, TPDF, spray, or sanitizer-large MPI selector was registered or run.

## Binary, configuration and external identities

| Artifact | SHA-256 |
| --- | --- |
| Debug `hundun` | `c78fe9fcd96daf241bd9919ed2af3a79acb3630c98f35de2154dcd8322d02cb5` |
| Release `hundun` | `962b39546c11e1e3cf9f26ad2e5ab79a320ad08f0bb019942901c67f5a42d59f` |
| Cantera backend test | `79ee397526338266e8d6ff2cb213501e854311e66eb42b2b1849459741a603a0` |
| Cantera thermo test | `c885929d9e446400ab739c5ec5221e46c3bb556ced1084561274f66a781a7124` |
| Cantera transport test | `6c8071276c96e98463fbe5e27f8b428051fd60b081df2be8208cbd15daae0a1b` |
| Cantera chemistry interval test | `8a0c217d113f09beaed5da3489c7df7c81fc4d7102e1a88ce77aaec42ceba98d` |
| 0D test | `ea546e441a83114c311ef5532b2033011eeb57eaf1900228d9f02a52a1b7f7a9` |
| PSR test | `10eb3b7a64e44ebc4e7f67d9eaee9b06de4534f241ef8984342eaaadce2e651e` |

The five CMake cache hashes (Debug, Release, ASan, UBSan, tests-off) are,
respectively,
`0bf44d0ebe32fb9bffeeb9ba88676dfe1e1da16fc1df754e3d71e6a7f31a68e7`,
`079b9ce1e56756749764beec2684530cd7d43ffd0395f0e62933ddb36b8ee0d9`,
`a36e3394fcfe900678cf0663eb2f54bf9eaca624a1ca0275db3dcb99f78eeca2`,
`a4b0e6d74a4b42cf573b8801d59817959e451dfe98c6e3cb1dc5bccd582a04f4`,
and `c7790fb3b6e8f4c02db83e2029661dcd616de0a9bdaa4c8522e85377bac1b26d`.

Cantera tests were rebuilt from C4 source and run with GCC 11 in the frozen
Ubuntu 22.04 rootfs. The package artifact SHA-256 remains
`093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760`;
the independently authored mechanism SHA-256 remains
`c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee`.
Host backend-neutral builds used Clang 15.0.6 with libc++; configured MPI
reports MPI standard 3.1 and OpenRTE 2.1.1.

## Scope and decision

The accepted evidence supports the Stage 4 reacting-flow scope only. It does
not claim a real-fuel mechanism, flame validation, TPDF/TCR, spray, AMR,
rank-changing restart, Python integration, or generic Linux portability. No
optional detached 24-cubed diagnostic was launched. The exact-head governance
seal remains Task 4A-4.


# S3-G1 acceptance inventory 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T23:02:43+08:00`

accepted parent：`e9b290a289fcd716d653d54bfad077a581a94dd2`

accepted commit：`cdbc7947ac33cf6f2cebb29d817f8430463980ef`

## Capability / inventory / manifest authority

capability ledger 精确列出九个 implemented-and-accepted profiles，以及 redistribution、
multigrid、AMR、GPU、moving bodies、rank-changing Restart 和 96-cubed 的 deferred/out-of-scope
边界。每行绑定 task/commit、acceptance test 与 final owner。

canonical TSV 固定 11 列和 57 个 rows，覆盖 low-cost、scientific、performance、sanitizer、
governance 五组；formal scientific/performance 全部绑定 Release producer，ASan/UBSan 各五个
冻结小规模 row。launcher 的 `--list` 不依赖 build root；group 模式要求 exact candidate 与
五个 absolute roots，按 build role 选择，逐 row 写 terminal manifest 并继续独立 rows。

manifest schema 强制 HEAD/tree/diff、root/cache/binary identity、compiler/libc++/MPI、argv/env/
cpuset/ranks、start/end/exit/duration/RSS、log/artifact SHA；missing exit/log/binary 与 duplicate
row mutation 均拒绝。projection contract 读取既有 product manifest，并拒绝 tests、
`.superpowers`、private/token path 与未登记 relation mutation。

## Verification

- `test_stage3_evidence_manifest` build：PASS；
- acceptance/ledger/manifest/projection contracts：4/4 PASS；
- `stage3_acceptance.sh --list`：57 data rows + header，PASS；
- shell syntax 与 `git diff --check`：PASS；
- 未调用 scientific、performance 或 sanitizer group。

内容 SHA-256：ledger `977e03373c4389d3fe0ddf0139fdd9d5b7d2c0202275b5a4fd6eb9763240b894`；
launcher `4212047658fb1da91dce89e346a0c7967fa79ad6f6b7030b22bf1e8df78b574c`；
inventory `8a7153e3bc0da5da240b0c23bad6cb1796d542ec8b2decd531f6dd74065f02bd`；
manifest header `d47d1bf7536702d9592263ddd26884ed79b9d98ed57cc9485cce883dde5edc17`；
manifest unit `6e9d95312205f6636a881753dd4950dd46806a11defb2a3aa7c6ba4f3c03a0d3`。

未运行 24/48/96-cubed；未访问私有源码、研究数据或研究进程；未 push/publish。

提交 subject：`test: define Stage 3 acceptance inventory`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

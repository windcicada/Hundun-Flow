# S3-P0 registration shards 实施回执

状态：`TASK_GATE_ACCEPTED`

accepted parent：`6d95fb83735383894c9848fd286ff424e6a18d4d`

activation commit：`2a6268d28666a65c176af0028e6c370ba12df85e`

accepted_at：`2026-08-09T12:42:06+08:00`

## 实现边界

- `tests/CMakeLists.txt` 只新增两条 P0 contract 注册和五个 registration fragment
  include；没有移动、重命名或改写任何既有测试；
- 五个 fragment 分别冻结 science、Checkpoint、diagnostics、framework 和 acceptance
  的独占写入者，初始仅含 SPDX、`include_guard(GLOBAL)` 和 ownership 注释；
- validator 要求五个 exact include 各出现一次，并验证 fragment 存在、SPDX 和 include
  guard；
- mutation fixture 只把 `tests/CMakeLists.txt` 与五个 fragment 复制到 build-tree
  sandbox，删除 acceptance include 后要求 validator nonzero；真实 source tree 不被修改；
- infrastructure lane 只冻结创建时点和命令，未提前创建；它仍必须等 S3-C1 签署接受。

## executable RED

先只注册 validator/fixture，不添加五个 fragment/include。`cmake --preset debug` 成功，
随后：

```text
ctest --test-dir build/debug -R '^test_stage3_registration_contract$' --output-on-failure
exit=8
missing stage3 registration include: stage3_science_registration.cmake
```

首次运行同时暴露 validator 的 script-mode `while(TRUE)` policy warning；该测试缺陷先被
修成 policy-independent 计数循环。第二次 RED 无 warning，仍以同一预期缺失 include
原因失败，之后才进入 GREEN。

## 授权的基线修复

P0 首轮 GREEN 中四项通过，既有 `source_layout_fixture` 因 accepted Task 12/13 文件
`les_wale.hpp/les_wale.cpp` 未进入旧 prefix allowlist 而失败。只读历史证明：

- fixture 最后修改于 `15001fbbf2acba3c61adae54f829bff01266ba74`；
- `les_wale.hpp` 新增于 accepted Task 12 commit
  `43d8446a4de354c40749b94d1f2143a8d744a4e9`；
- 全树不匹配项恰好只有两个 `les_*` 文件。

用户以“授权，请继续”明确授权把
`tests/cmake/source_layout_fixture.cmake` 窄增到 P0 白名单。实际修复只在 prefix alternation
中增加 `les`，未改变其他 layout 规则。

## GREEN 证据

命令：

```text
cmake --preset debug
cmake --build build/debug -j32 --target test_stage3_flow_header_contract
ctest --test-dir build/debug -R '^(test_stage3_registration_contract(_mutation)?|source_layout_fixture|cmake_include_authority_fixture|test_stage3_flow_header_contract)$' --output-on-failure
```

结果：configure/build exit 0；CTest `5/5 PASS`，0 failures。未运行 MPI、数值、Release、
ASan、UBSan 或长 selector。

证据身份：

- implementation/test diff SHA-256：
  `0f1873ac56e4f28759e28765f1c9008940c83f747cbd15dd3b4d51149713cbf9`；
- Debug `test_stage3_flow_header_contract` SHA-256：
  `896192904ce4ec25d77ef5c2335a50ad758b85fa00e5f0ac11e2dea3dcf4fac1`；
- final focused `LastTest.log` SHA-256：
  `c685fa39562f92ef50a0e1ddd2d56e88e0dbd94c1c5dfcfe56df9574858f1fb4`；
- CMake：`3.31.12`；
- codegraphf sync：exit 0；C++ product symbol/caller impact：none，build/test registration
  impact only。

## 主 agent review

Requirements：PASS。P0 registration ownership、RED semantics、mutation isolation、timeout/labels
和 C1 后建 lane 的时点均符合 v2 plan。

Quality/caller review：PASS。validator 计数输入非空且确定，fixture 删除范围受 build-root
检查约束，fragment include guard 不改变既有 build graph；没有 product source、public API、
schema、Restart、diagnostic enum、数值 authority 或 MPI 行为变化。

版权/运行边界：未采用外部算法、未访问私有源码或研究数据、未检查或干预研究进程、
未 push/publish。提交前未发现 HUNDUN-FLOW-owned build/test/MPI 后台进程。

提交 subject：`build: shard Stage 3 test registration`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

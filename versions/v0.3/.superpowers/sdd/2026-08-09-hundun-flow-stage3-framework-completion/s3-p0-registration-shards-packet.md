# S3-P0 registration shards 冻结执行包

状态：`TASK_GATE_ACCEPTED`

accepted parent：`6d95fb83735383894c9848fd286ff424e6a18d4d`

activation commit：`2a6268d28666a65c176af0028e6c370ba12df85e`

## 目标

只建立五个后续 Stage 3 测试 registration fragment、可执行 registration contract
和 build-tree-only mutation fixture。不得移动既有测试、修改产品代码、运行 MPI/数值
selector，或提前创建 C1 后才允许建立的 infrastructure worktree。

## 允许文件

- `tests/CMakeLists.txt`
- `tests/cmake/stage3_science_registration.cmake`
- `tests/cmake/stage3_checkpoint_registration.cmake`
- `tests/cmake/stage3_diagnostics_registration.cmake`
- `tests/cmake/stage3_framework_registration.cmake`
- `tests/cmake/stage3_acceptance_registration.cmake`
- `tests/cmake/stage3_registration_contract.cmake`
- `tests/cmake/stage3_registration_contract_fixture.cmake`
- `tests/cmake/source_layout_fixture.cmake`（用户窄授权：只补 accepted `les_*` prefix）
- 本 packet、对应 receipt 和 ledger status row

## mutation-sensitive contract

- 任一 exact include 缺失：真实 contract 失败并报告
  `missing stage3 registration include`；
- 任一 include 重复：contract 失败；
- fragment 缺失、SPDX 不匹配或 `include_guard(GLOBAL)` 不是恰好一次：contract 失败；
- mutation fixture 只复制 `tests/CMakeLists.txt` 和五个 fragments 到 build tree，删除
  一项 include 后要求 validator nonzero；它不得修改 source tree；
- 未 mutation 的真实 tree 必须通过。

## RED / GREEN

RED 先只注册 validator/fixture，不创建 fragments、不添加 includes；目标失败原因必须是
缺少 registration include，而不是 configure、syntax 或 target failure。

GREEN 仅执行 v2 P0 Step 3 的 registration/layout/header contract 矩阵，不运行 MPI、
数值、Release、sanitizer 或长测试。infrastructure lane 的创建命令只被冻结到 receipt，
必须等 S3-C1 签署接受后执行。

## 接受检查

- [x] executable RED 以预期原因失败；
- [x] 五个 include 各恰好一次；
- [x] 五个 fragment 存在、具有 SPDX 和 include guard；
- [x] mutation fixture 在 build-tree copy 上杀死缺失 include；
- [x] source layout、include authority 和 Stage 3 flow header contract 通过；
- [x] 完整 diff 仅含允许文件；
- [x] `git diff --check`、DCO、clean tree 和后台进程检查通过；
- [x] 不创建 infrastructure lane，不运行 MPI/数值，不访问私有路径，不 push/publish。

提交 subject：`build: shard Stage 3 test registration`

签署 identity：`WANG YUDONG <wangyudong@buaa.edu.cn>`

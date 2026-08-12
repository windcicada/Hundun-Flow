# S3-O1 WALE diagnostics 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T19:21:04+08:00`

## 委派、scope repair 与集成身份

- worker worktree：`/home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure`；
- worker baseline HEAD：`c2a3b71edf1fe5267e9a66642e191d4453bacb04`；
- signed handoff commit：`8f927a0afb9b941a7c96309c7166e57e81156ac6`；
- main integration parent：`efcdd25c80639816f14f25f1fa933adca1ee93a9`；
- main integration commit：`463d0486811d3120aea954006f6e71581c79fb98`；
- handoff staged diff SHA-256：
  `5a463bf2162d1938bcbdc7fb795753c1ccc762404d2b1dded28bd573232cff37`。

worker 严格只修改 O1 packet 的七个路径。Step 0 发现权威 Files list 的闭环缺口：新增 enum
18--22 后，schema-v1 的中央 `module_kind_name()` 仍只接受 0--17，而
`src/diag_structured.cpp` 不在 worker 白名单。绕开通用 validate/canonical 会产生残缺产品
路径，因此 worker 按 stop gate 保持 clean。主 agent 依据其跨模块组合职责，只在该中央
switch 追加五个名称 case；未改变 schema、旧 0--17 或 validator。worker 随后继续原七文件
TDD。最终八路径 handoff 中，中央文件是显式 main-owned integration repair，并有五名称精确
canonical regression。

## Enum、provider 与只读边界

`DiagnosticModuleKind::performance=17` 后只追加：

- `immersed_surface=18`；
- `ghost_stencil=19`；
- `local_flow_pattern=20`；
- `wall_force=21`；
- `les=22`。

旧 0--17 与 diagnostics schema v1 不重编号。O1 只实现 kind 22 的直接
`const les::WaleSummary&` provider；profile inventory 与 kinds 18--21 providers 仍归 O2。
不存在 optional/pointer/fake-summary overload；调用方没有 WaleSummary 时注册零 descriptor、
提交零 record。

descriptor 固定为 schema v1、kind `les`、module `hundun.les.wale`、instance `primary`，
只声明 summary/counters。summary 以严格字典序输出 `nu-t.l2/maximum/minimum` 三个
`m2/s` metric；counters 以严格字典序输出 `identity`、`nu-t.exact-zero-count`、
`owned-active-count` 三个精确 uint64 `count`。六字段共同进入 fingerprint；uint64 拆成
高低 32-bit 分量，避免 2^53 以上的 FP64 舍入失敏。六个逐字段 mutation 均改变 fingerprint。

nonfinite `nu_t` 直接复用 `DiagnosticValueStatus`，保留 positive/negative infinity 与 NaN
raw bits，不伪造 unavailable。重复 collection 的 canonical JSON 完全一致，且 WaleSummary
六字段 bitwise/exact 不变；不修改业务 counters、不 collective、不复制全场、不做第二次数值
计算。sink exception 只映射为既有 sink-failure contract。

公开参考只采用 AMReX pinned support/presence 分离和 OpenFOAM pinned WALE summary 责任边界；
未复制、翻译或改写 GPL 源码、命名、ABI、field class 或控制流。

## TDD 与 GREEN 证据

初次 build 因旧 Make graph 未 regeneration 而没有新 target，exit 2；该结果被明确拒绝为
RED 证据。正常 CMake regeneration 后，五个 target 均编译链接，CTest 3/4 PASS、exit 8；
唯一失败是 `test_wale_diagnostics.cpp` 对临时
`wale.diagnostics.capability` seam 的断言，属于预期可执行行为 RED。

- worker 最终 focused：4/4 PASS，real 0.06 s；
- worker tests-off Clang/libc++ Release `hundun` 与 provider symbols：PASS；
- 主 agent 隔离 worktree独立 focused：4/4 PASS，real 0.06 s；
- 主 agent隔离 worktree tests-off linkage/symbols：PASS；
- 集成 main R1/O1 组合 gate：8/8 PASS，real 30.61 s；O1 四行分别为
  0.01/0.01/0.00/0.03 s；
- 集成 main registration/layout：5/5 PASS，real 0.13 s；
- 集成 main Clang 15/libc++ Release tests-off `hundun` 与两个 WALE provider symbols：PASS。

Debug/tests-off binary SHA-256：

- `test_wale_diagnostics`：
  `5f182ee03e2f3eee609d3b5876942867f77c20997ecda5ce038ca49cbeb5024a`；
- `test_wale_diagnostics_header_contract`：
  `f4dd624d8add61c91e8d52b2789d7092899fca1600be1fb703d864a9b555b9fe`；
- `test_structured_diagnostics`：
  `14c4195f688912a686c48f8ba6fa10ca6b6e8820334a7173dc64960263cffa38`；
- Release tests-off `hundun`：
  `93bd599df0fb73ffd1b3555442ca094bbc297e11c3026ccca2b36964467e32ad`。

未运行 24/48/96-cubed、MPI 数值 case、sanitizer 或正式矩阵；未访问私有源码、研究数据
或研究进程，未 push/publish，未修改 flow orchestration、Checkpoint、density/closure、
Task 11 authority 或科学阈值。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: add WALE diagnostics`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

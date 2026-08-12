# S3-O2 Stage 3 providers 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T21:33:03+08:00`

accepted parent：`8f869aafd299617a36876c3380768b22fa8af2ce`

accepted commit：`041419480221f94af066172a6f7b2a697255d2e3`

## Inventory 与 provider authority

profile-driven inventory 已固定：body-fitted/WALE 仅 kind 22；IBM/none 为 kinds 18--21；
IBM/WALE 为 kinds 18--22。不存在的 module 不注册 descriptor、instance、counter 或伪造
zero record。Checkpoint report 与 presence overload 共享同一 profile truth table。

`ImmersedStaticDiagnosticSummary` 是不借用 numerical object 的 value-only snapshot：

- kind 18 从 accepted surface/query/domain immutable authority 输出 vertex/triangle/component、
  bbox、area、closed volume、area-vector closure、orientation 与 fingerprints；
- kind 19 从 Ghost/Wall public immutable plan 输出 links/donors、QR rank/condition、Halo reach、
  wall points、triangle coverage 与 fingerprints；合法 rank-local zero-work 分片保持 available；
- kind 20 不重算 LFP。`ImmersedOperatorAdapter` 在实际 row construction 中发布 row
  fingerprint、replacement groups/occurrences、dimensionless coefficient norm 与 limiting-case
  status；accepted `ImmersedFlowDiagnosticSource` 只复制该 report，再值化进 summary；
- kind 21 从 accepted report 与同 attempt `WallForceSample` 输出四组 pressure/viscous/total
  force vectors、moment、area closure、point count 与 lowest-rank evidence；
- kind 22 原样复用 O1 `WaleSummary` provider。

snapshot seal 覆盖 density variant、force、moment/area、LFP report 与 optional WALE；failed
attempt 的 kinds 20/21 source 为 unavailable。local 重复 collection canonical identical；
collective surface 检查 replicated authority，ghost/LFP 聚合 rank-local counts/norm/fingerprints，
wall force 独立验证全局 force/sample authority；sink failure 也 collective 收敛。

## Main-owned scope repair

冻结 field contract 要求 Halo reach=`cell`、moment=`N*m`，而旧 schema-v1 unit vocabulary
未包含这两个合法单位；主 agent仅向中央 unit validator additive 追加它们，未重编号 schema、
module kind 或旧单位。kind 20 又要求真实 operator row/coefficient authority，而原 public report
只有三个 counts；主 agent在既有 `ImmersedOperatorReport` additive 追加四个 value fields，并在
既有 row construction 处计算，避免用 ghost donor 近似产生第二 authority。对应修改
`src/diag_structured.cpp`、`include/hundun/fvm_immersed_operator.hpp` 与
`src/fvm_immersed_operator.cpp` 是完成冻结 contract 必需的跨模块 integration repair。

## TDD 与 GREEN 证据

- compile seam 后 combined RED：inventory 首 profile mismatch；MPI provider 稳定 capability
  failure；
- wall-force RED：accepted attempt 命中未实现 capability；
- empty local partition RED：合法 zero-work summary 被旧 preflight 拒绝；
- focused header/unit gate：9/9 PASS；
- Stage 3 diagnostics MPI：1/2-rank 2/2 PASS，real 16.31 s；
- impacted ideal/material/WALE/legacy immersed diagnostics：8/8 PASS，real 211.53 s；
- immersed operator header + uniform 1/2-rank：3/3 PASS，real 26.04 s；
- Clang 15/libc++ Release tests-off `hundun`：PASS；
- `test_stage3_diagnostics_mpi` SHA-256：
  `21fd253e2845b6fd826511a16f2d1e5ebf78c6e934e91bee644f0e0d5dfe8697`；
- `test_stage3_provider_inventory` SHA-256：
  `b1b3bb528144123f71f547ddfaa61e9241fe5927f99397eee9c06ebcf80ed5ab`；
- tests-off `hundun` SHA-256：
  `e2fe6e46d1376c2eca6aa80adf9b7d92254a59be3190064492f778d3febbf589`。

未运行 24/48/96-cubed 或正式矩阵；未访问私有源码、研究数据或研究进程，未 push/publish。
没有重算 pressure/force/WALE，没有修改 Task 11 科学阈值、两次 PISO、force sign、rollback、
Restart 或 MPI 一致性要求。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: complete Stage 3 diagnostics`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

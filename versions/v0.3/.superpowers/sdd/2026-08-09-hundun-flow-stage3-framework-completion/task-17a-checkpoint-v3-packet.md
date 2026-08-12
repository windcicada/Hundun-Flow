# Task 17A：Checkpoint v3 IBM-only 冻结证据包

状态：`TASK_GATE_ACCEPTED`

基线：`17b8434413a1dcc1c88070351f9669c6f6f7ed83`

## 科学目标

为 constant-density、static LFP-GCIBM、LES none 的已提交 19A 产品路径提供
相同分区 bitwise continuation。Checkpoint v3 必须同时恢复：

- `FlowState` history/committed 和 accepted metadata；
- fixed-time retry controller 的持久状态；
- Task 11 pressure authority 的 history/committed wall gradients 与 presence；
- IBM surface/domain/active layout/ghost/wall/LFP 静态身份。

只包装 `FlowState` 或只增加 geometry manifest 都不满足要求，因为重建出来的
pressure wall authority 会改变 Restart 后第一步的 BDF2/force 路径。

17A 冻结为 fixed-time constant IBM。adaptive time、WALE 和三密度 presence
组合进入 17B；17A 不伪造未经 time-controller 证明的 adaptive state。

## 架构

Checkpoint v3 复用 HUNDUN-FLOW 自有 Checkpoint v2 的字节编码、CRC-64、
verified temporary、no-overwrite publish、collective convergence 和
`FlowState` replacement transaction，但不改变 v2 格式。

目录权威：

```text
checkpoint-v3/
  rank-%06d.v3.bin        # 每 rank 的 accepted FlowState payload
  manifest.v3.bin         # presence、IBM identity、authority state
  COMPLETED               # 最后发布；缺失即不可读
```

读取顺序固定为：inventory/marker/manifest/CRC/presence/fingerprint/
partition/rank payload/authority-layout 全部通过 → 同时准备 FlowState 与
facade restore buffers → 所有 rank 收敛 → `noexcept` 发布 FlowState 与
pressure authority。任一 rank 在发布边界前失败时，两类状态均保持不变。

## presence 与 identity

- presence tag：IBM present、WALE absent、density constant；
- schema/config identity 继续由内层 v2 common payload认证；
- 数值 config identity 排除 case name、目标步数、Restart/diagnostics/
  performance 路径与频率；mesh、时间步参数、physics、boundaries、IBM/LES
  选择仍被认证。rank/process-grid 由 partition gate 独立认证，STL 内容由
  surface fingerprint 认证；
- outer identity 至少包含 surface、query、classification、surface coverage、
  active cell、active boundary、ghost plan、wall plan 和 LFP algorithm
  fingerprint；
- rank count、process grid 和 owned box 继续由 v2 manifest认证；
- authority rows 以稳定 link id 升序编码，保存 exact FP64 bits；pending/
  attempt-local force 不持久化。

## 文件白名单

产品：

- 新增 `include/hundun/flow_checkpoint_v3.hpp`
- 新增 `src/flow_checkpoint_v3.cpp`
- 新增 `src/flow_checkpoint_v3_detail.hpp`
- 修改 `include/hundun/flow_immersed.hpp`（仅 friend/opaque checkpoint access）
- 修改 `src/flow_immersed.cpp`（仅 snapshot/prepare/noexcept publish；不改算子）
- 修改 `src/flow_adaptive_time_control_detail.hpp` 或增加窄 factory，用于构造
  fixed-time、seal-valid 的持久 controller state
- 修改 `src/CMakeLists.txt`
- 17A codec 完成后再修改 `src/app_immersed_flow_driver.cpp` 接入 read/write；
  接入前先追加 driver RED。

测试/治理：

- 新增 `tests/mpi/test_checkpoint_v3.cpp`
- 新增 `tests/unit/test_checkpoint_v3_header_contract.cpp`
- 修改 `tests/CMakeLists.txt`
- 本 packet 与 receipt

若需要修改 v2 protocol，只允许增加格式无关的复用原语；禁止改变任何 v2
magic、schema、文件名、CRC 或现有读写行为。

## mutation-sensitive RED

1. 错 presence tag 被接受；
2. 漏掉任一 geometry/plan fingerprint；
3. authority link layout、presence 或 FP64 bits 未认证；
4. manifest/marker CRC 绕过、trailing bytes 或 outer inventory 放宽；
5. `COMPLETED` 在 manifest/rank payload 之前发布；
6. failed read 先改 `FlowState` 或 facade authority；
7. rank/process-grid/owned-box identity 被忽略；
8. pending gradient、attempt-local force 或 solver cache 被持久化；
9. continuous-vs-restart 下一步 state/force/report 不 bitwise；
10. v2 reader 接受 v3 outer directory，或 v3 reader把 v2误判为 v3。

## fast gate

- public header standalone contract；
- 1-rank encode/decode、CRC/truncation/trailing/inventory/presence/fingerprint；
- failed-read 对 state + authority bitwise neutral；
- 12³ fixed-time continuous-vs-restart 下一步 bitwise；
- 1/2-rank corruption/continuation；
- 1/2-rank driver write，1-rank driver restart continuation；
- focused Debug 与 tests-off；不运行 sanitizer、48³或完整矩阵。

## 版权与回退

只复用本项目 Apache-2.0 v2 protocol 和 Task 11 authority 数据模型；不采用
外部或私有源码，不新增依赖。允许回退范围仅限上述新增 v3 文件、非数值
checkpoint access 和直接 tests/driver wiring。

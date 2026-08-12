# Task 13 + 19B：body-fitted WALE flow 冻结证据包

状态：`TASK_GATE_ACCEPTED`

基线：`43d8446a4de354c40749b94d1f2143a8d744a4e9`

## 目标与科学合同

本节点只把 Task 12 的 WALE 系数接入 body-fitted constant-density PISO，并让同一
`hundun` 程序接受 schema v3 的 `immersed_boundary=none, les=wale` 组合。每次
attempt 的顺序冻结为：

```text
begin attempt
-> exchange committed/history
-> construct lagged velocity
-> compute lagged gradient
-> evaluate WALE exactly once
-> interpolate mu_sgs and form frozen mu_eff
-> predictor and exactly two PISO correctors
-> final momentum/transport residuals
-> collective commit or bitwise rollback
```

Backward Euler 使用 committed velocity；BDF2 使用
`u_lag = u_n + dt/dt_previous * (u_n - u_nm1)`。同一次 attempt 的 predictor、
history/committed spatial residual、两次 corrector 后的最终残差和 boundary
contribution 共享同一 face `mu_eff = mu_molecular + mu_sgs`。enthalpy 和 scalar
分别使用同一 `mu_sgs/Pr_t` 和 `mu_sgs/Sc_t`，不建立第二套 WALE authority。

## 硬不变量

- molecular-only 调用继续进入原有 uniform-viscosity overload，不构造或求值 WALE；
- WALE 只在 attempt 内求值一次；失败不发布 `WaleSummary`；
- accepted report 仍记录恰好两次 PISO corrector；
- retry 从 committed/history、dt/order、lagged gradient 和 density 重新派生稳定 identity；
- variable-viscosity FVM 在写 residual 前完整验证 face field；
- periodic face pair 使用一个确定性共享插值 authority；
- IBM+WALE 继续显式拒绝，留给 Task 16；
- WALE presence 的 Checkpoint v3 读写留给 Task 17B，不能声称已经支持 continuation。

## 文件白名单

- `include/hundun/flow_constant_density_piso.hpp`
- `include/hundun/flow_immersed.hpp`
- `include/hundun/fvm_cell_centered.hpp`
- `include/hundun/les_wale.hpp`
- `src/CMakeLists.txt`
- `src/app_immersed_flow_driver.cpp`
- `src/flow_constant_density_piso.cpp`
- `src/flow_immersed.cpp`
- `src/fvm_cell_centered.cpp`
- `src/les_wale.cpp`
- `tests/CMakeLists.txt`
- `tests/acceptance/task19a_immersed_flow_dispatch.sh`
- `tests/mpi/test_cell_centered_fvm.cpp`
- `tests/mpi/test_wale_body_fitted.cpp`
- `tests/unit/test_cell_centered_fvm_header_contract.cpp`
- 本 packet、receipt 和 ledger

## mutation-sensitive RED

1. WALE 在一个 attempt 中被重复求值；
2. BE/BDF2 lagged velocity 选择错误；
3. momentum 的任一消费者退回 molecular-only viscosity；
4. enthalpy/scalar 未使用同一 `mu_sgs`，或 Prandtl/Schmidt divisor 互换；
5. constant face viscosity 不再 bitwise 等价于旧 uniform-viscosity API；
6. periodic pair 各自插值，破坏守恒和 1/2-rank 一致性；
7. injected failure 后 trial、history、metadata 或 WALE summary 泄漏；
8. retry identity 依赖 attempt-local monotonic counter；
9. corrector count 不是 2；
10. driver 错误接受 IBM+WALE，或把 body-fitted WALE 当作 immersed force 路径；
11. WALE Checkpoint v3 presence 在 Task 17B 前被静默接受。

## task 与 milestone 门

- affected Debug unit/header/FVM/flow/driver；
- body-fitted WALE 1/2-rank fast；
- 同产品路径 12/24 screen；
- focused Release、ASan、UBSan；
- tests-off `hundun` build。

不运行 48-cubed、96-cubed、大型 MPI 或完整 Stage 3 矩阵。48-cubed WALE 与
24-cubed 1/2/4-rank 正式科学矩阵仍由 Task 21 冻结候选统一负责。

## 版权独立性

数学行为来自已批准 WALE 规格与有限体积守恒合同。实现、命名、控制流、periodic
authority 修复和测试场均为 HUNDUN-FLOW 独立代码；未复制、翻译或机械改写
AMReX/incflo、OpenFOAM、Basilisk 或任何私有源码，未新增运行时依赖。

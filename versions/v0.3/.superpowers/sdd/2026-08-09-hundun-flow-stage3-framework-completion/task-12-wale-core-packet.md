# Task 12：backend-neutral WALE core 冻结证据包

状态：`TASK_GATE_ACCEPTED`

基线：`da61bf9`

## 数学与尺度路径

实现 Nicoud--Ducros WALE：

```text
S = (g + g^T)/2
G = g g
Sd = (G + G^T)/2 - tr(G) I / 3
nu_t = (Cw Delta)^2 (Sd:Sd)^(3/2)
       / ((S:S)^(5/2) + (Sd:Sd)^(5/4))
Delta = cbrt(V)
mu_sgs = rho_attempt nu_t
```

为避免大/小梯度 overflow/underflow，先以 `m=max(abs(g_ij))` 归一化，使用
归一化不变量计算无量纲比值，再乘回唯一的 `m`。该变换严格齐次，不增加
epsilon；`m==0` 或两个数学不变量均为零时返回 bitwise `+0.0`。

## 文件白名单

- 新增 `include/hundun/les_wale.hpp`
- 新增 `src/les_wale.cpp`
- 修改 `src/CMakeLists.txt`
- 新增 `tests/support/les_wale_test_access.hpp`
- 新增 `tests/unit/test_wale.cpp`
- 新增 `tests/mpi/test_wale_mpi.cpp`
- 新增 `tests/unit/test_wale_header_contract.cpp`
- 修改 `tests/CMakeLists.txt`
- 本 packet、receipt、ledger

## 接口与所有权

保留 `WaleModel`、`WaleAttemptCoefficients`、`WaleCoefficientIdentity`、
`WaleSummary`。model 在构造期冻结 active-order local indices 与 cell volumes；
evaluate 只在 lexical callback 内取得 `KernelCellView`，不保存 checked/kernel
view。系数对象拥有一个 backend-neutral Buffer，按 active order 暴露 `nu_t` 和
`mu_sgs` 两个只读 view。identity 由配置、step/dt/order 和四个调用方提供的
fingerprint 派生，不是 committed revision。

## mutation-sensitive RED

1. `g=0` 返回负零或固定 epsilon；
2. `G` 误写为 `g g^T`；
3. 漏掉 `Sd` 的 symmetric/deviatoric 部分；
4. 分母指数错误；
5. Delta 不取 active background cell volume 的立方根；
6. `mu_sgs` 使用 stale/其他 density；
7. 读取 inactive cell 或重排 active IDs；
8. tensor 正交旋转、速度尺度、长度尺度不一致；
9. `y^3` 近壁斜率不在 `[2.9,3.1]`；
10. metadata/input fingerprint mutation 不改变 identity，或 model 持有 monotonic
    committed revision；
11. stale checked view 被接受；
12. device/non-host backend 被误报为可执行生产路径。

## task gate

只运行 tensor/exact-zero/y³/rotation/dimension/identity/stale-density 单元 RED，
1/2-rank active-order/failure smoke、standalone header、focused Debug 和 tests-off。
TGV 12/24 与 focused sanitizer 进入 WALE milestone，不在 Task 12 重复。

## 版权

公式来自已批准规格和公开论文数学定义；实现、命名、控制流和测试 oracle 均为
HUNDUN-FLOW 独立代码，不复制 OpenFOAM、AMReX/incflo、Basilisk 或私有源码，
不新增运行时依赖。

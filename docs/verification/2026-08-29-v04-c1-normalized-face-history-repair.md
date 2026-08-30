# HUNDUN-FLOW v0.4 C1 normalized-BDF 面历史修复回执

日期：2026-08-29

裁决：`P-H1 / C1 FACE-HISTORY AUTHORITY REPAIRED`。这不是完整的压力—焓
恢复结论；IBM 冷启动和同时间层压力—焓耦合仍是独立后续工作。

## 1. 源身份与隔离

- worktree：`/home/wyf/code_dev/hundun-flow-pressure-enthalpy-c1-production`
- branch：`codex/v04-pressure-enthalpy-c1-production`
- base commit：`9e013716b3525b59891f891ba2891d0bb1157337`
- base tree：`9001b2920e5a73940daeca60c0831e01ecd1da74`
- base parent：`2c234b14de71f3168f0bc4a7b5778666fbb688a8`
- `versions/v0.4` source/test patch SHA-256：
  `8b884d566f1eb028b562b8492330d2c5f750eae258f96819acd8ba1c906e6433`
- Release core SHA-256：
  `f0b06effd17a42285d3794aecbfabf293345c371f2d3e92d06a6ff202a5be635`
- production `hundun` SHA-256：
  `db17feb85815074ff6fbadcdd9b751aaed3575b71f5a0d2434fafc48c2e554c3`

原 dirty worktree 未承载 production 修改。

## 2. RED 与根因

旧 C1 在内部、周期和 pressure-correctable 面使用：

```text
phi_HbyA = current + a0*rho*rAU*(paired_EX2_flux - temporal_flux)
```

其中 `paired_EX2_flux` 与动量 predictor 的 normalized-BDF 历史不是同一时间
离散 authority。最小测试先要求 `PisoIntermediateInput` 和 certificate 显式携带
committed face history；旧代码因类型不存在而 RED。随后用真实
`FinalFaceFluxWriter` 发布 `phi_n=2, phi_nm1=-1`，对步长比
`r={0.5,1,2}` 分别要求 `{19/8,3,22/5}`，并确保与 paired EX2
`{7/2,5,8}` 可区分。

只读复核另发现跨轴 alias 漏检：旧预检只比较
`history[axis]` 与 `output/trial[axis]`。测试构造 x 输出指向 certified y 历史的
公开合法 bind seam；旧 refresh 成功，RED。修复后所有 3/6 个历史面与全部
output、trial、pressure 面做全配对检查，BDF2 accepted 与 previous 也做全配对
检查。该路径仅使用栈数组和定长循环。

最终复核又冻结了两个旧实现可触发的 RED：

- certified face-history storage 可以与 cell workspace `rAU/HbyA/grad(p)` 的连续
  内存区间重叠；旧代码在读取历史前先覆写 workspace，造成合法证书下的静默历史
  破坏。C1 预检现在对每个历史轴与三个可写 cell view 做完整 overlap 拒绝；
- periodic/internal C1 不消费 paired EX2 flux，却仍无条件读取并检查该 payload。
  测试把 paired payload 设为 NaN、保留有效 metadata 和 committed history，旧代码
  RED；现在只有 fixed physical boundary 才读取 paired payload，并且 dependency hash
  只在实际消费时纳入其 revision。

## 3. 实现

C1 内部、周期和 pressure-Dirichlet 面现在直接计算：

```text
phi_hist = (-a1*phi_n - a2*phi_nm1)/a0    # BDF2
phi_hist = phi_n                           # BE startup/retry
phi_HbyA = current + a0*rho*rAU*(phi_hist - temporal_flux)
```

固定物理边界仍复制 thermophysical paired flux；没有修改标准 variable-step BDF2
系数，也没有修改 IBM 无穿透约束。

新增 authority 契约包括：

- predictor certificate 发布 committed writer authority、storage 和 revision
  domain；
- C1 验证 accepted/previous certificate、revision、writer、storage、domain 和
  BDF2 两层相异性；
- BE 只接受 accepted 层，C2 必须清空两层 committed history；
- history lineage 进入 intermediate dependency/certificate；
- stale revision、writer、storage、revision domain、缺失/重复 previous、跨轴
  alias 及 history/workspace alias 均原子拒绝；
- fixed velocity inlet 继续复制 paired physical-boundary flux，pressure outlet 与内部面
  使用相应 authority，边界 oracle 对三类面分别冻结；以及
- retry 的 effective BE 会清空 previous history，C1 只消费 accepted committed flux。

## 4. 测试

构建：

```text
cmake --build build/c1-authority-release -j 4
```

结果：所有目标编译成功，包括 production `hundun`。

相关回归一次运行 19/19 通过，覆盖 thermophysical predictor、PISO authority /
mutation / temporal order、equations/PISO MPI 1/2/4、product freeze MPI 1/2、
I/O product path，以及 retry fallback。非零三轴 committed history 的 PISO seam
oracle 在 1/2/4 rank 均逐面得到 BDF2 期望值 `{4.375,6.5,8.625}`，包含全部
decomposition seam。

retry fallback 通过独立 `v04_core_temporal_fallback_mpi_{1,2,4}` 冻结：测试库在
下一次二阶 predictor 集体发布一次精确
`low_bdf_source_base_admissibility`，driver 必须在同一 attempt 内得到
`proposal order=2 / effective order=1 / predictor calls=2 / fallback=true`，随后完成
两次 PISO 并接受 step 2。该 hook 位于 `src` detail header，只由
`hundun_v04_test_core` 的 PRIVATE macro 编译；`nm -C` 已确认 production core 和
`hundun` 均不含 hook 的变量或 arm/clear 符号。

最终全套 CTest：134/136 通过。新增三个 fallback MPI 测试全部通过。两项失败已在
未修改 Route-A 源上独立复现：

- `v04_app_driver`：过期 warm-retry golden；当前稳定值为 C1
  `0.017762137667637876 / 0.013429488737914615`、19/19 iterations，C2
  `0.017762882167295341 / 0.013429901870250003`、3/3 iterations；
- `v04_core_product_freeze_mpi_4`：既有 forced implicit-enthalpy endpoint
  断言失败。

这两项测试及 app-driver production 源未为本修复改基线。`git diff --check` 通过，
最终静态复核未发现 Task 1 blocker。

## 5. exact 64-rank Re3900 门

Case：480×480×16、8×8×1 decomposition、`dt=0.006 s`、variable BDF2、IBM
cylinder。最终 production tracer：

- executable SHA-256：
  `93828ca53220901ce06bf27049e6a468bf7c8fbc876f57ba920a17fe0abbe43a`
- case JSON SHA-256：
  `cc502ee6478c05f53dec23c74400d5e996beb9dca1623d9de172c8e5e3ead0ee`
- 60-step trace SHA-256：
  `3404386edf9f84924efda5a5f9bd3f07505c48c063cc99d9a970b0c8c0b1790d`
- stderr SHA-256：
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

最终复验工件位于
`re3900_c1_face_history_production_01/run-exact-60-final-v2`；60-step trace 与上一轮
final trace byte-identical。

| 量 | step 60 | 门限 | 结果 |
|---|---:|---:|---|
| pressure span | 99.555263 Pa | < 150 Pa | PASS |
| `max|U|` | 1.964775 m/s | < 3 m/s | PASS |
| max cell CFL | 0.390118 | < 0.6 | PASS |
| outlet backflow faces | 0 | 0 | PASS |
| continuity residual | 4.558e-9 | 1e-6 | PASS |
| EOS residual | 4.441e-16 | 1e-6 | PASS |
| attempts | 1 | 无 retry | PASS |

全程最大速度 2.252217 m/s，最大 CFL 0.434087，均出现在 step 2；无回流、无
retry。BE startup 的 step 1 pressure span 为 167.193222 Pa，高于 step-60 参考
上限。这是已隔离的 IBM 冷初场不相容瞬态；若把 150 Pa 解释为“包括 startup 的
逐步硬上限”，则该更严格门仍未通过，不能由本 C1 修复宣称解决。

最终二进制的 no-IBM 64-rank 两步控制保持塞流不变量：`max|U|=1`、pressure
span=0、进出口质量流完全相等、EOS residual=0。trace SHA-256：
`4d1e7861d029dd3f8e5b74d5184a249cc8559169e7f79af3dc212d5e06b9d9be`；工件位于
`run-no-ibm-2-final-v2`，也与上一轮 final trace byte-identical。

## 6. 结论与下一门

原路径最后接受 step 34、attempted step 35 失败；本 production 修复通过 step 60，
消除了 C1 的快速放大链。历史 normalized-BDF 反事实在 step 160 仍出现较慢失稳，
因此本结果不能推导为压力—焓原型已经正确。

下一步应建立完整 `p -> Dp/Dt -> rate history -> predictor -> EOS/PISO` 闭环 RED，
验证完整 `p-h` block 及其 pressure-only exact Schur。IBM 冷初始化、IBM 双度量、
入口 final-flux authority 和出口禁回流继续保留为独立可回滚工作流。

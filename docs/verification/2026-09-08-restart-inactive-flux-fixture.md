# Restart 固体内部通量夹具纠正

基线 `95eb7d9a522098b618d6acd46175b01e51e5716b`。这是扩大 MG 观测回归时
发现的旧测试合同缺口，不是本轮发现的 Re3900 数值失败或新的恢复算法修复。

## 根因与独立证据

`initialize_restart()` 的 `canonicalize_ibm_flux` 复制输入通量后调用
`IbmEquationInterfacePlan::zero_interface_flux()`。当前静止占位方法除了界面，
还将所有 inactive cell 的六个相邻面置正零，退出无效 solid-solid 历史。
源文件、有效流体状态和质量目标合同未因此被放宽。

两个 Fresh/Restart 测试按 `low_solid != high_solid` 识别界面，把固体内部面
误列为必须逐位恢复的区域。夹具交替写入正零和负零，因而必然与现行方法冲突。
这是合法 legacy V1/BE-recovery 夹具；缺少当前 method signature 不是根因，
不能通过向 V1 单独加 signature 或绕过历史校验来修测试。

实际反馈链：

- 新工作树 Release 的相关组 12/15 通过，Fresh MPI 1/2/4 都在相同综合断言失败。
- 独立干净源码 `be721a7a581ae073cd703addffc692aab0a95617`，其 `versions/`
  和 `tools/` 与基线完全相同，在新建 `/tmp/hundun-mg-driver-baseline-20260908`
  构建也于 1 rank 失败（2.95 s）。没有增量 Ninja 恢复提示。
- 定向诊断 1 rank：Status=0/0、U 差异=0、通量位差异=672，全部672属于
  solid-solid；非零界面=0；跳过 Fresh 和派生量检查均通过。
- 独立串行 `v04_product_fresh_projection` 随后实际 RED（0.99 s），同样是
  solid-solid 负零预期过期。不是只凭搜索结果标记失败。

一个解析反例是 x 面 `(5,4,4)`：两侧 cell `(4,4,4)/(5,4,4)` 均在测试 cube
内；面索引1161为源负零，而现行 normalizer 输出正零。测试 STL 名为 cylinder，
该 fixture 的实际几何为 `[-1,1]^3` cube。

## 最小修改与不变项

只改 `tests/mpi/product_fresh_projection_mpi_test.cpp` 和
`tests/integration/product_fresh_projection_test.cpp`：

- 通量预期使用 `low_solid || high_solid`；fluid-fluid 的位比较仍保留。
- 串行夹具要求所有 solid-adjacent 为正零；MPI 重复分区面按既有合同允许任意零符号。
- 原192个负零切面计数保留。串行另核对672个负零 solid-solid 面，即
  `3×7×8×8/2`，将几何切面与内部面分开。
- 原字段、速度、未调用 Fresh 投影、派生量与速率历史检查不变。
  无生产源码、物性、阈值、方法签名或 checkpoint 格式修改。

临时 `[DEBUG-mg-driver-restart]` 探针已从源码移除，原始输出保存在
[数据目录](data/2026-09-08-product-mg-profile/)。
修正后 Release 的 Fresh串行/MPI1/2/4及MG重试MPI1/2/4共7/7通过（6.21 s）；
ASan+UBSan 的Fresh串行/MPI1/2/4共4/4通过（59.73 s）。保留了全部旧生产要求，
仅纠正区域预期；这不是新的物理方法接受。源码临时探针扫描为空、diff检查通过。
调试构建使用 `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1`，不宣称泄漏检查。

预防措施：维护方法历史语义时，同步审查测试中“哪些区域保持位相同”的合同。
把固体内部占位与有效流体控制体分开，不能用一个“非界面”分类替代二者。

# HUNDUN-FLOW Stage 3 紧凑科学框架设计

> **后续修订候选：**
> `2026-08-09-hundun-flow-stage3-parallel-completion-v2-design.md` 已按
> `7fc8c508...` 之后的真实代码状态重排剩余任务。用户批准前本文件仍为控制设计；
> 批准后，v2 只取代未完成任务顺序、并行边界和长测调度，本文的已接受能力与科学
> 约束继续有效。

**状态：** 用户批准，控制 Task 11 之后的 Stage 3。
**执行 profile：** `stage3-compact-scientific-v1`
**治理基线：** `ee4d2b18d0c68b3080edcc2e132045175961cfb8`
**Task 11 科学基线：** `66080e324089599711fdb26082af9b330bfdb5ce`

## 1. 决定

Stage 3 优先完成用户可运行的整体求解框架，再集中做不会阻断正确性和实际使用的性能、极端规模与维护性工作。每个 task 只承担能改变产品能力或阻断后续组合的最小闭环；昂贵科学矩阵集中到里程碑和最终候选。

本设计不降低以下门槛：守恒、收敛性、两次 PISO corrector、唯一 pressure/operator/final-flux/force authority、rollback、collective failure、MPI 分解一致性、Restart 完整性和用户可运行性。

永久禁止 `96^3`。正式网格最大为 `48^3`。

## 2. 仓库权威

- `/home/wyf/code_dev/hundun-flow-governance` 是唯一开发、测试、证据和验收仓库。
- `/home/wyf/code_dev/hundun-flow` 是无治理历史、无 remote 的产品投影仓库。
- Stage 3 开发工作树是 `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework`。
- product 只在仓库拆分完成时同步 `0.1.0`，在 Task 21 接受后同步 `0.2.0`。
- 旧 `/home/wyf/code_dev/.worktrees/hundun-flow-stage3` 保留全部脏修改和证据，但不再作为开发基线。

Task 11 数值 authority 未改变时，其哈希相同证据可复用。治理文档、索引和 seal 的变化不会单独使数值证据失效。

## 3. 数学与架构边界

HUNDUN-FLOW 继续拥有 topology、geometry、fields、execution、MPI、linear algebra、transaction、checkpoint 和 diagnostics。AMReX、IncFlo、OpenFOAM、Basilisk 等只提供公开数学行为和设计对照，不成为运行时依赖或代码模板。

每个 semantic port 必须按以下顺序实施：

```text
公开数学/行为来源
-> HUNDUN 接口与 mutation 清单
-> mutation-sensitive RED
-> 独立最小实现
-> fast
-> screen
-> acceptance
```

禁止复制、翻译或机械改写上游源码、注释、控制流、ABI、宏和错误文本。GPL 与无明确许可证项目只能用于数学和测试思想对照。AMReX/IncFlo 虽为 BSD-3，本阶段也只做语义复用。

## 4. 运行时组成

```text
hundun CLI / app_immersed_flow_driver
  -> immutable case + static geometry/plans
  -> density attempt adapter
  -> optional WaleAttemptCoefficients
  -> FixedStepImmersedFlow
  -> final FaceMassFlux + force + attempt report
  -> transaction commit / rollback
  -> Checkpoint v3 + diagnostics
```

固定规则：

- geometry、surface query、active layout、donor ownership、ghost stencil、wall quadrature 和通信计划只在构造期建立；
- attempt-local 的 donor values、wall density、pressure correction、`nu_t`、`mu_sgs`、`mu_eff` 和 force report 在 rollback 时全部丢弃；
- WALE 每个 trial 恰好求值一次；predictor、两次 corrector、最终残差和壁面力共享同一 `mu_eff`；
- IBM/WALE 都关闭时，molecular-only 路径保持 bitwise 不变；
- inactive cells 不作为物理状态读取；
- diagnostics 只读已发布 authority，不重算第二套数值结果。

## 5. Stage 3 任务结构

### 5.1 可运行 IBM MVP

1. `19A`：constant-density IBM driver，同一 `hundun` 支持 schema v3、静止 STL、两次 PISO、force、retry/rollback。
2. `17A`：Checkpoint v3 IBM-only，冻结 presence tag、CRC、publish-last、failed-read rollback 和相同分区 bitwise continuation。
3. `18A`：最小 diagnostics/counters，覆盖 continuity、final residual、wall penetration、四字段 force、失败分类及静态/每步计数。

### 5.2 WALE 与组合

4. `12`：WALE tensor core 与 attempt identity。
5. `13+19B`：body-fitted variable-viscosity WALE 和 `none/wale` driver 一次集成。
6. `14`：material-density IBM。
7. `15`：ideal-gas IBM。
8. `16`：三密度 IBM+WALE 固定 trial 顺序的不可分硬门。

Tasks 14、15、16 顺序执行，不并发修改 flow composition。

### 5.3 完整框架

9. `17B`：Checkpoint v3 presence combinations。
10. `18B`：完整 diagnostics/providers。
11. `19C`：全部合法 driver combinations。
12. `20`：capability ledger、exact counters、最终验收 launcher。
13. 公共文档终审。
14. `21`：单一 exact-HEAD Stage 3 接受。

## 6. 明确延期和禁止

延期到框架完成后，除非 RED 证明影响正确性、内存有界性或实际可运行性：

- ordinary-host 严格零 allocation；
- replicated-vector 消除和极端规模防御；
- 非关键性能精调；
- project-owned geometric multigrid。

Stage 3 不实现 AMR、移动/多物体边界、FSI、生产 GPU、rank-changing Restart、几何多重网格、化学、喷雾、颗粒，也不在没有 mutation RED 时加入守恒再分配。

## 7. 测试分层

### task gate

- mutation-sensitive RED；
- 直接受影响的 Debug unit/header/policy；
- 一个 `12^3` 或更小的同产品路径 fast case；
- collective 改动才补 1/2-rank；
- 公共头改动才补 standalone-header；
- build graph/test seam 改动才补 tests-off；
- 主 agent 一次完成 requirements、quality、caller-impact 和 task diff 审查。

task 内不机械重复完整 Debug、Release、ASan、UBSan、1/2/4-rank 和 `48^3`。

### milestone gate

- MVP：constant IBM `12^3`、1/2/4、Restart bitwise、rollback、collective failure、tests-off。
- WALE：tensor、`y^3`、12/24 TGV screen、1/2-rank、focused sanitizer。
- Combined：12/24 constant wake、三密度小型 1/2/4、守恒/closure/retry。
- Framework：presence、Checkpoint、diagnostics、driver 的小型 1/2/4，以及 `nm`、`ldd`。

### final gate

最终候选冻结后只运行批准矩阵：WALE 12/24/48；constant IBM+WALE 48 和 24 的 1/2/4；material/ideal 12/24 短程 1/2/4；Checkpoint v3 12；driver/diagnostics 小网格；完整受影响 Debug、focused Release、小型 ASan/UBSan、Stage 1 低成本回归、Stage 2 core whitelist 和治理检查。

## 8. 资源与长作业

- `L`：build/header/policy/non-MPI unit，`cmake --build -j32`、`ctest -j24`，最多两个低内存组。
- `M`：12/24 和 MPI 1/2/4，单作业 IBM thread budget 总计 96，最多两个 M 作业并绑定不同 NUMA CPU 集。
- `H`：48 和正式长 selector，任意时刻只运行一个；可与另一 NUMA 节点的 L 组重叠。

长作业必须记录 exact HEAD、dirty diff、binary SHA-256、命令、环境、CPU 集、日志、退出状态、时长、峰值 RSS 和日志 SHA-256。最终 H 作业只允许在冻结候选上运行。

## 9. Worker 与文档

主 agent 独占总体计划、数学推导、19A 构造顺序、Tasks 13--16 跨模块组合、完整 diff、版权独立性和最终验收。边界明确且证据矩阵冻结的实现、机械迁移或独立测试才可委派给配置有效的 `luna_worker`；返回后必须核验实际 turn context。

`humanizer-zh` 和 `shuorenhua` 只在 Task 20 完成后、Task 21 候选冻结前调用。法律原文、内部计划、日志和证据不润色；公共技术文档润色后由主 agent 逐项复核公式、术语、命令、JSON key、单位、版本、SHA 和法律文本。

## 10. 完成条件

Stage 3 只有在 constant IBM MVP 可运行/Restart/诊断、WALE 与三密度组合通过紧凑科学验收、Checkpoint/diagnostics/driver/counters 完整、product 投影为签署的 `0.2.0`、governance exact-HEAD 报告完成且无遗留后台作业时才接受。接受后停在 Stage 3 边界，不进入 Stage 4。

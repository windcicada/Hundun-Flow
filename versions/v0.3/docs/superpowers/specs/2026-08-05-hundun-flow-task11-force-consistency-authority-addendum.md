# HUNDUN-FLOW Task 11 Force-Consistency Authority Addendum

状态：2026-08-05 经用户确认；在 Task 11 acceptance 范围内立即生效。

本 addendum 只恢复并澄清 Task 11 的 force/consistency hard-gate authority。
它不预判尚未执行的 single-link 符号 RED，不在本文中静默改变现有正式网格、
制造解、阈值、PISO corrector 数量或其他 Stage 3 科学合同。若现有 signed-force
fixture 未通过下面的 candidate-independent non-degeneracy preflight，必须先建立
一份冻结具体场、分量、解析 reference 和 mutation 的 fixture addendum，之后才
能运行产品候选；不能从候选误差反推 fixture。

## 1. 被替代的窄范围结论

以下历史记录及其后继文件中“无收敛阈值、仅 diagnostic”的分类，在 Task 11
acceptance 上被本 addendum 替代：

- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-hybrid-lfp-and-force-measure-amendment.md`
  中将三项
  operator/surface difference 降为 diagnostic 的结论；
- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-viscous-traction-oracle-amendment.md`
  中将 signed viscous net force
  和三项 consistency 降为 diagnostic 的结论；
- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-a13-force-measure-oracle-reconciliation-brief.md`
  和同目录 `task-11-a13-force-measure-oracle-reconciliation-result.md` 的对应分类；
- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/task-11-root-cause-closure-sweep.md`
  和同目录 `ledger.md` 对该降级的传播或追认。

这些文件中的数值、日志、推导、失败案例和审查历史继续保留为诊断证据；只有
“可以不满足二阶仍接受 Task 11”这一结论失效。

## 2. 恢复的 hard gates

Task 11 的同一正式三层 refinement matrix 必须独立评估：

- true-surface pressure force；
- true-surface viscous force；
- cancellation-resistant viscous traction norm；
- true-surface total force；
- operator/surface pressure consistency；
- operator/surface viscous consistency；
- operator/surface total consistency。

每个适用量必须有限、严格为正、随 refinement 严格下降，并在两个相邻 segment
分别满足 frozen Stage 3 specification 的 `>= 1.8` 门槛。pressure、viscous 和
total 必须分别通过，不能只用 total cancellation；制造解和 normalization 必须
在运行候选前冻结，不能从结果反推。

`viscous traction norm` 不能替代 signed viscous force；它只提供抗闭合表面抵消
的额外局部分辨率证据。若现有制造解使某个 signed 分量退化到 roundoff，必须
停止 formal matrix，并按上一段建立事先冻结的非退化 fixture，而不是删除该
gate、运行后挑分量或把 traction norm 冒充 signed force。

## 3. 尚未冻结的符号语义

当前记录同时包含三种彼此不能全部成立的说法：

1. frozen specification 把 `operator_reaction_force` 定义为 immersed row
   residual 的负和；
2. frozen specification 又要求它与 fluid-on-solid 的 true-surface traction
   采用同向 difference；
3. 当前产品把 report 构造成负 wall residual，并用
   `operator_reaction + surface_traction` 形成 consistency。

因此本 addendum 不以文档投票选择加号或减号。最终语义必须由同一个可执行
single-link pressure/viscous RED 冻结，RED 必须明确：

- `n_s`（solid-to-fluid）和 fluid-domain outward normal 的关系；
- Cauchy stress、fluid-on-solid traction 和 momentum-residual wall term；
- raw wall-row contribution；
- momentum-budget reaction；
- 面向用户和 consistency 的 operator force；
- surface traction；
- 两者使用加法还是减法比较。

单-link RED 必须分别使用非零 pressure-only 和 viscous-only 场，且 normal
不得与坐标轴对齐。它必须能杀死 normal reversal、raw-row/report sign reversal、
pressure surface sign、viscous surface sign 和 consistency-operation mutation。
它只冻结方向、命名和独立 authority，不要求一个 full-background-cell operator
link 的离散 measure 与一个真实三角 quadrature weight 逐点相等；两种全局表达的
差值在正式多网格序列上按二阶合同判定。

若 RED 证明 frozen specification 的“negative residual”或“minus surface”措辞
错误，主 agent 必须在修改产品前提交一个只修正该矛盾的科学语义修订。不能在
测试 helper 中静默换号，也不能让 operator report 复制 surface result。

## 4. 独立性和 shared authority

operator 和 surface 是两条独立产品数值路径：

- operator 从唯一 immersed residual 及其共享 boundary-row authority 得到；
- surface 从真实 STL 三角面求积得到。

二者可以共享冻结的几何、normal、pressure boundary 和 reconstruction authority，
但不能共享最终 force 数值，也不能一方调用另一方作为 oracle。测试必须用解析
场独立计算预期量，并证明以下 mutation 可被检测：

- shared row 回退为 per-link reconstruction；
- 只取 multi-link row 的一个 link；
- background flux 未完整移除或被重复加入；
- pressure diagonal/neighbour defect 遗漏；
- viscous wall derivative 使用 stale 或不同 authority；
- operator force 直接复制 surface force。

## 5. Acceptance 边界

force/consistency 关闭顺序固定为：

```text
mathematical derivation
-> mutation-sensitive executable RED
-> explicit semantic amendment when needed
-> minimal implementation
-> 12³/24³ fast
-> 24³/48³ screen
-> nine-sequence 24³/48³/96³ acceptance
```

在 single-link RED、shared-row RED 和 exact-state consistency decomposition
全部成立前，任何已运行的 24/48/96 日志都只能作为诊断，Task 11 保持
unaccepted。

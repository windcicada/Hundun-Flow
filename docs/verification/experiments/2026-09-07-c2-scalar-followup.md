# C2 工作量、标量近抵消与 IBM 几何分离实验

本轮从 `e0fd326d3ed6cccd52f812c113712f34a08cec95` 继续核查。
本文件对应请求的第 4–7 项，与 [I/O 正确性修复](../2026-09-07-e0fd326-io-contract-audit.md)
分开记录。这里没有修改生产算子、32eps/128eps、时间步、refinement 上限或物性。
没有运行新的 Re3900 pilot，也没有把已有冻结验收标签授予当前构建。

以下本轮证据均位于：

```text
/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/module-io-audit-20260907/
```

## 4. C2 r1：先区分工作次数与单次成本

使用同一天 `module-audit-20260907/preopt-pilot` 的原始逐 rank 文件，
以当前 V3 observer 再次独立验证，输出 `loops-revalidated.json` 和
`loops-revalidated.jsonl`：128 ranks、5501–5510 全部十步、61 loops，
`complete=true`，无缺失、重复或来源问题。整个窗口 attempt/sweep 均为 1；
10 次 C1 spatial、51 次 C2 diagonal，没有 C2 spatial 样本。
这是冻结正确方法下的原始窗口，不是已撤回稀疏 IBM 实验的结果。

[只读分析脚本](2026-09-07-c2-loop-analysis.py)绑定 observer 的输入和 details 哈希，
记录自身及汇总哈希。复现命令中的 `ART` 为上述证据目录，`RUN` 为原 pilot：

```sh
python3 tools/v04_solver_observe.py "$RUN" \
  --output "$ART/loops-revalidated.json" --details-output "$ART/loops-revalidated.jsonl"
python3 docs/verification/experiments/2026-09-07-c2-loop-analysis.py \
  --run-root "$RUN" --observation "$ART/loops-revalidated.json" \
  --details "$ART/loops-revalidated.jsonl" --output "$ART/c2-cost-analysis-final.json"
```

输出文件应使用新路径，不覆盖原始记录。下表时间是窗口内各 loop 的 **rank 平均累计**，
不是进程总墙钟；M 已包含内部通信，不能再把等待时间相加。

| corrector / refinement | 求解次数 | Krylov 迭代 | solve / s | M / s | M 每次 / ms | A 每次 / ms |
|---|---:|---:|---:|---:|---:|---:|
| C1 r0 spatial | 10 | 219 | 8.97164 | 2.44495 | 11.1642 | 17.0091 |
| C2 r0 diagonal | 10 | 211 | 4.10159 | 2.35471 | 11.1598 | 3.57747 |
| C2 r1 diagonal | 10 | 329 | 6.25490 | 3.58890 | 10.9085 | 3.57970 |
| C2 r2 diagonal | 10 | 248 | 4.73904 | 2.69737 | 10.8765 | 3.57428 |
| C2 r3 diagonal | 10 | 235 | 4.48009 | 2.57172 | 10.9435 | 3.58847 |
| C2 r4 diagonal | 9 | 182 | 3.53483 | 1.99557 | 10.9646 | 3.59023 |
| C2 r5 diagonal | 2 | 44 | 0.82992 | 0.48080 | 10.9272 | 3.57294 |

r1 相比 r0 多做线性工作，但 A/M 单次成本接近。r1 的 true residual
收缩中位数为 7.87e-7，r0 为 4.08e-5。r1 初始 true residual 为
0.01521–0.01546，最终为 6.77e-9–1.41e-8；不是“初始残差较大”就能解释。

源码核查：`core_product_freeze_detail.hpp::product_pressure_inexact_forcing_control()`
保留接近终止带或不下降时回到基础线性容差的保护；
`core_product_freeze.cpp` 的 refinement 循环首次 previous merit 来自 C1。
将原始 CSV 中已选候选的 continuity/energy 代入未改变的控制规则，
本窗口 10 次 r1 都满足相对于 C1 的不下降分支，推得相对容差 1e-6；
C1/C2 的 mandatory r0 使用 1e-4。C2 r0 的 selected merit 是 C1 的
259–402 倍，但它们属于不同校正阶段，不能据此宣称同一残差序列发散。

这段容差说明是**源码规则重建**，并非 CSV 直接测得的实际容差字段。
`solver_krylov.cpp::convergence_limit()` 还取绝对/相对要求的最大值，
后续较小 RHS 会受绝对阈值影响。不能用 final/initial 单独反推全部停止条件，
也不能把这项观察直接变成放宽容差或提前 spatial 的理由。

现有 M 数据不能区分细层平滑、粗层求解和 MG 自有通信。只选择的下一项最小实验是
**低频 MG 层级计时**：按 level 记录平滑、限制/延拓、粗解、通信，明确嵌套时间归属，
并把 forcing 分支/实际线性阈值附到相同 loop。先以计时关闭/开启检查观测成本，
再在同起点、物理窗口与输出要求下选择一个优化。当前没有实施该计时或内核实验，
没有重做 multidot、稀疏候选或 MG storage swap。

### 最终动量：本窗口只支持单 dt 描述

保留已有归一化和最差单元/IBM 分区记录。窗口固定 dt=1.3808912271980336e-5 s、
BDF2，a0=108625.5 s⁻¹。原始 Rm 的 linf 范围分别为：
x 4.06e-7–1.29e-6 N，y 3.39e-7–7.02e-7 N，z 1.16e-7–2.11e-7 N。
Urms=2.96323–2.96334 m/s。归一化量范围仍在 JSON 中保留，但本轮没有
各项原始尺度的独立记录，也没有相同物理时间的 dt 细化运行，不能据此判断
缩小 dt 后动量误差更小。下一步需同时比原始 Rm、时间/对流/扩散/压力项尺度，
不增加经验动量接受门槛。

## 5. 有标量路径：补近抵消回归，未找到放宽阈值的反例

`product_scalar_contract_experiment.cpp` 增加 `--signed`、`--observe-cost`。
前者只在测试初场定义一个穿零且库存近抵消的被动标量；初始加权平移由解析
密度/单元宽度确定，不对推进结果做归一化或补账。参与 EOS 的组分仍非负。
独立库存审计增加 Qabs=Σrho|q|V，对近零净库存以 Qabs 衡量绝对误差；
这只改变测试的误差单位，不改 remap 残差、接受标准或任何生产字段。

实际 `--stretched --signed --observe-cost` 两 rank 回归通过，见
`scalar-signed-cost.log`。8×4×4 非均匀 COAST 网格，固定物理终点 1e-3 s，
dt 为 1.25e-4、6.25e-5、3.125e-5 s，分别推进 8/16/32 步；每种配置一轮。
1/4 ranks 的同场景 CTest 也通过，见 `scalar-signed-other-ranks.log`。

| 检查 | 被动标量 | EOS 组分 |
|---|---:|---:|
| 非均匀初始 | 穿零，\|Q\|/Qabs=3.08e-11 | 正组成，组分和/EOS 合法 |
| 最大独立库存误差 | 2.87e-15（以 Qabs 归一）；绝对值 1.03e-16 | 3.18e-13（相对库存） |
| 同终点时间阶数 | 1.93981 | 1.94580 |
| 库存、常量、有界、EOS、历史检查 | 通过 | 通过 |

这里只对被动标量施加了符号变化，不能声称组分允许负值。
所有矩阵/配对行残差原记录保留；本例未出现 32eps/128eps 导致的停滞。

粗时间步的非均匀场 8 步累计（rank 平均）：

| 路径 | advance / s | remap / s | Jacobi 总迭代 | 组成 sweep 总数 |
|---|---:|---:|---:|---:|
| signed passive | 0.038119 | 0.0004864 | 15 | 8 |
| EOS species | 0.137801 | 0.0011040 | 21 | 24 |

这是很小的回归夹具，不是目标规模热点证据。其结果不足以支持额外六份全场缓存、
缩 ghost 或省掉预测器。尚未将 remap 算术、halo 和控制归约进一步拆分；
若扩大有标量测试，冻结缓存必须同时绑定 rho、最终 phi、预测 Phi*、a0、
边界分支、几何及 workspace，sweep/retry 变化时失效。零标量生产路径没有新增工作。

## 6. IBM：独立分离几何，尚未设计生产曲面二阶替换

`solver_ibm_scalar_mms_test.cpp --geometry-separation` 新增独立实验入口，
不重复执行原来的平壁/斜壁/圆柱整组 MMS，也不改变其门槛。
固定圆柱半径 0.65、高度 0.9、1024 个周向平面，先用 24³/48³/96³ 网格；
再固定 48³ 网格，改变周向面数 48/96/192/1024。
解析圆柱体积为 1.1945906065；输出同时列解析体积、STL 多面体体积、二值活动体积，
以及解析解在多边形侧面中点的非零法向梯度，避免把几何边界数据差异当成离散失守恒。

实际一轮单 rank 实验通过，见 `ibm-geometry-separation.log`：

| 网格，固定 1024 面 | 二值活动体积 | 解 L2 | 细化阶数 |
|---|---:|---:|---:|
| 24³ | 1.3750000 | 0.00549450 | — |
| 48³ | 1.1347656 | 0.00174678 | 1.65329 |
| 96³ | 1.1723633 | 0.000961918 | 0.860711 |

STL 体积为 1.1945831106，与解析体积仅差 -7.50e-6；对应二值体积差为
+0.1804、-0.05982、-0.02222。固定 48³ 时，48/96/192/1024 面都得到
4648 个活动单元和完全相同的解误差，而 STL 体积差随面数明显下降。
本夹具的二值占据与配对零通量行不因进一步细分 STL 改变；结果说明提高 STL
分辨率本身不能恢复该序列的二阶精度。不能把体积差直接等同为全部解误差，
也尚未把二值几何误差与边界通量截断误差严格正交分解。

每个样本仍检查常量保持、配对行、净通量和 Helmholtz 真残差，均通过。
曲面高阶方法仍需整体设计体积/面几何和一致通量，再验守恒、有界、常量、EOS；
不得单换 ghost 插值或撤销零通量修正。未替换生产代码，方法/历史签名未改变。
零标量 Re3900 不走此标量扩散路径，结论不用于解释其热传导或长测故障。

## 7. 保留的有界扩展与发布限制

V3 observer 的独立完整性/来源绑定/增量哈希已有实现，本轮只复用。
上述窗口使用 details-output；默认 JSON 的全部 loop 数组仍随窗口增长。
64-loop 溢出仍拒绝完整归因。本窗口没有高 sweep/retry 溢出样本，
不以减少迭代或把 partial 标 complete 解决容量；有需求才设计 segment 协议。

显式初态和方法恢复不重做；自动非对流时间尺度仍需按具体显式/隐式算子设计。
完整产品峰值预算仍需补恢复物性暂存、arena 外数组、halo 和 I/O 活跃关系。
本轮 reader bulk 预算与 C++ 分配 profile 均不等于进程 RSS 上限。

干净 checkout/源码归档的构建与最小运行证据单列在
[I/O 与入口合同报告](../2026-09-07-e0fd326-io-contract-audit.md)，
不复用旧 build 的通过记录，也不为本轮构建授予原冻结长测的科学验收标签。

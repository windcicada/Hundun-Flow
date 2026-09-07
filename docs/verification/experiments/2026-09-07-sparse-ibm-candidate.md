# IBM 候选残差稀疏累加实验：不采纳

实验 diff 为同目录 `2026-09-07-sparse-ibm-candidate.patch`，**没有应用于默认源代码**。
只把两个全域 residual -= V*rate 循环改为既有唯一 interface_cells 列表；
没有改变方程、方法签名、物性、阈值、refinement、halo/consensus 或 checkpoint 合同。

## 冻结输入与实际结果

证据根目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/module-audit-20260907`。
`launch-pilot.sh`、各 frozen/FROZEN.sha256、build manifest、RUN.meta、performance.csv、
solver-rank-*.csv、health.csv、conservation.csv、Restart 与 Visit 均保留。

- 两程序均从原长测 generation-5500-64323211746098 开始，推进绝对步 5501–5510，各一轮。
- 128 ranks、456×256×52、D=0.02 m、Uc=2.89668 m/s、固定 dt=1.3808912271980336e-5 s。
- case/spec SHA256 相同，NASA/COAST-native 变物性不变；零输运标量。
- preopt binary SHA256：`387b52162b0749e1693e2cc4c45e1d11ff783e64c9b73f19db7904199bd205f6`。
- postopt binary SHA256：`8fcb3105b340f51a07e6d2f9de0089e4f54facbdb9ff17b768141262b1c22121`。
- 两份 V3 observer 均 complete=true、128 ranks、10 steps、61 loops、issues=[]。
- 128/128 rank checkpoint **完整文件逐字节相同**；128/128 Visit vtr 文件逐字节相同。
  conservation.csv 逐字节相同，health 只有 max_rank_step_ns 不同；
  solver-rank CSV 除 *_ns 和 source_meta_sha256 外所有列相同。

## 成本比较

除整段墙钟外，下表为 rank 均值，单位 s/步；子项已包含在父项中，不能重复相加。
full_step 不包含观测输出自身，整段墙钟还包含启动、恢复和末尾 checkpoint/Visit。

| 项目 | preopt | postopt |
|---|---:|---:|
| 整段进程墙钟，s/10步 | 98.94 | 99.59 |
| advance | 8.399503 | 8.416037 |
| full_step | 9.484090 | 9.563617 |
| 候选评估 | 1.973020 | 1.976085 |
| 候选 residual 子项 | 0.542895 | 0.541919 |
| 候选 boundary/derived | 0.465767 | 0.466266 |
| 候选 flux | 0.402617 | 0.412718 |
| A apply | 1.209042 | 1.203201 |
| M apply | 1.613401 | 1.606289 |
| structured halo wait | 0.249866 | 0.247200 |
| structured control consensus | 0.304520 | 0.309544 |
| 最终动量装配 | 0.207563 | 0.203143 |
| 最终 metrics | 0.005091 | 0.005005 |
| 边界守恒账本 | 0.033479 | 0.033731 |

没有可采纳的总耗时收益；residual 的微小差值也不足以据单轮实验认定提速。
因此撤回该源代码改动，保留支持域回归和本实验数据，不继续尝试第二项优化，
不增加重复性能轮数，也不以这组短窗口预测长期加速百分比。

## corrector / refinement 工作量

所有 loop 的 attempt=1、composition sweep=1，未发生 dt retry。
下表是 preopt 十步合计，两个程序的工作次数与数值收缩一致。
“候选”包括 baseline、外推和阶梯候选，收缩比为同一轮窗口内各 loop 的代表值
（中位数，不是多轮耗时中位数）。solve/A/M/refill/copy 是 rank 均值的窗口总秒数。

| 路径 | calls | iterations | A/M calls | 候选 | solve s | A s | M s | refill/copy s | 线性残差收缩比 |
|---|---:|---:|---|---:|---:|---:|---:|---|---:|
| C1 spatial | 10 | 219 | 323/219 | 20 | 8.971636 | 5.493947 | 2.444954 | 0.099043/0.014873 | 8.72e-5 |
| C2 diagonal r0 | 10 | 211 | 312/211 | 20 | 4.101591 | 1.116171 | 2.354707 | 0.107664/0.016430 | 4.08e-5 |
| C2 diagonal r1 | 10 | 329 | 486/329 | 20 | 6.254905 | 1.739734 | 3.588896 | 0.108437/0.015947 | 7.87e-7 |
| C2 diagonal r2 | 10 | 248 | 369/248 | 20 | 4.739040 | 1.318910 | 2.697368 | 0.106320/0.015813 | 5.13e-6 |
| C2 diagonal r3 | 10 | 235 | 344/235 | 26 | 4.480086 | 1.234433 | 2.571722 | 0.109080/0.016049 | 1.71e-5 |
| C2 diagonal r4 | 9 | 182 | 265/182 | 18 | 3.534827 | 0.951411 | 1.995566 | 0.095844/0.014273 | 6.44e-5 |
| C2 diagonal r5 | 2 | 44 | 66/44 | 4 | 0.829923 | 0.235814 | 0.480797 | 0.022160/0.003136 | 5.22e-5 |

C2 共 51 次 diagonal、1249 次迭代；没有 C2 spatial 样本。
MG refill/copy 计时说明复制不是这里的主要总成本，不支持冒险 swap 持久借用地址。
下一项应先区分 C2 refinement 方向带来的额外次数与 M apply 本身成本；
不能直接提前 spatial、删 collective 检查或放宽守恒/残差。multidot 不列入待办。

## 稳定性与最终动量观测

十步全部一次接受、BDF2、无恢复/重试；连续性最大 2.346e-7，能量最大 9.350e-7，
EOS=0，committed CFL(out) 最大 0.284829；两个残差门槛均保持 1e-6。
短窗口通过不等于长时间稳定性或统计收敛已验收。

step=5510，动量归一化为 `|R_m| / (a0*rho*V*U_rms)`，
`U_rms=sqrt(2*K/M)=2.96333876685 m/s`。这是时间惯性力尺度，不是新的接受门槛。
最差单元按归一化量定位，不能把它和全局绝对力残差的最大单元混同。

| 分量 | 最差 gid/rank | 全域归一化 Linf | IBM邻近 Linf/RMS | 内部 Linf/RMS |
|---|---|---:|---|---|
| x | 1929506/69 | 0.00813273 | 0.00656733 / 0.000185784 | 0.00813273 / 0.0000605813 |
| y | 2269140/50 | 0.00724030 | 0.00724030 / 0.000137205 | 0.00627236 / 0.0000466631 |
| z | 1798635/53 | 0.00409362 | 0.000556396 / 0.0000196582 | 0.00409362 / 0.0000219190 |

IBM 邻近区为触及控制 link 的 6500 个活跃单元，内部为其余 5,987,852 个活跃单元；
不是按距离随意划区。此表不声称动量已达到某个未定义的经验门槛。

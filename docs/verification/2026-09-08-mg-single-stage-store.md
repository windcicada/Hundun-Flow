# 单stage预平滑的中间方向写入实验

状态：**算法实验已撤回**。局部、sanitizer和干净构建验收通过；单轮产品数值等价，
但完整墙钟和advance没有收益，故恢复基线实现，保留新增回归及证据。
本轮只评估这一项性能候选，没有新增物理模块；两相继续暂挂。

## 测量依据与最小范围

[新版MG观测](2026-09-08-runner-mg-profile.md) 的9502–9510窗口中，
pre/post平滑占Native MG的51.32%，最细层pre为0.525153 s/步。
这只是父阶段成本，不是已测得的direction写入成本。
`solver_piso_detail.hpp::production_pressure_mg_policy()` 固定pre=1、post=2、Chebyshev F-cycle；
不能以公共policy默认值推测运行策略。

`solver_mg.cpp::chebyshev_point_smooth_streamed()` 的单stage pre具有以下生命周期：

1. 读取solution stencil，计算scaled、direction和next solution；StageZero不读取旧direction。
2. 将next写temporary，随后仍按原路径copyback到solution并交换solution halo。
3. 保留的最终stencil重新写整个residual内区，restriction才消费此最终缺陷。

第1步原来还将direction写入residual，但第2步没有读取它，第3步将其全部覆盖。
候选通过内部编译期`StoreDirection`仅省去`stages==1 && retain_final_defect`的这次写入。
不是改成原地Jacobi，也没有省略direction的计算或有限性检查。
更长recurrence、post平滑、line relaxation和终端求解仍走原路径。

| 保持不变的合同 | 核查 |
|---|---|
| 数值方法 | stencil算术顺序、scaled/direction/next、invalid_mask、Chebyshev参数、次数、F-cycle、最终投影均不改 |
| 存储与通信 | 不增加全场数组/动态分配/公共选项；原view地址、copyback、halo调用、revise和prepared关闭保持 |
| 失败 | 中间工作区不是已接受载荷；错误仍沿原汇总边界返回，最终correction不发布 |
| 历史 | 未改变存储状态、接受速率或演化语义；方法历史签名不变，仍需产品逐字节及重启验证 |

这是有源码数据流依据的性能实验，不是一个已发生的数值故障修复。
不存在“旧版算错、新版算对”的RED；不为制造RED而断言私有函数调用次数。
实际性能接受信号是同物理窗口的完整墙钟与已有MG成本，测量结果见后文。
从已测level visits和单元数可推算逻辑写入槽位，但没有测硬件store流量，
不能将缓存中的重复写入直接当成DDR带宽或峰值内存节省。

## 公开接口验证

复用已确认的`NativeCartesianMgPlan::compile/apply/update_coefficients`及ProductDriver接口，
数值oracle沿用已有独立degree解析、legacy生命周期和实际MPI测试，不暴露新的公共控制。

- 修改算法前，扩展degree=1的NaN/Inf失败不发布覆盖，Chebyshev测试通过，0.35 s。
  [基线](data/2026-09-08-mg-single-stage-store/baseline-chebyshev.txt)是正常基线，不标为失败重现。
- 初版候选Release **11/11通过，5.89 s**：Chebyshev、reuse、Krylov隔离、MG1/2/4、
  update合同2/4及ProductDriver retry1/2/4。
  [日志](data/2026-09-08-mg-single-stage-store/release-initial-LastTest.log)。
  覆盖degree1–4、非立方/小x、混合边界/周期/非均匀系数的correction、残差与Status逐字节oracle等价。
- 单rank Chebyshev内部非有限状态注入扩展到pre=1，并保留pre=2；Release MPI **3/3通过，1.32 s**。
  已有断言检查全rank一致错误、外部correction不变、无热路径C++分配及通信调度。
  [日志](data/2026-09-08-mg-single-stage-store/release-expanded-mpi-LastTest.log)。
- Clang Debug ASan+UBSan **5/5通过，16.21 s**：Chebyshev、reuse及MG MPI1/2/4。
  [日志](data/2026-09-08-mg-single-stage-store/sanitizer-LastTest.log)。
  `detect_leaks=0`；不声称LSan或进程RSS硬上限通过。

所有MPI串行，原128-rank长测维持SIGSTOP。开发Ninja目录再次出现日志恢复提示，
所以以上不能充当新可执行文件的发布身份；必须在独立干净checkout再次构建验收。
没有改变当前Re3900物性、网格、dt、残差门槛、refinement或checkpoint持久化合同。

## 干净验收与单轮产品结果

候选DCO提交 `ed9b09ab2a019a7957030844e5437421351fd044`，
tree `867be9a0697f71c1aed8953c5ad08f31fda917b3`。
独立干净checkout `/home/wyf/code_dev/.worktrees/hundun-flow-mg-single-stage-accept-20260908`，
新Ninja Release、Clang15/libc++、`-march=znver3 -mno-fma -ffp-contract=off`，`-j4`，
tests开启，ASan/UBSan/Hypre关闭；构建没有开发Ninja日志恢复问题。
全部MG测试、ProductDriver retry1/2/4及真实runner MG CLI1/2/4，**20/20通过，34.50 s**。
见[干净日志](data/2026-09-08-mg-single-stage-store/clean-acceptance-LastTest.log)与
[manifest](data/2026-09-08-mg-single-stage-store/clean-build-manifest.txt)。
runner SHA-256为 `a6f2f2d201bf1067e0bebc99bf2a560c85b852680dc48da18c0d16dbd8ea90fa`，
manifest为 `4269d8279ab316f60e93ec9e5ee79071f33b38fae1aa8c05d32cf2b6d8636ba4`；
head/tree前缀摘要已独立重算，源码清洁标志为true。

冻结配置、程序和检查器后，只运行一次128-rank、9500→9510，仍为原网格、变物性、dt、
阈值和MG观测开关。运行目录是
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/pilot-mg-single-stage-9500-9510-20260908`；
原始观测和收据在同级 `mg-single-stage-20260908`。基线为已经完成的MG观测窗口，未再跑基线。

10步BDF2、无retry，128 ranks/70 logical loops/6层账目完整；runtime validator返回0。
末次generation为 `generation-9510-161548994678342`，Visit和日志正常完成。
**128份rank checkpoint、manifest、128份Visit、statistics、accumulator与基线逐字节一致**；
每rank全部70行solver的非计时列（排除来源hash和`*_ns`）也完全一致。
complete附件仅合法generation名称不同。133份原checkpoint/统计文件与全部冻结输入哈希不变。
见[数值/载荷收据](data/2026-09-08-mg-single-stage-store/PILOT_ACCEPTED.json)。

| 相同口径 | 基线 | 候选 |
|---|---:|---:|
| 完整进程墙钟 / s | 110.86 | 111.97 |
| 9502–9510逐步max-rank advance均值 / s | 9.526839 | 9.601104 |
| rank-mean Native MG / s每步 | 2.053430 | 2.100171 |
| rank-mean MG最细层pre / s每步 | 0.525153 | 0.472403 |
| rank-mean MG全部层pre / s每步 | 0.703746 | 0.710865 |
| rank-mean MG嵌套halo wait / s每步 | 0.353292 | 0.448488 |
| 后9步iterations / A / M | 1609 / 2358 / 1609 | 1609 / 2358 / 1609 |

最细层pre局部变快，但全部层pre、MG、advance和完整墙钟没有改善。
未剔除较慢loop、未增加重复轮次，也不把单轮通信等待差值归因于某个已证实的系统或算法原因。
按预定规则撤回，而不是将该窗口外推为长期变慢百分比。
[性能决定及完整分层数值](data/2026-09-08-mg-single-stage-store/PERFORMANCE_DECISION.json)
绑定两侧来源与成本hash。

后处理首版`finalize-observation.py`曾因脚本生成时的字符串替换产生SyntaxError，
发生在Python执行前，没有生成收据或改写运行数据；原错误副本保存在audit目录，
修正后只重跑这一步读取核查，未重复CFD。来源和后处理hash分别由不变的FROZEN与新增ANALYSIS绑定。
复用的物理核查器输出含旧“observation-only”描述，本候选的实际范围以本报告和数值收据为准：
它是性能实验；固体只新增区域极值核查，未再做源9500到终9510的固体专用全载荷检查。
归档脚本保留当次audit路径合同，不是离开原始数据目录即可运行的通用程序。

## 撤回确认与后续

只撤回`solver_mg.cpp`中的算法改动，保留单stage异常覆盖。
该文件与实验前`5b993ea`逐字节相同，SHA-256
`f37cd9be2053715138fe91136a955bfcb60826c6e02128a96ded905150772555`。
恢复后Chebyshev及MG MPI1/2/4 **4/4通过，1.60 s**，见
[日志](data/2026-09-08-mg-single-stage-store/withdrawn-LastTest.log)。
没有再次运行相同128-rank配置。实验冻结程序与数据保留只读，不移作生产验收标签。

下一项先回到C2-r1工作次数：基线后9步记录58次`norm_breakdown_restarts`，
须核对该计数的实际分支和与真实残差重建的关系，再判断是否存在可减少的重复工作。
计数较多不等于已发生数值不稳定；不先放宽阈值、删保护或更换求解器。
原长测仍暂停，没有推送，也未宣布COAST替代完成。

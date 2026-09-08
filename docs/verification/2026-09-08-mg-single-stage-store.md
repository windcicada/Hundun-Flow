# 单stage预平滑的中间方向写入实验

状态：局部数值与sanitizer回归通过，**尚未完成干净构建和产品性能接受**。
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
实际性能接受信号是同物理窗口的完整墙钟与已有MG成本，尚待候选单轮测量。
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

## 待接受步骤

候选分项DCO提交后，以相同Clang/Release参数从干净checkout构建，复核相关公开测试及真实CLI。
通过后仅运行一次128-rank、9500→9510，保持原输入和MG观测开关；
使用已完成的MG观测窗口作为同配置基线，不再重复基线或三轮取中位数。
先比全部rank checkpoint、Visit、统计、守恒、工作次数，再看完整墙钟和分层成本。
无总成本收益则撤回算法实验；小幅单轮差值不外推为长期加速保证。
未完成这些步骤前不恢复长测、不推送、不标为可替代COAST。

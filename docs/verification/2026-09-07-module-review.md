# 分模块核查：第一批修复与待验收项

审阅/实现起点与本轮查询的远端 main 均为
`0127616f3716b7d7d6fcb4a00831129c2d6631b6`。
工作树起始只有未跟踪的 `.codegraphf/`，没有源文件修改。

本轮没有编译、启动 MPI 测试、暂停/替换冻结长测或推送未经验证的代码。
`hundun-re3900-scalar-repaired-20260906.service` 在检查时保持
`active/running`，MainPID 仍为 119112。冻结二进制、算例、checkpoint 与原始日志未改动。
Python 检查为低优先级离线后处理，不是 CFD 试算。

## 1. 标量与共享容量

静态确认：`ScalarMassRemap::prepare_passive_intervals()` 的一批上下界归约需要
`2*N_passive`；共享 ReductionEngine 原来仅考虑两个 Krylov 需求与 8。
这会返回容量错误，不是已经复现的内存越界；本轮没有运行 C++ 复现。

关键 diff：

- `solver_scalar_mass_remap_detail.hpp::reduction_requirement()` 由消费者申报槽数，
  检查空指针和乘法溢出。
- `core_product_freeze.cpp::ProductCompiler::compile()` 在 collective 汇总前取得需求，
  冻结容量为 `max(primary Krylov, auxiliary Krylov, 8, scalar requirement)`。
- 不改 `solver_linear_state.cpp` 的 Krylov 自身需求，不增大 restart、不减少合法标量数、
  不引入热路径扩容。
- `product_scalar_contract_experiment.cpp --capacity` 新夹具覆盖 4/5/12 个被动标量、
  中间插入 EOS 组分、FGMRES restart=2 / BiCGStab，以及各字段常量保持和非法时间提案回退。
  CTest 注册 1/2/4 ranks。**这些 C++ 测试已添加但未运行**。

容量交集核查：

| 资源 | 现有边界 | 本轮结论 |
|---|---|---|
| 产品字段登记 | 输运标量上限 253 | 不是最终可运行目录上限 |
| IoServicePlan / restart 状态目录 | 64 个字段；产品固定 U/pi/h 三个 | 当前产品交集最多 61 个输运标量；62 个在构建 I/O 计划时已拒绝，不是热路径溢出 |
| halo 字段表 | 按字段目录冷分配，FieldId 容量检查 | 未发现 4/5 个被动标量处另一个固定 8 槽限制 |
| Driver halo_views | `16 + N_scalar` | 覆盖现有各阶段拼装基址；混排映射沿用已有修复 |
| 数值失败质量分数诊断 | 64 | 对当前 I/O 交集不构成更紧约束 |
| 归约上下界 | `2*N_passive` | 本次纳入共享冻结容量 |

剩余：新增夹具不是非均匀全局极值回归，也不是多标量真实数值失败回退的完整验收。
已有 `--restart` 的真实数值失败回退记录属于旧的两标量夹具，不能冒充新增目录回归。
独占窗口还需核对非均匀范围、更多目录排列及输出/恢复交叉情形。

## 2. IBM 公共接口及边界

静态确认 6 处 `validate_bound()` 后、失败判断前调用 `kernels_->cells()`：
`constrain_pressure_predictor`、`constrain_momentum`、`correct_pressure_gradient`、
`correct_pressure_work`、`correct_zero_normal_diffusion`、
`correct_positive_bounded_zero_normal_diffusion`。已添加提前返回。
原 `zero_interface_flux` / `validate_interface_flux` 依靠短路求值安全，
`correct_velocity_gradient` / `correct_impermeable_scalar_diffusion` 已有提前返回，未重复改顺序。

三个 diffusion 方法使用现有 `field_views_overlap` 拒绝 rate 与 transported/diffusivity
的地址区间重叠，包含 ghost/stride 和不同 base 的部分重叠；两个只读输入可以重叠。
velocity gradient 同样拒绝写入区与输入 velocity 重叠。没有改 IBM 通量或曲面重构公式。

`solver_ibm_equation_interface_test.cpp` 增加默认对象、编译失败新对象的公共热接口调用，
要求 invalid_plan 且哨兵不变；增加三个 diffusion 的两种部分重叠夹具。
`v04_ibm.hpp` 明确本地错误返回与存储合同。**未执行 C++/sanitizer/MPI 回归**。

剩余：平壁/斜壁/曲面制造解及空间阶数仍未完成；不能用二值控制体零通量代替曲面二阶证据。
现有开边界/倒流/周期/IBM 标量回归保留，需在本轮候选版本上重新执行；未撤销守恒修复。

## 3. Remap / Picard

`solve()` 确实逐 Jacobi 扫描重算最终质量、六面增量与对角，读取一层六面邻居；
独立 halo 目前仍使用 `schemes.required_ghost_width()`。本轮没有缩减物理 ghost 或增加六份缓存。
`capture()` 每次 execute_attempt 都重新拷贝配对预测面通量、M*、Q*；
组成猜测在 capture 之后写入当前 species，不能仅靠相同步号作为跳过预测的证据。

`mass_pairing_residual` 没有独立接受门槛；上游
`pressure_energy_components_converged` 对有标量路径取
`min(public continuity tolerance, 16*epsilon)`，并保留配对预测通量。
这是现有保证链的一部分，不等于已证明所有分支上的配对行和；本轮没有构造出守恒反例，
也没有把“只有记录字段”写成已发生守恒错误。

组成 sweep 在同一 dt 下拒绝未提交事务并重新进入 attempt；真正 dt retry/恢复另走状态机。
未来若缓存，键至少覆盖已接受/历史状态与速率 revision、BDF/EX 系数、源项、limiter、
边界分支/权威、几何和配对通量证书。当前零标量 Re3900 的 remap 计时为 0，
不能用它估计通用组分路径的迭代、算术和 halo 收益。

## 4. 新版成本与下一项唯一优化候选

本轮重新读取 `scalar-recovery-pilot-20260906` 原始 CSV，给历史诊断显式指定 128 ranks，
逐步检查 1001–1010 的每个 loop 覆盖全部 ranks，结果为 63 loops：C1 spatial 10，
C2 diagonal 53，C2 spatial 0。由于旧 RUN.meta/CSV 没有本轮新增的来源绑定，
新版 observer 正确输出 `complete=false`，即使覆盖与账目通过，也不追认其为 V3 完整证据。

下表由已有同一窗口实测数据计算，单位 s/步、rank 均值，不是新程序加速测量：

| 项目 | s/步 |
|---|---:|
| advance | 9.013555 |
| 已记录 full_step（不含观测自身 gather/write） | 10.110506 |
| C1 spatial solve | 0.985524 |
| C2 diagonal solve | 2.814590 |
| C1 A / M apply（solve 的子项） | 0.580947 / 0.273753 |
| C2 A / M apply（solve 的子项） | 0.741720 / 1.630431 |
| 全部候选评估 | 2.074997 |
| 候选 residual / boundary_derived / flux（评估的子项） | 0.560043 / 0.493636 / 0.428576 |
| MG refill / copy（copy 是 refill 子项） | 0.066628 / 0.010011 |
| 最终动量装配 / terminal metrics / boundary ledger | 0.201297 / 0.005333 / 0.032035 |

不把子项重复相加，不套用旧方法百分比。MG coefficient 更新保留现有结构、在 inactive
区构建并经 consensus 后复制到固定 active 地址；这不是每次重建结构的证据，
也不支持直接 swap 已借用地址。prepared 结束前汇总与 MG 生命周期检查均未删除。

下一轮只选 **候选残差装配中的重复工作** 作为局部优化对象，不同时改后备方向、
MG 存储和 remap。先将重复计算缩到一个可证明输入不变的最小 diff，再做同方法/同物理
窗口的一轮对比。当前聚合计时还不能证明某个具体缓存有效，所以本轮没有提交算法优化，
也没有预报加速比例；multidot 不是待办。

## 5. 真实内存预算

`owned_payload_bytes` 已统计 remap 的 mass/quantity/next、每标量 padded storage、
配对 flux 和小容器容量，halo 则单独计数。但 remap 在 ProductDriver::create 分配，
晚于 CompiledCasePlan 的图工作区汇总；未找到把这些所有权及恢复 image / checkpoint
staging 的同时活跃集合统一用于产品硬预算的路径。不能称预算问题已修复。

后续预算需按唯一物理分配及阶段求峰值：arena（含 Krylov/MG）、arena 外 Driver/Picard
数组、独立 halo/flux、恢复 owned image、输出 staging。view 不另计，复用前检查事务回退、
revision、证书和未完成请求。现有 Visit 预建路径、局部阶段异常保护、载荷移动不重复改。
新增标量分配失败、恢复/输出、销毁后重建的动态检查留待独占窗口。RSS 高水位不是泄漏判据。

## 6. 应用与方法恢复

普通 ApplicationService 仍仅暴露 restart_storage_compatibility，传给
restart_expected/initialize_restart 的历史策略默认 require_compatible。
当前需要显式方法恢复的调用者应使用 ProductDriver 的 RestartHistoryPolicy，或专用 runner
的 `--restart-method-recovery`；普通入口没有此能力，不能自动绕过未知历史。
本轮只记录接口边界，没有在未测试时向普通应用加入新恢复策略。

V3 语义签名位于 `core_product_freeze.cpp::method_history_signature`，现有组件包括
BDF/EX、rho-h-p、scalar split、IBM thermal、momentum rates、fresh flux、C1 joint target、
open/periodic flux/metrics、AFC 算术、conditional boundary；有标量时再加入 paired remap、
Picard、mass roundoff closure、passive envelope、physical donor、inner accuracy 与 IBM
impermeable flux。此次容量/错误边界/观测修改不改变合法路径的离散历史含义，未绑定 Git SHA
或修改签名。独立已知迁移表仍待补充。

`runner_statistics_epoch_test.py` 已把链改为 source → 同方法 exact → 显式方法恢复 → 新 epoch
再次 exact；检查非零旧样本清空、新发展窗口、BE/BDF2、来源只读，并补 V3 观测来源/覆盖断言。
注册的 1/2/4-rank 真实 CLI 回归**尚未运行**；质量目标、速率来源的完整链式验收仍需补齐。

初场仍由最后遍历到的有效温度/入口覆盖；非对流时间尺度仍由调用者提供，只有对流约束接入
已接受状态。另见 `ProductDriver::initialize` 对所有标量统一限制 [0,1]，与已有有符号被动标量
源项/恢复能力之间存在接口范围差异，需在明确初场合同后处理，不能让容量夹具误触非法输入。
本轮容量夹具使用不同的合法 [0,1] 常量。当前固定 dt 的 Re3900 完全不变。

## 7. 观测完整性与长日志

实际 RED：`solver_observe_cli_test.py` 在 RUN.meta 声明 2 ranks、两份输入同时仅保留 rank 0
时，旧 observer 返回成功并给出 `ranks=1` 的归因；测试按要求拒绝而失败。

修复内容：

- runner 新 RUN.meta 记录 expected_ranks、requested_steps、observation_schema=3，
  加入运行启动时间和 launcher PID 区分重复运行；在纯本地阶段编码/计算完整元数据 SHA256，
  全体同意后一次广播固定 digest。
- loop/performance CSV 每行携带 `source_meta_sha256`；普通性能汇总器兼容这一非数值列。
- observer 以冻结 rank 集合和完整起止步为准；每个逻辑 loop 分开验证全 rank 覆盖、重复、
  调用次数、A/M 时间及完整步时间账目；拒绝负数、非有限值、非法判别量、截断尾行和混合来源。
- 失败/历史数据只能显式 `--allow-partial`；旧元数据另需 `--expected-ranks`，
  永不升级成 complete。严格失败不留下本次创建的结果/明细文件。
- 输入按步流式读取，SHA256 分块计算；`--details-output` 逐 loop 输出 JSONL，
  默认单个 JSON 的 loops 列表仍会随明细增长。源文件从不改写。
- 保留 64-loop 溢出拒绝。16 次组成 sweep 可能超过该容量；没有用 dropped 数据计算完整归因，
  也未擅自扩大运行时诊断数组。有界分段输出/诊断窗口仍待设计。

来源哈希增加每个观测行 65 个 ASCII 字节；按旧 pilot 平均 6.3 loops/步、128 ranks，
loop+performance 约新增 60,736 B/步（约 59 KiB），不是零 I/O 成本。
新 header 的一次性开销另计。已冻结长测不会生成这些新列。

实际 GREEN 命令：

```sh
nice -n 19 python3 versions/v0.4/tests/integration/solver_observe_cli_test.py tools/v04_solver_observe.py
python3 tools/v04_performance_observe.py --self-test
git diff --check
```

28 项 CLI 检查通过；后处理 8/2048 步各一轮、2-rank 合成 CSV，tracemalloc 峰值分别为
2,801,155 / 3,759,456 B。这是 Python 分配跟踪，不是进程 RSS 上限或 MPI 求解峰值。
性能账目 self-test 通过。新 runner 产出及全部 C++ 修改仍须编译和真实 CLI/MPI 验收。

## 下一安全边界

本地分项提交（均带 DCO，未推送）：`ffcba11` 为容量/IBM 候选修复，
`c8aab80` 为观测完整性与流式处理，`58eefc9` 为链式重启回归。
提交并不代表 C++/MPI 验收完成；没有算法实验或生产程序替换提交。

本轮是待验收的第一批修改，不是七项全部完成。长测继续占用 128 ranks 时，不执行以下命令；
等待独占窗口或另有明确资源授权后，再构建并按风险运行：

```sh
ctest --test-dir <candidate-release-build> --output-on-failure -j1 \
  -R 'v04_scalar_contract_capacity_mpi_[124]|v04_solver_ibm_equation_interface|v04_runner_statistics_epoch_mpi_[124]'
```

随后补非均匀多标量范围/数值回退、IBM 制造解、分配失败及方法历史完整链；
正确性验收通过后才对上述唯一算法优化候选计时。不得拿现有零标量长测代替通用组分验收。

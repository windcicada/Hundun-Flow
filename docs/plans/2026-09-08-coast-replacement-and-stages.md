# COAST 替代与 Stage 4–6 接入目标台账

更新日期：2026-09-08。状态：**进行中，尚未完成替代验收**。

当前范围调整（用户 2026-09-08 新指示）：**不中断正在进行的基础流动优化与验收，
暂时不加入两相模块。** Stage 6/S0–S2 暂挂，历史源码和方案保留只读；
即使燃烧验收结束，也须用户恢复该范围后才接入两相。Stage 4/5 燃烧仍安排在基础流动验收之后。

目标是在用户当前明确的业务范围内，完成 HUNDUN-FLOW 的 COAST 替代验收：先验收基础流动的正确性、稳定性、性能、模块合同与 I/O，再接入并完成已有 Stage 4 燃烧、Stage 5 ESF/TPDF/TCR。Stage 6 稀相喷雾仅保留后续方案，等待用户恢复范围。完成标准是实际 CLI 推进、共同事务、输出和 Restart 的组合行为通过约定案例；历史 seal、配置解析成功或纯核测试不自动成为当前产品验收。

本文件是持续更新的目标台账，不是发布收据。它集中记录当前工作、退出条件、源码来源和未关闭项。详细推导和原始失败留在链接的证据报告中；旧报告中的“正在运行”“尚未测试”和下一候选仅描述当时状态，以本文件的当前状态及后续验收收据为准。

## 1. 当前授权、基线与运行边界

| 项目 | 当前记录 |
|---|---|
| 工作区 | `/home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity` |
| 本次盘点基线 | `86542bb96678ec844ae5ac95d8f6391993da239e`；其生产源码继承 `48493ed2b8655440228e85396b26a070f4e5103f`，两者之间仅有验收资料 |
| 工作顺序 | 基础流动验收 → 当前接口下 Stage 4 → Stage 5；Stage 6 暂挂，当前不接入；每次冻结一个活动切片 |
| 实验次数 | **每个配置仅运行 1 轮**；不采用旧计划的三轮/九对重复、median 或置信区间门。未来若需要统计性能验收，应先单独确认协议，不作为默认执行 |
| 比较约束 | 保持所比较案例的物理、网格、固定 dt、容差、refinement 容量、工作量和输出策略；记录两侧不同的数值方法，不用网格缩放或迭代数差异冒充同工况速度 |
| 协作边界 | 保留其他任务的修改；历史仓库、冻结程序、源 checkpoint 和旧证据只读。当前授权不要求重复询问每个 Stage 的普通接入步骤 |
| 新增范围 | AMR、移动 IBM、稠密喷雾、液膜、多组分液滴、生产 GPU 等不由“替代”二字自动扩展进入本台账；需要业务需求再登记 |

当前长测状态（2026-09-08 交接后）：原 `48493ed` 作业已从暂停的 9951 内存状态推进到 durable 10000，保存并校验后正常行政停止；原 9000/9500 内容另存，未丢弃已有进度。新 `hundun-re3900-observed-long-20260908.service` 使用已验 `7c03a54` 二进制，在新目录以 128 ranks 同方法续算 10000→35000，无重新方法恢复或统计 reset。已验启动健康前缀为 10001–10016，V6 观测完整到 10017，但长测仍未完成。下一 checkpoint/Visit 为 10500；见[交接记录](../verification/2026-09-08-long-run-handoff.md)。以下各切片的旧暂停记录保留为历史状态，不代表当前进程仍暂停。

## 2. 活动切片与下一动作

| 项目 | 状态 / 完成条件 |
|---|---|
| 已完成切片 | RHS norm 与实际线性停止阈值观测；局部、干净候选和单轮 9500→9510 观测通过 |
| 当前活动切片 | FGMRES显式保留后续列实验已撤回，原路径及保护回归通过；稳定版本长测已安全交接并启动。当前保留128-rank独占运行，不再与测试/编译争抢算力 |
| 新观测字段 | `linear_criterion_valid`、`linear_rhs_norm`、`linear_atol`、`linear_rtol`、`linear_residual_limit` |
| 语义 | 记录 RHS norm 与实际线性停止阈值；初猜残差不能替代 RHS norm。它与已有补充物理审计的 `convergence_limit` 分开 |
| 兼容性 | 仅performance为 `observation_schema=4`；显式MG为5，新恢复观测为6（MG可选），读取器兼容3–5。旧数据不能补造缺失原因、阈值或MG层级成本 |
| 当前回归发现 | BiCGStab 测试夹具已按既有 `fixed_general` 合同纠正；新 observer 初版误拒纯绝对容差和误接受精确零阈值附近非零值，两项均 RED→GREEN。均未改变生产求解控制 |
| 下一切片退出条件 | 将现有norm_breakdown_restarts与实际分支、真实残差重建和额外A/M工作对齐；需要新观测时保持定长、默认关闭和失败回退；有局部证据后再选择一个最小优化 |
| 下一动作 | 核对10500的新checkpoint/Visit/统计附件，继续检查残差、RSS与新版分模块工作量；形成完整窗口证据后再选择一个最小优化。长测启动不等于35000稳定性或COAST替代验收，04badb1实验不进入生产 |

实际线性停止阈值为 `max(atol, rtol*||b||)`，其中 `||b||` 是本次求解真正采用的 RHS 范数。

本轮已闭合切片的身份与边界：

| 本轮验收字段 | 记录 |
|---|---|
| 候选代码 / tree / binary | 观测提交 `16fd8f085e845628a545dc7d4d371b3064e9ef38`，tree `8a53030af5f42f63dff6dac7fd3ce134aaa8e0c9`；干净 checkout `hundun-flow-criterion-accept-20260908`；runner hash `8cfc9fe3a613a84e677fe536a353779040b42bf525650ae10373616e23d33545` |
| 已运行测试及原始日志 | Release 7/7、ASan+UBSan Krylov MPI 3/3、干净 Release 7/7；最终 observer CLI 88 checks。[回归与单轮观测](../verification/2026-09-08-linear-criterion-observation.md)；不是全产品新验收 |
| 单配置单轮观测窗口 | 9500→9510，128 ranks，仅一轮；70 logical loops，全 rank/步来源完整，10 步 BDF2、无 retry；末次 Visit/checkpoint/关闭检查通过，133 个源文件不变 |
| 数值不变性、额外开销、是否接受 | 源码只增加观测，公开测试及该物理窗口通过；未做算法 A/B、COAST 比较或全载荷逐单元不变性检验。后 9 步 max-rank advance 均值 9.562760 s，rank-mean M 2.015330 s；不宣称加速 |
| 长测后续状态 | 原长测仍保留 128 个 SIGSTOP ranks，完整 health=9950，durable checkpoint=9500；新 pilot checkpoint=9510 只属于独立观测目录，不替换较前进的原进程 |

MG 观测的方案边界：复用 `NativeCartesianMgPlan::level_count()/level()` 与现有
32 层上限，plan-owned 定长 profile 默认关闭；开关和计时允许 rank-local，
不改变 collective fingerprint 或分支。各非递归阶段之和不超过 apply inclusive；
Halo/reduction/直接 MPI 是嵌套子项。保留 prepared epoch 入口共识、最终 checked-sum
汇总与发布，关闭时不调用新计时器。目标测试包括 `v04_solver_mg_mpi_[124]`、
`v04_solver_mg_update_contract_mpi_[24]` 及现有 reuse/line/coarse/Krylov 隔离入口。
Native 层的 [开发回归与计时合同](../verification/2026-09-08-mg-apply-profile.md)
已落地；产品逐层实测见下文，尚无新算法性能收益结论，不能把公开接口回归当成替代验收。

[ProductDriver 接线](../verification/2026-09-08-product-mg-profile.md) 使用一个4664字节
固定baseline和200字节每-loop摘要；关闭仍有固定存储/报告拷贝成本，不宣称零开销。
扩大回归发现并实测确认两个[旧Restart夹具](../verification/2026-09-08-restart-inactive-flux-fixture.md)
把solid-solid负零错误地当成权威流体通量；独立干净基线也失败。修正的是测试区域合同，
没有改变生产normalizer或方法签名。

[Runner/observer 接线](../verification/2026-09-08-runner-mg-profile.md) 已完成：
CLI与统计恢复链6/6，ASan+UBSan CLI 3/3，独立干净13/13；root及非零rank的新增日志关闭失败一致返回，
checkpoint字节不变。新observer有33项拒绝检查，旧88项保留。
Native开关仍可rank-local；runner因文件/通信分支要求两个观测开关冷入口一致。
候选`461d763`单轮128-rank、9500→9510完成，总墙钟110.86 s；10步/70 loops/6层来源与账目完整。
128份checkpoint、manifest、128份Visit及统计载荷与此前同起点criterion窗口逐字节一致。
新generation仅属于pilot，不替换原SIGSTOP进程。
后9步Native MG为2.053430 s/步，pre/post为0.703746/0.350072，terminal为0.200820；
51.32%的MG成本在平滑，copy只有refill内的0.011250 s/步。新观测不是加速或COAST替代证明。

[单stage中间写入实验](../verification/2026-09-08-mg-single-stage-store.md)已闭合并撤回：
候选`ed9b09a`干净20/20、ASan/UBSan5/5，通过128-rank同起点窗口的全部输出字节和
非计时solver列比较；总墙钟110.86→111.97 s，后9步advance均值9.526839→9.601104 s。
未显示总成本收益，不追加重复轮次。源MG文件恢复实验前内容，撤回后4/4通过；
保留新增单stage异常回归与全部证据。下一项使用基线r1的58次递推恢复记录作定位入口，
不能把该计数直接称为58次数值失败。

下一入口的静态核查：`solver_krylov.cpp::solve_fgmres()` 在
`unsafe_recurrence && column!=0`结束当前有效子空间后，以及`happy_breakdown`未达到最终
容差但仍有进展时，都会增加同一个`norm_breakdown_restarts`；普通restart和column=0的
显式重正交化不等价于该计数。现有58次不能区分这两个来源，尚不据此决定数值改动。

[Native恢复原因观测](../verification/2026-09-08-fgmres-recovery-observation.md)已完成开发回归：
公开`solve_fgmres()`默认null的96字节借用sink，区分负范数、丢弃列、显式恢复、
两类恢复restart、普通长度restart和五类A。RED捕获缺失观测，Release与ASan/UBSan
各1/2/4-rank通过；仅root/非零rank开启时解、归约和外部A/M工作不变。
非注入单位算子用例实际出现1次happy restart，不能把总计数统称数值失败。
代码提交`74e262d`独立干净生产构建及Krylov/PISO等14/14通过（4.84 s），
runner hash为`6ebc5559d54276cf17385e72eb02c7b0002925b0e0168c3334b58e9b15e0a78c`。
上述为Native切片的历史验收范围；后续接线见下文，仍未测目标原因比例或改变算法。

[ProductDriver/runner V6接线](../verification/2026-09-08-product-fgmres-recovery-observation.md)
已通过开发验收：Release 8/8、ASan/UBSan 4/4，旧reader 88/33 checks和新reader
28项拒绝检查。实际BiCGStab CLI独立验证恢复unavailable而非伪造零样本。
每个solve记录新增96字节，保持14-slot attempt和64-slot step上限；关闭仍有固定
初始化/拷贝成本。新列复用既有stream，root/非零rank关闭失败和checkpoint字节
等价均已验。DCO候选`7c03a54`独立干净31/31通过，191.10 s；含真实BiCGStab
的ASan/UBSan runner扩展也通过。新128-rank单轮窗口已接受，总进程113.58 s；
70 loops完整，末次128份Visit/checkpoint及统计与旧MG窗口逐字节一致，源133文件不变。
后九步347次unsafe丢列占1609次Arnoldi的21.57%，C2-r1有58次；没有happy restart。
这是工作量比例，不是可直接消除的时间或速度收益。r1单次M并不更贵，下一局部假设
针对保留已付A/M成本的方向；目前只是明确了下一实验，不将新算法混进已验观测提交。
原长测仍保留128个SIGSTOP ranks，pilot的9510不替换原9950之后的内存状态。

[后续列显式保留实验](../verification/2026-09-08-fgmres-salvage-trial.md)已闭合：
独立分支`04badb1`在无注入的小缩放非对称算例出现工作量/可接受性退化，
1e-8算例原44次收敛，候选80次仍拒绝。`33233d2`撤回数值源码，保留反例回归；
Release 3/3与主工作区ASan/UBSan 3/3通过。原路径并非在这些合法算例上失败，
不将保护计数多简单等同为错误。主工作区生产源码仍与7c03a54相同。
接下来按用户“满足长测需求”的优先级完成稳定版本交接；更深优化另立证据切片。

## 3. 基础流动模块台账

状态含义：**已修并验证**限定于已引用的历史候选和用例，本次仅只读复核；**实测限制**有数值证据；**静态缺口**由接口或调用链确认；**待测优化**尚无收益结论。后续代码变化只使受影响的证据重新待验，不抹去历史失败和通过记录。

| 模块 / 合同 | 已闭合项，避免重复实施 | 尚未关闭的边界 | 证据 |
|---|---|---|---|
| 动量、压力—焓与时间推进 | M1 新动量面通量、M5 热时间响应、M8 接受热速率、M9 C1 联合目标、M10 周期通量已修；最终动量与累计能量观测已接入 | `full_nonlinear_jacobian=false`；Schur 未含完整 EOS 密度对流、IBM 压力功/donor 导热及变物性响应。这是已确认近似，不能直接命名为当前故障根因；同物理终点的更完整动量 dt 细化待验 | [方法修复](../verification/2026-09-06-progressive-method-repairs.md) |
| 标量、组分与容量 | paired remap/Picard 已修有限 dt 库存；非均匀/混排、真实失败回退、`2*N_passive` 归约、有符号被动初态已验；61 个输运标量的 I/O 交集已测 | 目标规模有组分场成本、更多源模型与壁面物理待测；当前 I/O 为 64 个状态字段，3 个主字段后最多 61 个标量 | [标量守恒](../verification/2026-09-06-scalar-remap-repair.md)、[独占验收](../verification/2026-09-07-exclusive-module-acceptance.md) |
| IBM 与热物性 | 空对象/重叠 view 防护；NASA 相邻 double 反算；固体占位物性进入流体热扩散、固体占位状态演化均已修并局部/MPI/pilot 验证 | IBM 标量斜壁约一阶、圆柱未稳定二阶；不把封闭域零通量称为曲面二阶。当前零标量 Re3900 不使用该标量扩散路径 | [独占验收](../verification/2026-09-07-exclusive-module-acceptance.md)、[热恢复](../verification/2026-09-07-re3900-thermodynamic-recovery.md) |
| 应用初场、时间尺度与恢复 | 显式 `initial_state`、`restart_history_policy`、signed 初态、方法恢复后再次精确续算及统计 epoch 链已验 | 自动时间尺度目前只接入对流；其他尺度由调用者提供。新增燃烧/喷雾时须补相应合同 | [独占验收](../verification/2026-09-07-exclusive-module-acceptance.md) |
| MPI、事务与公共接口 | MG 可选本地 counters 不再控制 collective；应用七项冷控制一致性、分配失败与共同回退已有针对性验收 | 不能据局部失败矩阵声称所有产品路径已穷举；ESF/parcel/migration 还没有加入当前共同事务 | [I/O 与接口验收](../verification/2026-09-07-e0fd326-io-contract-audit.md) |
| I/O、Restart 与完成状态 | 日志 flush/close 统一决定全 rank 完成；读取前大小检查、reader bulk 预算、失败不发布、同方法与 rank relayout 恢复已验 | 普通 CSV close 不等于 fsync；reader bulk 不是产品总峰值/RSS。新增模型身份与持久状态仍需扩展 | [I/O 与接口验收](../verification/2026-09-07-e0fd326-io-contract-audit.md) |
| 内存、工作区与生命周期 | 多标量 C++ 唯一分配、恢复/输出交叉存活、销毁重建和 MPI owned 资源已有小型 profile | 全产品硬峰值预算仍缺；MPI/libc、allocator 余量、arena 外数组、halo、image/staging 的同时活跃集合需统一。局部零泄漏计数不等于全进程上限 | [独占验收](../verification/2026-09-07-exclusive-module-acceptance.md) |
| Krylov、MG 与性能观测 | 批量归约、融合 basis update / multidot、r1 guard已存在；V4阈值及MG层级观测已验；单stage写入实验无总收益已撤回 | C2-r1递推恢复来源待细分，guard触发率仍不可观测；64 loops/advance 溢出拒绝完整证据，尚无分段输出 | [阈值成本](../verification/2026-09-08-linear-criterion-observation.md)、[MG实测](../verification/2026-09-08-runner-mg-profile.md)、[撤回实验](../verification/2026-09-08-mg-single-stage-store.md) |
| 构建与可执行身份 | 有 Git 身份的干净 checkout 可独立构建，普通最小 CLI 已运行 | 无 `.git` 源归档可构建但 run 在 `invalid_plan/10505` 被身份合同拒绝；归档来源身份尚无运行合同 | [I/O 与接口验收](../verification/2026-09-07-e0fd326-io-contract-audit.md) |

两项关键量测限定：

- IBM 标量空间 MMS 的两次观测阶数：对齐壁 `1.974/1.994`，斜壁 `0.861/1.134`，圆柱 `1.653/0.861`。后续固定几何复核没有消除该限制；需要一致的几何/通量离散研究，不能只细分 STL。
- `48493ed` 的 7001–7250 单轮、128 ranks、1679 loops 的历史观测完整；纯 advance 平均 9.227 s/步，Krylov 3.780 s/步（M 1.866 s 是其子项）、候选装配 2.186 s/步。它记录了修正版方法的成本，未证明比 COAST 更快。CSV 没有 `||b||`，r1 的 `rtol=1e-6` 只是由冻结规则重建，完整停止阈值尚不可恢复。
- `16fd8f0` 的 9501–9510 新观测确认：r1 实际阈值约 `1.46e-8–1.56e-8`；后 9 步 r2–r6 共 36 loops 都由 atol=`1e-8` 主导。r1 M/apply 与其他 diagonal 相近，下一项选择 MG 内部归因，而不是再加初猜保护或放宽容差。

## 4. 替代验收的案例与比较身份

当前没有公开且冻结的完整替代业务案例集合。已有公共合同回归、Re3900 单例与历史 COAST 模板只能作为建表输入。基础流动总验收前，根代理需按用户业务范围冻结案例、观察量、接受容差和两侧工作量；已有软件短测不能替代该集合。

| 案例层次 | 现有证据 / 当前缺口 | 退出条件 |
|---|---|---|
| 公共基础合同 | 小型流动/标量/物性/IBM/MPI/失败/Restart 覆盖已存在 | 在冻结候选上执行受影响回归，保留数值与故障目标 |
| Re3900 稳定性与科学量 | 7000→7250 已验；本次长测暂停于 9950 完整 health，durable 点 9500 | 按已确认运行协议完成目标窗口；网格、统计收敛和文献量分别验收 |
| 普通 COAST 对照 | 历史 ordinary 模板为 native-air/NASA、可压缩、两次 SIMPLE；现有目录为 104 层，未找到当前 52 层同窗口数据 | 先冻结可比物理/网格/初始与边界状态、工作量、输出及两侧算法差异，再做每配置单轮比较 |
| 反应与喷雾 | 当前主线未接入；外部 mechanism/fuel、科学 oracle 与安全比较案例身份仍需确定 | 完成下文 Stage 接入及其真实 CLI/事务/Restart 验收，再对已确认业务案例给出替代结论 |

历史 ordinary 可执行来源为 `/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/bin/coast-ordinary-compressible-D0p02-native`，本轮只读核实 SHA-256 为 `ab8f5985456960a662334f19aed3902c2ee6785202fd5f72200740cccd60d878`。它是历史比较资产，不是当前 52 层的新结果。旧 constant-cp/480 网格/PISO 收据，以及旧 `19.40/20.16/1.75 s` 计时，均不用于当前速度排名；104 层耗时不除以 2 代替 52 层实验。

本轮进一步确认：普通 COAST fresh 内场速度为零，现有 HUNDUN 圆柱 runner
默认从入口设置全场速度；因此旧 fresh 窗口也不能称为相同初态。
新启动比较可通过 HUNDUN 普通 CLI 已有的显式初态使用静止内场，单独建立新案例；
若使用 runner，先补相同的公开初态接口。它不替换或重置当前 Re3900 长测。
COAST 的 `dtim/tim` 为默认 real，输入相同字面量不保证与 HUNDUN double
相同的实际时间；比较必须列名义/实际 dt、累计时间和步数，严格时间合同仍待冻结。
当前没有已验证的 HUNDUN→COAST checkpoint 转换器；热窗口须从共同初态各自形成原生历史。
52 层网格使用独立生成入口产生新 Decomp/runtime_mesh/IBM 身份，生成后逐点/分区核实；
现有独立生成器尚缺完整 build receipt，生成峰值也待测。旧缓存继续只读。

## 5. Stage 来源与实际完成范围

历史 governance 对象库为 `/home/wyf/code_dev/hundun-flow-governance`。表内 SHA 与 DCO trailer 已只读核实；测试数字属于历史收据，本轮未重跑。Stage 4 工作树存在；Stage 5 `/tmp` 工作树和 portable 工作树已删除，但分支/对象仍可读。路径消失不等于实现不存在。

| 来源 | 可复用身份 | 当前判定 |
|---|---|---|
| Stage 4 | tested code `6407cd7c591ce088db7f1dd7e296d77acd18da1c`；seal `033a685c90c1f9c674e93a4b82db10db4c381abe`；`coast/stage4-reacting-flow` | DCO；历史 `STAGE4_ACCEPT` 限定于其代码/合同；当前 Cartesian 接入待做 |
| Stage 5 框架 | tested code `41b2aac97d28da3949a6a4bd629c079d7aa6a8b7`；seal `02b57cce45311bfd6c1507b6e6b32eb0617a59cd`；`coast/stage5-esf-tpdf-tcr` | DCO；框架验收有效，实际 field/performance 门后来重新打开 |
| Stage 5 场执行后续 | `coast/stage5-field-validation`，`8ffdf2b5673374fb14639fc4dce09a1b586ee5db` | DCO tip；`STAGE5_FIELD_VALIDATION_ACTIVE`；没有 field/performance seal |
| Portable Stage 5 | `25287ececc4c15a19d2bf248ffd41496ce60ab09` | DCO；纯闭包/候选，不发布产品状态 |
| Portable Stage 6 | `a3cdd925d1575fd0e58e43ece3a0ba02fbfd501c` | DCO；纯 parcel/property/exchange/correlation/TAB 内核 |
| Portable handoff | 当前仓库 `codex/stage5-stage6-portable` tip `a20b85b399e96e537ddc3d241da3a615d1153216`；`docs/handoff/2026-08-30-stage5-stage6-portable.md` | DCO；`hundun-flow-four-branches-precleanup-20260902.bundle` 另有该 tip 备份 |

按阶段读取历史资料时，使用对应对象库的 `git show <SHA>:<path>`；不要把 governance 的 branch 名直接交给当前 Cartesian 对象库。Stage 4 收据位于其现存工作树 `.superpowers/sdd/stage4-final-acceptance-report.md`；Stage 5 的关键资料是 field 分支中的 `.superpowers/sdd/stage5-field-validation-reopen.md`、`docs/numerics/stage5-capability-ledger.md` 和 `docs/validation/stage5-1d-premixed-final-receipt.md`。

| Stage | 已有可复用行为 | 真正未接入或未完成项 |
|---|---|---|
| 4 | 完整 composition/mechanism 身份、全 species `rhoY`、总热化学焓；Cantera 3.2.0 backend/workspace、thermo/transport/interval、反应源事务、两次 PISO、schema/Checkpoint/diagnostics 合同 | 当前主线没有 backend target/source admission，仍拒绝 `reacting=true`。旧 schema-v4 `run_reacting_flow_case` 只验证计划和 MPI 一致性后返回成功；真实 mean reacting CLI 仍须完成 |
| 5 框架与场执行 | N=2/4、Philox/Wiener、输运、IEM、逐场 chemistry、ensemble、元素一致性、PSR shadow、TCR algebra/history；后续分支有真实 ESF/off executor、cellwise 状态和重启测试 | `execute_tpdf_tcr_field_case` 明确拒绝 non-off TCR；field candidate mapper 未完成。schema-v5 field 路径依赖 `HUNDUN_PRODUCT_HAS_CANTERA`，无 backend 时旧分支仍成功退出。当前产品应对未实现路径显式失败 |
| 5 portable | `finite_rate_mean_shadow`、`pasr_algebraic_v1`、chemistry-only ESF ensemble、LES 混合时间、命名化学时间尺度、candidate-only 报告 | 没有 direct Cantera adapter、完整 ESF 输运/TCR 或源项发布。两次纯积分调用不等于流动中 C–T–C 顺序通过 |
| 6 portable | 128-bit ID、无状态 RNG、液体物性状态、成对守恒 delta、Schiller–Naumann、Ranz–Marshall、Spalding 单组分蒸发、TAB 振子/首次事件 | 无容器/轨迹/注入/MPI migration/IBM rebound/TAB 子滴/共同事务/持久化。该蒸发核不是计划中的 Abramzon–Sirignano；TAB 仅请求 breakup |

Stage 5 一维收据的最终科学状态为 `SOFTWARE_PASS_ATTRIBUTION_COMPLETE_MULTI_SOURCE_ORACLE_PENDING`：多来源统计 oracle 和合格 COAST 三维比较尚欠。它不是已发现的反应速率实现缺陷，也不是 Stage 5 总接受。旧 field 计划中的九对重复/median/置信区间要求按第 1 节由当前单轮规则覆盖；剩余科学身份和业务验收问题保留，统计性能门不默认启动。

## 6. 当前公开接口的接入合同

| 接口 / 权威 | 保留并扩展的责任 |
|---|---|
| `ThermodynamicsPlan` / `ThermophysicalCompiler` | 对齐完整 mechanism 与 Ns−1 独立 species 的映射、element/phase/文件身份及焓基准。native-air 的 273.15 K 参考平移不直接当 reacting 总热化学焓；heat release 不重复加入能量 |
| `ThermophysicalPredictorPlan` / rate history | 明确区间积分 chemistry delta 与非对流 RHS 的单位、时间层和采样点，保持已接受/历史状态权威 |
| `ContributionRegistry` / graph | 当前仅准入 `inert_source`；graph 12/45 是空贡献阶段。先冻结实际化学—输运—PISO 顺序、读写集、单位和资源，再扩准入 |
| `AttemptTransaction` | 使用既有 collective prepare / accept / reject；把 chemistry、ESF、TCR、parcel、injector、owners 和迁移暂存纳入共同回退 |
| field/revision/halo / `FinalFaceFluxAuthority` | ESF 与标量消费同一最终 flux、几何/边界权威；增加 replica/source/parcel 采样声明，禁止无证书的陈旧缓存 |
| `TurbulencePlan` | 消费已有 WALE/Vreman 与 `Sc_t`。PaSR 使用运动学 `nu_t`；ESF 随机项不重复计入分子扩散 |
| Restart / diagnostics | 扩展 mechanism、ESF、TCR、parcel、injector、RNG 及算法版本身份；旧 reader 和 validate-then-publish 行为保持，新增状态须真实续算验收 |
| chemistry/parcel 热路径 | portable 值类型可复用；批量、借用视图与冷分配 workspace 按当前资源合同适配。SplitMix 与 Philox 的算法身份分开保存，不静默换序列 |

## 7. 依赖顺序与阶段退出条件

以下是接入顺序，具体实现文件与 focused selectors 在进入每一包时冻结。旧 driver/schema/checkpoint 只作行为和测试参考，在当前接口上重表达；历史 seal 不继承给新二进制。

| 包 | 依赖与工作 | 退出条件 |
|---|---|---|
| B0 基础流动 | 完成当前观测切片，单一有依据优化，关闭业务案例中的稳定性/模块/I/O 项 | 冻结源码/二进制/案例/证据；明确通过范围及剩余限制；基础流动达到用户当前业务验收要求 |
| C0 身份与纯接口 | B0；引入 portable 值类型/纯核，定义 mechanism、焓、源项区间和时序；Cantera/rate-query 批量 adapter | 纯核原有测试、真实 backend thermo/transport/interval、身份错配、异常/守恒、无 backend 显式失败；公开头不暴露外部 backend 类型 |
| C1 Stage 4 产品 | C0；chemistry 准入、共同事务、真实 mean CLI | 0D 独立 Cantera 对照；两次 half chemistry/PISO 与共享 flux；各阶段故障回退；fresh 字段/步数/时间；连续/重启及 1/2/4-rank 全局 cell 比较 |
| C2 mean closures | C1；finite-rate/PaSR 按产品时序接入 | `kappa=0/1` 极限、species/热报告一致缩放、timescale unavailable 失败、shadow 不发布源项；流动中时序验收 |
| C3 Stage 5 产品 | C1；ESF N=2/4 → IEM/逐场 chemistry → TCR field mapper | RNG retry 不变、2N 调用、13 阶段、共享 WALE/flux、逐 cell root history、off/shadow/validated、共同回退、Checkpoint 与实际 CLI |
| C4 Stage 5 业务验证 | C2/C3；0D/MMS、小场/IBM/失败，再进入已冻结业务案例和可比性能窗口 | 科学 oracle 与案例身份明确；所有要求的 field/continuation/性能观察量有证据；采用当前单轮协议，不以纯核速度替代场执行 |
| S0 Stage 6 力学（暂挂） | 用户恢复范围 + C4；SoA/ID、轨迹、stencil、注入、物性、migration、IBM rebound | ballistic/Stokes/Galilean、插值、ownership/rollback、Restart、small 1/2/4-rank 合同 |
| S1 Stage 6 交换（暂挂） | S0 + Stage 4；film sampling、A–S、事件终止、gas/parcel 事务 | 单滴加热、d² oracle、A–S 亚步收敛；单/多 parcel 质量/动量/总热化学焓闭合；第二次 PISO 源项账本 |
| S2 Stage 6 组合（暂挂） | S1 + Stage 5；TAB 子滴、ESF common-source、driver/schema/Restart/diagnostics | 子滴守恒；每 parcel 只算一次交换；N=2/4 共用气源；injector/TAB/RNG/迁移持久化；两个 surrogate 的有界 smoke、失败矩阵和实际 CLI |
| A0 当前范围替代接受 | B0、C1–C4；S0–S2暂挂，不计入本阶段退出条件 | 对冻结流动/燃烧业务案例逐项给出接受/拒绝和准确证据身份；新代码引入的问题已关闭，业务限制已说明；不得只用 development seal 宣称全部替代，也不据此声称两相已完成 |

## 8. 待确认的外部输入与维护规则

普通 Stage 进入、接口实现、当前 direct Cantera / ESF N=2/4 / dilute single-component / A–S / none–TAB 选择已有任务上下文，不重复申请。`fgm_table` 只有 enum 预留，未进入现有已完成能力清单。

真正缺少的输入按依赖出现时处理：最终 mechanism/phase/composition 与液体物性包身份；私有 COAST 可安全比较的 manifest/executable/case/restart 身份；一维多来源 oracle 的机器可读曲线或预注册重建方案；如需公开分发，核实机制/物性资料的分发权。公开 fixture 足以支撑的接口与软件验证可继续，不为这些后续输入暂停无依赖工作。

更新本台账时：

1. 先更新第 2 节当前切片，只保留一个活动优化；同时更新实际运行/暂停状态。
2. 状态从“待测”变“通过”时附确切候选、命令、原始日志和作用范围；失败保留并区分生产缺陷、夹具问题、外部身份和方法限制。
3. 新候选只重新验证其影响的合同；已通过的独立工作不重新实施。每配置单轮，额外重复或统计门先确认。
4. 新字段、源项、RNG 或模型改变持久语义时同步检查 Restart、history identity、诊断和共同回退；不能只更新 kernel 测试。
5. 完成各包后回填退出证据；总目标仅在 A0 满足时完成，剩余阶段不可由历史 seal 或“能运行一步”省略。

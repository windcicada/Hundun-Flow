# 递进的方法检查与修复

工作基线：b779bff4bbff7067e691952b25c90cbc86be7795 加此前未提交修复，实际工作区为 `/home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity`。保留现有修改。上一轮方法清单见 [审查报告](2026-09-06-method-first-numerical-audit.md)。

本文件按实施顺序保留中间失败和冻结记录；下文“未提交/未推送”描述的是各阶段当时的状态，不代替当前 Git 状态。本次按用户授权整理源码、测试及证据提交 GitHub main，不改写原冻结程序、源码归档、checkpoint 或历史 receipt。提交前已用 `tar --compare` 核对当前 CMake、cmake、v0.4、tools 和 AGENTS.md 与 `method-frozen-periodic-20260906/source-snapshot.tar.gz` 完全一致；本次仅补充本报告的状态说明，不重建或中断正在运行的冻结程序。

## 顺序与提交边界

| 项目 | 当前状态 | 完成条件 |
|---|---|---|
| M1：SIMPLE 新动量面通量 | 本轮集成回归通过 | 权限、时间阶数、checkerboard 与产品 1/2/4 rank 未出现回退 |
| M4/M6：最终方程与累计收支 | 已实现并通过独立库存核对 | 最终动量重装配；区分内能/动能、方程误差和物理边界能量收支；仍是观测，未擅自新增动量容差 |
| M2/M3/M5：压力—焓线性化 | M5 已修复；M2/M3 近似范围已验证 | 冻结块方向导数通过；完整密度对流/IBM donor 响应仍不在 Schur 中，不能称为完整 Newton |
| M7/M8：标量分裂、IBM 历史闭合 | M8 已修复；M7 已定量记录限制 | 被动标量解与库存漂移均为二阶，但有限 dt 非严格守恒；未实现耦合组分的严格守恒校正 |
| M9：SIMPLE C1 联合目标 | 已修复，时间细化通过 | C1 也求解与实际候选验算一致的 p/h 联合方向；不改变 SIMPLE 的第二次动量 sweep |
| M10：开放/周期混合边界 | 周期端面通量误清零已修复；1/2/4 rank 回归通过 | 新 C1/C2 线性目标与候选零修正目标一致，保留解析非零周期通量 |
| 全体集成与 Re3900 | 最终 156/156 聚焦集成通过；修复后干净初场长测已启动 | 不能以局部通过推断长期稳定性；原恢复步失败保留，当前长测仍待验收 |

不改变 Re3900 的网格、物性、时间步、12 次 refinement 容量或收敛门槛。测试每种配置执行一轮，不取三轮中位数，不设置耗时门禁。未知问题不可能穷尽；报告将明确已检查范围、实际验证和剩余风险，不将“能跑过一个时间步”当作完整正确性证据。

## M4 第一项：最终动量方程

确认的缺口：程序在最终 p/h 校正之后，检查了连续性和能量，却没有用同一最终状态与面通量重装配动量残差。只看最后一次动量线性求解的残差不能补上这个检查。

实现约束：复用现有动量工作区；放在 pending final flux 发布之后、时间步事务提交之前，以免覆盖仍参与候选状态校验的系数。报告有量纲的动量 L1/L∞、质量/能量方程误差及内能/动能存量；不把组合方程误差冒充物理边界能量收支，不擅自套用新的动量容差。

### M4/M6 已完成的观测与验证

`DriverStepReport` 增加最终方程和边界收支两个报告。时间步提交前，用最终状态、物性、梯度、IBM 约束和最终面通量重新装配动量。终端方程报告提供动量 L1/L∞（N）、质量误差（kg/s）、能量误差（W）、质量（kg）、内能与动能（J）。组合方程误差 `RE + U·Rm − K RC` 不是物理总能量守恒账目。

独立账目计算 BDF 总能量变化与净边界焓输运、动能输运、导热和黏性功。压力功包含在边界焓通量中，不再重复相加。当前产品无体源，IBM 为固定、绝热壁面；不能将该报告无修改地用于运动壁面、化学反应或喷雾源。IBM 笛卡尔重构带来的总能量离散缺陷归入收支误差，不假装成真实壁面热源。

累计账目按实际 BE/变步长 BDF2 系数推进独立库存，仅在时间步接受后提交。不能用 `sum(dt*BDF残差)` 冒充 BDF2 的累计物理漂移。重启后以加载的库存/历史建立新 epoch；`epoch_start_step` 明确表示它不是原始零时刻开始的全程累计值。运行器新增 `conservation.csv`，与 checkpoint 一起刷新，记录单位和最终通量身份，不改变原残差门槛。

实际验证：

- 公开 ProductDriver 的 SIMPLE 均匀场报告先 RED（报告未填写），实现后 GREEN。解析单位体积内能为 253312.5 J。
- 1/2/4 rank 产品回归通过：`method-m4-boundary-balance-final-ctest.log`。连续推进、开边界、质量流量、暖启动均保留。
- 非均匀周期波：从公开重启快照和独立理想气体关系积分库存，核对每步质量/总能量 BDF 变化及累计误差。原 PISO 时间细化 1/2/4 三档通过（是三个不同时间离散，不是三轮性能取中位数）。记录 `method-m4-independent-balance-ctest.log`。
- 实际 runner CLI 的 1/2/4 rank checkpoint 单 rank 分配故障检查及新增 CSV 解析均通过，记录 `method-m4-m5-green-ctest.log` 的前四项。该批次整体不是全通过，见下节。

实现过程发现并修正了新增观测的工作区复用错误：候选焓的物理边界角点仍会被下一步热力学闭合读取，不能将整个 padded 字段临时覆盖为 K。改用已结束本步生命周期的 HbyA 动量工作区，不增加全场数组。初次错误日志保留为 `method-m4-boundary-balance-ctest.log`；这是本轮新增代码的问题，不是旧圆柱故障的证据。

## M5：廉价 Schur 的局部热力学时间响应

原 diagonal Eh 直接采用普通焓对角 `a0 V rho + diffusion_proxy`，确实没有 `a0 V h rho_h`。新增公共数值入口 `PressureEnergyEnthalpyOperator::bind_diagonal` 保留该精确时间响应，同时继续省略空间耦合、保留原扩散代理。使用原有 compiled-local 工作区，无新增全场分配或 halo；不改变 COAST 原生物性。

公开算子测试的独立解析场为定压理想气体、h=cpT、rho=1：rho*h 对 h 的导数为零，所以方向 dh=3、扩散代理=0.7 的响应必须是 2.1，而不是旧路径的 8.1。旧 diagonal 委托的 RED 记录在 `method-m5-temporal-diagonal-red.log`，修正后通过；同时检查焓参考平移（期望 7.7）、工作区别名拒绝及 IBM solid 行不读取无效 EOS。

这不是完整 Newton，也不是已证实的圆柱加速。小 IBM 产品在本修改后用了 7 次 refinement，成功达到原残差要求，却触发旧测试“所有 refinement 必须 diagonal”的断言。源码原本在后半预算切换 spatial；该断言已改为核对实际已存在的 diagonal→spatial 契约。不得为使测试通过而修改终端阈值或把整个失败批次记为通过。

## M9：新增发现——SIMPLE C1 求解与候选目标冲突

将原公开 ProductDriver 周期压力—焓波时间阶数检查扩展到 SIMPLE 后，三个 dt 均在第二步出现 `6/5792`（stage 44），随后重试减半 dt，并从 BDF2 降为 BE。见 `method-simple-temporal.log`。这不是 Re3900 第 1297 步重放，而是对 SIMPLE 时间离散契约的独立检查。

源码依据：C1 特殊分支只解连续性压力方程，强制 dh=0，并跳过 C1 的压力—能量块准备；后续却以联合连续性/能量残差的 Armijo 下降筛选该方向。因此压力方程解得再精确，也不能保证它是联合目标的下降方向。

已验证的修复：移除压力单独求解的 C1 快捷分支，使 C1 也冻结目标焓、准备块并求解同一联合目标的 p/h 方向。继续复用现有 Schur/FGMRES/MG 与真实候选验算；SIMPLE 仍有第二次新动量 sweep，与 PISO 不同。没有将 C1 能量残差伪装成零，也没有放宽联合筛选或最终能量门槛。

`2026-09-06-1297-fix-evidence/method-m9-c1-complete-target-ctest.log` 的五项通过；后续 M8 集成日志再次验证。SIMPLE 在 dt=1.25e-4、6.25e-5、3.125e-5 三档均只执行一个 BE 恢复步，随后分别 7/15/31 个 BDF2 步，无重试或降阶。p/h/rho/T 时间自收敛阶为 1.965187，U 为 1.991381，最终面通量为 1.984975。PISO 同一检查仍通过。它们是小型光滑周期波结果，不是圆柱速度或长期稳定性结论。

## M8：接受步的 IBM 热扩散速率历史

源码确认：最终目标能量方程使用未截断的零法向热流 donor 重构，而 `evaluate_thermophysical_rates` 先前使用有界温度 donor 重构。后者存入下一步 EX2 历史，两者不是同一空间算子。

新增公开速率入口回归使用具有负二次 donor 权重的温度场，材料系数冻结，独立计算目标扩散。旧路径局部速率差为 493.128，修复后为 0。现在接受步的速率历史与目标方程统一使用 `correct_zero_normal_diffusion`。预测器的正性检查、低阶端点和保守混合保持不变。

RED/GREEN 见 `2026-09-06-1297-fix-evidence/method-m8-accepted-rate-red.log`、`method-m8-target-rate-green-ctest.log`；后者八项通过，包括 IBM 接口、热物性 IDP、产品 1/2/4 rank 以及 PISO/SIMPLE 时间细化。

## M2/M3：线性化范围的验证，而非完整 Newton 声明

`method-m2-m3-scope-ctest.log`：公开四块算子的 p-only、h-only、混合方向在三档差分扰动下，与冻结密度的目标一致；最大能量方向误差约 1.60e-8。另将面密度对流恢复为 EOS 随动时，连续性与能量的一阶响应分别多出 1.61e-5 和 4.8325，证实这部分没有进入当前块系统。该检查显式要求 `full_nonlinear_jacobian=false`，不能用冻结目标的导数检查冒充完整候选导数检查。

IBM 原始 donor 导热算子的中心差分方向误差为 2.50e-9。它证明目标 donor 响应的符号/离散一致，不证明当前 masked Schur 已包含它。IBM 模式下被省略的压力功导数、donor 导热响应以及物性随动响应仍属于准 Newton 近似；候选筛选与最终连续性/能量门槛使用实际目标残差。完整补齐需要处理非局部 Ch 消元与 donor 通信，不能只往 Cp 随意加项后仍逐单元除 Ch。

## 方法变更时的 checkpoint 恢复

新增 runner 参数 `--restart-method-recovery`。严格验证原 checkpoint 的布局、清单、哈希和统计绑定后，只丢弃旧时间/派生速率历史，不改变检查点物理流场，也不改写源 checkpoint。通过既有恢复路径重建速率，首步使用 BE，此后再使用 BDF2；RUN.meta 记录显式选项，health/Evidence 保留恢复标志，累计收支建立新的 checkpoint epoch。

此选项不能与要求精确历史的 MG 布局迁移混用。它不等于程序已经能自动辨别所有跨版本算法历史。当前恢复由提交脚本明确选择。

CLI 验收包含无该选项的旧程序 RED，增加选项后的 1/2/4 rank GREEN；核对 step=2、BDF=1、attempts=1、恢复标志与 epoch=1。见 `2026-09-06-1297-fix-evidence/method-restart-recovery-red.log` 和 `method-restart-recovery-green.log`。源 checkpoint 保持不变。

## M7：标量分裂的实际后果和未解决边界

使用同一公开 ProductDriver 周期波，加入常值 0.2 和 `0.2+0.05 sin(x)` 两个被动标量。独立从最终密度和标量积分 `sum(rho*q*V)`，周期边界无净标量外流。当前方法对常值的最大误差为 4.44e-16；非均匀解的时间自收敛阶为 2.0101。

dt=1.25e-4、6.25e-5、3.125e-5 时，1e-3 s 后的非均匀标量质量漂移分别为 -1.9267e-7、-4.6335e-8、-1.1286e-8 kg，误差阶 2.0560、2.0376。见 `2026-09-06-1297-fix-evidence/method-m7-scalar-valid-input.log`。`v04_product_scalar_splitting_observation` 仅证明此光滑场的常量保持/时间阶数，标签为 method_observation，**不是严格守恒验收**。

确认的限制仍在：预测器保守更新的 rho-star*q-star 在后续 p/h 改变密度时没有配套通量校正，有限 dt 存在库存误差；对会影响 EOS 的组分尚不能据此宣称同样的长期误差控制。本轮没有用全局归一化或简单乘 rho-star/rho-final 掩盖缺陷，后者会破坏常量保持。严格修复应同时校正密度、标量面通量及组分—EOS 耦合，需要独立任务包。当前 Re3900 的 transported_scalars 为空，所以此限制不阻止其诊断性长测；并不等于软件全部已知问题都已修完。

## 最终集成与冻结候选

重新编译 Release 全部目标后，单轮选定 153 项方法/产品/重启/输出/MPI 回归，153/153 通过，总墙钟 139.81 s；`git diff --check` 通过。日志为 `2026-09-06-1297-fix-evidence/method-final-integrated-ctest.log`。这不是仓库全部 289 项测试、sanitizer 验证或 Re3900 长测验收；其中 scalar_splitting 是明确标记的近似范围观测。

候选位于 `/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/method-frozen-20260906`。包含完整 v0.4 源码归档、构建清单/编译命令、原样复制的 case/STL/网格/物性/统计计划，以及独占输出目录的直接 mpirun 脚本。未提交或推送 GitHub。

- executable SHA-256：`1c6046a488f493814b6f0100e40de252cacbe303c1ec41faab7c1cecf7add640`
- source archive SHA-256：`1fa7a7c55e75ac76ddd5f2b7f6aef3f95802c292d7469b57b855e70914bbf72a`
- target build manifest SHA-256：`73ffbcc80dd822343b87d266adb4d8b665afd3ae9f40b40e9101d492ff5c731b`
- 原 checkpoint manifest SHA-256：`6fa7b17f2e00145bebded67822aeeee793327a1812ff6f18d3bf2dd4b27f657e`

### 首次长测提交的真实结果

2026-09-06 10:18:15 CST 通过 systemd 用户服务 `hundun-re3900-method-20260906` 启动直接 mpirun 128 ranks，原 case、dt、容限、12 次 refinement 不变，显式方法恢复，目标 35000 步。程序 28.36 s 后受控返回 7：首次 BE 恢复步 1001 用完 12 次 refinement，连续性 4.80875e-6、能量 2.59397e-6，未达到 1e-6。**没有接受新步，不存在新的长测通过或运行中状态。**

日志 `/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/long-method-1000to35000-20260906.log`；同名目录中的 RUN.meta 验证 checkpoint=1000、方法恢复=1、上述 executable/build SHA 一致。原 checkpoint 未改写。

后续先测试“线性 RHS 是否等于被候选验算的 alpha-zero 残差”及实际线性方向预测，尤其核查开放边界目标冻结的时序。不直接认定省略项就是唯一根因，不扩大迭代预算或提交相同配置反复碰运气。新增手动 `v04_product_method_direction_probe` 使用公开 ProductDriver/RestartReader 和既有测试观测接口，源 checkpoint 只读，不生成生产结果或 checkpoint，其退出成功仅表示得到诊断，不表示流场计算通过。

## M10：开放进出口与周期端面混合时错误清零通量

这是后续定位并修复的确定性方法缺陷，不是调参或放宽残差。

`PressureCorrectionFaceRule::physical` 表示几何上的全局计算域端面，**包括周期端面**。`stage_frozen_momentum_flux()` 在开放域候选路径中仅按该布尔值判断外边界所有权，把周期端面与入口/壁面一起设为零占位。独立的物理边界 finalizer 正确跳过周期面，却因此保留了这个错误的零通量。其结果是：新 C1/新 SIMPLE C2 的线性系统包含周期输运，而候选验算将其截断；后续 refinement 从已清零的候选开始，所以两套残差看似又一致了。

独立诊断排除了速度梯度/有效黏度未刷新的猜测。校正诊断工作区复用后，128 ranks、旧 checkpoint 1000 的 C1/C2 目标差如下：

| 状态 | C1 连续性目标差 | C2 连续性目标差 | C1/C2 能量目标差 |
|---|---:|---:|---:|
| 修复前 | 1.9917886e-3 | 1.0723078e-3 | 5.5460845e-4 / 3.0063914e-4 |
| 修复后 | 1.1542255e-17 | 1.1543120e-17 | 0 / 0 |

修复前最差全局单元分别为 (135,146,51)、(93,152,0)，均紧邻 Z 周期端面，而不是 X 入口/出口。日志见 `method-linear-target-location-probe.log` 与 `method-open-periodic-checkpoint-green.log`。最初 `method-linear-target-gap-probe.log` 的连续性差值错误复用了被候选能量覆盖的工作区，已明确废弃；不得引用该日志的巨大 C 差值作为程序缺陷证据。

修复在 `PressureCorrectionFaceRule` 中明确非周期外边界判定；候选构造保留周期面通量，并将周期面纳入内部/交换面的数值来源校验。没有通过在 RHS 末端补偿差值来掩盖错误，也没有增加永久场、改物性/时间步/容差。检查了该标志的其他使用点：压力系数与初始化投影均已有显式周期分支，保持不变。

新公共 ProductDriver 回归使用 8×8×8 均匀斜向流 U=(3,0,0.25)，X 开放、Y 对称、Z 周期，无 IBM、无 checkpoint。修复前 PISO/SIMPLE 均在第一步返回 6/5792，C1 目标差 C=1.9723866e-3、E=1.5555556e-3；修复后 1/2/4 rank 均通过。验收同时检查两套零修正目标、一次接受、解析非零 Z 质量流率以及速度常量保持，不只看退出码。见 `method-open-periodic-red.log`、`method-open-periodic-green.log`。

真实 128-rank checkpoint 的只读诊断在修复后接受恢复步 1001，使用 6 次 refinement，未改变原 12 次上限。它仅证明该步以及本次缺陷修复，**不是长时间稳定性或旧流场正确性的证明**。

旧 checkpoint 曾在被错误截断的周期输运下演化，不能自动当作正确物理初场继续累计正式 Re3900 统计。旧输出、checkpoint、冻结程序和失败日志均保留。已向用户说明：无进一步指定时，从干净初场提交新长测；若选择旧 checkpoint，只标记为稳定性诊断，不复用旧统计。

### M10 集成结果与新冻结版本

完整 Release 构建成功。前述 153 项加上 3 项开放/周期混合边界回归，单轮 **156/156 通过**，墙钟 135.56 s，见 `method-m10-integrated-ctest.log`。这是选定的相关回归，不是全部 292 项或 sanitizer 验收。`git diff --check` 通过。

新冻结目录 `method-frozen-periodic-20260906`；程序 SHA-256 为 `9c29bb403bfc3a4e85b29a1604c2990858c76de939af7f9d5b1ac14e12214d8b`，源码归档 SHA-256 为 `9c403356db6fc298b6650dce5b42ab1ef9a7f67ea2a7b791fd3d2c18a779ea87`，目标构建清单 SHA-256 为 `b4fb675bb3335b537a6b082f54709e6488fee3e198a41bf54303e012bb8baec2`。仍为含已保留修改的工作区冻结，不是新的 Git commit；未推送 GitHub。

准备提交的干净初场长测目录是 `long-periodic-fresh-35000-20260906`：同一 Re3900/456×256×52/128 ranks/固定 dt，前 10000 步发展、收集至 35000 步，checkpoint 与 Visit 周期均为 500。未恢复旧统计。启动及实时接受状态以该目录 RUN.meta、health.csv 和外部同名日志为准；目录计划本身不证明进程已启动或接受了时间步。

本轮不能宣称穷尽“全部可能问题”：M2/M3 的准 Newton 近似边界、M7 的有限步长标量守恒缺口仍按上文记录。新长测验收重点是周期通量、最终动量残差、质量/总能量收支、重试/修正次数与 RSS 的趋势；不能用一个短时通过代替这些长期检查。

### 新长测已实际启动

2026-09-06 11:21:31 CST，`hundun-re3900-periodic-fresh-20260906.service` 启动成功。RUN.meta 的 restarted=0、requested_steps=35000、程序/清单哈希与冻结版本一致。启动检查已接受至 51 步；完整 health 前 43 步均为一次接受，refinement 最大 4，最大连续性残差 2.8892e-7、最大能量残差 9.8162e-7。36–43 步推进最慢 rank 的平均推进为 8.372146 s；只是启动短窗，不作与旧错误程序或不同网格 COAST 的加速比。

运行目录 `long-periodic-fresh-35000-20260906/submission.md` 保留提交信息。测试列表 `test-results.md` 的独立标记区域由轻量观测服务每 60 秒更新，仅采完整 MPI 步，明确保留普通可压缩 COAST 和 PISO 的未测行，进程结束后保存最终状态。观测器不修改流场、checkpoint 或收敛参数。

# HUNDUN-FLOW Stage 4--6 Linux CPU v1 架构设计

日期：2026-08-09

状态：架构方向、阶段边界和任务图已经用户逐节批准。本文件是书面复核稿；
用户审阅并批准本文件及配套实施计划前，不授权修改 Stage 4--6 产品代码。

规划基线：governance commit
`ee4d2b18d0c68b3080edcc2e132045175961cfb8`。该哈希仅表示编写本规格时
使用的治理快照，不预先声明正在开发的 Stage 3 已完成。Stage 4 实施必须以
未来正式接受的 Stage 3 product/code head 为 parent，并在任务 `4F-0` 重新登记。

本设计保留 Stage 0--3 的已接受合同和历史记录，取代
`docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md` 中尚未实施的
Stage 4--6 粗粒度路线。旧文件不删除、不回写失败证据；新计划通过 superseding
notice 说明差异。

## 1. 产品目标与非目标

Stage 4--6 交付 HUNDUN-FLOW 第一个完整的 Linux CPU-reference 版本：

- Stage 4：多组分低马赫反应流、Cantera C++ 后端和无 Python 安装包；
- Stage 5：低雷诺一致 ESF/TPDF、IEM 和用户提出的 TCR；
- Stage 6：稀相 Lagrangian point-parcel 喷雾、单组分液滴蒸发、双向守恒、
  static IBM impact 和 TAB breakup；
- 最终版本：`1.0.0`，能力声明严格受实际证据限制。

唯一生产平台为：

```text
Linux x86_64
Ubuntu 22.04 / glibc 2.35+
GCC 11 / libstdc++
C++17
_GLIBCXX_USE_CXX11_ABI=1
x86-64 generic ISA（禁止 -march=native 进入发行包）
MPI-3 CPU reference
```

Clang 15/libc++ 只用于独立源码可移植性、headers 和 sanitizer 测试，禁止与
GCC/libstdc++ 构建的 Cantera 二进制混链。

Stage 4--6 明确不实现：

- AMR、非结构网格或 rank-changing Restart；
- 移动/多部件 IBM、FSI、cut-cell 守恒再分配；
- 生产 GPU、GPU-aware MPI、混合精度；
- 全可压缩、声学、激波、WENO 或 DG；
- MISDC/SDC、动态增加 PISO corrector 或经验 pressure damping；
- NativeChemistryBackend；
- ESF field-index MPI decomposition、differential-diffusion ESF、EMST/MMC；
- 多组分液滴、碰撞/聚并、液膜、一次雾化、KH--RT、LISA、稠密喷雾；
- 96^3 及更大正式网格。

NativeChemistryBackend 是 post-v1 独立路线，首版不同时重造 Cantera。

## 2. 共同架构原则

### 2.1 权威和依赖方向

HUNDUN-FLOW 始终拥有：

- 网格、全局编号、Halo、MPI 和空间分解；
- 质量、动量、组分和总热化学焓守恒输运；
- PISO、机械压力、热力学压力 `p0` 和最终面质量通量；
- IBM、WALE、边界和表面语义；
- operator schedule、时间步、retry、rollback 和 collective failure；
- Checkpoint、Restart、diagnostics、driver 和 capability ledger；
- Stage 5 TPDF/TCR 与 Stage 6 spray 的组合顺序。

外部后端不得返回可直接提交的网格场、面通量或 pressure correction。

### 2.2 事务

一个 macrostep 只有 `history`、`committed` 和 `trial` 三类物理状态。任何
chemistry、PISO、ESF、TCR、parcel、migration、breakup 或 collective 失败均使
所有 rank 拒绝整步。拒绝后：

- 物理时间、BDF2 history、`p0`、RNG accepted clock 不推进；
- Chemistry workspace 重新初始化；
- 随机场值、TCR root history、parcel owner 和 injector 余量恢复；
- 不复用失败 trial 的 chemistry half-step 或 parcel predictor；
- 不增加 corrector、不放宽阈值、不对非物理解 clipping 后继续提交。

### 2.3 状态、服务与产品 API

Stage 4--6 不扩展当前 metadata-only plugin ABI v1。公开头只暴露稳定 value、
configuration 和 report 类型。下列服务保持 in-tree internal interface：

```text
ThermodynamicsService
TransportPropertyService
ChemistryBackend
ReactingSourceTransaction
TpdfTcrCombustionModel
LiquidPropertyService
SprayCouplingTransaction
```

公共 C++ 头不得出现 Cantera、SUNDIALS、COAST、OpenFOAM、AMReX 或 Pele 类型。

### 2.4 构建与运行

正常 configure、build、install、runtime 和正式验收均不得要求 Python、Conda
或在线 fetch。所谓 network-independent build 只表示构建不会发起网络访问，
不得关闭宿主机网络、修改 firewall/routing/NetworkManager 或影响 Codex 服务。

### 2.5 Post-v1 原生化学后端兼容性

NativeChemistryBackend 不进入 v1 实现或验收，但 v1 internal services 不得把未来
替换路线锁死。`CompositionIdentity`、`ThermochemicalPoint`、
`ThermodynamicsService`、`TransportPropertyService`、`ChemistryBackend` 和
`ChemistryIntervalReport` 必须只表达 HUNDUN 的物理量、单位、积分区间、积分增量、
失败状态和 composition fingerprint，不暴露 Cantera mechanism、Solution、Reactor、
Integrator 或 error 类型。

post-v1 原生路线至少预留：

- NASA-7/9 组分热力学；
- Arrhenius、可逆反应、third-body 和 falloff；
- 独立 mechanism parser 或 AOT mechanism compiler；
- 解析/自动生成 Jacobian 及其 sparsity identity；
- 自适应刚性积分、dense/sparse linear solve；
- chemistry batch scheduling 和 MPI/thread 负载均衡。

这些能力依据公开报告、标准方程和论文独立实现。Cantera 作为黑盒交叉验证和兼容后端；
不得依据 Cantera 源码机械翻译。Native backend 正式接受前，Cantera backend 继续是 v1
生产后端和科学对照。

## 3. Stage 4：反应流和 Cantera

### 3.1 Cantera 分工

Cantera v1 后端只负责：

- 理想气体混合物热力学；
- 组分热容、焓、熵和状态反演；
- 反应速率、生成率和单元局部刚性化学积分；
- mixture-averaged 生产级气相输运性质；
- 0D/PSR；
- mechanism 解析和明确失败诊断。

首选版本固定为 Cantera `v3.2.0`：

```text
tag: v3.2.0
commit: 4a8358eb80cfeb50474386b5f9ec0b3a83519889
source archive SHA-256:
a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b
```

全部传递依赖版本、archive hash、二进制 hash、ABI、许可证和 patch 在
`4P-1` 冻结；设计候选包括 SUNDIALS 5.3.0、yaml-cpp 0.7.0、fmt 9.1.0 和
Eigen 3.4.0，但在 provenance receipt 生成前不得写成已接受事实。

### 3.2 打包

方案优先级固定为：

1. **A：官方发行流水线预构建共享库。** HUNDUN 与 Cantera/依赖使用同一
   GCC 11/libstdc++ ABI；安装包用相对 RPATH 携带共享库和 Cantera data。
2. **B：可复现源码 escrow。** `third_party` 保存固定 source archive 或
   maintainer rebuild metadata；configure 不联网；若上游构建强制 Python，
   B 不能成为普通用户唯一构建方式。
3. **C：受控 CMake port。** 只有 A/B 无法满足无 Python 和可复现安装时另行
   批准。所有移植文件仍属于 Cantera third-party，不进入 HUNDUN `src/`。

普通安装包不分发 glibc，不要求用户安装 compiler、Python 或 Cantera。

官方包只携带冻结的 Release Cantera/依赖共享库。GCC 11 的普通 HUNDUN Debug
构建可链接同一共享库，但必须保持 libstdc++、exceptions、RTTI、
`_GLIBCXX_USE_CXX11_ABI=1` 一致，并禁止 `_GLIBCXX_DEBUG` 等改变 STL ABI 的选项。
Clang/libc++ 不得与 GCC/libstdc++ Cantera 二进制混链。若维护者另建 Debug 或
sanitizer Cantera，必须使用独立 artifact root、manifest 和 build profile，不能与
Release `.so` 混装，也不进入默认发行包。

安装包升级按完整 bundle 原子替换；不得只替换 Cantera 或某个传递依赖 `.so` 后沿用
旧 ABI manifest。包大小只记录，不作为科学 gate。

### 3.3 反应状态

Stage 4 保存全部 `N_s` 个保守组分密度和总热化学混合物焓：

```text
rho Y_1 ... rho Y_Ns
rho h_tc
h_tc(T,Y) = sum_k Y_k h_k(T)
```

`h_tc` 含 formation enthalpy。禁止把 Stage 2 的 `h=cp*T` 原地重新解释为
新变量。温度只由 `(p0,h_tc,Y) -> T` 反演；`T`、`cp`、`mu`、`lambda`、
`D_k` 和 reaction rate 是 attempt-local 派生量。

所有组分一起通过共享 convex-simplex/element consistency 算子处理。不得逐组分
clipping 后归一化；仅允许有固定 machine-precision budget 的 roundoff
canonicalization，超出预算即失败。

### 3.4 输运

Cantera 提供 molecular mixture-averaged diffusion coefficient。HUNDUN 形成
组分扩散通量并施加 correction velocity，使：

```text
sum_k F_k = 0
```

WALE 只提供 SGS/turbulent transport。multicomponent diffusion、Soret、Dufour
和 barodiffusion 延期。

壁面 v1 只支持 non-catalytic、impermeable species、adiabatic 或 isothermal。

### 3.5 时间推进

正式名称为 `partitioned second-order C-T-C coupling`，不把完整多步流动推进
误称为经典一步 Strang flow：

- 组分和 `h_tc` 使用 `C(dt/2)-T(dt)-C(dt/2)`；
- momentum 保留 BDF2；
- 每个成功 step 恰好两次 PISO corrector；
- BDF2 history 只保存完整 accepted endpoint；
- 第二个 PISO 必须看到第二 chemistry half 后的最终热化学状态；
- chemistry 由两个 half-step 的积分增量形成全步 source，不能用 endpoint
  instantaneous rate 代替；
- chemistry 不向 `rho h_tc` 添加独立 heat-release source；
- open/fixed-`p0` 先实现，closed/partially-closed 通过独立 `p0` proof 后加入。

`4R-0` 必须在固定两次 PISO 内选择：预测通量已足够，或 PISO #2 内需要一个
守恒 predictor-to-final delta-flux correction。不得加入第三 corrector、滤波或
经验 damping。若两个候选均无法证明二阶和守恒，MISDC 只能通过新设计进入
post-v1，不在当前 task 内临时添加。

### 3.6 workspace

每个 MPI rank 一个 `CanteraBackendRuntime`；每个 worker thread 一整套独立
Solution/Thermo/Kinetics/Transport/Reactor/Integrator workspace。禁止跨线程共享
可变 Cantera 对象。workspace 不进入 Checkpoint；Restart 后按 mechanism identity
重建。

### 3.7 Stage 4 出口

Stage 4 接受后必须具备：

- 无 Python 的 bundled Cantera package；
- open 和已证明的 closed reacting-flow；
- all-species/total-enthalpy conservative transport；
- IBM+WALE reacting path；
- Checkpoint v4 和 diagnostics provider registry；
- 0D/PSR、rollback、Restart、MPI 和短程 reacting smoke；
- governance version receipt `0.3.0`。

## 4. Stage 5：ESF/TPDF 与 TCR

### 4.1 状态

```text
MeanState
StochasticFieldSet[N]
PsrShadowState
TcrRootState
```

- `MeanState` 是守恒、PISO、密度、面通量和 spray force 的唯一流动权威；
- 每个 spatial rank 保存其 owned/ghost 空间上的全部随机场；
- MPI 只做空间分解，不建立 field-owner communicator；
- 随机场跨时间步连续，只在每个 accepted stochastic interval 生成新 Wiener
  increments；
- `PsrShadowState` 非输运、每步从 committed mean state 确定性构造；
- `TcrRootState` 保存 signed branch 和 versioned 九项跨步历史；
- 持久布局为 component-major SoA，chemistry pack 只存在于 scratch。

N 必须为偶数：

```text
minimum N = 2
recommended transient N = 4
N = 2: debug 和依赖时间平均的稳态算例
N = 8/16/32: 可选采样敏感性 screen
```

`N=2` 诊断标记 `minimal_pair_sampling`，不得宣称瞬时高阶统计已收敛。

### 4.2 低雷诺一致 ESF

Stage 5 采用低雷诺一致、unity-Lewis ESF。分子扩散只进入确定性扩散项，
Wiener 系数只包含 SGS/turbulent diffusivity；`mu_t=+0.0` 时随机项必须精确
退化为 `+0.0`。TPDF 模式使用单一 composition diffusivity，要求
`Pr_t=Sc_t`；Stage 4 非 TPDF 模式仍可使用独立配置。

同一随机场内所有 species 和 `h_tc` 共享同一 Wiener 向量。不得按 cell 或
species 生成空间白噪声。

### 4.3 counter RNG

HUNDUN 使用 versioned counter RNG，首选依据公开 Philox 算法独立实现并用
公开 golden vectors 验证。ESF key 固定包含：

```text
algorithm_version
user_seed
accepted_step
stochastic_stage
field_pair
spatial_direction
purpose_tag
```

禁止包含 `global_cell_id`、`species_id`、MPI rank、PISO iteration 或 failed
attempt ordinal。场成对满足 `DeltaW_(2m+1)=-DeltaW_(2m)`。retry 保留符号，
按新 `dt` 重算幅值；拒绝步不推进 accepted clock。Brownian bridge/tree 和
stochastic-error adaptive stepping 延期。

### 4.4 IEM、reaction ensemble 和 mean closure

唯一 IEM 更新为：

```text
zeta_new = mean + exp(-dt/(2*tau)) * (zeta_old - mean)
```

均值不变，方差按 `exp(-dt/tau)` 衰减。`kappa=1` 是同一实现的严格 IEM
极限，不维护第二套 mixer。

每个随机场只通过 Stage 4 `ChemistryBackend` 反应。两个 chemistry half 的积分
增量由 `ReactionEnsemble` 求平均后，通过 Stage 4 source transaction 更新
MeanState chemistry closure；不得用 field 0、不得复制 sample mean 覆盖
MeanState、不得使用 endpoint rate。

元素和 moment consistency 依据公开修正 ESF 方程独立推导。禁止逐字段 clipping。

### 4.5 TCR

TCR 是用户提出的方法，不是开源第三方组件。HUNDUN 独立实现并在产品代码和
文档引用：

- DOI `10.1016/j.proci.2026.106128`；
- DOI `10.1016/j.cja.2026.104123`。

核心合同：

```text
D = 1 - 4 eta (1-eta) R
s_kappa = +/- sqrt(D)
kappa = (1+s_kappa)/(2(1-eta))
```

- old accepted state 决定 branch sign，current `D` 决定 magnitude；
- 无 event certificate 不跨 fold 猜测分支；
- source control 只调整 `R`，不调整 `eta`；
- `D<0` 不 clipping；
- weak denominator、algebraic invalid、statistical unresolved、branch crossing、
  fold unresolved 和 fallback 分别报告；
- 模式为 `off/shadow/experimental/validated`，默认 `shadow`；
- 只有 COAST point/state replay 通过后才允许 `validated` feedback。

### 4.6 COAST 参考边界

COAST 可用于两类私有 oracle：

1. **TCR exact oracle：** 绑定用户确认的准确 source root、commit、manifest、
   theory status 和文件 hash；只允许 allowlisted pure modules。
2. **ESF semantic/oracle：** 可参考任意成熟 COAST 版本，但每次仍登记 realpath、
   commit 和 hash；用于核对 field persistence、noise staging、antithetic pairing、
   IEM、边界、Restart 和单步结果。

正式读取/截取前必须请用户确认实际 COAST 根目录。allowlisted source 只进入
generated temporary build tree；HUNDUN 自有 standalone driver 通过单独进程和
合成输入调用。COAST source、case、research data 和 executable 不进入 Git、
安装包或产品 ABI。产品实现不得翻译、机械重写或模仿 COAST 控制流。

HUNDUN TCR differential 主配置使用 `N=4`；`N=2` 使用解析、守恒和 steady
time-average 证据，不宣称 finite-field statistics 与 COAST 完全等价。

### 4.7 Stage 5 出口

- ESF/IEM、reaction ensemble、PSR shadow 和 TCR；
- N=2/4、small 1/2/4-rank、Restart 和 rollback；
- Checkpoint v5、diagnostics 和 driver combinations；
- COAST ESF/TCR oracle receipt；
- 无 Vblowoff、Flame D 或其他长科学矩阵；
- governance version receipt `0.4.0`。

## 5. Stage 6：稀相喷雾

### 5.1 物理范围

v1 使用 dilute Eulerian--Lagrangian point parcels：

- 球形液滴、parcel multiplicity、内部均温；
- 单组分挥发性液体；
- 一个 liquid material 映射到一个 gas species；
- Schiller--Naumann drag，含 Stokes 极限；
- Ranz--Marshall heat/mass transfer；
- Abramzon--Sirignano 单组分 Stefan-flow 修正蒸发；
- mass、momentum、species、total thermochemical enthalpy two-way coupling；
- deterministic injector、MPI migration、static IBM rebound；
- `none` 和 TAB breakup。

`d^2` law 只作解析 oracle，不是产品蒸发路径。

### 5.2 状态和内存

`SprayParcelState` 至少保存：

```text
global id
position, velocity
single-droplet mass/diameter
parcel multiplicity
droplet temperature
liquid material identity
owner cell/rank
age
injector identity
TAB deformation/history
```

container 使用 pure SoA。一个 parcel 同时只有一个 owner rank。migration 的
outgoing/incoming buffer 属于 trial，collective accept 后才发布。

### 5.3 液体物性

Cantera 不承担液体数据库职责。`LiquidPropertyService` 支持 constant、temperature
polynomial、Antoine 和 Clausius--Clapeyron 形式，并返回单位、适用范围和明确
状态。每套 property pack 独立登记来源、许可证和 hash。

低成本通用性验收使用两个单组分 surrogate：

```text
kerosene surrogate: n-dodecane
gasoline surrogate: iso-octane
```

真实航空煤油和汽油是多组分混合物；文档不得把上述 surrogate 写成真实燃料
预测。具体 reduced mechanisms 优先从用户确认的 COAST `EXEC/Fuels` 筛选，
各自记录 mechanism identity、species mapping、license 和 SHA。无再分发许可的
机制只能作为本机外部验收资产或用户输入。

### 5.4 交换和事件

气相 sampling 只通过 MeanState、ThermodynamicsService 和
TransportPropertyService。parcel 不直接访问 Cantera。采样和源项沉积共享
`ParcelGridCouplingStencil`，权重和为 1。

液滴损失与气相增加通过一个 `SprayCouplingTransaction` 成对登记。drag、vapor
momentum、sensible transfer、latent contribution 和 vapor thermochemical
enthalpy 必须由 `6F-1` 的总量推导防止漏项或双计。

完全蒸发使用 event-to-zero，不允许负 mass/diameter 后 clipping。trajectory
跨 cell 时分段沉积；outlet removal 进入流出 ledger。

### 5.5 与 Stage 5 组合

parcel trajectory 和交换量只由 MeanState 计算一次。相同 vapor/enthalpy
source 通过确定性共同源算子施加到所有 stochastic fields，随后调用 Stage 5
element/moment consistency。禁止每个随机场创建或推进一套 parcels。

spray exchange 属于中央 transport/source block；第二 chemistry half 和 PISO #2
必须看到最终蒸发源。每个 gas macrostep 固定一次 parcel predictor 和一次
corrector；parcel 内部可按 particle CFL/thermal/evaporation timescale subcycle，
但只能形成一份积分交换报告。不得增加第三次 PISO。

### 5.6 IBM 与 TAB

static IBM 使用 trajectory segment earliest-hit 和 Stage 3 几何 authority。
v1 只支持有明确 restitution 的 rebound；沉积、film 和 splash 延期。

TAB 在 spray evaporation MVP 后实现，但属于 Stage 6 最终接受条件。child mass、
momentum、temperature、multiplicity 和 ID 必须守恒/确定；breakup RNG 使用与
ESF、injector 隔离的 counter domain。

### 5.7 Stage 6 出口

- parcel container、injection、migration、Restart 和 rollback；
- drag、heating、single-component evaporation 和 exact two-way budgets；
- IBM rebound 和 TAB；
- n-dodecane/iso-octane 双机制低成本接口证据；
- MeanState + N=2/4 ESF common-source 组合；
- Checkpoint v6、diagnostics 和 driver；
- `0.5.0-rc.1` development-complete candidate。

## 6. 跨阶段执行

### 6.1 默认串行

任务图保留并行边，但执行策略默认串行：

```text
Stage 4 accepted
-> 主 agent 在阶段节点提出串行/并行建议并等待用户指示
-> 未明确批准并行时执行 Stage 5
-> Stage 5 accepted
-> 再执行 Stage 6
```

只有用户在阶段节点明确批准，Stage 5/6 的纯模块才可从 Stage 4 service-contract
freeze commit 并行开发。计划中“可以并行”从不等同于自动授权并行。

无论实际串行或并行，accepted history 固定为 Stage 4 -> Stage 5 -> Stage 6。
并行 development commits 必须更新到前一 stage accepted head 后才进入正式验收；
accepted commits 不重写。

### 6.2 版本

```text
Stage 3 accepted       0.2.0
Stage 4 receipt        0.3.0
Stage 5 receipt        0.4.0
Stage 6 dev complete   0.5.0-rc.1
cross-stage accepted   1.0.0
```

governance 是唯一开发和验收仓库。默认不把 `0.3.0`、`0.4.0` 中间状态投影到
product repo；`V1_ACCEPT` 后只同步一次 `1.0.0`。任何中间投影需用户单独授权。

### 6.3 共享文件

顶层 CMake、root schema/driver dispatch、Checkpoint registry、diagnostics kind
registry、capability root table 和 version banner 只由 integration lane 修改。
stage task 只提供 module-local descriptor/registration function，避免并发修改 central
switch。

## 7. 测试和资源

### 7.1 task gate

每个 task 只要求：

- mutation-sensitive RED；
- 直接受影响的 unit/header/policy；
- 最多一个 12^3 或更小产品路径 smoke；
- collective 改动才补 small 1/2-rank；
- public header 改动才补 standalone header；
- build graph 改动才补 tests-off/linkage；
- task-focused Debug；
- 主 agent 一次 requirements/quality/caller-impact/complete-task-diff review。

不为每个 task 机械重复 Release、ASan、UBSan、1/2/4-rank 或长数值矩阵。

### 7.2 milestone gate

milestone 最大使用 24^3，覆盖 Stage 4 reacting、Stage 5 N=2/4 和 Stage 6 parcel
迁移/双向守恒。长作业不得阻塞后续开发。

### 7.3 最终冻结矩阵

全部产品、测试和公共文档完成后冻结 `0.5.0-rc.1`：

- M1：Stage 4 reacting-flow 24^3，1/2/4-rank；
- M2：Stage 5 TPDF/TCR 24^3，N=4，1/2/4-rank；
- H1：唯一 48^3 n-dodecane integrated run，包含 IBM、WALE、reacting、
  N=4 ESF/TCR 和 two-way spray；
- iso-octane 只运行 8^3/12^3 接口验收；
- Debug full affected、focused Release、small focused ASan/UBSan；
- Checkpoint、diagnostics、driver、headers、tests-off、package、RPATH、ABI、DCO、
  provenance 和 capability ledger。

永久不运行 96^3、Vblowoff、两套燃料各一遍 48^3、高负载 sanitizer MPI 或
为“更放心”添加的重复矩阵。

H1 correctness failure 阻塞 `V1_ACCEPT`；wall time/RSS 只记录，不设跨机器阈值。

### 7.4 detached runner

M/H 作业记录 candidate HEAD、tree、dirty diff、binary hash、command、environment、
CPU binding、MPI、log、exit、elapsed、peak RSS 和 log hash。同时最多一个 H。
冻结 H 运行期间不修改产品或测试源码。

## 8. 版权、来源和机制

- HUNDUN 原创代码 Apache-2.0，提交遵守 DCO；
- Cantera 和全部传递依赖保留原许可证、版权、SPDX、upstream revision 和 patch；
- `THIRD_PARTY_NOTICES` 记录 binary 和 transitive dependency；
- 不使用 Cantera、Caltech、Sandia 或贡献者名称背书；
- 机制文件是独立版权对象；Cantera BSD 不覆盖 mechanism data；
- 不复制 GPL OpenFOAM/Basilisk 源码；GPL 软件只提供公开数学/架构参考；
- AMReX/Pele/Code_Saturne 设计只作算法与架构参考，是否复用代码必须单独通过
  license/provenance task；当前计划不复制其源文件；
- BOFFIN 和私有研究数据不访问；
- COAST 只按本设计的受控私有 oracle 边界使用，不是产品 source ancestor；
- 不发布、不 push，除非用户以后明确授权。

## 9. 能力声明和接受结论

Stage 4、5、6 各自产生 exact-HEAD stage receipt。v1 最终分两步：

```text
V1_DEVELOPMENT_COMPLETE
V1_ACCEPT
```

`V1_DEVELOPMENT_COMPLETE` 表示 `0.5.0-rc.1` 功能、低成本门和文档完成，可启动
最终冻结矩阵；不等于发布接受。

`V1_ACCEPT` 要求 M1/M2/H1、package/ABI/RPATH、provenance、DCO、exact-HEAD
manifest 和 product projection 全部完成，最终 product 版本为 `1.0.0`。

即使 `V1_ACCEPT`，仍不得宣称：

- 通用 Linux 或非冻结 toolchain 兼容；
- NativeChemistryBackend 已存在；
- Vblowoff、Flame D、真实 Jet-A/汽油实验级预测已经重新验证；
- multicomponent droplets、dense spray、liquid film、KH--RT、GPU、AMR 或
  rank-changing Restart；
- 48^3 以外的大规模性能和强扩展性。

## 10. 配套文档

本规格与以下文件共同构成 Stage 4--6 权威规划：

- `docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md`；
- `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-reacting-flow.md`；
- `docs/superpowers/plans/2026-08-09-hundun-flow-stage5-esf-tpdf-tcr.md`；
- `docs/superpowers/plans/2026-08-09-hundun-flow-stage6-spray.md`；
- `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md`。

当规格与计划冲突时，本规格控制科学、状态、API、版权和能力边界；各 stage
计划控制 task 文件白名单、RED、依赖、执行和局部验收；跨阶段计划控制串行/
并行选择、集成顺序、版本和最终 exact-HEAD 接受。

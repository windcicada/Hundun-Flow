# HUNDUN-FLOW Stage 4--6 公开参考与独立实现目录

日期：2026-08-09

本文件记录 Stage 4--6 可以复用的公开数学、物理模型和软件架构思想。它不是
源代码搬运清单。除单独登记的 third-party component 外，HUNDUN 根据公开方程
独立推导、使用项目命名和代码风格，并通过 mutation-sensitive RED 验证。

## 1. 参考等级

| 等级 | 用途 | 允许 | 禁止 |
|---|---|---|---|
| P | peer-reviewed/public equation | 推导方程、量纲、极限和 oracle | 复制论文附带代码或大段文字 |
| D | 官方软件文档 | 理解接口和职责分层 | 假设其控制流适合 HUNDUN |
| S | 公开源码，只读架构参考 | 记录模块边界和设计得失 | 复制、翻译、机械改写 |
| T | 合法 third-party component | 保留原文件/许可证/SPDX/patch | 改名后混入 HUNDUN 原创源码 |
| C | 私有 COAST oracle | 受控合成输入、进程外差分 | 进入产品、Git、安装包或公共 ABI |

GPL 软件只能作为 P/D/S 参考，不能提供产品源码。

## 1.1 开源项目活跃性快照

快照日期为 2026-08-09。“活跃”只表示近期仍有官方 release、development 或
issue/PR 活动，不承诺 API/ABI 稳定，也不改变本文件的源码复用边界。

| Project | Snapshot evidence | Assessment and use |
|---|---|---|
| [Cantera](https://github.com/Cantera/cantera) | stable 3.2.0，development 4.0.0a2；官方仓库仍有 issue/PR 和 CI | 活跃；v1 固定 3.2.0，不跟随 development |
| [AMReX](https://github.com/AMReX-Codes/amrex) | release 26.06（2026-06），持续 PR/issue | 高度活跃；只参考 SoA、particle ownership 和 parallel checkpoint，不引入 AMR runtime |
| [PeleLMeX](https://github.com/AMReX-Combustion/PeleLMeX) | 2026 年发布 1.2.0/1.3.0，持续 submodule、CI 和 bug-fix release | 活跃；参考 low-Mach reacting architecture，不复制 SDC/AMR/GPU 控制流 |
| [PelePhysics](https://github.com/AMReX-Combustion/PelePhysics) | development branch、持续 issue/PR；已承接 spray、soot、radiation modules | 活跃；当前 spray/chemistry architecture reference |
| [PeleMP archive](https://github.com/AMReX-Combustion/PeleMP) | 2024-06-27 archived；官方说明功能已迁往 PelePhysics | 不再活跃；只作历史论文/文档 reference，开发参考改用 PelePhysics |
| [SUNDIALS](https://github.com/LLNL/sundials) | version 7.8.0（2026-06），持续 issue/PR | 活跃；v1 只消费 Cantera 实际冻结的兼容版本，不能混用最新 ABI |
| [Code_Saturne](https://github.com/code-saturne/code_saturne) | EDF public mirror，当前仍有 issue/PR | 维护中；GPL，仅参考 Lagrangian 职责分层和 two-way bookkeeping |
| [OpenFOAM Foundation v11](https://github.com/OpenFOAM/OpenFOAM-11) | 固定版本仓库；GPL model catalog | 版本化只读参考；不作为依赖，不以其活跃性决定 HUNDUN 架构 |

论文、NASA/CHEMKIN 报告和 Random123 算法按稳定公开数学资料处理；它们不需要
持续 release 才具有参考价值。实际 third-party 锁定仍由 task provenance receipt 的
URL、revision、archive/binary SHA-256 和许可证决定。

## 2. 低马赫反应流和 operator coupling

### Strang、低马赫和二阶约束

1. G. Strang, “On the Construction and Comparison of Difference Schemes,”
   DOI [10.1137/0705041](https://doi.org/10.1137/0705041)。
   - 参考点：对称 `C(dt/2)-T(dt)-C(dt/2)` 的一步算子条件。
   - 避免：把多步 BDF2 momentum map 直接称为经典 Strang 子流。
2. M. S. Day and J. B. Bell, low-Mach reacting-flow algorithm，DOI
   [10.1088/1364-7830/4/4/309](https://doi.org/10.1088/1364-7830/4/4/309)。
   - 参考点：reaction half / advection-diffusion / reaction half；积分 reaction
     increment；总焓闭合。
3. A. Nonaka et al., DOI
   [10.1080/13647830.2012.701019](https://doi.org/10.1080/13647830.2012.701019)。
   - 参考点：low-Mach reaction/transport splitting、constraint source timing、
     MISDC 对照。
   - 避免：将其 SDC/MAC projection 控制流照搬到固定两次 PISO。
4. A. Nonaka et al., DOI
   [10.1080/13647830.2017.1390610](https://doi.org/10.1080/13647830.2017.1390610)。
   - 参考点：EOS drift 和 thermodynamic consistency 的可测误差边界。
5. J. Sportisse, DOI
   [10.1006/jcph.2000.6495](https://doi.org/10.1006/jcph.2000.6495)。
   - 参考点：stiff reaction splitting order reduction。

### PeleLMeX/PeleC

- PeleLMeX model：
  [official model documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html)。
- PeleLMeX source/architecture：
  [official implementation documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/Implementation.html)。
- PeleC：
  [official repository](https://github.com/AMReX-Combustion/PeleC)。

可复用设计：

- thermodynamic pressure 与 mechanical pressure 分离；
- integrated species source 和 divergence constraint 使用同一数据；
- chemistry/transport/property 职责分层；
- final projection/synchronization 的必要性；
- mechanism identity、transport 和 chemistry 的一致输入。

不复用：

- AMR hierarchy；
- SDC iteration；
- MAC/nodal projection 控制流；
- GPU portability layer；
- case-local compile-time user code。

## 3. Cantera

- 官方项目：[Cantera](https://cantera.org/)。
- C++ API：[stable C++ documentation](https://cantera.org/stable/cxx/)。
- mechanism YAML：
  [input reference](https://www.cantera.org/stable/yaml/index.html)。
- ReactorNet：
  [C++ ReactorNet](https://www.cantera.org/3.2/cxx/d7/def/classCantera_1_1ReactorNet.html)。

采用：

- ideal-gas mixture thermo；
- species thermo、kinetics 和 mixture-averaged gas transport；
- stiff local integration；
- mechanism parsing 和 failure diagnostics。

不采用：

- Python runtime/interface；
- Cantera mesh、flow、PISO、Checkpoint 或 driver；
- 将 constant-pressure/constant-volume reactor 名称当成 HUNDUN conservation
  contract 的替代证明；
- 从 Cantera 零散复制函数进 HUNDUN `src/`。

固定候选：Cantera v3.2.0，commit
`4a8358eb80cfeb50474386b5f9ec0b3a83519889`，archive SHA-256
`a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b`。

### Post-v1 NativeChemistryBackend 的公开数学路线

1. B. J. McBride, M. J. Zehe, S. Gordon，*NASA Glenn Coefficients for
   Calculating Thermodynamic Properties of Individual Species*，NASA/TP-2002-211556，
   [NASA NTRS](https://ntrs.nasa.gov/citations/20020085330)。
   - 参考点：NASA polynomial 的 `cp/R`、`h/RT`、`s/R` 形式和温区数据身份；
   - v1 只冻结接口兼容性，不实现 parser 或 polynomial evaluator。
2. R. J. Kee, F. M. Rupley, J. A. Miller，*Chemkin-II: A Fortran Chemical
   Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics*，
   DOI [10.2172/5681118](https://doi.org/10.2172/5681118)。
   - 参考点：elementary reaction、reversible/third-body/falloff 和 mechanism
     interpreter 的标准问题分解；
   - 只采用公开方程和格式概念，不复制 Fortran package 或兼容 ABI。
3. A. C. Hindmarsh et al.，*SUNDIALS: Suite of Nonlinear and
   Differential/Algebraic Equation Solvers*，DOI
   [10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020)，以及
   [LLNL CVODE documentation](https://computing.llnl.gov/projects/sundials/cvode)。
   - 参考点：自适应 stiff BDF/Newton、Jacobian、dense/sparse/Krylov linear solve
     的职责边界；
   - 当前 v1 由 bundled Cantera/SUNDIALS 提供，未来 native backend 独立规划。

未来原生实现以 Cantera value-level outputs 作黑盒 oracle。不得从 Cantera、CHEMKIN
或 SUNDIALS 源码翻译函数，也不得在 v1 同时启动该路线。

## 4. ESF/TPDF

1. L. Valiño, “A Field Monte Carlo Formulation...”，DOI
   [10.1023/A:1009968902446](https://doi.org/10.1023/A:1009968902446)。
   - 参考点：Eulerian stochastic fields 和 field-wise Wiener process。
2. L. Valiño, R. Mustata, K. B. Letaief, DOI
   [10.1007/s10494-015-9687-0](https://doi.org/10.1007/s10494-015-9687-0)。
   - 参考点：低雷诺一致形式；从 Wiener coefficient 移除 molecular diffusion；
     laminar limit。
3. S. Xu, S. Zhong, F. Zhang, X.-S. Bai, DOI
   [10.1016/j.combustflame.2021.111577](https://doi.org/10.1016/j.combustflame.2021.111577)。
   - 参考点：有限 stochastic fields 的 element-mass consistency。
4. Mustata et al., DOI
   [10.1016/j.combustflame.2005.12.002](https://doi.org/10.1016/j.combustflame.2005.12.002)。
   - 参考点：LES/PDF Eulerian Monte Carlo field 应用和统计解释。

HUNDUN 独立实现：

- spatially shared Wiener increment per field/direction；
- antithetic field pairs；
- exact IEM；
- all fields per spatial rank；
- element/moment consistency；
- N=2 steady/debug 和 N=4 transient evidence split。

## 5. Counter RNG

J. K. Salmon et al., “Parallel Random Numbers: As Easy as 1, 2, 3”，DOI
[10.1145/2063384.2063405](https://doi.org/10.1145/2063384.2063405)。

采用：

- Philox counter/key/domain separation；
- 无 mutable cursor；
- published golden vectors；
- decomposition/restart identity。

ESF semantic key 不含 cell/species/rank。spray injector、TAB breakup 和其他 RNG
用途使用不同 `purpose_tag`，不能共享 stream。

## 6. TCR

TCR 是用户方法，产品必须引用：

- [10.1016/j.proci.2026.106128](https://doi.org/10.1016/j.proci.2026.106128)；
- [10.1016/j.cja.2026.104123](https://doi.org/10.1016/j.cja.2026.104123)。

相关 ESF/IEM 和守恒参考：

- [10.1023/A:1009968902446](https://doi.org/10.1023/A:1009968902446)；
- [10.1007/s10494-015-9687-0](https://doi.org/10.1007/s10494-015-9687-0)；
- [10.1016/j.combustflame.2021.111577](https://doi.org/10.1016/j.combustflame.2021.111577)。

COAST 只提供 C 级 point/state oracle，不是公开上游。正式取证前记录用户确认的
source root、commit、manifest、status 和 allowlist hashes。

## 7. Lagrangian parcels 和 MPI

### AMReX ParticleContainer

[AMReX Particle documentation](https://amrex-codes.github.io/amrex/docs_html/Particle.html)。

可复用设计：pure SoA、stable identity、owner redistribution、parallel
checkpoint、application-owned component schema。

不复用：AMR level、DistributionMap API、AMReX checkpoint 格式或 GPU layer。

### Code_Saturne

[Lagrangian particle module](https://www.code-saturne.org/cms/web/documentation/overview/modules/lagrangian)。

可复用设计：frozen/one-way/two-way capability separation；动力、热量和质量耦合
分别登记；boundary interaction 与 statistics 独立。

Code_Saturne 的 GPL 源码不复制。

## 8. 液滴传热、蒸发和 breakup

1. B. Abramzon and W. A. Sirignano, DOI
   [10.1016/0017-9310(89)90043-4](https://doi.org/10.1016/0017-9310(89)90043-4)。
   - 参考点：single/multicomponent droplet vaporization、Stefan-flow 修正、
     film properties；v1 只采用 single-component 子集。
2. W. E. Ranz and W. R. Marshall, “Evaporation from Drops”。
   - 参考点：sphere Nusselt/Sherwood correlation 和 `Re->0` 极限。
3. Schiller--Naumann drag correlation。
   - 参考点：Stokes limit 和 finite-Re sphere drag；产品先验证极限、符号和
     range status。
4. P. J. O'Rourke and A. A. Amsden, TAB，DOI
   [10.4271/872089](https://doi.org/10.4271/872089)。
   - 参考点：deformation ODE、threshold 和 child generation budget。

### PelePhysics spray 与 PeleMP 历史资料

- [spray equations](https://amrex-combustion.github.io/PelePhysics/Spray.html)；
- [spray inputs](https://amrex-combustion.github.io/PeleMP/SprayInputs.html)。

PeleMP 仓库已归档，继续开发的 spray modules 已迁入 PelePhysics。前者只保留历史
论文/输入文档价值；代码架构判断以当前 PelePhysics 为准。

可复用：dilute point parcel、infinite-conductivity droplet、one-third film rule、
liquid property service、species mapping、single-drop fixtures。

不复用：Python validation runtime、AMR/GPU container、compile-time fuel count、
multicomponent liquid path 或源码实现。

### OpenFOAM spray

[official source guide](https://cpp.openfoam.org/v11/dir_001528_001553.html) 显示 ETAB、
Pilch--Erdman、Reitz--Diwakar、KH--RT 等成熟模型目录。它只用于比较模型职责和
识别应延期的经验模型；GPL 源码不得复制或翻译。

## 9. Surrogate mechanisms

- n-dodecane：
  [Cantera example](https://cantera.org/stable/examples/python/reactors/ic_engine.html)、
  [LLNL mechanism page](https://combustion.llnl.gov/mechanisms/alkanes/n-dodecane)。
- iso-octane：
  [LLNL mechanism page](https://combustion.llnl.gov/mechanisms/alkanes/iso-octane-version-3)。

这些链接只证明公开候选存在，不自动授予再分发许可。Stage 6 `6F-4` 优先审计
用户确认的 COAST `EXEC/Fuels` reduced mechanisms；每个 mechanism 单独记录
来源、license、revision、conversion command 和 source/output SHA。

## 10. 设计得失摘要

| 参考 | 复用 | 避开 |
|---|---|---|
| Cantera | thermo/kinetics/transport/local integration | Python、mesh、flow authority |
| PeleLMeX | low-Mach constraint 和 integrated source 思路 | SDC/AMR 控制流 |
| Valiño ESF | SPDE、field ensemble、laminar limit | per-cell noise |
| Random123 | stateless counter RNG | shared mutable engine |
| AMReX particles | SoA、identity、migration、checkpoint | AMR/runtime dependency |
| PeleMP spray | A--S、film properties、liquid service | Python validation、compile-time fuels |
| OpenFOAM | 模型目录和职责比较 | GPL 源码/模板架构 |
| COAST | 私有 ESF/TCR oracle | 产品 source/control-flow ancestor |

## 11. 引用和变更记录要求

每个 semantic-port task receipt 至少记录：

```text
reference URL / DOI
upstream software release or revision（若适用）
adopted mathematical behavior
rejected architecture/code behavior
independent RED and mutations
license/provenance decision
```

若实现发现公开参考之间存在冲突，由主 agent先更新规格和数学推导，再允许修改
产品代码。不得以“某开源软件这样做”为唯一正确性依据。

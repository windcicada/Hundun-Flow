# HUNDUN-FLOW Stage 3 公开算法与软件参考

**状态：** `REFERENCE_BASELINE_FOR_REVIEW`

**核对日期：** 2026-08-09

**用途：** 为 Stage 3 剩余任务提供数学、数据流和测试思想；不授权复制实现。

## 1. 使用规则

本文件区分三件事：

1. 可复用的数学关系和公开行为；
2. 可借鉴的责任划分、数据流和验证方法；
3. 不进入 HUNDUN-FLOW 的源码、ABI、宏、数据结构和运行时依赖。

任何实现任务都要先把参考行为改写为 HUNDUN-owned equation/interface note，再写
mutation-sensitive RED。worker 只能看到本文件列出的 reference point，不得把上游
函数当作翻译模板。不得访问 BOFFIN、COAST、COAST-2 或研究数据。

## 2. 固定公开 revision

| 项目 | Revision | 许可证 | 活跃性背景 | 本项目使用方式 |
| --- | --- | --- | --- | --- |
| AMReX | `5273b558c13011573e1b6bf71860db4fcc1f0cfb` | BSD-3-Clause | development 持续更新 | EB 数据分层、regular/irregular dispatch、flux accounting 语义 |
| AMReX-Hydro | `ae9a742c17aa064d299826c60585c1e61b1d6117` | BSD-3-Clause | development 持续更新 | projection 阶段职责、EB face extrapolation 与测试组织 |
| incflo | `7307d8725c2a538f09cafbeacbfeb63e0fb11d22` | BSD-3-Clause | 官方文档说明 development 活跃 | variable-density predictor/corrector、density update、EB/no-EB 分离 |
| OpenFOAM-dev | `627e2f909cb6a08dbb3a74e9a34aa632a975650e` | GPL-3.0-or-later | master 持续更新 | WALE 公式责任边界、effective viscosity 消费关系 |
| Basilisk | 2026-08-09 public-page snapshot；下列逐文件 SHA-256 | GPL-3.0 | 公开源码持续维护 | topology guard、embedded gradient/flux/force 数学对照 |
| Nek5000/gslib | `95acf5b42301d6cb48fda88d662f1d784b863089` | BSD-3-Clause | master 持续更新 | global donor owner/local element/weights 的查找语义 |
| PETSc | `9183a15b9d20dbb91d91a365783fab658bad1796` | BSD-2-Clause | main 持续更新 | 未来 preconditioner capability 评估；Stage 3 不依赖 |
| Trilinos | `38df4d5b4da2233f9150d3644589803f45942a40` | BSD-3-Clause | develop 持续更新 | 未来 solver package 边界评估；Stage 3 不依赖 |
| Bell--Marcus projection paper | DOI `10.1016/0021-9991(92)90011-M` | 论文，仅引用数学/阶段思想 | 1992 fixed publication | variable-density predictor/projection 的独立背景 |

固定 revision 只用于可复现比较，不表示 HUNDUN-FLOW 包含这些 revision 的代码。

### 2.1 项目问题、架构和依赖对照

| 项目 | 主要解决的问题 | 公开架构/关键技术依赖 | 2026-08-09 活跃性判断 | HUNDUN 结论 |
| --- | --- | --- | --- | --- |
| AMReX | block-structured AMR、并行网格、EB 和 coarse/fine 守恒 | C++、MPI、可选 GPU backend；MultiFab/EB2/FluxRegister 分层 | development branch 和文档持续更新 | 复用几何事实分层、work partition 和守恒登记思想；拒绝 AMR 容器与 redistribution 默认化 |
| AMReX-Hydro | 可复用的 MAC/nodal projection 与 EB face treatment | AMReX、MLMG、Godunov/projection components | AMReX-Fluids 仍维护 | 复用 provisional/constraint/final flux 的阶段责任；拒绝替换 HUNDUN PISO/布局 |
| incflo | incompressible/low-Mach、variable-density、EB 应用级推进 | AMReX/AMReX-Hydro、MPI、可选 GPU；advance/predictor/corrector 分层 | 官方 development 文档和仓库仍更新 | 复用 coarse-grained composition 与 density stage order；拒绝输入系统、AMR 和 small-cell 路径 |
| OpenFOAM-dev | 通用有限体积 CFD 与 turbulence model integration | C++ field/object registry、runtime selection、MPI domain decomposition | `OpenFOAM-dev` 持续开发 | 仅复用 WALE model/output responsibility；GPL 实现、registry 和 wall-function 架构不进入产品 |
| Basilisk | 自适应 Cartesian flow、embedded boundaries 和紧凑离散 | C-like DSL/macros、tree/grid modules、MPI/OpenMP options | public source pages 持续维护 | 仅把 gradient/flux/force/topology 公式作为独立 oracle；GPL 宏和 cut-cell 布局不复用 |
| Nek5000/gslib | 分布式高阶网格中的全局 point location/interpolation | C、MPI、owner/local-element/weights、crystal-router family | gslib master 仍维护 | 复用 donor 是全局对象的语义；拒绝运行时依赖和每步全局 allgather |
| PETSc | scalable linear/nonlinear solvers 与 preconditioner lifecycle | C/MPI、KSP/PC/DM，可选外部 solver packages | main 和 release 文档持续更新 | 只用于未来 solver capability 边界；Stage 3 不引入依赖 |
| Trilinos | composable solver/preconditioner packages | C++/MPI、Tpetra/Belos/Ifpack2/MueLu | develop 仍维护 | 只用于未来 package/ownership 对照；Stage 3 不引入依赖 |

“活跃”只说明这些公开项目适合做设计对照，不提高其代码在本项目中的许可，也不把
live branch 当作固定算法基线；真正可复现的实现参考仍是上表的 commit 或页面 hash。

Basilisk public-page snapshot hashes：

| Public URL | SHA-256 on 2026-08-09 |
| --- | --- |
| `https://basilisk.fr/src/embed.h` | `27bdb129ec0a88f5055e940acf11d502f2edb5083216edb5eced0bbadb55c923` |
| `https://basilisk.fr/src/navier-stokes/centered.h` | `10be4cccbceb65d6d2c93fb7d6d1d0dc2aadd30b52b1d82b0e1d9f90b88b72d3` |
| `https://basilisk.fr/src/viscosity-embed.h` | `ce92bfb9fa11accfa746a1284a497fa009d390a8e67d12c47f274a306b7e614c` |
| `https://basilisk.fr/src/poisson.h` | `2d8e58452a9e466bee531c1ce2f05dfc110939b498bcb5ca958dfd108ceba38f` |
| `https://basilisk.fr/src/COPYING` | `0dcc9a456af8180079f663dea66faed97a983eec0fec51c54147e57179ddd5f8` |

若以后页面内容 hash 不同，只能显式建立新 reference receipt；不能静默把 live page 当成
本基线。

## 3. AMReX：几何事实、稀疏工作和守恒登记

公开入口：

- [AMReX repository at fixed revision](https://github.com/AMReX-Codes/amrex/tree/5273b558c13011573e1b6bf71860db4fcc1f0cfb)
- [Embedded Boundary documentation](https://amrex-codes.github.io/amrex/docs_html/EB.html)
- [`Src/EB/AMReX_EBCellFlag.H`](https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/Src/EB/AMReX_EBCellFlag.H)
- [`Src/EB/AMReX_EBFabFactory.H`](https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/Src/EB/AMReX_EBFabFactory.H)
- [`Src/EB/AMReX_EBFluxRegister.H`](https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/Src/EB/AMReX_EBFluxRegister.H)
- [`Src/EB/AMReX_EB_Redistribution.H`](https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/Src/EB/AMReX_EB_Redistribution.H)
- [AMReX BSD-3-Clause license](https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/LICENSE)

### 3.1 可复用点

`EBSupport::basic/volume/full` 表达“消费者需要哪些几何事实”，值得映射为 HUNDUN
现有 `ImmersedDomain`、active layout、Ghost/Wall plan 的只读 view。它提醒我们让
geometry 决定合法工作域、operator plan 决定离散、flow composition 决定时间顺序。

cell flag 的 regular/single-valued/covered 分类值得复用于 deterministic work
partition：

```text
regular active -> background kernels
interface active -> sparse immersed authority
inactive -> never read as physical state
```

`EBFluxRegister` 的价值不是其 AMR 容器，而是“每个修正都要有明确来源、目标和闭合
登记”。HUNDUN 的 final `FaceMassFlux`、wall zero-mass-flux 和 conservation diagnostic
应共享一个 provenance，而不是在 density、transport、PISO 中各自重算。

### 3.2 不复用点

- 不引入 MultiFab、EB2、FabFactory、FluxRegister、ParmParse 或 AMR level；
- 不把 LFP-GCIBM 改成 volume-fraction cut-cell；
- 不实现 redistribution，除非独立 RED 证明存在实际守恒/稳定缺陷；
- 不复制 AMReX kernel、lambda、macro、index layout 或错误文本。

### 3.3 对应任务

- S3-C1/C2/C3：regular/interface/inactive partition 和单一 flux provenance；
- S3-D1/D2：active-cell mass/closure accounting；
- S3-O2：presence-driven provider inventory；
- S3-E1：初始化/每步 exact counters 分离；
- S3-G1：capability ledger 中把 redistribution 明确列为 deferred。

## 4. AMReX-Hydro：projection 阶段职责

公开入口：

- [AMReX-Hydro fixed revision](https://github.com/AMReX-Fluids/AMReX-Hydro/tree/ae9a742c17aa064d299826c60585c1e61b1d6117)
- [AMReX-Hydro documentation](https://amrex-fluids.github.io/amrex-hydro/docs_html/)
- [Projection methods](https://amrex-fluids.github.io/amrex-hydro/docs_html/Projections.html)
- [`Projections/hydro_MacProjector.cpp`](https://github.com/AMReX-Fluids/AMReX-Hydro/blob/ae9a742c17aa064d299826c60585c1e61b1d6117/Projections/hydro_MacProjector.cpp)
- [`Projections/hydro_NodalProjector.cpp`](https://github.com/AMReX-Fluids/AMReX-Hydro/blob/ae9a742c17aa064d299826c60585c1e61b1d6117/Projections/hydro_NodalProjector.cpp)
- [`EBGodunov/hydro_ebgodunov_extrap_vel_to_faces.cpp`](https://github.com/AMReX-Fluids/AMReX-Hydro/blob/ae9a742c17aa064d299826c60585c1e61b1d6117/EBGodunov/hydro_ebgodunov_extrap_vel_to_faces.cpp)
- fixed-revision `Tests/MAC_Projection_EB`、`Tests/Nodal_Projection_EB`、`Tests/Slopes`

### 4.1 可复用点

- projection 是清楚的阶段边界：构造 provisional flux、施加约束、发布 corrected
  flux；
- EB face extrapolation 测试把 regular 与 irregular case 分开，同时使用同一高层
  product path；
- projection tests 独立检查 boundary solvability、nullspace 和 divergence，而不是
  只看最终流场。

HUNDUN 对应规则：两次 PISO corrector 数量不变，pressure wall condition 每次由当前
`rho_wall/D_wall` 更新，最终 transport 只消费第二次 corrector 后的共享
`FaceMassFlux`。

### 4.2 不复用点

- 不采用 MAC/nodal AMReX field layout 或 MLMG；
- 不更换 HUNDUN collocated PISO、Rhie--Chow、matrix-free operators；
- 不从 AMReX-Hydro 引入 limiter、small-cell correction 或 AMR coarse/fine BC。

### 4.3 对应任务

- S3-D1/D2：density prediction 与 PISO 的阶段接口；
- S3-C1--C3：WALE coefficient 在 predictor 前冻结；
- S3-V1：divergence/nullspace/1/2/4-rank formal evidence。

## 5. incflo：variable-density flow composition

公开入口：

- [incflo repository at fixed revision](https://github.com/AMReX-Fluids/incflo/tree/7307d8725c2a538f09cafbeacbfeb63e0fb11d22)
- [incflo documentation](https://amrex-fluids.github.io/incflo/)
- [`src/incflo_advance.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/incflo_advance.cpp)
- [`src/incflo_apply_predictor.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/incflo_apply_predictor.cpp)
- [`src/incflo_apply_corrector.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/incflo_apply_corrector.cpp)
- [`src/incflo_update_density.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/incflo_update_density.cpp)
- [`src/convection/incflo_compute_MAC_projected_velocities.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/convection/incflo_compute_MAC_projected_velocities.cpp)
- [`src/projection/incflo_apply_cc_projection.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/projection/incflo_apply_cc_projection.cpp)
- [`src/diffusion/incflo_diffusion.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/diffusion/incflo_diffusion.cpp)
- `Docs/sphinx_documentation/source/EB.rst`
- `Docs/sphinx_documentation/source/InputsCheckpoint.rst`
- [incflo BSD-3-Clause license](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/LICENSE)

### 5.1 可复用点

incflo 把 advance、predictor、corrector、density update、diffusion 和 projection 分成
粗粒度函数，这支持 HUNDUN 的 `ImmersedDensityAttemptAdapter`：密度模型参与指定阶段，
但不拥有 pressure solver。EB/no-EB branch 在高层选择，底层 regular kernel 不承担
字符串或动态对象查找。

可复用的数据流思想：

```text
known committed/history
-> predict density/transport
-> construct momentum coefficients
-> predictor
-> projection/correction
-> final transport from corrected flux
-> diagnostics/output
```

### 5.2 不复用点

- 不引入 AMR、subcycling、small-cell correction 或 incflo input system；
- 不照搬 incflo class hierarchy、file decomposition、runtime parameters；
- 不把 incflo 的一次 projection 数量映射成 HUNDUN 的 corrector 数量；HUNDUN 保持
  accepted two-corrector PISO；
- 不用 incflo regression result 代替 HUNDUN manufactured/rollback tests。

### 5.3 对应任务

- S3-D1：material density conservative update；
- S3-D2：ideal-gas closure 插入点；
- S3-C2/C3：`rho_attempt` 在 WALE 前冻结；
- S3-A1：same-executable construction order。

## 6. OpenFOAM：WALE 模型责任边界

公开入口：

- [OpenFOAM-dev repository at fixed revision](https://github.com/OpenFOAM/OpenFOAM-dev/tree/627e2f909cb6a08dbb3a74e9a34aa632a975650e)
- [WALE API documentation](https://cpp.openfoam.org/v13/classFoam_1_1LESModels_1_1WALE.html)
- [`WALE.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/627e2f909cb6a08dbb3a74e9a34aa632a975650e/src/MomentumTransportModels/momentumTransportModels/LES/WALE/WALE.C)
- [`WALE.H`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/627e2f909cb6a08dbb3a74e9a34aa632a975650e/src/MomentumTransportModels/momentumTransportModels/LES/WALE/WALE.H)
- fixed-revision `src/MomentumTransportModels/momentumTransportModels/LES/LESeddyViscosity/`
- [OpenFOAM GPL-3.0-or-later license](https://github.com/OpenFOAM/OpenFOAM-dev/blob/627e2f909cb6a08dbb3a74e9a34aa632a975650e/COPYING)

### 6.1 可复用点

- WALE 的输入是 velocity gradient、filter width 和 density/viscosity context；输出是
  eddy viscosity，不拥有 momentum/pressure time integration；
- effective viscosity 是下游 diffusion/stress 消费的单一 coefficient；
- model summary 与 field publication 可以分开。

数学仍以 Nicoud--Ducros 原始定义为 authority：

```text
G = grad(u) grad(u)
Sd = sym(G) - trace(G) I / 3
nu_t = (Cw Delta)^2 (Sd:Sd)^(3/2)
       / ((S:S)^(5/2) + (Sd:Sd)^(5/4))
```

论文：F. Nicoud and F. Ducros, “Subgrid-Scale Stress Modelling Based on the
Square of the Velocity Gradient Tensor”, DOI
[10.1023/A:1009995426001](https://doi.org/10.1023/A:1009995426001)。

### 6.2 不复用点

- OpenFOAM 为 GPL；禁止复制、翻译或机械改写实现；
- 不采用 run-time selection table、objectRegistry、field class 或 wall function；
- 不用固定 epsilon 让静止场产生非零 `nu_t`；
- 不在 corrector 中重复计算 WALE。

### 6.3 对应任务

- 已完成 Task 12/13：tensor core、body-fitted WALE；
- S3-C1--C3：IBM-aware gradient、一次求值、wall `mu_eff`；
- S3-O1/O2：WALE summary/provider；
- S3-V1：TGV、near-wall y-cubed 和 finite-body matrix。

## 7. Basilisk：拓扑、embedded gradient、flux 和 force 对照

公开入口：

- [Basilisk `embed.h`](https://basilisk.fr/src/embed.h)
- [Basilisk `centered.h`](https://basilisk.fr/src/navier-stokes/centered.h)
- [Basilisk `viscosity-embed.h`](https://basilisk.fr/src/viscosity-embed.h)
- [Basilisk `poisson.h`](https://basilisk.fr/src/poisson.h)
- [Basilisk GPL-3.0 license](https://basilisk.fr/src/COPYING)

### 7.1 可复用点

- `face_condition` 表达“插值前先证明拓扑连通”，与 HUNDUN donor/interface link
  authority 同类；
- `dirichlet_gradient` 提供 boundary value 到法向梯度的独立数学 oracle；
- `embed_flux` 展示 homogeneous Neumann 为零 wall flux、Dirichlet gradient 进入
  diffusive flux 的闭合关系；
- `embed_force` 将 pressure 与 viscous traction 分开，可作为二维切片交叉 oracle。

HUNDUN 继续使用自己的 full-strain low-Mach stress、真实三角面 quadrature、Task 11
signed-force 语义和 operator/surface consistency，不能改成 Basilisk 的数据布局或
符号约定。

### 7.2 不复用点

- Basilisk 为 GPL；不复制宏、operator override、fallback、tree/AMR 代码；
- 不引入 volume/face fractions；
- 不采用固定 fallback 或 limiter 掩盖 HUNDUN reconstruction RED；
- 不用 Basilisk force 方向替换 Task 11 已接受方向。

### 7.3 对应任务

- S3-C1：IBM velocity gradient、wall viscosity 和 force 共用 authority；
- S3-D1/D2：zero wall mass flux、density/h one-sided reconstruction；
- S3-C2/C3：same `mu_eff` 的 operator/surface force consistency；
- S3-V1：pressure/viscous/total force 独立验收。

## 8. Nek5000/gslib：全局 donor owner 语义

公开入口：

- [gslib repository at fixed revision](https://github.com/Nek5000/gslib/tree/95acf5b42301d6cb48fda88d662f1d784b863089)
- [`src/findpts.c`](https://github.com/Nek5000/gslib/blob/95acf5b42301d6cb48fda88d662f1d784b863089/src/findpts.c)
- [`src/findpts_local.c`](https://github.com/Nek5000/gslib/blob/95acf5b42301d6cb48fda88d662f1d784b863089/src/findpts_local.c)
- [`src/findpts_el_3.c`](https://github.com/Nek5000/gslib/blob/95acf5b42301d6cb48fda88d662f1d784b863089/src/findpts_el_3.c)
- [`tests/findpts_test.c`](https://github.com/Nek5000/gslib/blob/95acf5b42301d6cb48fda88d662f1d784b863089/tests/findpts_test.c)
- [gslib BSD-3-Clause license](https://github.com/Nek5000/gslib/blob/95acf5b42301d6cb48fda88d662f1d784b863089/LICENSE)

### 8.1 可复用点

findpts 将结果表示为 owner rank、owner-local element 和 interpolation coordinates/
weights。可复用的核心是“donor 是全局对象，halo reach 只是通信实现，不是数学合法性”。
这与 Task 11 已完成的 global donor authority 一致。

### 8.2 不复用点

- 不引入 gslib runtime、crystal router、element maps 或高阶谱元数据；
- 不把每步 global allgather 当作正式 donor exchange；
- 不重做已经接受的 Task 11 donor 实现，除非新 density/WALE RED 证明调用方违反现有
  authority。

### 8.3 对应任务

- S3-C1/D1/D2：只消费现有 global donor plan；
- S3-E1：记录 donor payload/message counters；
- S3-V1：1/2/4-rank shared-row consistency。

## 9. PETSc 与 Trilinos：为什么 Stage 3 不引入

公开入口：

- [PETSc repository at fixed revision](https://github.com/petsc/petsc/tree/9183a15b9d20dbb91d91a365783fab658bad1796)
- [PETSc PCMG manual](https://petsc.org/release/manualpages/PC/PCMG/)
- [Trilinos repository at fixed revision](https://github.com/trilinos/Trilinos/tree/38df4d5b4da2233f9150d3644589803f45942a40)
- Trilinos MueLu/Ifpack2 packages

它们证明成熟 solver stack 需要清楚的 operator、nullspace、preconditioner lifecycle
和 communicator ownership。HUNDUN 可以复用这些接口思想做未来 project-owned
geometric multigrid 设计，但 Stage 3 不增加 PETSc/Trilinos/HYPRE 依赖，不改变当前
matrix-free solver ABI，也不把 solver 性能问题混入科学组合任务。

## 10. Bell--Marcus：variable-density projection 的方程背景

固定引用：J. B. Bell and D. L. Marcus, “A second-order projection method for
variable-density flows”, *Journal of Computational Physics* 101(2), 334--348 (1992),
[DOI 10.1016/0021-9991(92)90011-M](https://doi.org/10.1016/0021-9991(92)90011-M)。

只复用“先得到 trial density/transport，再形成 momentum coefficient，并由 projection
发布 divergence-constrained flux”的阶段思想。HUNDUN 保留自己的 collocated FVM、
Rhie--Chow、两次 PISO、final-flux transport 和 rollback；不复制论文文字、伪代码、
离散式排版或测试数据，也不把该文的一次 projection 解释为改变 corrector 数量。

## 11. Checkpoint、diagnostics 与 transaction 的来源

Checkpoint v3 的 little-endian、binary64、CRC-64/ECMA-182、publish-last 和 failed-read
rollback 是 HUNDUN-owned protocol；上游 checkpoint 目录只能提供测试思想。关键
authority 是已接受 17A fixture 和当前 byte hash，不是 AMReX 或 OpenFOAM 文件格式。

Diagnostics/counters 同样由 HUNDUN schema v1、stable enum 和 read-only provider
合同控制。上游 profiler 只说明初始化成本、每步成本和 wall time 应分开，不能复制
其 artifact schema。

## 12. 任务参考索引

| Task | 必读参考点 | 只复用 | 明确拒绝 |
| --- | --- | --- | --- |
| S3-A0 | Git/DCO/exact-hash governance | immutable candidate activation、single authority switch | product/test edits、implicit approval、scope expansion |
| S3-P0 | Git/CMake/CTest project contracts | isolated worktree、registration ownership | numerical or source-layout refactor |
| S3-C1 | WALE paper；OpenFOAM WALE；Basilisk embed gradient/force；AMReX-Hydro projection | tensor责任、topology guard、同 coefficient 多消费者 | GPL源码、wall function、额外 corrector |
| S3-D1 | incflo density update；AMReX flux accounting；Basilisk zero wall flux | conservative final flux、active mass、one-sided wall rho | cut-cell、redistribution、clipping |
| S3-C2 | S3-C1 + incflo variable density | `rho_attempt` before WALE、same final flux | stale rho、第二 WALE |
| S3-D2 | Bell--Marcus DOI；incflo stage order；AMReX active geometry facts | active-volume p0、closure transaction | species/cp iteration、second density state |
| S3-C3 | S3-C1/C2 + ideal closure | one frozen attempt authority | final-velocity refresh、third solve |
| S3-S1 | WALE paper；smooth TGV self-convergence；HUNDUN Task 11 | cell-average restriction、two-segment observed order、fixed physical time | external result as pass oracle、alternate PDE solver、filtering |
| S3-R1/R2 | accepted 17A protocol | additive profiles、publish-last、rollback | protocol rewrite、cached transient WALE |
| S3-O1/O2 | AMReX support/presence separation | presence inventory、read-only static/per-attempt providers | second numerical calculation、fake absent record |
| S3-A1 | incflo high-level construction separation | optional module construction order | runtime registry/macros、implicit fallback |
| S3-E1 | AMReX initialization/step accounting；HUNDUN artifact v1 | checked exact work counters、compatibility metadata | wall-time threshold、new JSON authority |
| S3-G1/V1 | all above | traceability and independent oracles | fast case labeled formal acceptance |

## 13. Provenance receipt template

每个 semantic-port task receipt 记录：

```text
reference project and fixed revision
exact public URL/path
mathematical or architectural behavior adopted
upstream behavior explicitly rejected
independent HUNDUN RED names
implementation/test diff SHA-256
statement: no source/comment/control-flow/ABI copied
```

仅写“参考了 AMReX/OpenFOAM”不合格；必须指出具体算法责任和拒绝边界。

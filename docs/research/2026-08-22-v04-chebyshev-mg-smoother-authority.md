# v0.4 Chebyshev/多项式 MG smoother 公开方法边界

日期：2026-08-22
状态：**公开方法调查完成；实现合同由 candidate ledger 单独冻结**

## 1. 问题与结论

Re=3900 的 64-rank profile 把下一轮压力架构调查指向 native MG 的重复
point-row/smoother 工作。Chebyshev 或其他多项式 smoother 的吸引力在于其主要
kernel 是 matrix-vector product，避免红黑更新的颜色内依赖，并且其效果不随并行
分区、进程数或未知量排序而改变。公开实现也把 Chebyshev 与对角缩放组合用于 MG
smoothing。

但是，HUNDUN 当前的公共 `NativeCartesianMG` 合同只要求 cell diagonal 为正、face
coefficient 非负；它不要求

```text
diagonal >= x_low + x_high + y_low + y_high + z_low + z_high
```

也没有一项可消费的 SPD、对称界面系数或谱区间证书。因此不能由现有合同推出
`rho(D^-1 A) <= 2`，也不能把 `lambda_max = 2` 固定写入通用 solver。对行 `i`，
Gershgorin 只能给出以 1 为中心、半径
`sum_j(abs(a_ij))/diagonal_i` 的圆盘；当前合法输入允许该半径任意大。公开研究还明确
指出低估最大特征值会使 Chebyshev 多项式在目标区间之外迅速增长。

此外，`MgHierarchyPolicy::{pre_sweeps,post_sweeps}` 是公共策略，当前生产语义是每个
sweep 执行一次有序的两色 point relaxation（强各向异性层则执行 line relaxation）。
把一个 sweep 静默改成某个次数的多项式既改变数学，也改变通信与工作量含义。它必须
作为显式 smoother policy 和 setup identity 进入规范、fingerprint、复用/刷新规则、
证据和 COAST scientific-work receipt，不能伪装成 exact-math 内部优化。

因此本调查的结论是 **REVISE**：多项式路线在原理上适合继续研究，但在冻结实现矩阵
前必须先补齐第 4 节的合同。当前保留的 replicated-coarsest/fused-point candidate 仍是
性能基线；本文件不授权 candidate freeze、COAST pairing 或 full20。

## 2. 固定公开来源与许可边界

| authority | 固定身份与范围 | 本项目只复用的公开思想 | 许可/禁止项 |
|---|---|---|---|
| PETSc | [`petsc/petsc@4f687e3f`](https://gitlab.com/petsc/petsc/-/tree/4f687e3ff0ff5a1ff31f706e200fa1eff85633bb), `KSPCHEBYSHEV`, `PCJACOBI`, `PCGAMGSetUseSAEstEig` 文档 | Chebyshev 需要目标谱区间；MG smoother 通常瞄准高端谱；Jacobi 是 diagonal scaling；已有 setup 的特征值估计可以在 smoothing 生命周期复用 | PETSc license 为 BSD 风格。只采用数学和生命周期结论；不复制源码、递推实现、默认参数控制流、类型或命名布局 |
| hypre | [`hypre-space/hypre@667da829`](https://github.com/hypre-space/hypre/tree/667da8293887d7680964e8ac0310d8e43754f092), `HYPRE_parcsr_ls.h`, PFMG manual | Chebyshev 是公开 relaxation 选项；公开 API 的 order 默认 2、有效 1--4，默认处理高端 30% 谱；pointwise smoother 每 cycle 便宜但可弱于更稳健方法 | Apache-2.0 OR MIT。参数只是研究起点，不是 HUNDUN 默认值；不复制 hypre 源码、对象/API 或控制流 |
| Baker et al. | *Multigrid Smoothers for Ultra-Parallel Computing: Additional Theory and Discussion*, LLNL-TR-489114；[公开 PDF](https://github.com/hypre-space/hypre/wiki/pubs/smoothers-2011-tr.pdf) | polynomial smoother 主要依赖 matvec、对分区/进程数/排序不敏感；Chebyshev 是目标区间上的 min-max polynomial；不得低估 `lambda_max` | 只独立推导数学与测试要求；不复制论文伪代码或公式排版。下载 receipt SHA-256 `35ccf17357962da84ea24b7e2ee5c1d267c37ef340415f133eb35a33edec2df1` |
| Chow et al. | *A Survey of Parallelization Techniques for Multigrid Solvers*；[公开 PDF](https://github.com/hypre-space/hypre/wiki/pubs/parmg-survey.pdf) | 多项式 kernel 是 matvec；最大特征值可用 Gershgorin 或少量 Lanczos 估计；作为 smoother 时低端界可取高端界的一个比例；通信/缓存可使其优于并行 GS | 只采用方法比较和生命周期思想。下载 receipt SHA-256 `f758b38b070abc4f054d1bf403ff0a2a5bac20fca78453499ebd0c48eb08ae79` |

检索当日上游 branch heads 分别是 PETSc
`4f687e3ff0ff5a1ff31f706e200fa1eff85633bb` 和 hypre
`667da8293887d7680964e8ac0310d8e43754f092`。权威引用使用这些不可变 commit，网页
“latest”内容不作为后续实现漂移的许可。

## 3. 对 HUNDUN 当前算子的证明边界

当前七点 level operator 的 owner-local 行为是

```text
(A x)_i = d_i x_i - sum_j(f_ij x_j),  d_i > 0, f_ij >= 0.
```

对 inactive cell，numeric refresh 把 `d_i` 置为 1 并屏蔽连接 face。粗层先限制
`max(0, fine_diagonal - fine_face_sum)` 这一非负的非 face 部分，再加粗层六个 face
coefficient。因此，只要 finest operator 原本对角占优，active masking 和现有粗化会保持
该性质；但 **compile 并不把 finest 对角占优列为公共前置条件**。正 diagonal 和非负
face 也不足以单独证明 SPD。跨 rank 的对应 interface face 数值一致性同样没有由现有
MG certificate 明示。

所以存在三种合规选择：

1. 增加显式、集体一致的 operator-class certificate，并在 numeric refresh 时验证其
   全部前提；
2. 每层从当前 numeric rows 形成保守 Gershgorin 上界并以已登记的 collective 合并，
   同时仍为需要 SPD 的路线提供独立对称/正半定证书；
3. 未获证书的所有输入继续使用现有 red/black 或 line fallback。

“压力生成器在本算例中通常对角占优”不能替代上述任何一项，因为目标要求算例通用。

## 4. 进入实现前必须冻结的合同

### 4.1 公共数学与 policy

- 新增显式 point-smoother kind；默认值和已有构造保持现有 red/black 行为，不能静默
  改变 `pre_sweeps/post_sweeps`。
- 明确定义 polynomial degree 与 sweep 的关系，并进入 policy validation、public/setup
  fingerprint、tests-on/off ABI/符号检查和 COAST scientific-work 描述。
- line-relaxation 层及 replicated coarsest solve 首轮保持不变；多项式只可作用于经证书
  允许的 distributed point levels。
- exact outer operator、FGMRES policy、FP64 terminal true residual、null-space projection、
  failure consensus 和最终接受准则保持不变。

### 4.2 谱界与 numeric lifecycle

- 每个受影响 level 的 `lambda_max` 必须是当前 numeric operator 的安全上界，不得从
  Re=3900 常数或历史 solve 猜测。
- 若用 Gershgorin，界至少为所有 owner rows 的
  `1 + sum(abs(offdiag))/diagonal` 的集体最大值，并检查 finite/positive/overflow；若需要
  SPD/对称性，还必须由独立 certificate 保证，不能由这个标量推断。
- 谱界属于 `ExactNumericState` 或明确的 `PreconditionerSetupState`：系数 refresh 后旧界
  不得继续认证新 operator。任何复用阈值不能跳过当前 exact/coarse numeric refresh。
- 所有额外 collective、payload、stage ordinal 和 failure path 必须在 compile/seal 时
  登记；不得在 hot apply 中临时估谱、分配或创建 communicator。
- 上界无效、证书缺失或适用性失败必须集体、确定地回退到注册好的 red/black 路线，或
  在 seal 前拒绝；不能 rank-local 分叉。

### 4.3 workspace 与 hot schedule

- recurrence 只能复用四个持久 level slots（solution/RHS/residual/temporary）或在
  `make_mg_workspace_requirements` 中登记新增最大容量；不得产生 hot allocation。
- 每个 polynomial stage 的 halo begin/finish、matvec、scale/update 和 revision publish
  必须形成可编号的固定 schedule。pre-smoother 最终 defect retention 与 restriction 的
  现有 seam 必须重新证明，不能额外再算一个未登记的完整 `A*x`。
- reverse post-smoothing 的语义要显式规定；多项式本身无遍历方向，不能让 `reverse`
  形参成为无审计的死策略。

## 5. 冻结前证据矩阵草案

在主 agent 完成接口/数学选择前不得委派实现。选择后至少需要：

- serial row-spectrum oracle：严格占优、仅弱占优、constant-nullspace、inactive mask、
  Dirichlet/Neumann/periodic、非占优但仍为当前公共合法输入，以及 NaN/Inf/overflow；
- 1/2/4-rank 谱界/证书一致性和 rank-interface mismatch/failure 注入；
- polynomial/reference smoother 的独立 residual-reduction authority，覆盖 constant/
  variable coefficient、非均匀 RHS、null space、odd partition、odd-periodic、activity
  mask、各向同性与触发 line fallback 的强各向异性；
- exact/coarse numeric refresh、setup reuse/forced refresh、workspace move/rebind、borrowed
  service failure 的事务测试；
- complete MG/Krylov/update/odd/replicated/line/convergence/hot-allocation focused matrix 和
  UBSan；tests-on/off production binary identity 仍须成立；
- 先做小型 K/R 数值工作量 screen，再做相邻交替 retained-R/new-N full-grid screen。
  只有两个方向都表现出 material total/pressure 改善，且 terminal true residual、迭代、
  collective/halo/IBM/heap/resource 合同均通过，才允许 immutable candidate freeze 和
  COAST reseal。

## 6. 当前决定

本 receipt 自身不冻结实现。主 agent 随后选择了第一条路线：扩展公共 MG policy 和
operator certificate，承担完整的谱界、事务验证、fallback 和 COAST 工作量变化。精确
数学、生命周期、文件边界和验证矩阵以同日 candidate ledger 的
“Certified Chebyshev-Jacobi pressure-smoother experiment freeze”为准。在 ledger 封印
前不启动 worker；无论实现结果如何，本 receipt 都不单独授权正式 pairing、full20 或
长统计。

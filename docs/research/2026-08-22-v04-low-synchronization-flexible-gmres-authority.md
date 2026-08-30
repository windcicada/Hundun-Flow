# v0.4 低同步/通信规避 Flexible GMRES 公开方法证据

日期：2026-08-22
用途：为 v0.4 Cartesian/低速可压缩路线提供公开方法边界与待证明项。
状态：研究 receipt；不构成候选路线选择、数值接受或性能接受。

## 1. 问题边界

本 receipt 只回答以下公开方法是否能为一个 restarted flexible GMRES 提供低同步依据：

- 右预条件、预条件器可变/非线性/不精确；
- 每次 restart 的 Krylov 长度为 (m)（HUNDUN 公开输入允许有界地选择
  (m)，而不是把某个算例中的数值写进算法）；
- 当前基线包含 CGS2 两轮正交化、每若干步的 FP64 true-residual audit，以及 MPI fail-together 语义。

这里的“可用”必须逐项核验，而不能把普通固定预条件 GMRES 的结论外推到 FGMRES：

1. 算法是否明确允许 variable/inexact/nonlinear right preconditioner；
2. 归约数量是每个逻辑 Krylov step、每个 block，还是仅把多个 dot product 合并为一次 MPI 调用；
3. 有限精度下的正交性、残差可信度、归一化/平方根 breakdown 和 restart 语义；
4. 额外向量、三角小矩阵、block basis、pipeline 中间量的存储；
5. 是否保持 Krylov 子空间和预条件器/矩阵应用工作，或者改变了迭代关系和科学工作；
6. 公开数学能否支持 HUNDUN 的独立实现，而不复制 GPL/COAST 旧 Fortran 或其他项目源代码。

本文没有读取 HUNDUN 算例结果、性能结果、candidate 或 COAST receipt；仅核对公开方法材料。

## 2. 权威来源与可核验定位

下表只列原论文、作者公开稿/技术报告、DOI 出版页和官方软件文档/源码。页码是 PDF 页码；官方源码定位使用生成的源码页行号，随 release 页面可核验。

| ID | 来源 | 本 receipt 使用的定位 |
|---|---|---|
| SSM16 | Sanan, Schnepp, May, *Pipelined, Flexible Krylov Subspace Methods*, SIAM J. Sci. Comput. 38 (2016), DOI [10.1137/15M1049130](https://doi.org/10.1137/15M1049130)；作者公开稿 [arXiv:1511.07226](https://arxiv.org/abs/1511.07226)，[PDF](https://arxiv.org/pdf/1511.07226) | §4.2--§4.4，PDF pp.16--19，Algorithms 13--15；§5 的 norm breakdown、storage 与实现说明。2026-08-22 获取的 PDF SHA-256 为 `2a853f6f09ada329b020d65be4803febe62dd06e81063c837d87404e81288f5b` |
| PETSc-PFG | 官方 PETSc [KSPPIPEFGMRES 文档](https://petsc.org/release/manualpages/KSP/KSPPIPEFGMRES/) 与 [pipefgmres.c 源码](https://petsc.org/release/src/ksp/ksp/impls/gmres/pipefgmres/pipefgmres.c.html) | 文档说明 right variable/nonlinear/inexact preconditioning 与 shift；源码 lines 135--143（合并 dot/norm 与非阻塞 reduction）、181--195（square-root breakdown/restart）、211--247（递推、happy breakdown 与 Hessenberg caveat） |
| PETSc-FGMRES | 官方 PETSc [KSPFGMRES 文档](https://petsc.org/release/manualpages/KSP/KSPFGMRES/) 与 [KSP 用户手册](https://petsc.org/release/manual/ksp/) | `KSPFGMRES` 使用的 `-ksp_gmres_restart` 是需要正交化的 Krylov directions 数；手册给出通用 GMRES restart 默认值 30，并说明 FGMRES 默认采用右预条件。这里只采纳公开参数语义与 GMRES-family 通用起点，不采纳或复制实现。 |
| G13 | Ghysels, Ashby, Meerbergen, Vanroose, *Hiding Global Communication Latency in the GMRES Algorithm on Massively Parallel Machines*, SIAM J. Sci. Comput. 35 (2013), DOI [10.1137/12086563X](https://doi.org/10.1137/12086563X) | 出版页摘要：延迟使用 reduction 结果、用非阻塞 reduction 与 operator work 重叠；这是普通 GMRES 的历史管线依据，不作为 flexible 适用性的证明 |
| SLAT21 | Świrydowicz, Langou, Ananthan, Yang, Thomas, *Low synchronization Gram--Schmidt and generalized minimal residual algorithms*, Numer. Linear Algebra Appl. 28 (2021), DOI [10.1002/nla.2343](https://doi.org/10.1002/nla.2343)；作者公开稿 [arXiv:1809.05805](https://arxiv.org/abs/1809.05805)，[PDF](https://arxiv.org/pdf/1809.05805)；DOE 作者稿 [OSTI](https://www.osti.gov/servlets/purl/1710196) | PDF pp.1--5，Algorithms 4、6、7；一同步 MGS 的 mass inner-product/normalization 合并、backward normalization lag、(O(epsilon)kappa) 结论 |
| B16 | Bielich, Langou, Thomas, Świrydowicz, Yamazaki, Boman, *Low-Synch Gram--Schmidt with Delayed Reorthogonalization for Krylov Solvers*, Parallel Computing 112 (2022), DOI [10.1016/j.parco.2022.102940](https://doi.org/10.1016/j.parco.2022.102940)；作者公开稿 [arXiv:2104.01253](https://arxiv.org/abs/2104.01253)，[PDF](https://arxiv.org/pdf/2104.01253)；DOE 作者稿 [OSTI](https://www.osti.gov/servlets/purl/1872059) | PDF pp.1--6，Algorithms 1--4，Tables 1--3，§8；DCGS2 一 reduction/column、延迟 reorthogonalization、有限精度证据与“conjectured”边界 |
| H10 | Hoemmen, *Communication-avoiding Krylov subspace methods*, UC Berkeley EECS Technical Report 2010-37，[作者/学校 PDF](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-37.pdf) | §2.2（preconditioned matrix-powers kernel，PDF pp.79--87）、§2.3（TSQR，pp.88--97）、§3.4（CA-GMRES，pp.143--147）；明确固定 (M) 的 left/right/split 形式与 preconditioner-compatible kernel |
| D12 | Demmel, Grigori, Hoemmen, Langou, *Communication-Optimal Parallel and Sequential QR and LU Factorizations*, SIAM J. Sci. Comput. 34 (2012), DOI [10.1137/080731992](https://doi.org/10.1137/080731992)；作者 PDF [paper84.pdf](https://people.eecs.berkeley.edu/~demmel/Demmel_Pubs_07_11/paper84.pdf) | TSQR 的 parallel reduction tree、通信下界和 Householder 级稳定性 |
| XD24 | Xu, Alonso, Darve, *A Numerically Stable Communication-Avoiding (s)-Step GMRES Algorithm*, SIAM J. Matrix Anal. Appl. (2024), DOI [10.1137/23M1577109](https://doi.org/10.1137/23M1577109)；作者公开稿 [arXiv:2303.08953](https://arxiv.org/abs/2303.08953)，[PDF](https://arxiv.org/pdf/2303.08953) | §2（Algorithms 2.1/2.2，BCGS2/CholQR 通信与稳定性）、§3（partial CholQR、adaptive (s)、scaled Newton）、§5（TSQR/ILU 实验与 O(1) reductions/block） |
| PETSc-L | 官方 [PETSc 2-clause BSD license](https://petsc.org/release/install/license/) | 只用于软件许可证边界；不复制 PETSc 源码 |
| Belos | 官方 [Trilinos/Belos 主页](https://trilinos.github.io/belos.html) 与 [Belos TSQR orthogonalization 文档](https://docs.trilinos.org/r12.16/packages/belos/doc/html/classBelos_1_1GCRODRSolMgr_3_01ScalarType_00_01MV_00_01OP_00_01true_01_4.html) | 证明公开软件栈提供 block/flexible GMRES 接口及 TSQR 选项；不把接口存在当作低同步 FGMRES 稳定性证明 |

## 3. 事实矩阵

### 3.1 直接针对 flexible GMRES 的 single-reduction / pipelined 路线

| 项目 | 公开事实与边界 |
|---|---|
| Flexible 适用性 | SSM16 明确把 FGMRES 与 variable right preconditioning 放在同一推导中；PETSc 文档明确 `KSPPIPEFGMRES` 只支持 right preconditioning，但该 preconditioner 可以 nonlinear、variable、inexact。此项是四类路线中唯一被官方接口直接确认匹配当前 flexible 语义的路线。 |
| 每步 reduction | SSM16 Algorithm 14 是 single-reduction FGMRES；Algorithm 15 是 single-stage PIPEFGMRES。两者每个逻辑 step 各有一个合并的全局 reduction phase；只有 Algorithm 15 将该 phase 与后续 preconditioner/operator 工作重叠。PETSc 源码先把旧 basis 与 shifted `z` 组织到 reduction workspace，再发起 split/nonblocking reduction（源码 lines 135--143）。这不是“零 reduction”，也不等同于删除 true-residual audit。 |
| Algorithm 14 的迭代关系 | Algorithm 14 直接计算 `u_i=B(p_i)` 和 `z_i=A u_i-sigma_i p_i`，没有 Algorithm 15 的近似 `u_tilde_i` 递推。由其算法式可推得：若 Pythagorean 范数恒等式精确、没有 breakdown，shift 只改变 basis 组织而不改变 exact-arithmetic FGMRES Arnoldi 空间；这是根据 Algorithms 13--14 的逐式比较作出的推论，不是有限精度等价声明。 |
| Algorithm 15 的迭代关系 | SSM16 明确指出 PIPEFGMRES 在 exact arithmetic 中也不等同于普通 FGMRES，因为 line 21 的 `u_tilde_i` 是对 `B(p_i)` 的近似；只有在线性 `B` 下该近似才算术等价。其稳定性依赖 shift，不能继承普通 FGMRES 的 bitwise 或残差轨迹。 |
| Breakdown/failure | Algorithms 14--15 以 `t=||z_bar||^2-sum(h_bar^2)` 构造新范数，`t<0` 即 norm/square-root breakdown；论文说明需要 restart，Algorithm 15 还要 refill pipeline，且 breakdown 实际上常在收敛附近。PETSc 源码对小的负舍入量作零截断，仍对实质负值退出当前 cycle 并重启（lines 181--195）。还要处理 happy breakdown 与 FGMRES Hessenberg singularity；PETSc 源码明确指出 FGMRES 不保证 Hessenberg nonsingular（lines 225--247）。多 rank 运行时必须把该事件提升为 fail-together。 |
| 额外存储 | SSM16 §4.4 明确给出 PIPEFGMRES 为 `4m+2` 个向量，而 FGMRES 与 single-reduction CGFGMRES 均为 `2m+2`；`m=12` 时分别为 50 与 26 个向量槽位（公式代入，不是 HUNDUN 测量）。PETSc 源码还分配 preconditioned vectors、shifted `z` vectors、pipeline vectors 与 reduction workspace（lines 31--44、64--70）。 |
| 科学工作/调用图 | Algorithm 14 与基线 FGMRES 都在每个有效 Arnoldi column 直接做一次 `B` 和一次 `A`，主要变化是正交化、归约与有限精度轨迹。Algorithm 15 用 pipeline `B(z_bar)`/`A(q_bar)` 取代稳态的 `B(p)`/`A(u)`，并有额外的启动/flush 状态；不能笼统声称它增加“每一步”调用，也不能声称与当前解轨迹等价。两条路线都必须单独冻结 operator/preconditioner call graph、true-residual audit、breakdown/restart 与终止规则。 |
| 独立实现依据 | 公开论文给出变量预条件器的算法式，PETSc 官方源码给出 reduction overlap 和 failure 处理的可核验实现行为；可作为独立实现的强依据。不能复制 PETSc 代码；PETSc 为 2-clause BSD，但本项目应只复用公开数学/生命周期思想并保留独立实现。 |

### 3.1.1 Restart 长度的公开参数边界

PETSc 官方 `KSPFGMRES` 文档把 restart 定义为要正交化的 Krylov
directions 数；官方 KSP 手册给出通用 GMRES 默认值 30，并说明
FGMRES 默认使用右预条件。由此只能得到一个算例无关的参数语义和可审计
GMRES-family 通用起点，不能推出 restart 30 对 HUNDUN 必然更快、必然更稳定，也不能把
PETSc 的停止准则或实现细节带入 HUNDUN。HUNDUN 若评估 30，仍须保持自身
的 variable Native-MG apply、single-reduction recurrence、每 4 次迭代及终止
时的 FP64 true-residual audit、fail-together、最大迭代数和两次 PISO 修正；
其 workspace、收敛轨迹、通信工作与性能必须由 HUNDUN 自己的冻结矩阵验证。

### 3.2 单 reduction / lagged-normalization / delayed reorthogonalization

| 项目 | 公开事实与边界 |
|---|---|
| SLAT21 一同步 MGS | Algorithm 4 将 mass inner-product 与 normalization 合并进一个 `MPI_AllReduce`，Algorithm 7 将 normalization lag 到下一列；论文摘要和正文给出稳定性目标为 (O(epsilon)kappa([r_0,AV_m]))。它是 one-column MGS-GMRES 的正交化重排，不是 flexible preconditioned GMRES 的推导。 |
| B16 DCGS2 | 普通 CGS2 每列 3 次 global reduction（两次 projection 加一次 normalization）；DCGS2 将第二次 projection 与 normalization 延后，用一次批量 inner-product，达到每个 Arnoldi column 1 次 reduction。论文给出与 CGS2 相同的经验 LOO/representation error，但 Table 3 将 DCGS2-Arnoldi 的 (O(epsilon)) 记为 conjectured，而非一个可直接迁移到 variable-preconditioned FGMRES 的定理。 |
| 有限精度 | SLAT21 的一同步 MGS 结论围绕 compact-WY/lower-triangular correction；B16 的 Pythagorean norm 与 delayed reorthogonalization 用来避免 cancellation。两篇工作都以输入列满秩/正交化假设为前提；B16 明确当 Krylov vectors 失去线性独立时 GMRES 会停滞或失败。不得把“同样 LOO”解释成当前 FGMRES 的 true residual 或 fail-together 保证。 |
| 每步 reduction | 一同步 MGS：每个 column 1 次；DCGS2：每个 column 1 次。归约只是把多个局部 dot/norm 合并，仍需所有 rank 完成一致的 delayed state。对于 restart 边界，首列、末列以及 delayed normalization 尚未消费的状态要单独定义。 |
| 额外工作/存储 | compact-WY correction 需要下三角/三角小矩阵 (T) 或 (L)，并增加 dense triangular solve/matrix work；B16 Table 1/2 对 DCGS2-Arnoldi 记录额外 (j^2) 表示修正项，且 one-reduce 并不自动减少浮点工作。论文没有给出可直接换算为当前 FGMRES 的固定向量槽位数。 |
| Flexible 适用性 | 原文的 Arnoldi relation 是 (A Q=QH)，没有变量 (B_i) 或 (u_i=B_i(p_i)) 的 FGMRES 递推；B16 的实验也不提供 variable/nonlinear/inexact right preconditioner 证据。若要用于 FGMRES，必须重新推导 (A u_i)、Hessenberg 表示、delayed normalization 与 restart/breakdown，而不是替换 CGS2 函数后宣称等价。 |
| 科学工作/迭代关系 | 在固定 operator 的 exact arithmetic 中，它们是正交化/Arnoldi 重排，通常保持同一列空间；有限精度下顺序、三角修正和归一化时刻改变。DCGS2-Arnoldi 还增加表示修正项。它们不能证明当前 variable MG preconditioner 的 Krylov work、预条件器调用数或收敛轨迹不变。 |
| 独立实现依据 | 可作为一个-column 正交化内核的公开数学依据，尤其适合研究 reduction 合并；对当前 flexible 语义仍“不足”，直到完成变量预条件器推导、rank-wide breakdown/restart 规则和真实残差审计证明。 |

### 3.3 (s)-step / communication-avoiding GMRES

| 项目 | 公开事实与边界 |
|---|---|
| 基本结构 | H10 的 CA-GMRES 每个 outer block 用 matrix-powers kernel 产生 (s) 个 Krylov basis vectors，再做 block QR/TSQR/BGS；exact arithmetic 下可与标准 GMRES 产生同一近似解（正确 basis、QR 和 first-vector scaling 条件成立）。XD24 的 Algorithm 2.1 同样先以 MPK 生成 (s+1) 列，再 block-orthogonalize 和更新 Hessenberg。 |
| 每块 reduction | 通信复杂度按 block 计：MPK 将 (s) 次 SpMV 的通信降到约一次 SpMV 级别；block QR 使用 O(1) global reductions per (s) vectors，但具体数由 BCGS2/CholQR/TSQR 选择决定。XD24 明确 BCGS/CholQR 的 reorthogonalization 会把一次 block 的 synchronization 加倍（BCGS2 + 两次 CholQR 语义可达 4 次），其 adaptive 方案目标为 (O(1)) 次 per block。不能把它写成每个 FGMRES step 都只有一次 reduction。 |
| Flexible 适用性 | H10 的 preconditioned CA-GMRES 通过固定 (M) 替换为 (M^{-1}A)、(AM^{-1}) 或 split operator，并要求 matrix-powers kernel 与 preconditioner compatible。公式是一个固定 preconditioned operator 的多项式；没有 (M_i^{-1}) 随迭代变化的 flexible 推导。XD24 的预条件数值例子是固定 ILU/ILUTP 或 equilibration，也没有 flexible/nonlinear/inexact preconditioner 结论。 |
| 有限精度 | 风险来自两层：MPK 产生的 polynomial/Krylov block 可能病态，且 block QR 可能丢失正交性。H10 说明 monomial basis 可能病态，TSQR 与 Householder 同级稳定；XD24 说明 BCGS2/CholQR 仍受 (kappa(V)) 限制，Gram matrix 在有限精度下可能非 SPD。adaptive (s)、scaled Newton 和 partial CholQR 可监视/缩小 block，但其整体 backward-stability proof 在论文中仍被保留为后续工作。 |
| Breakdown/failure | CholQR 的非 SPD/Cholesky breakdown 由 partial CholQR 在发生前保存可用 leading block，并丢弃/另行生成剩余高次列；这会改变 SpMV/MPK 工作。任意 rank deficiency、basis scaling/condition failure、restart 边界仍需显式策略。TSQR 不依赖 Gram matrix 的 Cholesky，但 rank-deficient block 的 happy-breakdown/deflation 语义仍要由 solver 定义。 |
| 额外存储 | 至少需要当前 block 的 (s+1) 个向量、block (Q/R)（或 TSQR tree reflectors）以及 MPK 临时量；(R/H) 的小矩阵随 block/restart 增长。H10 可保持 restart 有界，但 block staging 的峰值存储高于一列一列的 FGMRES。 |
| 科学工作/迭代关系 | 固定 (M) 时，exact arithmetic 的 Krylov subspace/GMRES least-squares 可保持；实际执行把 (s) 次 operator/preconditioner 应用成 block，并改变归一化和舍入顺序。对 variable MG，(A M_i^{-1}) 不是一个固定矩阵的幂，直接套用 MPK 会生成错误的子空间，不能宣称与 FGMRES scientific work 等价。 |
| 独立实现依据 | 可作为固定预条件 CA-GMRES、MPK、block QR 和自适应 block size 的公开方法依据；对当前 variable-preconditioned FGMRES 只提供“需要重新推导”的边界，不足以直接实现或接受。 |

### 3.4 TSQR / block orthogonalization 作为内核

| 项目 | 公开事实与边界 |
|---|---|
| 角色与 Flexible 适用性 | TSQR 是 tall-skinny QR/block orthogonalization kernel，不是一个 flexible GMRES solver。只要 solver 已明确产生一个 block (Z=[z_1,ldots,z_s])，TSQR 可以正交化这组显式向量；这不等于已经解决如何从 variable (B_i) 生成该 block，也不证明 block FGMRES 的 Arnoldi relation。Belos 官方文档只证明其存在 TSQR orthogonalization 选项及 block flexible GMRES 接口，不提供两者组合的稳定性/科学等价证明。 |
| 每块 reduction | D12/H10：对 (P) 个 MPI 进程，parallel TSQR 使用 reduction tree，约 (Theta(log P)) messages；相对逐列 Householder/MGS 少一个 (Theta(s)) 因子，且达到 tall-skinny QR 的通信下界。这里的计数是一个 block QR 的 tree message 复杂度，不是每个 Krylov step 的一次 `MPI_Allreduce`。 |
| 有限精度/正交性 | TSQR 由局部 Householder QR 与树上小 QR 组成，D12/H10 给出与 Householder QR 同级的 (O(epsilon)) 正交性、无 (kappa(A)) 放大项的 bound（对数值满秩 QR 输入）。这比单次 CholeskyQR 的 (O(epsilon)kappa(A)^2) 更稳；但 TSQR 不会自动修复 MPK basis 的病态，也不规定 FGMRES 的 residual/restart 语义。 |
| Breakdown/failure | TSQR 本身避免 Cholesky Gram matrix breakdown，但列相关/数值秩亏仍会表现为小/零 (R_{jj}) 或 happy breakdown。block FGMRES 必须定义 deflation、缩小 block、restart 或 fail-together；这些不是 TSQR paper 的 solver policy。 |
| 额外存储 | 需要本地 block 的 Householder reflectors、树上 (R) 因子以及回代/显式 (Q) 所需工作；若显式 materialize (Q)，峰值 storage 随 block size 增加。它通常替换逐列正交化的 reduction pattern，但不消除 Krylov history/Hessenberg storage。 |
| 科学工作/迭代关系 | 对已经给定的 block，TSQR 的 exact-arithmetic QR 只改变正交化组织；与 block Arnoldi/CA-GMRES 结合时，basis 生成和 least-squares 仍由上层算法决定。它没有公开证据表明可在不改变 variable-preconditioner call graph 的情况下替换当前 FGMRES CGS2。 |
| 独立实现依据 | 可独立实现 TSQR 的公开数学/树布局/稳定性思想；可作为 block orthogonalization candidate 的基础。不能把 TSQR 单独冻结为 flexible GMRES、不能用 TSQR 的 QR 稳定性替代 variable-preconditioner 或 true-residual 证明。 |

## 4. 对当前 flexible/FMGRES 问题的证据分级

这里的分级是“证据是否足以进入独立实现的研究队列”，不是路线优先级，也不是 ACCEPT/REJECT。

### A：公开方法直接覆盖 variable/inexact/nonlinear right-preconditioned FGMRES

- SSM16 的 single-reduction/pipelined FGMRES，加上官方 PETSc `KSPPIPEFGMRES` 文档和源码，是直接覆盖 flexible 语义的最强公开证据。
- 仍然缺少 HUNDUN-specific 证据：shift 的可选取值、variable MG 下的 norm breakdown 频率、与 true-residual audit 的兼容性、fail-together restart、每次 restart 的最终残差和性能工作矩阵。论文的“flexible”不等于当前 two-PISO scientific-work equivalence。

### B：公开正交化内核，可在 flexible solver 中重新推导后使用

- SLAT21 的 one-sync compact-WY MGS 和 B16 的 one-reduce DCGS2 给出 one-column 归约合并/延迟归一化/延迟 reorthogonalization 的公开数学。
- D12/H10 的 TSQR、H10/XD24 的 BGS/BCGS2/CholQR 给出 block orthogonalization 的公开数学和通信模型。
- 这些材料没有把 (u_i=B_i(p_i)) 的 variable preconditioner 带入相应 Arnoldi relation；因此只能作为 kernel/re-derivation 依据，不能直接把固定 (A M^{-1}) 或 (A Q=QH) 代码替换到 FGMRES。

### C：目前不足以作为当前 FGMRES 的直接实现依据

- 任何只证明 fixed (M) 的 (s)-step/CA-GMRES 结果；
- 任何只证明 one-sync MGS/CGS2 的 unpreconditioned 或 fixed-operator Arnoldi 结果；
- 只统计 MPI reduction 次数、没有定义 delayed state、negative norm、rank loss、restart、true-residual 和 fail-together 的实现；
- 只给出 TSQR 正交性、没有给出 variable preconditioner 生成 block 的 Krylov/scientific-work 关系。

## 5. “是否改变 Krylov 迭代/科学工作”的可审计判据

为避免把通信优化误写成科学等价，任何后续独立实现都至少应记录以下不变量/变化：

| 审计项 | 需要回答的问题 |
|---|---|
| Krylov relation | 是 (A U=P H)（variable right FGMRES）还是固定 (A M^{-1}) 的多项式/MPK？若后者，是否有严格理由 (M_i=M)？ |
| operator/preconditioner work | 每个 logical step/block 的 (A)、(B_i)、halo exchange、smoother、coarse solve 次数是否改变？pipeline 是否引入额外 (B/A) 应用？ |
| orthogonalization work | dot product、mass GEMM、triangular solve、TSQR tree、reorthogonalization 是否只是重排，还是新增浮点工作？ |
| residual semantics | 递推残差是否只是监视量？true residual audit 是否仍按原定频率、精度和全局归约执行？何时允许递推值驱动停止？ |
| failure semantics | 任一 rank 检测到 negative norm、零 (R_{jj})、Cholesky non-SPD、Hessenberg singular 或 happy breakdown 时，是否所有 rank 一致 restart/abort？ |
| restart state | delayed normalization/reorthogonalization、pipeline 中间向量和 block 未消费列在 restart 前如何清空/重建？ |
| finite-precision claim | 证明/实验的是正交性、backward error、true residual、最终解，还是仅 runtime？这些指标不能互换。 |

在这些答案冻结前，“one reduction”只能描述通信形状，不能描述与当前 FGMRES 的 scientific-work equivalence。

## 6. 许可与公开方法边界

- 本 receipt 只复述公开论文中的数学关系、算法结构、通信模型和失败条件；不复制 PETSc、Trilinos、Hypre 或 COAST 源代码，不复制 GPL/旧 Fortran 实现。
- PETSc 官方许可证页确认 PETSc 为 2-clause BSD。若仅依据文档/源码行为做独立实现，不会自动引入 PETSc 代码许可证；若未来复制代码，必须另行做文件级版权/许可证审计。
- DOI/出版社页面和作者公开稿的论文文字/图表仍受其版权约束；本文件使用短定位和改写，不替代原文，也不把论文中的 benchmark 数字当作 HUNDUN 证据。
- 本文不为 HUNDUN 授权任何第三方代码、COAST scientific-work equivalence、candidate freeze 或 gate 判断；这些仍由主 agent 按项目交接规则亲自完成。

## 7. 结论（证据，不是路线选择）

1. 对 variable/inexact/nonlinear right preconditioner，公开证据最直接的是 Sanan--Schnepp--May 的 single-reduction FGMRES 与 pipelined flexible GMRES，以及 PETSc `KSPPIPEFGMRES`。前者保留直接 `B(p_i)` 和 `2m+2` 向量但承担 Pythagorean norm breakdown/有限精度风险；后者再引入近似递推、`4m+2` 向量、pipeline refill 和非算术等价。两者都必须保留真实残差与 fail-together 审计，不能合并为同一个候选结论。
2. 一同步 MGS、DCGS2、TSQR、BGS/CA-GMRES 能把 reduction 从列内多次或逐列变成每列一次/每 block (O(1))，但其中绝大多数公开稳定性和通信结论属于 fixed/unpreconditioned Arnoldi、固定 (M) 或独立 QR kernel。它们不能自动外推到 variable MG FGMRES。
3. (s)-step/CA-GMRES 的核心限制不是“有没有 block reduction”，而是 MPK 需要固定 preconditioned operator 的 powers；variable (M_i^{-1}) 破坏这一前提。对当前问题，若研究这一路线，必须先提供独立的 flexible block Arnoldi/FGMRES 推导和 call-graph 证明，再谈性能。
4. 因此，当前可作为 HUNDUN 独立实现依据的是：single-reduction FGMRES 与 PIPEFGMRES 各自的 flexible 算法式，以及 one-sync/DCGS2/TSQR 的独立正交化数学内核。当前尚不足的是：任何未经重新推导和 failure/residual 证据封印的“低同步 FGMRES 等价替换”。

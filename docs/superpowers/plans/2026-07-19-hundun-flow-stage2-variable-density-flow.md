# HUNDUN-FLOW Stage 2 变密度流动实施计划

## 1. 总结与全局约束

**目标：** 在已验收的 C++17/MPI-3 Stage 1 运行时上，实现瞬态低马赫变密度有限体积求解器，依次通过定密度、材料变密度和最小单组分焓反馈三个门。

**架构：** 保留同一个 `hundun` 可执行程序。现有 schema v1、被动标量输出、Restart v1、primitive VTK 和插件 metadata ABI 原样冻结；新增 schema v2 和独立 Stage 2 组合路径。依赖固定为：

```text
application
  -> flow / boundary / finite-volume / density closure
  -> backend-neutral linear algebra
  -> mesh / fields / execution / Halo / MPI
```

执行前先完成文档冻结：

- 将本计划保存到 `docs/superpowers/plans/2026-07-19-hundun-flow-stage2-variable-density-flow.md`。
- 精确提交当前总体设计文档 SHA-256 `0581043ba2b0bfabb222e584219dc99d7910c936c76a70faed931ffbba801516`，不再改写内容。
- 更新 `AGENTS.md`，将 Stage 2 已批准规格和计划加入 required reading，并将实现范围改为 Stage 2；保留全部独立性、私有目录、Python、发布和 push 禁令。
- 从 HEAD `28d539f33bfc4de253a538688e2c979c970ec404` 开始；核对 `stage1-runtime` 仍指向 `934931652e0102e2cac4e9d797a9c0c86f6297f2`。
- 只使用已批准的公开规格和脱敏能力评估；不进入 COAST、BOFFIN 或私有基线源码目录。

每个任务严格执行：

```text
RED test
-> fresh implementation worker
-> 主 agent 检查完整 diff 并重跑测试
-> fresh requirements reviewer
-> 必要时 fresh repair worker 并重新 requirements review
-> fresh code-quality reviewer
-> 必要时 fresh repair worker 并重新 code-quality review
-> 主 agent 验收并终止全部相关 worker
```

同时最多保留 5 个 worker，只允许 1 个 implementation worker 活跃。每个提交带 DCO；不发布、不 push；公开构建和运行不依赖 Python。

明确排除 Stage 3 及以后能力：LES、IBM、化学、物种、多组分热力学、TPDF-TCR、喷雾、颗粒、生产 GPU、PETSc/HYPRE/AMG、WENO/DG、运动壁面、复杂热化学边界、rank-changing Restart、全可压缩、声学和激波。

## 2. 冻结的公共接口与数值合同

### 配置与应用入口

新增：

```cpp
enum class SimulationType { passive_scalar, variable_density_flow };
struct FlowCaseConfig;
using ResolvedCase = std::variant<CaseConfig, FlowCaseConfig>;

ResolvedCase load_resolved_case(const std::filesystem::path&);
std::string to_resolved_json(const ResolvedCase&);
ResolvedCase broadcast_resolved_case(
    MPI_Comm, int root, const ResolvedCase* root_case);
```

现有 `CaseConfig`、`load_case_config()`、`broadcast_case_config()` 及 schema v1 行为保持不变。

Schema v2 固定包含：

- `simulation.type = "variable_density_flow"`；
- `density_model = constant | material | ideal_gas`；
- `mesh.mapping = uniform_box | analytic_warped_box`；
- `time`：固定/自适应、初始/最小/最大 `dt`、CFL、扩散数和重试参数；
- `physics`：`rho_ref`、`mu`，ideal-gas 路径另含常 `cp`、常 `R` 和 `p0`；
- 零个或多个命名通用标量及其扩散率；
- 六个固定逻辑 patch：`x_min/x_max/y_min/y_max/z_min/z_max`；
- `periodic/no_slip_wall/symmetry/velocity_inlet/pressure_outlet`；
- Restart v2、诊断输出和性能记录配置。

每个开放算例最多一个 velocity inlet 和一个 pressure outlet；闭域不含二者。周期边界必须按同轴 min/max 成对。未知键、重复 patch、非法组合和路径逃逸在 rank 0 解析阶段拒绝。

无滑移壁面的焓和通用标量固定使用零法向扩散通量，但不提供“绝热壁/等温壁”公共配置能力；可配置热壁仍属于 Stage 4。

### 字段接口

现有 `FieldView<T>` 的公开方法保持不变，但内部增加存活控制块和 generation 检查。

新增：

```cpp
using ActorId = std::uint32_t;
using PhaseId = std::uint32_t;
enum class AccessMode { read, write, read_write };

class FieldAccessPlan;
template<class T> class FaceFieldView;
template<class T> class KernelCellView;
template<class T> class KernelFaceView;
```

`FieldAccessPlan` 在 freeze 前声明 phase/actor/field/mode；每个 phase/field 最多一个 writer。`FieldStorage` 新增 capability-aware `acquire_read()`、`acquire_write()`，以及支持 cell/face 数量的 `FieldLayoutSet` 构造重载。Stage 2 产品代码不得使用 legacy unrestricted `view()`。

`with_kernel_cell_view()` 和 `with_kernel_face_view()` 是 kernel view 的唯一构造入口。Kernel view 必须 trivially copyable，无控制块、权限、共享所有权、虚调用或逐元素检查。

### 网格接口

新增：

```cpp
class MeshTopology;
class BoundaryPatch;
class MeshGeometry;
class UniformBoxMapping;
class AnalyticWarpedBoxMapping;
```

`MeshTopology` 使用公式化 cell/face 编号，拥有 owner/neighbour、owned/ghost、周期配对、patch membership 和分区无关全局 ID。六个 patch stable ID 固定为 x−/x+/y−/y+/z−/z+ 的 0–5。

解析曲线映射固定为：

```text
X = x0 + Lx [ξ + ax sin(2πξ) sin(πη) sin(πζ)]
Y = y0 + Ly [η + ay sin(πξ) sin(2πη) sin(πζ)]
Z = z0 + Lz [ζ + az sin(πξ) sin(πη) sin(2πζ)]
```

其中 `ξ,η,ζ ∈ [0,1]`，每个 `|a| <= 0.02`。运行时仍须计算并验证正 Jacobian、正 cell volume、共享面互反和 cell closure，不能仅靠幅值范围推断合法。

### 执行与线性代数

新增项目自有类型：

```cpp
enum class ExecutionSpace { host, device };
class ExecutionContext;
class CpuReferenceContext;
class Buffer;
template<class T> class VectorView;
class ExecutionEvent;

ExecutionEvent transfer(
    VectorView<const double>, VectorView<double>, ExecutionContext&);
```

生产只注册 `CpuReferenceContext`。Device context、direct/staged Halo 和 event 路径只由不可注册的 test double 验证；公共头不出现厂商类型。

线性接口固定为：

```cpp
class LinearOperator {
 public:
  virtual VectorLayout domain_layout() const = 0;
  virtual VectorLayout range_layout() const = 0;
  virtual const ExecutionContext& context() const = 0;
  virtual std::uint64_t revision() const = 0;
  virtual ExecutionEvent apply(
      VectorView<const double> x, VectorView<double> y) const = 0;
  virtual bool has_diagonal() const = 0;
  virtual ExecutionEvent diagonal(VectorView<double> d) const = 0;
};

class Preconditioner {
 public:
  virtual void update(const LinearOperator&, std::uint64_t revision) = 0;
  virtual ExecutionEvent apply(
      VectorView<const double> r, VectorView<double> z) const = 0;
};

class LinearSolver {
 public:
  virtual SolveReport solve(
      const LinearOperator&, Preconditioner&,
      VectorView<const double> b, VectorView<double> x,
      const SolveControl&) const = 0;
};
```

`SolveControl` 默认 `atol=1e-12`、`rtol=1e-10`、`max_iterations=500`、独立残差重算周期 20。停止条件固定为：

```text
||b-Ax||₂ <= max(atol, rtol*||b||₂)
```

成功必须由项目自有 FP64 路径独立重算。`SolveReport` 记录结束原因、迭代、初始/递推/最终残差，以及 matvec、preconditioner 和 reduction 计数。

### 时间推进和失败分类

可变步 BDF2 系数固定为：

```text
r = dt_n / dt_(n-1)
alpha0 = (1 + 2r)/(1 + r)
alpha1 = -(1 + r)
alpha2 = r²/(1 + r)
```

首步使用 backward Euler。自适应控制固定为：

- `CFL_target = 0.5`；
- `diffusion_number_target = 0.25`；
- 成功且上一步所有线性求解迭代不超过各自上限 50% 时，增长最多 `1.25`；
- 否则增长因子为 `1.0`；
- BDF2 步长比限制 `[0.5, 2.0]`；
- 可恢复失败令 `dt *= 0.5`；
- 最多重试 8 次，到达 `min_dt` 后再次失败即终止。

可恢复失败：trial 状态非有限/非正、线性不收敛或数值 breakdown、出口回流、最终物理残差超限。配置、layout、capability、已提交输入、MPI 操作和文件完整性错误不可重试。

## 3. 按八个硬门执行的任务

### Gate 1：规格、配置与测量

1. **数值合同和性能 artifact**
   - 增加 `docs/numerics/stage2-contracts.md` 与 `hundun_diagnostics`。
   - 冻结残差、守恒、面方向、曲线映射、测试参数和 JSON artifact 字段。
   - 实现 fake-clock 聚合、raw per-rank sample、median/max-rank 公式和 incomparable 判定。
   - 登记 Thomas–Lombard GCL、Bell–Marcus 变密度方法、Choi/Zhang 时间一致 Rhie–Chow 等公开来源。[Thomas–Lombard](https://doi.org/10.2514/3.61273)、[Bell–Marcus](https://doi.org/10.1016/0021-9991(92)90011-M)、[Choi](https://doi.org/10.1080/104077899274679)、[Zhang–Zhao–Bayyuk](https://doi.org/10.1016/j.jcp.2013.11.006)。

2. **Schema v2 与同一 executable 分派**
   - 新增 `FlowCaseConfig/ResolvedCase` parser、canonical resolved JSON 和 MPI broadcast。
   - 将现有 `main.cpp` 拆成无行为变化的 Stage 1 driver 与新的 dispatch shell。
   - schema v1 的 banner、stdout、Restart、VTK 和 85 项 gate 必须保持原样。
   - 增加 UBSan preset，但不改变现有 preset 语义。

### Gate 2：字段安全

3. **Epoch/liveness**
   - `FieldStorage` 创建控制块和单调 generation。
   - 析构、替换、rebuild/repartition 入口及 Restart v2 读取事务入口先失效旧 view。
   - generation 回绕通过 test seam 明确拒绝。

4. **Capability 与 face layout**
   - 实现 `FieldAccessPlan`、actor/phase、单 writer、多 reader。
   - 增加 `FaceFieldView` 和 cell/face `FieldLayoutSet`；其余 FunctionSpace 继续明确拒绝分配。
   - Stage 1 unrestricted acquisition 仅供既有路径保留。

5. **Kernel views**
   - 增加受限 cell/face kernel callback。
   - 迁移 Stage 2 测试 kernel；复查 Halo、I/O、solver 无长期保存 view。
   - 运行 ASan/UBSan stale checked-view 测试，但禁止解引用 stale/OOB kernel view。

### Gate 3：网格分离

6. **MeshTopology 与 BoundaryPatch**
   - 实现全局 cell/face ID、owner/neighbour、周期配对和六个稳定 patch。
   - 验证 1/2/4 rank 排序后的拓扑一致，分区面不成为物理 patch。

7. **MeshGeometry 与解析曲线映射**
   - 实现 uniform adapter 和 `analytic_warped_box`。
   - 面面积向量由共享 face 顶点顺序一次计算；相邻 cell 使用相反方向。
   - 实现 cell centre/volume、face centre/area、Jacobian、skewness 和 non-orthogonality。
   - Uniform adapter 必须精确复现 Stage 1 spacing、centre、volume、extent 和 owned box。

### Gate 4：执行与线性系统

8. **ExecutionContext/Buffer/View/Event/transfer**
   - 实现同步 `cpu_reference` context、allocation identity、epoch、显式 transfer 和 inline-complete event。
   - Device test double 只测试 capability、lifetime 和拒绝行为。

9. **GhostedVector 与 Buffer Halo**
   - 增加线性 workspace 所需的连续 owned+ghost 布局。
   - 生产路径使用 host-direct；test double 覆盖 runtime-confirmed device-direct、host-staged fallback 和无路径拒绝。
   - 保留现有 `HaloExchange` Stage 1 接口。

10. **VectorOps**
    - 实现 fill/copy/scale/axpy/linear-combination/norm/`dot_batch`。
    - 所有 global reductions 通过 `MpiContext`，并准确计数。

11. **线性接口、Identity 与 Jacobi**
    - 冻结 `LinearOperator/Preconditioner/LinearSolver`。
    - Jacobi 只在 revision/layout 匹配时复用 diagonal；无 diagonal 能力时明确拒绝。

12. **CG**
    - 实现 SPD 检查、zero RHS、breakdown、最大迭代和独立 FP64 残差。
    - 使用连续 workspace，不在迭代中分配 FieldStorage 或临时字段。

13. **BiCGStab 与 collective failure**
    - 覆盖非对称系统、breakdown denominator、非有限输入和任一 rank 失败的一致结果。

14. **Matrix-free Poisson**
    - Uniform/orthogonal路径使用 CG；曲线非正交路径使用 BiCGStab。
    - 非正交面通量分解固定为：

```text
S_orth = (|S|²/(S·d)) d
S_nonorth = S - S_orth
flux = Gamma[(phi_N-phi_P)|S|²/(S·d)
             + grad(phi)_f · S_nonorth]
```

   - `apply()` 固定执行 `Halo begin -> interior -> wait -> partition boundary`。
   - 周期/全 Neumann 做 RHS 投影和零均值规范化；存在 pressure outlet 时禁用重复 nullspace 约束。

### Gate 5：边界与有限体积算子

15. **五类边界**
    - periodic：互反 pairing；
    - no-slip wall：`u_f=0`、零质量通量、压力零法向梯度、h/标量固定零法向扩散通量；
    - symmetry：反射法向速度、复制切向速度；
    - velocity inlet：常速度及密度/热状态/h/标量；
    - pressure outlet：`pi` Dirichlet，其余量纯出流，最终负通量只触发整步拒绝，不截断。

16. **共享面通量有限体积算子**
    - density/h/scalars 使用 MUSCL+MC；动量使用中心形式。
    - 当 face skewness `>0.25` 或 non-orthogonality `>70°` 时，动量对流回退到 MC-limited upwind reconstruction。
    - 实现 Green–Gauss/least-squares gradient、黏性、扩散和共享 `FaceMassFlux`。
    - 所有守恒方程只引用同一个 face-mass-flux field。

### Gate 6：定密度 PISO

17. **Momentum predictor 与时间一致 Rhie–Chow**
    - 三个动量分量使用实际 `a_P`。
    - 面速度由完整离散动量方程推导，并保留 BE/BDF2 历史 face/cell discrepancy；不得用只含当前压力项的简化插值。
    - 单元测试覆盖时间步缩小后 checkerboard 抑制不消失。

18. **固定步长 PISO 事务**
    - 实现 `FlowState` committed/history/trial 分层、`PisoCoupler` 和 `StepAttemptReport`。
    - 每次成功 trial 恰好两次 corrector。
    - 第一次 corrector 后的输运只是 provisional；第二次后从步初状态按最终通量重做 transport finalization。
    - finalization 后若连续性或压力残差失败，只回滚/减小 `dt`，不增加第三次 corrector。
    - 完成 Taylor–Green、pressure checkerboard、nullspace、1/2/4-rank 分解不变性验收。

### Gate 7：变密度与焓反馈

19. **材料密度输运**
    - 保守推进 `rho`、`rho*h` 和全部 `rho*phi`。
    - 密度波验证二阶、正性、质量守恒和 final-flux provenance。

20. **Variable-density vortex**
    - 使用周期域：

```text
psi = sin(x) sin(y)
rho = 1 + 0.1 psi
u = sin(x) cos(y)
v = -cos(x) sin(y)
w = 0
pi = 0
```

   - 测试专用解析体源维持稳态，不进入公共模型配置。
   - 验证质量、三分量动量、共享通量和 MPI 不变性。

21. **最小 ideal-gas closure**
    - `h=T*cp`、`rho=p0/(R*T)`。
    - 闭域：

```text
p0 = M_target*R / sum_cells(V_i/T_i)
```

   - 开放域 `p0` 为正、有限、固定配置量；pressure outlet 只约束 `pi`。
   - inlet 的 h/T 恰有一个权威值，冗余 h/T/rho 只允许按配置容差交叉校验。
   - 覆盖均匀闭域加热、开放 plug-flow、非正 T/rho 和配置矛盾拒绝。

### Gate 8：集成出口

22. **Adaptive BDF2 与 collective retry**
    - 固定步长 BDF2 通过后才启用批准的 CFL controller。
    - 覆盖跨 rank 单点失败、完整 rollback、`dt` 减半、8 次上限、min-dt 终止和 Restart 后下一次 `dt/order` 一致。

23. **Checkpoint v2**
    - 独立格式：`manifest.v2.bin`、每-rank `rank-XXXXXX.v2.bin`、最后发布 `COMPLETED`。
    - 固定 little-endian、binary64、CRC-64/ECMA-182。
    - 持久化 committed n/n−1、dt history、startup/order、controller state、最终 face flux、动态闭域 p0、field schema 和 fingerprints。
    - 只支持相同 rank/process-grid/owned-box；读取事务进入即使旧 checked view 失效，失败不改变字段值或 committed step/time。
    - Restart v1 字节和 API 不修改。

24. **Flow driver 与诊断输出**
    - 同一 `hundun` 按 `ResolvedCase` 分派 Stage 1 或 Stage 2 driver。
    - Stage 2 组合 decomposition、topology、geometry、field access、execution、linear、boundary、closures 和 time driver。
    - 新增 `meshdiag.v2.rank-XXXXXX.bin`，包含版本、全局 IDs、patch、vertices、cell/face geometry、守恒状态和 CRC；提供纯 C++ reader 测试。
    - primitive VTK 保持冻结。

25. **曲线集成与性能 artifact**
    - 在 `analytic_warped_box` 上完成 free-stream、扩散制造解和流动制造解。
    - 接入 allocation/Halo/collective/I/O/SolveReport 计数。
    - Release 手工性能矩阵：ranks `1/2/4`；strong global `64³`；weak 每 rank `32³`；warmup 5 步、measured 20 步、5 次重复。
    - Portable CI 只硬判精确 bytes/message/collective/matvec 等计数；wall-clock、RSS、带宽和吞吐只保存兼容基线，不作硬阈值。

26. **最终 Stage 2 验收**
    - 从 accepted HEAD 全新构建 Debug、Release、ASan、UBSan。
    - 运行全部 Stage 2 unit/MPI/numerical/acceptance tests。
    - 单独重新运行完整 Stage 1 85 项 gate、source-policy、provenance、离线构建和 `ldd` 检查。
    - 核验每项任务的 accepted commits、requirements/code-quality verdict、DCO、无跳过测试、无残留 worker/process、无发布或 push。
    - 通过后停止，不进入 Stage 3。

## 4. 测试与验收阈值

| 验收项 | 固定案例与阈值 |
|---|---|
| Geometry | warped amplitude `[0.02,-0.015,0.01]`；closure `<=256*eps*sum(|S|)`；共享面面积向量互反；所有 Jacobian/volume 正 |
| Linear | 63×63 SPD/非对称 manufactured systems；`atol=1e-12`、`rtol=1e-10`；报告残差与独立重算差 `<=64*eps*max(1,||r||)` |
| Poisson | periodic sine solution，`16³/32³/64³`；L2 收敛阶 `>=1.8`；零均值 `<=1e-12` |
| Taylor–Green | `[0,2π]²` 周期、z 向 4 cells；空间和时间收敛阶 `>=1.8`；相对总质量误差 `<=5e-12` |
| Checkerboard | 初始 parity pressure amplitude 1；两次 corrector 后 parity amplitude `<=1e-8`，规范化 continuity L2 `<=1e-10`，corrector count 必须等于 2 |
| Density wave | `rho=1+0.2 sin(2π(x-t))`，grids `32/64/128`；L1 阶 `>=1.8`；质量误差 `<=5e-12`；rho 保持正 |
| Variable-density vortex | 上述稳态解析场；L2 阶 `>=1.8`；1/2/4-rank max field difference `<=5e-12*max(1,||q||∞)` |
| Ideal-gas | `h/(cp*T)`、`rho*R*T/p0` 相对误差 `<=1e-12`；闭域质量误差 `<=5e-12`；非正状态 collective 拒绝 |
| Final residuals | continuity 规范化 L2 `<=1e-10`；momentum/h/scalars 规范化 L2 `<=1e-9`；全局守恒相对误差 `<=5e-11` |
| Restart v2 | 同分区连续/续算 fields、history、next `dt`、order 和 final flux bitwise 相同；所有损坏/不匹配案例明确拒绝 |
| MPI | 所有适用案例覆盖 1/2/4 ranks；无挂起；失败类别和最低失败 rank 一致 |
| Performance | exact counters 与独立公式完全一致；所有 raw timing 为正且有限；元数据不兼容只标记 incomparable |

最终命令至少包括：

```bash
cmake --preset debug
cmake --build --preset debug -j2
ctest --preset debug --output-on-failure

cmake --preset release
cmake --build --preset release -j2
ctest --preset release --output-on-failure

cmake --preset asan
cmake --build --preset asan -j2
ctest --preset asan --output-on-failure

cmake --preset ubsan
cmake --build --preset ubsan -j2
ctest --preset ubsan --output-on-failure

bash tests/acceptance/stage1_acceptance.sh
bash tests/acceptance/stage2_acceptance.sh
```

## 5. 已选默认与能力声明

- 使用同一个 `hundun`，schema v1 不变，schema v2 显式选择 Stage 2。
- 曲线网格首版只支持 `uniform_box` 和受限 `analytic_warped_box`，不增加坐标文件、AMR 或一般网格导入。
- 自适应时间步采用确定性 CFL/扩散控制，不实现嵌入误差 PI controller。
- wall-clock/RSS/带宽/I/O 吞吐先记录而不硬判；exact counters 是 CI 门禁。
- 无滑移壁面上的 h/标量固定零法向通量不是可配置热壁能力。
- `cpu_optimized` 只有在已接受热点证据后才能另立 Stage 2 内任务；它不是本计划出口条件。
- Device、GPU-aware MPI、mixed precision 和 vendor backend 只冻结可替换接口及 test-double 合同，不宣称已经实现。
- 计划覆盖 Stage 2 低马赫流动、五类基础边界、曲线网格、线性求解、Restart 和性能证据；COAST 黑盒矩阵中的化学、LES、IBM、TPDF-TCR、喷雾和复杂边界继续明确标记为后续阶段，Stage 2 不宣称完全替代旧工作流。

# HUNDUN-FLOW Stage 3 Static LFP-GCIBM and WALE Implementation Plan

> **2026-08-09 compact execution notice:** Historical task definitions and
> scientific contracts below remain evidence, but the active post-Task-11
> order, repository topology, test cadence and milestone ownership are
> superseded by `2026-08-09-hundun-flow-stage3-framework-completion.md` and
> `../specs/2026-08-09-hundun-flow-stage3-compact-scientific-design.md`.

> **2026-08-05 execution notice:** Task numbering and scientific boundaries in
> this plan remain authoritative. Post-Task-11 governance order, test layering,
> resource scheduling and force-consistency gate authority are amended by
> `2026-08-05-hundun-flow-stage3-science-first-execution-amendment.md` and
> `../specs/2026-08-05-hundun-flow-task11-force-consistency-authority-addendum.md`.
> Historical task and review text below is preserved.

> **2026-08-08 no-96 execution amendment:** The active Task 11 and Stage 3
> numerical schedules permanently use `12^3/24^3/48^3`; no command, CTest
> entry or detached runner may launch a `96^3` flow solve. Earlier `96^3`
> wording is retained only as historical specification context and is
> superseded for execution by
> `docs/superpowers/specs/2026-08-08-hundun-flow-task11-no96-extrema-disambiguation.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` to implement this plan
> task-by-task. Only tasks marked **bounded worker eligible** may be delegated;
> all cross-module integration, complex review, and final acceptance remain
> with the main agent. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在已验收的 Stage 2 C++17/MPI-3 低马赫求解器上，实现一个静态闭合
STL 的二阶 Local-Flow-Pattern Ghost-Cell IBM、真实表面受力、WALE LES，以及
constant/material/ideal-gas 三种密度路径下的事务、Checkpoint、诊断和 MPI
验收。

**Architecture:** 新增独立 `hundun_immersed` 和 `hundun_les` 模块，保持
`application -> flow composition -> immersed/LES/FVM/PISO -> mesh/runtime`
单向依赖。Stage 1/2 的公开类型和运行路径不改；schema v3、Stage 3 flow
composition、Checkpoint v3 和 Stage 3 diagnostics 全部采用加法式接口。
LFP-GCIBM 在算子内部执行 row-local replacement-group 变换和一次 simultaneous
Ghost 消元，任何近壁梯度都不能读取 inactive storage slot。

**Tech Stack:** C++17、MPI-3、CMake 3.21+、项目自有 CPU-reference execution、
CG/BiCGStab、PISO、可变步 BDF2、yyjson、纯 C++ STL reader、CRC-64/ECMA-182；
公开构建、运行和测试生成器均不依赖 Python。

## Global Constraints

- 实施基线的规格提交是
  `d682be331cace1db57e2c11b43ab0482062a5928`。
- 已批准 Stage 3 规格是
  `docs/superpowers/specs/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale-design.md`，
  SHA-256 固定为
  `9220bc03da8df13d55f9b8a0a311879af1191e052f8aae2dbd443eb4be3f6de3`。
- Stage 2 accepted product base 是
  `7b9dbcfedfe93918bdf6cf635b71a563103c412b`；任何 Stage 3 重构都必须证明
  schema v1/v2、Restart v1、Checkpoint v2、primitive VTK、banner、stdout 和
  Stage 1/2 数值路径未变。
- 继承
  `.superpowers/sdd/2026-07-21-stage2-coordinator-execution-protocol.md`
  （SHA-256
  `a6b03fe1a2fad24c23597bdfffb074ba197621dea2599143aa16bc151fba8c2c`）
  和
  `.superpowers/sdd/2026-07-22-coordinator-acceptance-acceleration-protocol-v2.md`
  （SHA-256
  `ee3d30242dcb4088fb3a25a3592891ac35f4a370cba62b5bfd5d961588a68f42`）。
- 用户最新执行修订覆盖上述协议中的 reviewer 调度方式：主 agent 亲自完成完整
  task diff 的 requirements review、code-quality review 和 exact-HEAD
  acceptance；不把需要完整上下文的复杂审查派给 worker。只有边界窄、输入短、
  结果可独立验证的实现或 repair 才可交给 fresh worker。
- 同时最多保留 5 个 worker，但默认只保留当前唯一 implementation/repair
  worker；任何时刻最多 1 个 implementation worker 活跃。完成或不用后立即
  结束，不允许 worker 联系用户、扩大范围或请求批准。
- 每个 task 开始前由主 agent 冻结 evidence matrix、allowed-file list、RED、
  comparison semantics、1/2/4-rank matrix、sanitizer matrix、commit subject、
  exact DCO 和 brief SHA-256。复杂 task 的 brief 可以有多个 acceptance
  clusters，但 task 只能整体 accepted/rejected。
- 每个 tracked commit 使用 DCO
  `Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`；不发布、不 push。
- 仅使用批准规格、用户论文和公开文献。不得进入或查看 BOFFIN、COAST、
  COAST-2、固定私有比较基线、私有研究目录或旧实现源码；不得修改、停止、
  清理或打包 `/home/wyf/code_dev/Coast_software` 及其研究数据。
- 优先用 `codegraphf` 做符号、调用方和影响范围导航，用 `rg` 做精确和穷举
  搜索；手工修改使用 `apply_patch`。
- 保持 C++17 和 MPI-3；公共头不得出现 CUDA、HIP、SYCL、PETSc、HYPRE 或
  其他 vendor 类型。Stage 3 生产只注册 CPU reference backend。
- Stage 3 仅支持一个静态、闭合、可定向、单 connected-component STL；
  静止 no-slip、不渗透壁；h/标量零法向扩散通量。
- Stage 3 不实现 cut-cell、moving/multi-part IBM、wall function、热壁、化学、
  物种、TPDF-TCR、喷雾、颗粒、生产 GPU、GPU-aware MPI、AMR、WENO/DG、
  PETSc/HYPRE/AMG、rank-changing Checkpoint、全可压缩、声学或激波。
- 所有 Stage 3 新产品接口必须双向映射：
  `approved requirement -> implementation -> positive test ->
  failure/rollback test -> MPI/numerical acceptance`。

---

## 1. Main-agent-led Execution Protocol

### 1.1 Ownership

下列任务允许使用一个 fresh bounded implementation worker：

```text
Task 2  ImmersedSurface / SurfaceQuery
Task 3  classification / active layouts
Task 4  cell-average moments / deterministic QR
Task 5  GhostStencilPlan / WallQuadraturePlan
Task 6  LocalFlowPatternTransform / replacement groups
Task 7  IBM-aware reconstruction provider
Task 10 WallForceIntegrator
Task 12 WALE core
Task 17 Checkpoint v3 codec/transaction
Task 18 Stage 3 diagnostic adapters
```

这些 worker 只收到该 task 的已冻结 brief、必要接口摘录、allowed files 和测试
命令，不承担整阶段规格审查或跨任务架构判断。

下列复杂任务由主 agent 直接实现、审查和验收：

```text
Task 1  authorization/contracts/schema-v3 compatibility
Task 8  immersed residual and shared-flux integration
Task 9  PISO pressure Ghost and trial transaction
Task 11 full laminar second-order hard gate
Task 13 body-fitted WALE flow composition
Task 14 material-density IBM composition
Task 15 ideal-gas IBM composition
Task 16 IBM+WALE combined gate
Task 19 Stage 3 application driver
Task 20 performance/acceptance/capability ledger
Task 21 final Stage 3 exit
```

若 bounded worker 的 finding 需要跨模块重构，worker 立即停止并返回证据；主
agent 重新冻结范围并亲自处理，不能把长上下文继续堆入 worker prompt。

### 1.2 Per-task sequence

```text
main freezes evidence matrix and brief SHA-256
-> RED is recorded
-> bounded worker implements only when task ownership allows
-> main inspects complete diff from previous accepted task
-> main independently reruns task-focused matrix
-> main requirements review
-> bounded repair worker only for a narrow enumerated closure set
-> main repeats requirements review
-> main code-quality review
-> main exact-HEAD complete Debug
-> focused Release/ASan/UBSan/tests-off as frozen by task matrix
-> main verifies parent/subject/DCO/workers/processes
-> main records accepted HEAD and evidence SHA-256
```

Requirements 和 code-quality 是主 agent 的两个独立 review pass。每一 pass
都必须：

- 从 previous accepted HEAD 到 candidate HEAD 审查完整 task diff；
- 用 codegraphf 穷举 changed symbols、callers、callees、affected tests；
- 用 `rg` 搜索同类实现、断言、test seam、failure classification 和 public
  header；
- 按 Product correctness、Contract violation、Evidence gap、Test quality、
  Nonblocking maintainability 分类，并标 Blocker/Important/Minor；
- 连续两轮同类 finding 时先做 read-only closure sweep，再冻结一次性 repair。

### 1.3 Test evidence identity

每条可复用证据必须记录：

```text
candidate HEAD
preset and build type
compiler / standard library / MPI versions
test names, selectors, rank count and parameters
test-binary SHA-256 where available
exit status, duration and log SHA-256
```

同一角色、同一 exact HEAD、同一 artifact 和同一参数不机械重复。不能复用旧
HEAD、不同 rank、不同 selector 或 repair 前证据。每个 task 最终 candidate
仍必须由主 agent 运行一次 complete Debug。

### 1.4 Bounded-worker handoff packet

派出任何 bounded implementation/repair worker 前，主 agent 生成一份短且闭合
的 task-local packet，固定包含：

```text
previous accepted HEAD and expected candidate parent
task brief path and SHA-256
required-reading files limited to the task
exact allowed files
requirement-to-test evidence matrix
accepted invariants that cannot be reopened
all current findings and objective closure conditions
RED, GREEN, affected regression and rank commands
comparison semantics for each oracle
commit subject
Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>
forbidden adjacent task, Stage 4, private paths, publication and push
```

Packet 不附整份阶段历史。worker 不做架构审批、requirements verdict、
code-quality verdict 或最终 acceptance；worker 返回 commit、完整 diff 摘要、
命令/结果和仍存在的风险后结束。主 agent 核验 worker 已关闭，再进入自己的
完整 diff 双 review。

### 1.5 Stage 3 start gate

Task 1 RED 前，主 agent 只读核验：

```text
current HEAD is the accepted plan commit whose parent chain contains
  d682be331cace1db57e2c11b43ab0482062a5928
approved specification SHA-256 is
  9220bc03da8df13d55f9b8a0a311879af1191e052f8aae2dbd443eb4be3f6de3
working tree has no unapproved tracked change
the unrelated untracked Stage 7 plan remains untouched
no worker or project-owned build/test/MPI process from an earlier task remains
no publication or push target is configured as an execution step
```

记录 plan SHA-256、accepted plan HEAD、compiler/MPI versions 和 start-gate
report SHA-256。若 HEAD/规格/计划不匹配则停止 Task 1；不得通过改写已批准
规格、读取私有目录或清理用户文件来“修复”起点。

---

## 2. File and Target Map

### 2.1 New production targets

`hundun_immersed`：

```text
immersed/include/hundun/immersed/immersed_surface.hpp
immersed/include/hundun/immersed/surface_query.hpp
immersed/include/hundun/immersed/immersed_domain.hpp
immersed/include/hundun/immersed/quadratic_reconstruction.hpp
immersed/include/hundun/immersed/ghost_stencil_plan.hpp
immersed/include/hundun/immersed/local_flow_pattern.hpp
immersed/include/hundun/immersed/wall_force.hpp
immersed/src/immersed_surface.cpp
immersed/src/surface_query.cpp
immersed/src/immersed_domain.cpp
immersed/src/quadratic_reconstruction.cpp
immersed/src/ghost_stencil_plan.cpp
immersed/src/local_flow_pattern.cpp
immersed/src/wall_force.cpp
immersed/src/stl_reader_detail.hpp
immersed/src/surface_bvh_detail.hpp
immersed/src/deterministic_qr_detail.hpp
immersed/src/immersed_test_access.hpp
```

职责：不可变三角表面、查询/分类、active layouts、cell-average 二次权重、
Ghost/Wall plans、LFP coefficient/replacement groups 和真实表面受力。该 target
依赖 `hundun_boundary`、`hundun_mesh`、`hundun_runtime`、MPI 和
`hundun_options`，不依赖 `hundun_flow`。

`hundun_les`：

```text
les/include/hundun/les/wale.hpp
les/src/wale.cpp
les/src/wale_test_access.hpp
```

职责：WALE tensor、attempt-local coefficients/identity 和只读 summary。它
依赖 mesh/runtime/execution，不拥有 communicator 或 `FlowState`。

Stage 3 finite-volume/flow composition：

```text
finite_volume/include/hundun/finite_volume/immersed_reconstruction.hpp
finite_volume/include/hundun/finite_volume/immersed_operator.hpp
finite_volume/src/immersed_reconstruction.cpp
finite_volume/src/immersed_operator.cpp
flow/include/hundun/flow/stage3_flow.hpp
flow/include/hundun/flow/checkpoint_v3.hpp
flow/src/stage3_flow.cpp
flow/src/stage3_piso_detail.hpp
flow/src/checkpoint_v3.cpp
flow/src/checkpoint_v3_detail.hpp
```

职责：只在粗粒度调用边界组合现有 PISO/transport/closure 与新 IBM/WALE。
Stage 2 public create/attempt APIs 保持原样；共享 private helpers 的重构必须由
Stage 2 bitwise/regression tests 证明行为不变。

Stage 3 config/application/diagnostics：

```text
config/include/hundun/config/resolved_case_v3.hpp
config/include/hundun/config/resolved_case_v3_loader.hpp
config/src/resolved_case_v3_loader.cpp
config/src/resolved_case_v3_loader_detail.hpp
applications/hundun/resolved_case_v3_broadcast.hpp
applications/hundun/resolved_case_v3_broadcast.cpp
applications/hundun/stage3_driver.hpp
applications/hundun/stage3_driver.cpp
diagnostics/include/hundun/diagnostics/stage3_module_diagnostics.hpp
diagnostics/include/hundun/diagnostics/checkpoint_v3_diagnostics.hpp
diagnostics/src/stage3_module_diagnostics.cpp
diagnostics/src/checkpoint_v3_diagnostics.cpp
```

### 2.2 Test support

```text
tests/support/stage3_test_contracts.hpp
tests/support/stage3_stl_fixture.hpp
tests/support/stage3_stl_fixture.cpp
tests/support/stage3_mms.hpp
tests/support/stage3_mms.cpp
tests/support/stage3_dual3.hpp
tests/support/stage3_state_equality.hpp
tests/support/stage3_case_generator.cpp
tests/cmake/stage3_source_policy.cmake
tests/cmake/stage3_header_preprocess.cpp
tests/acceptance/stage3_acceptance.sh
```

所有 STL、case JSON 和 corruption fixtures 由这些纯 C++ test tools 在 build
tree 临时目录生成；不提交研究算例、不依赖 Python、不写 source tree。

---

## 3. Frozen Additive Interfaces

### 3.1 Schema v3

```cpp
namespace hundun::config {

enum class ImmersedBoundaryModel : std::uint8_t {
  none,
  local_flow_pattern_ghost_cell
};
enum class ImmersedFluidSide : std::uint8_t { inside, outside };
enum class LesModel : std::uint8_t { none, wale };

struct StlGeometryConfig final {
  std::filesystem::path file;
  double length_scale_to_m{};
  ImmersedFluidSide fluid_side{ImmersedFluidSide::outside};
};
struct StaticImmersedWallConfig final {
  runtime::Real3 velocity_m_per_s{};
};
struct ImmersedBoundaryConfig final {
  ImmersedBoundaryModel model{ImmersedBoundaryModel::none};
  std::optional<StlGeometryConfig> geometry;
  std::optional<StaticImmersedWallConfig> wall;
};
struct WaleConfig final {
  double coefficient{};
  double turbulent_prandtl{};
  double turbulent_schmidt{};
};
struct LesConfig final {
  LesModel model{LesModel::none};
  std::optional<WaleConfig> wale;
};
struct ImmersedFlowCaseConfig final {
  int schema_version{3};
  FlowCaseConfig common_flow;
  ImmersedBoundaryConfig immersed_boundary;
  LesConfig les;
};

using ResolvedCaseV3 =
    std::variant<CaseConfig, FlowCaseConfig, ImmersedFlowCaseConfig>;

ResolvedCaseV3 load_resolved_case_v3(const std::filesystem::path&);
std::string to_resolved_json_v3(const ResolvedCaseV3&);
ResolvedCaseV3 broadcast_resolved_case_v3(
    MPI_Comm, int root, const ResolvedCaseV3* root_case);

}  // namespace hundun::config
```

`ImmersedFlowCaseConfig::common_flow` 保存与 schema v2 相同的已验证 common
values；其内部 `schema_version` 规范化为 2，但 canonical schema-v3 JSON 只输出
顶层 `schema_version:3`。v1/v2 输入必须委托原 loader，并由 byte-for-byte
canonical snapshot 证明未改变。

### 3.2 Surface and domain plan

```cpp
namespace hundun::immersed {

using TriangleId = std::uint64_t;
using ImmersedLinkId = std::uint64_t;

struct SurfaceTriangle final {
  TriangleId id{};
  std::array<runtime::Real3, 3> vertices_m{};
  runtime::Real3 geometric_outward_normal{};
  double area_m2{};
};
struct ClosestPointResult final {
  TriangleId triangle{};
  runtime::Real3 point_m{};
  runtime::Real3 geometric_outward_normal{};
  double squared_distance_m2{};
};
struct SegmentIntersection final {
  TriangleId triangle{};
  runtime::Real3 point_m{};
  double segment_fraction{};
};
enum class CellRegion : std::uint8_t { fluid, solid };

class ImmersedSurface final {
 public:
  static ImmersedSurface load_collective(
      const std::filesystem::path&, double length_scale_to_m,
      const runtime::MpiContext&, int root);
  std::size_t vertex_count() const noexcept;
  std::size_t triangle_count() const noexcept;
  const SurfaceTriangle& triangle(TriangleId) const;
  runtime::Real3 bounding_box_min_m() const noexcept;
  runtime::Real3 bounding_box_max_m() const noexcept;
  double closed_volume_m3() const noexcept;
  std::uint64_t fingerprint() const noexcept;
};

class SurfaceQuery final {
 public:
  static SurfaceQuery create(const ImmersedSurface&);
  ClosestPointResult closest_point(runtime::Real3) const;
  std::vector<SegmentIntersection>
  segment_intersections(runtime::Real3 a_m, runtime::Real3 b_m) const;
  CellRegion classify(runtime::Real3, config::ImmersedFluidSide) const;
  std::uint64_t fingerprint() const noexcept;
};

struct ImmersedLink final {
  ImmersedLinkId id{};
  mesh::GlobalCellId fluid_cell{};
  mesh::GlobalCellId solid_cell{};
  TriangleId triangle{};
  runtime::Real3 wall_intercept_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double fluid_to_wall_fraction{};
};

class ActiveCellLayout final {
 public:
  std::size_t owned_active_count() const noexcept;
  std::size_t local_active_count() const noexcept;
  bool active(mesh::LocalCellId) const;
  std::optional<std::size_t> active_index(mesh::LocalCellId) const;
  const std::vector<mesh::GlobalCellId>& ordered_global_ids() const noexcept;
  std::uint64_t fingerprint() const noexcept;
};

class ActiveBoundaryLayout final {
 public:
  const std::vector<mesh::GlobalFaceId>&
  patch_faces(std::uint32_t stable_patch_id) const;
  bool open_domain() const noexcept;
  bool has_pressure_reference() const noexcept;
  std::uint64_t fingerprint() const noexcept;
};

class ImmersedDomain final {
 public:
  static ImmersedDomain create(
      const ImmersedSurface&, const SurfaceQuery&,
      config::ImmersedFluidSide, const mesh::MeshTopology&,
      const mesh::MeshGeometry&, const boundary::BoundaryRegistry&,
      const runtime::MpiContext&);
  CellRegion region(mesh::LocalCellId) const;
  const std::vector<ImmersedLink>& links() const noexcept;
  const ActiveCellLayout& active_cells() const noexcept;
  const ActiveBoundaryLayout& active_boundaries() const noexcept;
  std::uint64_t classification_fingerprint() const noexcept;
  std::uint64_t surface_coverage_fingerprint() const noexcept;
};

}  // namespace hundun::immersed
```

### 3.3 Reconstruction, Ghost and LFP

```cpp
namespace hundun::immersed {

struct WeightedDonor final {
  mesh::GlobalCellId global_cell{};
  double weight{};
};
struct AffineGhostConstraint final {
  ImmersedLinkId link{};
  std::vector<WeightedDonor> donors;
  double wall_value_weight{};
  double wall_normal_gradient_weight_m{};
};
struct FluidExtrapolation final {
  ImmersedLinkId link{};
  std::vector<WeightedDonor> donors;
};
struct ReconstructionQuality final {
  std::uint32_t rank{};
  double condition_estimate{};
  std::uint32_t halo_reach{};
  std::uint64_t pivot_fingerprint{};
};

class QuadraticReconstruction final {
 public:
  double value(runtime::Real3 point_m,
               const runtime::FieldView<const double>&,
               std::size_t component) const;
  runtime::Real3 gradient(
      runtime::Real3 point_m,
      const runtime::FieldView<const double>&,
      std::size_t component) const;
  const ReconstructionQuality& quality() const noexcept;
};

class GhostStencilPlan final {
 public:
  static GhostStencilPlan create(
      const ImmersedSurface&, const SurfaceQuery&, const ImmersedDomain&,
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const runtime::StructuredDecomposition&, const runtime::MpiContext&);
  const AffineGhostConstraint&
  velocity_constraint(ImmersedLinkId, std::size_t component) const;
  const AffineGhostConstraint& zero_normal_constraint(ImmersedLinkId) const;
  const FluidExtrapolation& density_extrapolation(ImmersedLinkId) const;
  const QuadraticReconstruction& reconstruction(ImmersedLinkId) const;
  std::uint32_t maximum_halo_reach() const noexcept;
  std::uint64_t fingerprint() const noexcept;
};

struct WallQuadraturePoint final {
  TriangleId triangle{};
  std::uint32_t point_index{};
  runtime::Real3 position_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double weight_m2{};
  int owner_rank{};
  QuadraticReconstruction reconstruction;
};
class WallQuadraturePlan final {
 public:
  static WallQuadraturePlan create(
      const ImmersedSurface&, const SurfaceQuery&, const ImmersedDomain&,
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const runtime::MpiContext&);
  const std::vector<WallQuadraturePoint>& local_points() const noexcept;
  std::uint64_t fingerprint() const noexcept;
};

struct LocalCoefficientRow final {
  std::array<double, 6> neighbour{};
  double diagonal{};
  double source{};
};
struct ReplacementGroup final {
  std::uint64_t stable_group_id{};
  std::vector<ImmersedLinkId> links;
  std::vector<std::uint32_t> algebraic_occurrences;
  std::uint64_t evaluator_fingerprint{};
};
struct RowReplacementPlan final {
  mesh::GlobalCellId active_cell{};
  std::vector<ReplacementGroup> groups;
  std::uint64_t fingerprint{};
};

class LocalFlowPatternTransform final {
 public:
  LocalCoefficientRow transform_full(
      const LocalCoefficientRow&, double normal_scale,
      runtime::Real3 solid_to_fluid_normal) const;
  RowReplacementPlan plan_row(
      mesh::GlobalCellId, const std::vector<ImmersedLinkId>&,
      const LocalCoefficientRow&) const;
  double evaluate_wall_replacement(
      const RowReplacementPlan&, const LocalCoefficientRow& immutable_snapshot,
      const std::vector<double>& link_local_symbols) const;
  std::uint64_t algorithm_fingerprint() const noexcept;
};

}  // namespace hundun::immersed
```

`evaluate_wall_replacement()` 必须先构造完整 `W_P`，然后由调用方对所有
`q_G,l` 执行一次 row-wide simultaneous substitution。任何 singleton group
退化为批准规格中的 marginal formula；multi-link group 使用 pure joint
evaluator，不按 link 顺序修改 row。

### 3.4 WALE and Stage 3 flow

```cpp
namespace hundun::les {

enum class WaleTimeOrder : std::uint8_t {
  backward_euler = 1,
  bdf2 = 2
};
struct WaleControl final {
  double coefficient{};
  double turbulent_prandtl{};
  double turbulent_schmidt{};
};
struct WaleCoefficientIdentity final {
  std::uint64_t value{};
};
struct WaleAttemptInput final {
  std::uint64_t step{};
  double attempted_dt_s{};
  WaleTimeOrder order{WaleTimeOrder::backward_euler};
  std::uint64_t committed_state_fingerprint{};
  std::uint64_t history_state_fingerprint{};
  std::uint64_t lagged_gradient_fingerprint{};
  std::uint64_t density_fingerprint{};
  runtime::FieldView<const double> lagged_velocity_gradient;
  runtime::FieldView<const double> rho_attempt;
};
struct WaleSummary final {
  WaleCoefficientIdentity identity{};
  double minimum_nu_t_m2_per_s{};
  double maximum_nu_t_m2_per_s{};
  double l2_nu_t_m2_per_s{};
  std::uint64_t exact_zero_count{};
  std::uint64_t owned_active_count{};
};

class WaleAttemptCoefficients final {
 public:
  ~WaleAttemptCoefficients() noexcept;
  WaleAttemptCoefficients(WaleAttemptCoefficients&&) noexcept;
  WaleAttemptCoefficients& operator=(WaleAttemptCoefficients&&) = delete;
  WaleAttemptCoefficients(const WaleAttemptCoefficients&) = delete;
  WaleAttemptCoefficients& operator=(const WaleAttemptCoefficients&) = delete;

  WaleCoefficientIdentity identity() const noexcept;
  const WaleSummary& summary() const noexcept;
  std::size_t owned_active_count() const noexcept;
  std::size_t local_active_count() const noexcept;
  execution::VectorView<const double> nu_t_m2_per_s() const;
  execution::VectorView<const double> mu_sgs_pa_s() const;
};
class WaleModel final {
 public:
  static WaleModel create(
      WaleControl, const mesh::MeshTopology&, const mesh::MeshGeometry&,
      std::size_t owned_active_count,
      const std::vector<mesh::GlobalCellId>& ordered_local_active_global_cells,
      execution::ExecutionContext&);
  WaleAttemptCoefficients evaluate(const WaleAttemptInput&) const;
};

}  // namespace hundun::les

namespace hundun::immersed {

struct ForceComponents final {
  runtime::Real3 pressure_N{};
  runtime::Real3 viscous_N{};
  runtime::Real3 total_N{};
};
struct MomentComponents final {
  runtime::Real3 pressure_N_m{};
  runtime::Real3 viscous_N_m{};
  runtime::Real3 total_N_m{};
};

}  // namespace hundun::immersed

namespace hundun::flow {

struct ForceAttemptReport final {
  immersed::ForceComponents operator_reaction;
  immersed::ForceComponents surface_traction;
  immersed::ForceComponents consistency;
};
struct Stage3StepAttemptReport final {
  DensityStepAttemptReport base;
  std::optional<ForceAttemptReport> force;
  std::optional<les::WaleCoefficientIdentity> wale_identity;
};
struct Stage3AttemptSummary final {
  TimeAttemptSummary base;
  std::optional<immersed::ForceComponents> force_consistency;
  std::optional<les::WaleCoefficientIdentity> wale_identity;
};
class Stage3TimeAdvanceReport final {
 public:
  TimeAdvanceDisposition disposition() const noexcept;
  StepFailureReason reason() const noexcept;
  int lowest_failing_rank() const noexcept;
  std::size_t attempt_count() const noexcept;
  const Stage3AttemptSummary& attempt(std::size_t) const;
  bool final_attempt_available() const noexcept;
  const Stage3StepAttemptReport& final_attempt() const;
  double accepted_dt_s() const noexcept;
  double proposed_next_dt_s() const noexcept;
};
struct Stage3Physics final {
  config::DensityModel density_model{config::DensityModel::constant};
  double rho_ref_kg_per_m3{};
  double dynamic_viscosity_pa_s{};
  std::optional<double> cp_J_per_kg_K;
  std::optional<double> gas_constant_J_per_kg_K;
  std::optional<double> thermodynamic_pressure_pa;
};

class FixedStepStage3Flow final {
 public:
  static FixedStepStage3Flow create(
      const runtime::StructuredDecomposition&,
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const boundary::BoundaryRegistry&,
      const immersed::ImmersedDomain*,
      const immersed::GhostStencilPlan*,
      const immersed::WallQuadraturePlan*,
      const immersed::LocalFlowPatternTransform*,
      const les::WaleModel*,
      const runtime::MpiContext&, execution::ExecutionContext&,
      runtime::HaloExchange&, const linear::LinearSolver&,
      std::array<linear::Preconditioner*, 3>,
      const linear::LinearSolver&, linear::Preconditioner&);
  Stage3StepAttemptReport attempt(
      FlowState&, const Stage3Physics&, const MomentumTimeStencil&,
      const linear::SolveControl&, const linear::SolveControl&) const;
};

class Stage3RetryController final {
 public:
  static Stage3RetryController create(
      const config::FlowTimeConfig&, const Stage3Physics&,
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const runtime::MpiContext&, const FlowState&);
  Stage3TimeAdvanceReport advance(
      FlowState&, FixedStepStage3Flow&,
      const linear::SolveControl& momentum,
      const linear::SolveControl& pressure);
  TimeControlState state() const noexcept;
};

}  // namespace hundun::flow
```

Stage 3 flow 先按规格中的 BE/BDF2 规则构造 attempt-local lagged velocity，
再通过 body-fitted 或 IBM-aware gradient provider 写入九分量 transient
gradient field；`hundun_les` 只读该 gradient、同一 trial density 和 active
global-cell order，因此不依赖 `FlowState`、IBM 或 finite-volume target。
`WaleAttemptInput` 的 checked views 只在 `evaluate()` 调用期间借用，模型不得
保存。`WaleAttemptCoefficients` 自己拥有 backend-neutral `Buffer`，两个
`VectorView` accessor 按 owned-first local-active order 暴露
`nu_t/mu_sgs`；flow 必须在对象
存活期间将它们插值为本 attempt 冻结的 face `mu_eff`、h diffusivity 和 scalar
diffusivity，不能长期保存 view。
两个 vector 的 size 等于 `local_active_count()`；`WaleSummary` 的 min/max/
norm/zero count 只统计 owned active cells，collective diagnostics 再做显式
rank aggregation，避免 ghost double count。

Material 和 ideal-gas 分支保持批准的 Favre-filtered 语义：
`mu_sgs=rho_attempt*nu_t` 使用同一 trial density，质量、动量、h 和所有标量
继续引用同一最终 `FaceMassFlux`；不得引入第二份速度状态、额外密度闭合或
WALE 对 `p0` 的修改。

`Stage3TimeAdvanceReport` 在 `stage3_flow.hpp` 中用固定容量
`std::array<Stage3AttemptSummary, 9>` 保存最多 9 次 value-only summaries，
并只拥有最终 attempt 的完整 report；不保存 checked/kernel views。失败
attempt 的 `ForceAttemptReport` 和 WALE identity 都不进入 committed state。

---

## 4. Frozen Numerical and Test Contracts

### 4.1 Geometry constants

令

```text
L_ref = max(surface_bbox_diagonal, background_domain_diagonal)
h_i   = cbrt(cell_volume_i)
h_max = collective max over active cells of h_i
eps   = numeric_limits<double>::epsilon()
```

固定：

```text
weld tolerance                    = 128 * eps * L_ref
minimum triangle area             = 1024 * eps * L_ref^2
intersection coincidence bound    = 512 * eps * L_ref
minimum surface/domain separation = 2 * h_max
coverage point/link distance      <= 2 * h_local
fluid/solid witness distance      <= 2 * sqrt(3) * h_local
```

STL reader 支持 standard binary little-endian STL 和严格 ASCII STL。ASCII
数字只接受 C locale finite decimal/scientific tokens；binary record count 必须
与 exact file size 相符。三角法向从顶点重算，文件法向只作诊断。大于
`2^31-1` triangles、size overflow、trailing binary bytes、ASCII trailing
tokens、NaN/Inf、open/nonmanifold/self-intersecting/multi-component surface
全部 collective 拒绝。

Surface quadrature 固定为三点 degree-2 triangle rule：

```text
barycentric points:
  (2/3, 1/6, 1/6)
  (1/6, 2/3, 1/6)
  (1/6, 1/6, 2/3)
weight per point = triangle_area / 3
```

### 4.2 Reconstruction constants

二次 basis 顺序固定为：

```text
1, n, t1, t2, n^2, n*t1, n*t2, t1^2, t1*t2, t2^2
```

donor selection：

```text
minimum donors       = 14
maximum donors       = 32
maximum halo reach   = 4 logical cells
required matrix rank = 10
condition limit      = 1.0e8
normal layers        >= 3 distinct positive-normal coordinate bands
tangential coverage  >= 4 deterministic quadrants
```

候选按 physical distance、global cell ID 排序。使用 deterministic
column-pivoted Householder QR；pivot magnitude tie 在
`64*eps*max(1,max_column_norm)` 内时选择较小 basis index。condition estimate
固定为 accepted `max(abs(R_ii))/min(abs(R_ii))`。uniform/warped cell moments
使用与 `MeshGeometry` 相同的 12-tetra polyhedral representation，二次 monomial
在 tetra 上解析积分；不能把 cell average 当 point value。

权重 reproduction tolerance：

```text
abs(error) <= 512 * eps * max(1, analytic_scale, weight_l1 * data_scale)
```

### 4.3 LFP oracle

standalone paper oracle 使用固定 row：

```text
[A_N,A_S,A_W,A_E,A_L,A_R,A_P,S]
  = [2,3,5,7,11,13,41,17]
k0 = 0.37
k1 = 0.40
k2 = 0.65
k3 = 0.25
```

scale transform 按论文 Eq. 13，三次 coefficient rotations 按 Eqs. 14--16；
test oracle 直接构造 7x7 FP64 matrices 并乘原 row，产品不得调用 oracle。
容差：

```text
abs(product - oracle)
  <= 256 * eps * max(1, row_l1_norm)
```

replacement tests 使用 1、2、3 solid-neighbour rows。group masks 必须对
algebraic occurrence 构成 exact partition；随机置换 group/link 输入后
fingerprint 和 FP64 result bitwise 相同。漏 group、重复 group、复制
cross-link occurrence、sequential Ghost substitution 四种 mutation 必须失败。

### 4.4 Manufactured laminar IBM cases

背景域固定为 `[0,1]^3`，六个 body-fitted patch 均为 no-slip wall，fluid 在
closed body 外部。三种 body：

```text
sphere:
  centre = (0.5,0.5,0.5), radius = 0.18

finite cylinder:
  centre = (0.5,0.5,0.5)
  axis = normalize(1,1,0.5)
  radius = 0.12, length = 0.36, closed planar caps

oblique rectangular prism:
  centre = (0.5,0.5,0.5)
  half lengths = (0.14,0.11,0.09)
  intrinsic rotations = (17 deg,23 deg,31 deg) about x,y,z
```

STL refinement 保证 triangle maximum edge `<=0.45*h_max`；sphere/cylinder
chord error 必须以独立 analytic surface 计算并达到两段 `>=1.8`。prism 面由
exact planar triangles 构成。

令 `L_ref=1 m`、`U_ref=1 m/s`、`rho_ref=1 kg/m^3`，
`xhat=x/L_ref`、`chat=c/L_ref`。对 cylinder 定义
`s=a_hat dot (x-c)`、
`r_perp_squared=|x-c|^2-s^2`；对 prism 定义
`q=transpose(Rz(31 deg)*Ry(23 deg)*Rx(17 deg))*(x-c)`。三个固定的
dimensionless implicit factors 是：

```text
F_sphere =
  (|x-c|^2 - radius^2) / L_ref^2

F_cylinder =
  (r_perp_squared - radius^2)
  * (s^2 - (length/2)^2) / L_ref^4

F_prism =
  (q_x^2 - half_x^2)
  * (q_y^2 - half_y^2)
  * (q_z^2 - half_z^2) / L_ref^6
```

每个 `F_b` 在该 body 的全部真实 faces 上为零；测试只在 fluid domain 评价
该平滑 extension，不用其符号替代产品 `SurfaceQuery` 分类。令：

```text
E(x,y,z) =
  [xhat(1-xhat)yhat(1-yhat)zhat(1-zhat)]^2
asym(x,y,z) = 1 + 0.17*xhat + 0.11*yhat - 0.07*zhat
psi(x,t) =
  U_ref*L_ref
  * [1 + 0.1*sin(2*pi*t*U_ref/L_ref)]
  * E * asym * F_b^2
A = (0.7*psi, -0.4*psi, psi)
u_exact = curl(A)
pi_exact =
  rho_ref*U_ref^2
  * [0.6*xhat + 0.4*yhat - 0.3*zhat
     + 0.15*sin(2*pi*xhat + 0.2)
     + 0.10*sin(2*pi*yhat + 0.4)
     - 0.08*sin(2*pi*zhat + 0.6)]
  * cos(pi*t*U_ref/L_ref)
mu = 0.01 Pa*s
```

`u_exact` 在 outer walls 和 immersed wall 都精确为零且
`div(u_exact)=0`。test-only source：

```text
f = partial_t(u_exact)
    + div(u_exact tensor u_exact)
    + grad(pi_exact)
    - mu * laplacian(u_exact)
```

由 `tests/support/stage3_dual3.hpp` 的独立二阶 forward automatic
differentiation oracle 计算；每个 polyhedral cell 用 tensor 3-point
Gauss-Legendre rule 积分 source cell average。产品路径不包含 source schema。

空间 acceptance 使用 exact BDF2 history：

```text
t_(n-1) = -dt
t_n     = 0
t_(n+1) = dt
dt      = 0.05 * h_max^2 / (U_ref*L_ref)
grids   = 12^3, 24^3, 48^3
```

每个 body 在 uniform 和批准的 warped amplitude
`[0.02,-0.015,0.01]` 上运行。每层 error 必须 finite、strictly positive、
strictly decreasing，且大于
`4096*eps*reference_scale`。两段观测阶均 `>=1.8`：

```text
velocity volume L2 / Linf
pressure gauge-normalized volume L2 / Linf
near-wall-band velocity / pressure L2
wall-normal penetration surface L2 / Linf
pressure force vector error
viscous force vector error
total force vector error
operator/surface pressure consistency
operator/surface viscous consistency
operator/surface total consistency
```

这里的 `near-wall-band` 按控制规格使用固定物理厚度。对批准的
`12^3/24^3/48^3` 序列，在运行前冻结为 coarse 两单元宽度
`2*(L_ref/12)=L_ref/6`；不得用随 refinement 收缩的 `2*h_max` 支持替代
正式误差行。后者只允许作为定位网格相位与单元族的诊断量。

除上述六个 shape/mapping sequences，还要把 sphere centre 固定平移
`(0.013,-0.009,0.007) m`，在 uniform 和同一 warped mapping 上各运行完整
三层序列；解析 `F_b`、source 和 force oracle 使用平移后的 centre。两个
translated sequences 必须通过同一全部误差和两段阶数门，不能只比较单层
结果。

`fluid_side=inside` cavity 使用 centre `(0.5,0.5,0.5)`、radius `0.32 m`
的 sphere，uniform `12^3/24^3/48^3`、相同 exact BDF2/source/oracle 和
`F_sphere`。全部六个 background patches 必须 inactive，pressure operator
使用 RHS projection 和 active-volume zero-mean gauge，不存在 outlet
reference。该序列执行同一 velocity、pressure、penetration、force 和
consistency order gates，并额外要求每层 gauge mean `<=1e-12`。

若三层 penetration 都小于
`8192*eps*max(1,U_ref)`，按 exact-enforcement branch 通过，不计算 `0/0`
order。force reference 使用 analytic body surface 上独立 C++ Gauss
quadrature；节点数加倍后的 reference 变化必须
`<=1e-13*max(1,force_scale)`。

力的归一化尺度在运行前固定为：

```text
A_ref                 = analytic body surface area
pressure_force_scale  = rho_ref * U_ref^2 * A_ref
viscous_force_scale   = mu * U_ref * A_ref / L_ref
total_force_scale     = max(pressure_force_scale, viscous_force_scale)
```

独立 analytic oracle 必须在启动产品求解前证明 pressure、viscous 和 total
reference vector norm 分别大于对应 scale 的 `1e-6`；否则 fixture
`non_discriminating_reference` 失败，不能修改容差或从候选结果反推尺度。

MPI decomposition：

```text
ranks 1: process grid (1,1,1)
ranks 2: process grids (2,1,1) and (1,2,1)
ranks 4: process grids (4,1,1) and (2,2,1)
max field difference
  <= 5e-12 * max(1, global reference Linf)
```

### 4.5 Laminar engineering cases

这些案例验证物理趋势，不替代 §4.4 的形式二阶证明。所有 reference 参数在
RED 前固定；fast 与 acceptance 使用同一产品路径。

低 Reynolds 有限圆柱：

```text
D = 1 m
cylinder centre = (4D,0,0)
cylinder axis = z
cylinder length = 4D
domain = [0,16D] x [-6D,6D] x [-3D,3D]
x_min velocity inlet: (U_ref,0,0)
x_max pressure outlet: pi=0
y/z patches: symmetry
rho = 1 kg/m^3
Re_D = 20 and 40, set mu = rho*U_ref*D/Re_D
fast grid = 96 x 72 x 36, final time = 8D/U_ref
acceptance grid = 192 x 144 x 72, final time = 40D/U_ref
CFL <= 0.25
statistics window = final 10D/U_ref
```

For this case,
`C_D=F_x/(0.5*rho*U_ref^2*D*length)` and
`C_L=F_y/(0.5*rho*U_ref^2*D*length)`. `L_r` is the distance from the
downstream cylinder surface on the `z=0,y=0` centreline to the first
downstream zero of linearly reconstructed `u_x` after a contiguous
negative-velocity interval. Upper/lower separation angles are the first
wall-arc zero of reconstructed tangential traction, measured from the
upstream stagnation direction; ambiguous or multiply alternating zeros are a
failed diagnostic.

在 acceptance window 内，drag 和 recirculation-length 的
coefficient of variation 均 `<=1e-2`；两种 Reynolds 数的 mean drag 都正，
`C_D(Re=20) > C_D(Re=40)`，并且
`L_r(Re=40) >= L_r(Re=20)+0.2D`。上下分离点角度差 `<=5 deg`，
`abs(mean C_L) <= 0.02*mean C_D`。若任一 case 无可解析的闭合 recirculation
region，明确失败，不把零长度当作通过。

低 Reynolds 球：

```text
D = 1 m
sphere centre = (3D,0,0)
domain = [0,12D] x [-3D,3D] x [-3D,3D]
x_min velocity inlet: (U_ref,0,0)
x_max pressure outlet: pi=0
remaining patches: symmetry
rho = 1 kg/m^3
Re_D = 1, mu = rho*U_ref*D
fast grid = 96 x 48 x 48, final time = 8D/U_ref
acceptance grid = 192 x 96 x 96, final time = 24D/U_ref
CFL <= 0.25
statistics window = final 6D/U_ref
```

For the sphere,
`C_D=F_x/(0.5*rho*U_ref^2*pi*D^2/4)`; pressure and viscous fractions use
their separately reduced positive streamwise components divided by total
streamwise force.

令 `C_D,Stokes=24/Re_D`。验收要求 mean total drag 正且
`0.90 <= C_D/C_D,Stokes <= 1.60`，pressure drag fraction 在
`[0.20,0.48]`，viscous drag fraction 在 `[0.52,0.80]`，两者之和与 total
drag 的相对差 `<=5e-11`；两个横向 mean force magnitude 各自
`<=0.02*mean drag`。窗口内 total、pressure、viscous drag 的 coefficient
of variation 均 `<=1e-2`。

封闭 immersed transient：

```text
domain = [0,1]^3, all outer patches no-slip
sphere centre = (0.5,0.5,0.5), radius = 0.18
rho = 1, mu = 0.01
initial velocity = §4.4 u_exact at t=0
body source = zero after initialization
fast grid = 48^3, final time = 0.02 s
acceptance grid = 48^3, final time = 0.05 s
CFL <= 0.25
```

每个 accepted step 的 kinetic energy 必须 finite 且 non-increasing
（允许 `64*eps*max(1,E_previous)` roundoff），active-fluid mass 相对误差
`<=5e-12`，最终 residual 和 wall-penetration 使用 Stage 2/§4.4 阈值，force
components finite，且每步 corrector count 恰为 2。

圆柱和球的 full acceptance 各在 4 ranks 运行一次；相同产品路径的 fast
case 在 1/2/4 ranks 运行并使用 §4.4 decomposition threshold。封闭 transient
在 1/2/4 ranks 运行 full acceptance。工程阈值只声明 finite-body、有限域
baseline 的趋势，不宣称二维无限圆柱、无限域 Stokes 精确值或实验统计精度。

### 4.6 WALE cases

Schema-v3 validation range is frozen as:

```text
1.0e-6 <= Cw   <= 1.0
0.1    <= Pr_t <= 10.0
0.1    <= Sc_t <= 10.0
```

All endpoints are finite; values outside the closed intervals are rank-0
configuration failures. These are input-safety ranges, not claims that every
allowed combination is a validated turbulence calibration.

所有批准案例显式使用：

```text
Cw   = 0.50
Pr_t = 0.90
Sc_t = 0.70
```

unit tensor cases 覆盖 exact zero、rotation invariance、scale homogeneity、
non-finite rejection 和 `y^3`：

```text
y = 2^-k, k = 4..12
2.9 <= log(nu_t(y)/nu_t(y/2))/log(2) <= 3.1
```

body-fitted channel：

```text
domain = [0,2*pi] x [-1,1] x [0,pi]
periodic x/z, no-slip y
rho = 1
nu = 1/180
grid fast       = 32 x 33 x 16, 20 steps
grid acceptance = 64 x 65 x 32, 100 steps
CFL <= 0.25
test-only streamwise pressure-gradient source = 2/180
```

验收：无非有限值；`mu_eff>=mu`；质量相对误差 `<=5e-12`；最终 Stage 2
residual thresholds；上下壁 mean shear 相对差 `<=2e-2`；mean velocity
镜像差 `<=2e-3*U_ref`；第一至第三 cell layers 的 `nu_t/y^3` positive finite；
1/2/4-rank decomposition threshold 同上。不声称 DNS/statistical accuracy。

periodic decaying Taylor--Green：

```text
Re = 1600
domain = [0,2*pi]^3
grid fast = 32^3, grid acceptance = 64^3
dt = 1e-3, final time = 0.1
```

要求 kinetic energy finite and non-increasing、enstrophy finite、
`mu_eff>=mu`、mass/residual/MPI contracts。

IBM+WALE wake：

```text
finite sphere, diameter D=1
domain = [0,12D] x [-3D,3D] x [-3D,3D]
sphere centre = (3D,0,0)
velocity inlet x-min, pressure outlet x-max
remaining patches symmetry
Re_D = 500
grid fast = 96 x 48 x 48, 10 steps
grid acceptance = 192 x 96 x 96, 50 steps
CFL <= 0.25
```

要求 wall penetration、continuity/momentum residual、outlet backflow、
mass、force decomposition、`mu_eff`、rollback 和 1/2/4-rank contracts。该
case 只证明合法 finite-body IBM+WALE composition，不宣称复现论文跨展向
`Re=3900` cylinder。

### 4.7 Equality and failure semantics

- rollback、failed trial、Checkpoint continuation 和要求 bitwise 的 metadata
  使用 FP64 bit patterns；
- nested fields 同时比较 outer size、inner size 和每个 element；
- numerical error 使用上文容差，不允许容器 `operator==` 替代；
- inactive owned cell/face slots 在初始化、trial、rollback、commit、Halo 和
  restore 后都必须是 bitwise `+0.0`；
- test helper 必须证明 exact copy passes、普通字段 mutation fails、nested
  field mutation fails、negative-zero inactive mutation fails；
- rank-local failure collective 汇总相同 classification 和 lowest failing rank；
- config/layout/capability/MPI/file-integrity 失败不可 retry；trial 非正/非有限、
  linear breakdown/non-convergence、backflow 和 final residual 超限可 retry。

---

## Gate 1 — Contracts, Schema and Static Geometry

### Task 1: Freeze Stage 3 contracts and add schema-v3 APIs

**Ownership:** main agent only.

**Files:**

- Modify: `AGENTS.md`
- Create: `docs/numerics/stage3-contracts.md`
- Create: `config/include/hundun/config/resolved_case_v3.hpp`
- Create: `config/include/hundun/config/resolved_case_v3_loader.hpp`
- Create: `config/src/resolved_case_v3_loader.cpp`
- Create: `config/src/resolved_case_v3_loader_detail.hpp`
- Create: `applications/hundun/resolved_case_v3_broadcast.hpp`
- Create: `applications/hundun/resolved_case_v3_broadcast.cpp`
- Create: `tests/unit/test_resolved_case_v3.cpp`
- Create: `tests/mpi/test_resolved_case_v3_broadcast.cpp`
- Create: `tests/unit/test_resolved_case_v3_header_contract.cpp`
- Create: `tests/cmake/stage3_source_policy.cmake`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: frozen `FlowCaseConfig`, v1/v2 loaders, `MpiContext`.
- Produces: exact interfaces in §3.1; Stage 3 numerical constants in
  `stage3-contracts.md`; updated required reading and scope in `AGENTS.md`.

- [ ] **Step 1: Freeze the task evidence matrix**

Record every schema branch:

```text
none/none                         reject
none/wale                         accept
LFP-GCIBM/none                    accept
LFP-GCIBM/wale                    accept
ibm none with geometry/wall       reject exact JSON Pointer
les none with wale/Pr_t/Sc_t      reject exact JSON Pointer
missing/unknown STL format         reject exact JSON Pointer
non-positive/non-finite scale      reject exact JSON Pointer
missing fluid_side                 reject exact JSON Pointer
nonzero h/scalar wall mode         reject exact JSON Pointer
WALE coefficient outside [1e-6,1] reject exact JSON Pointer
Pr_t or Sc_t outside [0.1,10]      reject exact JSON Pointer
nonzero wall velocity             reject
path absolute/.. / escape         reject
v1/v2 through v3 loader           canonical bytes unchanged
v3 through old loader             reject
```

- [ ] **Step 2: Write RED config and header tests**

The header contract must compile:

```cpp
static_assert(std::variant_size_v<hundun::config::ResolvedCase> == 2);
static_assert(std::variant_size_v<hundun::config::ResolvedCaseV3> == 3);
static_assert(
    std::is_same_v<
        decltype(hundun::config::load_resolved_case_v3(
            std::declval<const std::filesystem::path&>())),
        hundun::config::ResolvedCaseV3>);
```

Run:

```bash
cmake --preset debug
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(resolved_case_v3|resolved_case_v3_header_contract)' \
  --output-on-failure
```

Expected RED: missing v3 headers/functions.

- [ ] **Step 3: Implement strict rank-0 parser and canonical JSON**

Use the existing yyjson helpers without changing v1/v2 code paths. Parse all
common values into `common_flow`, set its internal schema marker to 2, and
emit one flat schema-v3 object. Reject all unknown/forbidden keys before
constructing the variant.

- [ ] **Step 4: Implement typed collective broadcast**

Broadcast a versioned byte payload with explicit variant tag and bounded
length chunks no larger than `INT_MAX`. Root-only parse failures must become
one collective failure with the same message and lowest rank.

- [ ] **Step 5: Add compatibility and failure oracles**

Tests must compare v1/v2 canonical strings byte-for-byte against the
Stage 2 loader and cover root/non-root null pointer misuse, tag mismatch,
truncation, rank disagreement and MPI operation failure.

- [ ] **Step 6: Update authorization documents and source policy**

`AGENTS.md` required reading becomes:

```text
Stage 0/1 handoff and plan
overall design
Stage 2 spec and plan
Stage 3 approved spec
Stage 3 approved implementation plan
```

Scope becomes exactly Stage 3 Tasks 1--21. Preserve every private-directory,
no-Python, no-publication and no-push prohibition. Source-policy test rejects
Python dependencies, vendor headers, Stage 4 schema keys and private-path
strings in new production files.

- [ ] **Step 7: Run focused GREEN and MPI**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(resolved_case_v3|resolved_case_v3_header_contract|stage3_source_policy)' \
  --output-on-failure
mpiexec -n 1 build/debug/test_resolved_case_v3_broadcast
mpiexec -n 2 build/debug/test_resolved_case_v3_broadcast
mpiexec -n 4 build/debug/test_resolved_case_v3_broadcast
```

- [ ] **Step 8: Run task-focused sanitizers and tests-off**

```bash
cmake --preset asan
cmake --build --preset asan -j2
ctest --test-dir build/asan -R 'resolved_case_v3' --output-on-failure
cmake --preset ubsan
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R 'resolved_case_v3' --output-on-failure
cmake -S . -B build/tests-off-stage3-task1 \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/tests-off-stage3-task1 -j2
```

- [ ] **Step 9: Main-agent reviews and exact-HEAD acceptance**

Inspect all config/public-header callers with codegraphf; run complete Debug;
verify only v3 APIs are additive and v1/v2 snapshots are unchanged.

- [ ] **Step 10: Commit**

```bash
git add AGENTS.md docs/numerics/stage3-contracts.md \
  config/include/hundun/config/resolved_case_v3.hpp \
  config/include/hundun/config/resolved_case_v3_loader.hpp \
  config/src/resolved_case_v3_loader.cpp \
  config/src/resolved_case_v3_loader_detail.hpp \
  applications/hundun/resolved_case_v3_broadcast.hpp \
  applications/hundun/resolved_case_v3_broadcast.cpp \
  tests/unit/test_resolved_case_v3.cpp \
  tests/unit/test_resolved_case_v3_header_contract.cpp \
  tests/mpi/test_resolved_case_v3_broadcast.cpp \
  tests/cmake/stage3_source_policy.cmake CMakeLists.txt
git commit -s -m "feat: add Stage 3 resolved configuration"
```

### Task 2: Implement deterministic ImmersedSurface and SurfaceQuery

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/immersed_surface.hpp`
- Create: `immersed/include/hundun/immersed/surface_query.hpp`
- Create: `immersed/src/immersed_surface.cpp`
- Create: `immersed/src/surface_query.cpp`
- Create: `immersed/src/stl_reader_detail.hpp`
- Create: `immersed/src/surface_bvh_detail.hpp`
- Create: `immersed/src/immersed_test_access.hpp`
- Create: `tests/support/stage3_stl_fixture.hpp`
- Create: `tests/support/stage3_stl_fixture.cpp`
- Create: `tests/unit/test_immersed_surface.cpp`
- Create: `tests/unit/test_surface_query.cpp`
- Create: `tests/mpi/test_immersed_surface_mpi.cpp`
- Create: `tests/unit/test_immersed_surface_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: `MeshGeometry` units, `MpiContext`, schema-v3 STL path/scale/side.
- Produces: `ImmersedSurface` and `SurfaceQuery` from §3.2; immutable stable
  triangle IDs and fingerprints.

- [ ] **Step 1: Write RED parser/topology/query tests**

Generate strict ASCII and binary tetrahedron/cube fixtures. Add rejection
fixtures for exact-size mismatch, trailing bytes/tokens, NaN, zero area,
open edge, three-triangle edge, duplicate triangle, disconnected component,
inconsistent orientation and self-intersection.

- [ ] **Step 2: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(immersed_surface|surface_query)' --output-on-failure
```

Expected RED: `hundun_immersed` and public headers absent.

- [ ] **Step 3: Implement normalized collective STL load**

Root reads exact bytes, validates arithmetic/size, applies
`length_scale_to_m`, welds only within the frozen tolerance, assigns stable
vertex/triangle IDs, checks closed oriented manifold and one component, then
broadcasts normalized vertices/indices in `INT_MAX` chunks. Every rank rebuilds
the same immutable surface and verifies the same fingerprint.

- [ ] **Step 4: Implement deterministic BVH and queries**

BVH split axis is largest bbox extent with x/y/z tie order; split key is
triangle-centroid coordinate then stable triangle ID. Closest-point ties use
triangle ID. Segment intersections use half-open edge/vertex ownership and
return sorted `(segment_fraction,triangle_id)` records.

- [ ] **Step 5: Implement ray-parity classification**

Use fixed normalized rays:

```text
(1, sqrt(2), sqrt(3))
(sqrt(5), 1, sqrt(7))
(sqrt(11), sqrt(13), 1)
```

All three must agree after half-open intersection handling; disagreement is a
typed geometry error, never random perturbation.

- [ ] **Step 6: Add mutation-sensitive determinism tests**

Permute input triangle order before canonical IDs, alter a coordinate by one
ULP, reverse one face and change MPI partition. Prove intended fingerprint
changes/rejections and identical 1/2/4-rank normalized output.

- [ ] **Step 7: Run GREEN**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(immersed_surface|surface_query)' --output-on-failure
mpiexec -n 1 build/debug/test_immersed_surface_mpi
mpiexec -n 2 build/debug/test_immersed_surface_mpi
mpiexec -n 4 build/debug/test_immersed_surface_mpi
```

- [ ] **Step 8: Run parser sanitizers and tests-off**

```bash
cmake --build --preset asan -j2
ctest --test-dir build/asan -R \
  'test_(immersed_surface|surface_query)' --output-on-failure
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R \
  'test_(immersed_surface|surface_query)' --output-on-failure
cmake -S . -B build/tests-off-stage3-task2 \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/tests-off-stage3-task2 -j2
```

- [ ] **Step 9: Candidate commit**

```bash
git add immersed/include/hundun/immersed/immersed_surface.hpp \
  immersed/include/hundun/immersed/surface_query.hpp \
  immersed/src/immersed_surface.cpp immersed/src/surface_query.cpp \
  immersed/src/stl_reader_detail.hpp immersed/src/surface_bvh_detail.hpp \
  immersed/src/immersed_test_access.hpp \
  tests/support/stage3_stl_fixture.hpp \
  tests/support/stage3_stl_fixture.cpp \
  tests/unit/test_immersed_surface.cpp \
  tests/unit/test_surface_query.cpp \
  tests/unit/test_immersed_surface_header_contract.cpp \
  tests/mpi/test_immersed_surface_mpi.cpp CMakeLists.txt
git commit -s -m "feat: add deterministic immersed surface queries"
```

主 agent 随后完成完整 diff、requirements、quality、complete Debug 和 exact-HEAD
acceptance。

### Task 3: Build classification, ActiveCellLayout and ActiveBoundaryLayout

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/immersed_domain.hpp`
- Create: `immersed/src/immersed_domain.cpp`
- Create: `tests/mpi/test_immersed_domain.cpp`
- Create: `tests/mpi/test_active_boundary_layout.cpp`
- Create: `tests/unit/test_immersed_domain_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Task 2 surface/query, `MeshTopology`, `MeshGeometry`,
  `BoundaryRegistry`.
- Produces: `ImmersedLink`, `ActiveCellLayout`, `ActiveBoundaryLayout`,
  `ImmersedDomain` and SurfaceCoverage fingerprints from §3.2.

- [ ] **Step 1: Write RED classification/link tests**

Cover fluid-side outside and inside, uniform/warped geometry, all-active
identity only when IBM is absent, stable link IDs, exactly one centre-segment
intersection, partition link ownership and sorted global IDs.

- [ ] **Step 2: Write RED boundary/coverage failures**

Include:

```text
zero active inlet                   reject
zero active pressure outlet         reject before pressure operator
one-sided active periodic pair       reject
both periodic sides inactive         accept
inside cavity / six inactive patches use pressure nullspace
small body in one cell               reject
surface with no solid centre         reject
coverage point beyond 2*h_local      reject
missing fluid/solid witness          reject
multi-intersection cell-centre link  reject
```

- [ ] **Step 3: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(immersed_domain|active_boundary_layout)' --output-on-failure
```

- [ ] **Step 4: Implement collective classification and link IDs**

Classify every local owned+ghost cell centre through `SurfaceQuery`. Build
links from topology internal faces with opposite regions. Canonical
`ImmersedLinkId` is the rank-independent sorted index of
`(min(global_cell_ids),max(global_cell_ids),triangle_id,intercept_bits)`.

- [ ] **Step 5: Implement active layouts and validation order**

Construct active cell mapping, then active boundary view, then validate open
patch counts, periodic pairs and pressure reference. Never call
`BoundaryRegistry::open_domain()` directly in Stage 3 flow code after this
task; use `ActiveBoundaryLayout`. `ordered_global_ids()` is owned-first:
owned active cells sorted by global ID, followed by required ghost active
cells sorted by global ID; `owned_active_count()` is the split point and
`local_active_count()` equals the complete vector size.

- [ ] **Step 6: Implement SurfaceCoverage**

Use the frozen three triangle quadrature points, require each point to map to
one or more legal links/active rows within the frozen distances, and count each
triangle contribution once. Record witnesses and max distance in the
fingerprint.

- [ ] **Step 7: Run 1/2/4-rank GREEN and process-grid variants**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_immersed_domain
  mpiexec -n "$n" build/debug/test_active_boundary_layout
done
ctest --test-dir build/debug -R \
  'test_immersed_domain_header_contract' --output-on-failure
```

- [ ] **Step 8: Run sanitizer focus**

ASan/UBSan run header/unit plus 1-rank and 2-rank failure matrices; Release
runs all 1/2/4-rank classification cases.

- [ ] **Step 9: Candidate commit**

```bash
git add immersed/include/hundun/immersed/immersed_domain.hpp \
  immersed/src/immersed_domain.cpp tests/mpi/test_immersed_domain.cpp \
  tests/mpi/test_active_boundary_layout.cpp \
  tests/unit/test_immersed_domain_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: add immersed classification and active layouts"
```

主 agent 验证 Gate 1：schema、surface、query、classification、coverage、
inside/outside、active-boundary 和 1/2/4-rank evidence 全部闭合后才能开始
Task 4。

---

## Gate 2 — Quadratic Ghost Plans and Local Flow Pattern

### Task 4: Implement cell-average moments and deterministic QR

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/quadratic_reconstruction.hpp`
- Create: `immersed/src/quadratic_reconstruction.cpp`
- Create: `immersed/src/deterministic_qr_detail.hpp`
- Create: `tests/unit/test_quadratic_reconstruction.cpp`
- Create: `tests/mpi/test_quadratic_reconstruction_mpi.cpp`
- Create: `tests/unit/test_quadratic_reconstruction_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: active fluid donors, topology/geometry vertices and volumes.
- Produces: `QuadraticReconstruction`, `WeightedDonor`,
  `ReconstructionQuality`; no Ghost boundary semantics yet.

- [ ] **Step 1: Write RED exact polynomial tests**

Use:

```text
q(n,t1,t2) =
  1 + 2*n - 3*t1 + 0.5*t2
  + 0.7*n^2 - 0.2*n*t1 + 0.4*n*t2
  + 0.3*t1^2 - 0.6*t1*t2 + 0.9*t2^2
```

Generate donor **cell averages** analytically from the same 12-tetra cell
representation. Test wall value, normal/tangential gradients and arbitrary
ghost-centre value on uniform and warped cells.

- [ ] **Step 2: Write RED rank/condition/tie tests**

Test rank 9, condition just below/above `1e8`, pivot ties, reordered donors,
14/32 donor limits, missing normal layers, missing tangential quadrants and
halo reach 5 rejection.

- [ ] **Step 3: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_quadratic_reconstruction' --output-on-failure
```

- [ ] **Step 4: Implement analytic polyhedral moments**

For each of the 12 tetrahedra `(reference,face_triangle)`, integrate every
degree-0/1/2 monomial exactly using tetra barycentric moment identities; sum
signed moments and divide by the exact `MeshGeometry::cell_volume_m3`.
Reject any volume/moment mismatch exceeding
`256*eps*max(1,volume_scale)`.

- [ ] **Step 5: Implement deterministic column-pivoted Householder QR**

Use the fixed pivot tie rule and form value/gradient weights by solving
`M^T w = functional`. Do not form normal equations. Reject rank/condition
before publishing any plan.

- [ ] **Step 6: Prove decomposition and ordering identity**

Run 1/2/4 ranks with remote donors and two legal process grids; sorted donor
IDs, pivots, weights and fingerprint must match within the approved FP64
weight tolerance, while final evaluated polynomial values meet §4.2.

- [ ] **Step 7: Run GREEN and sanitizers**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_quadratic_reconstruction' --output-on-failure
mpiexec -n 1 build/debug/test_quadratic_reconstruction_mpi
mpiexec -n 2 build/debug/test_quadratic_reconstruction_mpi
mpiexec -n 4 build/debug/test_quadratic_reconstruction_mpi
cmake --build --preset asan -j2
ctest --test-dir build/asan -R \
  'test_quadratic_reconstruction' --output-on-failure
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R \
  'test_quadratic_reconstruction' --output-on-failure
```

- [ ] **Step 8: Candidate commit**

```bash
git add immersed/include/hundun/immersed/quadratic_reconstruction.hpp \
  immersed/src/quadratic_reconstruction.cpp \
  immersed/src/deterministic_qr_detail.hpp \
  tests/unit/test_quadratic_reconstruction.cpp \
  tests/unit/test_quadratic_reconstruction_header_contract.cpp \
  tests/mpi/test_quadratic_reconstruction_mpi.cpp CMakeLists.txt
git commit -s -m "feat: add quadratic cell-average reconstruction"
```

### Task 5: Build GhostStencilPlan and WallQuadraturePlan

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/ghost_stencil_plan.hpp`
- Create: `immersed/src/ghost_stencil_plan.cpp`
- Create: `tests/mpi/test_ghost_stencil_plan.cpp`
- Create: `tests/mpi/test_wall_quadrature_plan.cpp`
- Create: `tests/unit/test_ghost_stencil_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Task 3 domain/links, Task 4 reconstruction core,
  `StructuredDecomposition`.
- Produces: `GhostStencilPlan`, affine velocity/Neumann constraints,
  material-density fluid extrapolation and `WallQuadraturePlan` from §3.3.

- [ ] **Step 1: Write RED Ghost constraints**

For each link verify:

```text
ghost coordinate = actual solid-neighbour cell centre
velocity q_wall = +0.0
h/scalar g_wall = +0.0
material rho has no wall derivative coefficient
ideal-gas rho derives from h/T closure only
```

Use the Task 4 quadratic polynomial and prove value/gradient/ghost-centre
reproduction.

- [ ] **Step 2: Write RED donor/Halo/failure tests**

Cover remote donors, exact halo reach, rank disagreement, invalid donor side,
rank deficiency, condition limit, no fallback, no per-rank link deletion and
generation of a non-positive material `rho_wall` sample for later retry
classification.

- [ ] **Step 3: Write RED wall quadrature tests**

Verify three points per triangle, sum of weights equals surface area, stable
single-rank triangle ownership, complete coverage, pressure/gradient
reproduction and missing quadrature stencil rejection.

- [ ] **Step 4: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_(ghost_stencil|wall_quadrature)' --output-on-failure
```

- [ ] **Step 5: Implement deterministic donor selection and plans**

Build candidates once during initialization, perform declared structured Halo
reach planning, publish immutable weights/quality/fingerprints only after
collective success. No STL/BVH query or allocation is allowed in later
operator apply.

- [ ] **Step 6: Implement WallQuadraturePlan**

Reuse the exact Task 3 surface quadrature points and Task 4 moment/QR engine.
Assign each triangle to the lowest rank owning an associated active row; tie
by active global cell ID. Every triangle must be integrated once globally.

- [ ] **Step 7: Run GREEN 1/2/4 and sanitizers**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_ghost_stencil_plan
  mpiexec -n "$n" build/debug/test_wall_quadrature_plan
done
cmake --build --preset asan -j2
ctest --test-dir build/asan -R \
  'test_(ghost_stencil|wall_quadrature)' --output-on-failure
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R \
  'test_(ghost_stencil|wall_quadrature)' --output-on-failure
```

- [ ] **Step 8: Candidate commit**

```bash
git add immersed/include/hundun/immersed/ghost_stencil_plan.hpp \
  immersed/src/ghost_stencil_plan.cpp \
  tests/mpi/test_ghost_stencil_plan.cpp \
  tests/mpi/test_wall_quadrature_plan.cpp \
  tests/unit/test_ghost_stencil_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: add ghost and wall quadrature plans"
```

### Task 6: Implement LocalFlowPatternTransform and replacement groups

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/local_flow_pattern.hpp`
- Create: `immersed/src/local_flow_pattern.cpp`
- Create: `tests/unit/test_local_flow_pattern.cpp`
- Create: `tests/mpi/test_local_flow_pattern_mpi.cpp`
- Create: `tests/unit/test_local_flow_pattern_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Task 3 links/geometry, Task 5 Ghost symbols.
- Produces: `LocalFlowPatternTransform`, `ReplacementGroup`,
  `RowReplacementPlan`; no PISO/FVM assembly yet.

- [ ] **Step 1: Write RED paper coefficient-map oracle**

Encode Eqs. 13--16 as independent matrix multiplication using §4.3 row and
`k0..k3`. Tests separately label:

```text
paper full-row coefficient map
HUNDUN singleton marginal replacement
HUNDUN multi-link joint replacement
affine Ghost/quadratic enhancement
```

No test may call product `transform_full()` to generate expected values.

- [ ] **Step 2: Write RED replacement-group tests**

Construct rows with 1, 2 and 3 solid neighbours, including one deferred term
depending on two link symbols. Verify complete/disjoint occurrence partition,
pure joint evaluation, stable sorted accumulation and one simultaneous
substitution map.

- [ ] **Step 3: Write RED mutation tests**

The following altered evaluators must fail:

```text
delete LFP transform
apply full transform to active-active shared contribution
omit group
duplicate group
copy cross-link occurrence into two groups
sequentially substitute Ghost symbols
add WallQuadrature source to row
retain ordinary non-orthogonal P-G term
```

- [ ] **Step 4: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_local_flow_pattern' --output-on-failure
```

- [ ] **Step 5: Implement full transform and marginal/joint evaluators**

`transform_full()` is the standalone paper map. `plan_row()` partitions only
immersed algebraic occurrences, leaves active-active shared coefficients
unchanged, and creates singleton/multi-link groups. Product code evaluates
from one immutable snapshot and never mutates row state while traversing
groups.

- [ ] **Step 6: Prove limiting cases**

Constant and linear fields on uniform planar wall are exact. Product row may
claim equality to the paper whole-row transform only for the tested special
case where `T_LFP(R_shared)=R_shared`; all other reports identify the
masked/marginal HUNDUN extension.

- [ ] **Step 7: Run GREEN, MPI and sanitizers**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_local_flow_pattern' --output-on-failure
mpiexec -n 1 build/debug/test_local_flow_pattern_mpi
mpiexec -n 2 build/debug/test_local_flow_pattern_mpi
mpiexec -n 4 build/debug/test_local_flow_pattern_mpi
cmake --build --preset asan -j2
ctest --test-dir build/asan -R \
  'test_local_flow_pattern' --output-on-failure
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R \
  'test_local_flow_pattern' --output-on-failure
```

- [ ] **Step 8: Candidate commit**

```bash
git add immersed/include/hundun/immersed/local_flow_pattern.hpp \
  immersed/src/local_flow_pattern.cpp \
  tests/unit/test_local_flow_pattern.cpp \
  tests/unit/test_local_flow_pattern_header_contract.cpp \
  tests/mpi/test_local_flow_pattern_mpi.cpp CMakeLists.txt
git commit -s -m "feat: add local flow pattern transform"
```

主 agent 验证 Gate 2：cell-average moments、QR、Ghost、WallQuadrature、paper
oracle、replacement groups、no-fallback 和 MPI plan identity 全部 accepted 后
才能进入有限体积/PISO。

---

## Gate 3 — Constant-density Laminar LFP-GCIBM

### Task 7: Add the unified IBM-aware reconstruction provider

**Ownership:** bounded worker eligible.

**Files:**

- Create: `finite_volume/include/hundun/finite_volume/immersed_reconstruction.hpp`
- Create: `finite_volume/src/immersed_reconstruction.cpp`
- Create: `finite_volume/src/immersed_reconstruction_test_access.hpp`
- Create: `tests/mpi/test_immersed_reconstruction.cpp`
- Create: `tests/unit/test_immersed_reconstruction_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: `ImmersedDomain`, `GhostStencilPlan`, existing
  `CellCenteredFvmOperators`, field access plans and Halo.
- Produces:

```cpp
namespace hundun::finite_volume {

class ImmersedReconstruction final {
 public:
  static ImmersedReconstruction create(
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const boundary::BoundaryRegistry&,
      const immersed::ImmersedDomain&,
      const immersed::GhostStencilPlan&,
      const runtime::MpiContext&, runtime::HaloExchange&);

  void compute_gradient(
      GradientScheme, FiniteVolumeQuantity,
      const runtime::FieldView<const double>& cell_values,
      const runtime::FieldView<double>& cell_gradients) const;
  void reconstruct_transport_faces(
      FiniteVolumeQuantity, const FaceMassFlux&,
      const runtime::FieldView<const double>& cell_values,
      const runtime::FaceFieldView<double>& face_values) const;
  void reconstruct_momentum_faces(
      const FaceMassFlux&,
      const runtime::FieldView<const double>& velocity,
      const runtime::FaceFieldView<double>& face_velocity) const;
  std::uint64_t dependency_fingerprint() const noexcept;
};

}  // namespace hundun::finite_volume
```

- [ ] **Step 1: Write RED inactive-read tracing tests**

Use a test-only read adapter that fails immediately on any owned/ghosted
inactive cell access. Exercise pressure/Rhie--Chow gradients, momentum
convection/viscous/non-orthogonal terms, density/h/scalar MUSCL/diffusion,
WALE velocity gradients and wall-traction reconstruction.

- [ ] **Step 2: Write RED polynomial and shared-face tests**

For all quantities, reproduce Task 4 quadratic values/gradients on uniform and
warped interface bands. An active--active face adjacent to the interface band
must be computed once and added with bitwise opposite signs to owner/neighbour.

- [ ] **Step 3: Add the isolated inactive-payload mutation oracle**

At kernel-fixture level only, run with two different inactive backing payloads
without constructing a product `FlowState`. Active gradients and shared-face
fluxes must be bitwise identical. End-to-end product tests still reject any
non-`+0.0` inactive state.

- [ ] **Step 4: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_immersed_reconstruction' --output-on-failure
```

- [ ] **Step 5: Implement explicit interface-band dependency sets**

Build mutually exclusive active interior, partition boundary and immersed
interface-band row/face sets at creation. Every interface-band derivative
uses active donors plus link-local Ghost symbols; never call Stage 2
full-topology `compute_gradient()` for those rows.

- [ ] **Step 6: Implement fixed apply order**

```text
Halo begin
-> active interior reconstruction
-> Halo wait
-> update remote-donor Ghost symbols
-> partition-boundary reconstruction
-> immersed-interface rows/faces
```

No STL query, BVH traversal, field allocation, plan mutation or virtual
per-element dispatch is permitted during these calls.

- [ ] **Step 7: Run GREEN 1/2/4**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_immersed_reconstruction
done
ctest --test-dir build/debug -R \
  'test_immersed_reconstruction_header_contract' --output-on-failure
```

- [ ] **Step 8: Run focused sanitizers and tests-off**

ASan/UBSan execute the 1-rank and 2-rank interface-band matrix, including
access tracing; tests-off builds `hundun_finite_volume` and a standalone
consumer of the new header.

- [ ] **Step 9: Candidate commit**

```bash
git add finite_volume/include/hundun/finite_volume/immersed_reconstruction.hpp \
  finite_volume/src/immersed_reconstruction.cpp \
  finite_volume/src/immersed_reconstruction_test_access.hpp \
  tests/mpi/test_immersed_reconstruction.cpp \
  tests/unit/test_immersed_reconstruction_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: add IBM-aware finite-volume reconstruction"
```

### Task 8: Assemble the unique immersed residual and shared flux

**Ownership:** main agent only.

**Files:**

- Create: `finite_volume/include/hundun/finite_volume/immersed_operator.hpp`
- Create: `finite_volume/src/immersed_operator.cpp`
- Create: `finite_volume/src/immersed_operator_test_access.hpp`
- Modify: `finite_volume/include/hundun/finite_volume/cell_centered_fvm.hpp`
- Modify: `finite_volume/src/cell_centered_fvm.cpp`
- Create: `tests/mpi/test_immersed_operator.cpp`
- Create: `tests/unit/test_immersed_operator_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Tasks 6/7, one existing `FaceMassFlux` field, active layout.
- Produces:

```cpp
namespace hundun::finite_volume {

struct ImmersedResidualParts final {
  std::array<double, 3> pressure{};
  std::array<double, 3> viscous{};
  std::array<double, 3> convective{};
};
struct ImmersedOperatorReport final {
  std::uint64_t active_row_count{};
  std::uint64_t replacement_group_count{};
  std::uint64_t simultaneous_substitution_count{};
  ImmersedResidualParts operator_reaction_N{};
};

class ImmersedOperatorAdapter final {
 public:
  static ImmersedOperatorAdapter create(
      const mesh::MeshTopology&, const mesh::MeshGeometry&,
      const immersed::ImmersedDomain&,
      const immersed::GhostStencilPlan&,
      const immersed::LocalFlowPatternTransform&,
      const ImmersedReconstruction&);
  void accumulate_momentum(
      const runtime::FieldView<const double>& velocity,
      const runtime::FieldView<const double>& pressure,
      const runtime::FieldView<const double>& velocity_gradient,
      const runtime::FaceFieldView<const double>& dynamic_viscosity_by_face,
      const runtime::FieldView<double>& residual) const;
  void accumulate_transport(
      FiniteVolumeQuantity,
      const runtime::FieldView<const double>& values,
      const runtime::FieldView<const double>& gradients,
      const runtime::FaceFieldView<const double>& gamma_by_face,
      const runtime::FieldView<double>& residual) const;
  ImmersedOperatorReport report() const;
};

}  // namespace hundun::finite_volume
```

- [ ] **Step 1: Freeze the row-algebra RED fixture**

For each active row:

```text
R_final(P)
  = R_bg_full(P)
  - R_immersed_bg(P)
  + simultaneous_substitute(W_P, Q_G(P))
```

The fixture must expose direct, non-orthogonal, deferred and joint multi-link
occurrences separately and compare against an independent direct quadratic
full-row evaluator.

- [ ] **Step 2: Write RED conservation and no-double-count tests**

Verify one shared active--active face record, canonical-zero fluid--solid face
mass-flux slot, exact-zero wall mass flux, exact-zero h/scalar diffusive wall
flux, one immersed momentum contribution and no `WallQuadraturePlan` write
path.

- [ ] **Step 3: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_immersed_operator' --output-on-failure
```

- [ ] **Step 4: Refactor only the necessary Stage 2 private helpers**

Extract coefficient-row construction and shared-face accumulation into
private functions callable by both accepted Stage 2 and Stage 3 code. Keep all
Stage 2 public signatures and arithmetic order unchanged; v2 calls must still
select the original body-fitted branch bitwise.

- [ ] **Step 5: Implement immutable snapshot and group replacement**

Create one snapshot per active row, remove the complete group-mask partition
once, build `W_P` once, apply the row-wide substitution map once, then add
`W_P`. Sort row/group/link IDs only during plan construction, not each apply.

- [ ] **Step 6: Add mutation-sensitive tests**

Run every Task 6 mutation through the actual finite-volume adapter. Also
mutate shared-face sign, leave one P--G non-orthogonal term, inject a second
wall source and vary inactive backing payload. Each mutation must fail a
specific oracle.

- [ ] **Step 7: Run task matrix**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_immersed_operator
done
ctest --test-dir build/debug -R \
  'test_(cell_centered_fvm|immersed_operator)' --output-on-failure
```

- [ ] **Step 8: Main requirements and quality passes**

Use codegraphf to inspect every `compute_gradient`,
`accumulate_viscous_residual`, `accumulate_scalar_diffusive_residual`,
Rhie--Chow and transport caller. Confirm Stage 2 callers do not accidentally
select IBM code and no test-only access macro leaks into tests-off.

- [ ] **Step 9: Run exact-HEAD Debug and focused Release/ASan/UBSan**

Release runs the full 1/2/4 operator matrix. Sanitizers run unit/header and
1/2-rank operator cases. Complete Debug is mandatory after any private-helper
repair.

- [ ] **Step 10: Commit**

```bash
git add finite_volume/include/hundun/finite_volume/immersed_operator.hpp \
  finite_volume/src/immersed_operator.cpp \
  finite_volume/src/immersed_operator_test_access.hpp \
  finite_volume/include/hundun/finite_volume/cell_centered_fvm.hpp \
  finite_volume/src/cell_centered_fvm.cpp \
  tests/mpi/test_immersed_operator.cpp \
  tests/unit/test_immersed_operator_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: add immersed residual operators"
```

### Task 9: Couple pressure Ghost conditions to two-corrector PISO

**Ownership:** main agent only.

**Files:**

- Create: `flow/include/hundun/flow/stage3_flow.hpp`
- Create: `flow/src/stage3_flow.cpp`
- Create: `flow/src/stage3_piso_detail.hpp`
- Create: `flow/src/stage3_flow_test_access.hpp`
- Modify: `flow/include/hundun/flow/flow_state.hpp`
- Modify: `flow/include/hundun/flow/constant_density_piso.hpp`
- Modify: `flow/include/hundun/flow/momentum_predictor.hpp`
- Modify: `flow/src/constant_density_piso.cpp`
- Modify: `flow/src/momentum_predictor.cpp`
- Modify: `finite_volume/src/matrix_free_poisson.cpp`
- Modify: `finite_volume/src/poisson_boundary_adapter.cpp`
- Create: `tests/mpi/test_immersed_piso.cpp`
- Create: `tests/mpi/test_immersed_pressure_operator.cpp`
- Create: `tests/unit/test_stage3_flow_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: `ImmersedOperatorAdapter`, active boundaries, Stage 2 solvers,
  `FlowState`.
- Produces: `FixedStepStage3Flow` and constant-density subset of
  `Stage3StepAttemptReport` from §3.4.

- [ ] **Step 1: Write RED pressure-wall coefficient tests**

For every link verify:

```text
mdot_wall^(k+1)
  = mdot_wall^* - D_wall * normal_gradient(pi')_wall
  = 0

D_wall =
  rho_wall * effective_transformed_measure
  * actual_momentum_velocity_correction
```

`mdot_wall^*` must include complete predictor and BE/BDF2 face-history
discrepancy. `D_wall` must change when actual momentum diagonal, density,
geometry or attempt coefficients change.

- [ ] **Step 2: Write RED failure and revision tests**

Cover non-positive/non-finite `rho_wall`, non-positive/non-finite `D_wall`,
stale Jacobi revision, rank-local failure, nullspace for inside cavity,
active pressure outlet reference, exactly two correctors and no final velocity
overwrite.

- [ ] **Step 3: Write RED transaction tests**

Force failure after corrector 1, after corrector 2 and after transport
finalization. Compare all fields/history/controller/final flux with the
authoritative bitwise helper. Corrector count stays 2; failure never adds a
third correction.

- [ ] **Step 4: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_immersed_(piso|pressure_operator)' --output-on-failure
```

- [ ] **Step 5: Implement Stage 3 pressure operator and revision**

Build RHS projection/zero-mean only for active closed domain. Dynamic wall
constants update per corrector; operator revision changes for density,
momentum diagonal, geometry, active layout or wall coefficients. Apply order
remains `Halo begin -> interior -> wait -> partition/interface`.

- [ ] **Step 6: Implement Stage 3 fixed-step transaction**

Use existing `FlowState::begin_attempt/rollback/commit`. The final transport
is recomputed from step-start state with second-corrector final shared flux.
Final wall penetration, continuity, pressure and momentum residuals are
independently recomputed before commit.
Add only the private Stage 3 friend/access hooks required by the new facade;
do not change any existing public `FlowState` method or layout.

- [ ] **Step 7: Preserve Stage 2 arithmetic**

Any shared helper extraction must leave accepted v2 byte/output and numerical
tests unchanged. Add a source-level dispatch assertion proving null
immersed/LES pointers are never consulted by Stage 2 create/attempt calls.

- [ ] **Step 8: Run GREEN 1/2/4 and checkerboard**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_immersed_pressure_operator
  mpiexec -n "$n" build/debug/test_immersed_piso
done
ctest --test-dir build/debug -R \
  'test_(fixed_step_piso|time_consistent_face_velocity|immersed_piso)' \
  --output-on-failure
```

Checkerboard parity amplitude remains `<=1e-8`; continuity normalized L2
`<=1e-10`; corrector count exactly 2.

- [ ] **Step 9: Main review and focused configurations**

Run complete Debug. Release runs 1/2/4 PISO and pressure operator. ASan/UBSan
run header, rollback, one- and two-rank cases. Tests-off full build and
standalone `stage3_flow.hpp` consumer are mandatory.

- [ ] **Step 10: Commit**

```bash
git add flow/include/hundun/flow/stage3_flow.hpp \
  flow/src/stage3_flow.cpp flow/src/stage3_piso_detail.hpp \
  flow/src/stage3_flow_test_access.hpp \
  flow/include/hundun/flow/flow_state.hpp \
  flow/include/hundun/flow/constant_density_piso.hpp \
  flow/include/hundun/flow/momentum_predictor.hpp \
  flow/src/constant_density_piso.cpp flow/src/momentum_predictor.cpp \
  finite_volume/src/matrix_free_poisson.cpp \
  finite_volume/src/poisson_boundary_adapter.cpp \
  tests/mpi/test_immersed_piso.cpp \
  tests/mpi/test_immersed_pressure_operator.cpp \
  tests/unit/test_stage3_flow_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: couple immersed walls to PISO"
```

### Task 10: Implement real-surface pressure, viscous and total force

**Ownership:** bounded worker eligible.

**Files:**

- Create: `immersed/include/hundun/immersed/wall_force.hpp`
- Create: `immersed/src/wall_force.cpp`
- Create: `tests/mpi/test_wall_force.cpp`
- Create: `tests/unit/test_wall_force_header_contract.cpp`
- Modify: `flow/include/hundun/flow/stage3_flow.hpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: `WallQuadraturePlan`, pressure, velocity gradients,
  attempt-frozen `mu_eff`.
- Produces:

```cpp
namespace hundun::immersed {

struct WallForceSample final {
  ForceComponents surface_traction;
  MomentComponents moment_about_global_origin;
  runtime::Real3 area_vector_closure_m2{};
  std::uint64_t quadrature_point_count{};
  int lowest_failing_rank{-1};
};

class WallForceIntegrator final {
 public:
  static WallForceIntegrator create(
      const WallQuadraturePlan&, const runtime::MpiContext&);
  WallForceSample integrate(
      const runtime::FieldView<const double>& mechanical_pressure,
      const runtime::FieldView<const double>& velocity,
      const runtime::FieldView<const double>& velocity_gradient,
      const runtime::FieldView<const double>& mu_eff_by_cell) const;
};

}  // namespace hundun::immersed
```

- [ ] **Step 1: Write RED sign/unit/quadrature tests**

Use constant pressure on a closed surface for zero net force and nonzero
linear pressure for analytic force. Verify:

```text
sigma = -pi*I + tau
F_pressure = integral(-pi*n_s)dA
F_viscous  = integral(tau*n_s)dA
F_total    = F_pressure + F_viscous
M_origin   = integral(x cross (sigma*n_s))dA
```

`p0` never enters force. Force units are N; moment units are N m and the
reference point is the fixed global coordinate origin `(0,0,0) m`.

- [ ] **Step 2: Write RED ownership and determinism tests**

Every triangle is owned/integrated once globally. Repartition and triangle
input permutation produce the same stable contribution ordering and the
approved decomposition tolerance.

- [ ] **Step 3: Write RED operator/surface separation tests**

Changing `WallQuadraturePlan` must change only surface traction, never solver
row. Changing operator replacement changes operator reaction. The three
pressure/viscous/total consistency residuals remain separate.

- [ ] **Step 4: Implement read-only integration**

Reconstruct `pi,u,grad(u),mu_eff` at each real quadrature point from the
cell fields and the quadrature point's accepted quadratic reconstruction. Use
deviatoric stress for molecular+WALE viscosity. Accumulate locally in stable
triangle/point order, then three explicit force and three explicit moment MPI
reductions with collective failure classification.

- [ ] **Step 5: Prove diagnostic/transaction neutrality**

Repeated collection from committed state is deterministic and does not change
fields, cache, revision, generation, allocation identity or business
counters. Failed attempt report is discarded and not persisted.

- [ ] **Step 6: Run task matrix**

```bash
cmake --build --preset debug -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_wall_force
done
ctest --test-dir build/debug -R \
  'test_wall_force_header_contract' --output-on-failure
```

Run ASan/UBSan one/two-rank force cases and tests-off header/build.

- [ ] **Step 7: Candidate commit**

```bash
git add immersed/include/hundun/immersed/wall_force.hpp \
  immersed/src/wall_force.cpp tests/mpi/test_wall_force.cpp \
  tests/unit/test_wall_force_header_contract.cpp \
  flow/include/hundun/flow/stage3_flow.hpp CMakeLists.txt
git commit -s -m "feat: add immersed wall force integration"
```

### Task 11: Close the full constant-density second-order IBM hard gate

**Ownership:** main agent only.

**Files:**

- Create: `tests/support/stage3_test_contracts.hpp`
- Create: `tests/support/stage3_dual3.hpp`
- Create: `tests/support/stage3_mms.hpp`
- Create: `tests/support/stage3_mms.cpp`
- Create: `tests/support/stage3_state_equality.hpp`
- Create: `tests/numerical/test_laminar_ibm_order.cpp`
- Create: `tests/numerical/test_laminar_ibm_engineering.cpp`
- Create: `tests/mpi/test_laminar_ibm_decomposition.cpp`
- Create: `tests/mpi/test_immersed_closed_transient.cpp`
- Create: `tests/mpi/test_immersed_transaction.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: complete Tasks 2--10 product path.
- Produces: no new public product API; authoritative Stage 3 equality helpers
  and the full Gate 3 scientific verdict.

- [ ] **Step 1: Implement independent C++ MMS and STL fixtures**

Implement §4.4 shapes, vector potential, pressure, automatic derivatives,
cell-average source quadrature, exact history and analytic force reference.
The test oracle may use product `MeshGeometry` coordinates/volumes but may not
call product IBM reconstruction, LFP, force or residual functions.

- [ ] **Step 2: Prove equality helper mutations**

```text
exact FlowState copy                        pass
one ordinary velocity element changed      fail
one nested transported element changed     fail
one metadata double changed by 1 ULP        fail
one inactive +0.0 changed to -0.0           fail
one final face flux changed                 fail
```

- [ ] **Step 3: Record fast RED**

Fast selector uses `12^3/24^3` on the sphere and intentionally disables one
quadratic term to prove the order oracle fails. Restore product code before
GREEN; fast case cannot claim acceptance.

- [ ] **Step 4: Run full uniform and warped acceptance sequences**

```bash
cmake --build --preset release -j2
ctest --test-dir build/release -R \
  '^test_laminar_ibm_order_.*_acceptance$' \
  --output-on-failure -j1
```

Run sphere, finite cylinder and oblique prism on `12^3/24^3/48^3`, uniform
and warped, plus the two translated-sphere sequences and the inside-sphere
cavity sequence. Enforce every §4.4 error/order/force/consistency/nullspace
contract.

- [ ] **Step 5: Run 1/2/4-rank decomposition and process grids**

```bash
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_laminar_ibm_decomposition acceptance
done
```

For ranks 2/4, run both frozen legal process grids; compare fields, forces,
link/group fingerprints and lowest failure rank.

- [ ] **Step 6: Run the laminar engineering matrix**

```bash
cmake --build --preset release -j2
mpiexec -n 4 build/release/test_laminar_ibm_engineering \
  cylinder_acceptance
mpiexec -n 4 build/release/test_laminar_ibm_engineering sphere_acceptance
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_laminar_ibm_engineering \
    decomposition_fast
  mpiexec -n "$n" build/release/test_immersed_closed_transient acceptance
done
```

Enforce all §4.5 run lengths, steady-window, drag decomposition,
recirculation, symmetry, kinetic-energy, residual and MPI contracts. A fast
case supplies only decomposition/development evidence and cannot satisfy the
engineering acceptance row.

- [ ] **Step 7: Run rollback and failure classification**

Inject one rank-local failure at:

```text
rho_wall validation
D_wall validation
momentum solve
pressure correction 1
pressure correction 2
final wall penetration
final force reconstruction
```

Verify collective classification, lowest rank, complete bitwise rollback and
no third corrector.

- [ ] **Step 8: Main full-diff requirements review**

Review from Gate 2 accepted HEAD through candidate. Map every Stage 3 spec
section 6--12 and 18.1--18.3 to implementation and tests. Explicitly search
for inactive reads, post-solve overwrite, wall-function terms, duplicate wall
sources and full-row paper overclaim.

- [ ] **Step 9: Main code-quality review**

Audit allocation-free apply, deterministic ordering, exception/collective
boundaries, integer/count overflow, public-header ownership, PIMPL moves,
tests-off, test seam macros and Stage 2 call paths.

- [ ] **Step 10: Run exact-HEAD complete Debug and focused sanitizers**

```bash
cmake --preset debug
cmake --build --preset debug -j2
ctest --preset debug --output-on-failure
cmake --build --preset asan -j2
ctest --test-dir build/asan -R \
  'test_(immersed_|wall_force|laminar_ibm_.*fast)' --output-on-failure -j1
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R \
  'test_(immersed_|wall_force|laminar_ibm_.*fast)' --output-on-failure -j1
```

- [ ] **Step 11: Commit**

```bash
git add tests/support/stage3_test_contracts.hpp \
  tests/support/stage3_dual3.hpp \
  tests/support/stage3_mms.hpp tests/support/stage3_mms.cpp \
  tests/support/stage3_state_equality.hpp \
  tests/numerical/test_laminar_ibm_order.cpp \
  tests/numerical/test_laminar_ibm_engineering.cpp \
  tests/mpi/test_laminar_ibm_decomposition.cpp \
  tests/mpi/test_immersed_closed_transient.cpp \
  tests/mpi/test_immersed_transaction.cpp CMakeLists.txt
git commit -s -m "test: accept second-order laminar IBM"
```

Gate 3 只有在本 task 的完整 constant-density scientific matrix accepted 后
才关闭；Tasks 7--10 的组件 GREEN 不能单独宣称二阶 IBM。

---

## Gate 4 — Standalone WALE and Body-fitted LES

### Task 12: Implement the backend-neutral WALE core

**Ownership:** bounded worker eligible.

**Files:**

- Create: `les/include/hundun/les/wale.hpp`
- Create: `les/src/wale.cpp`
- Create: `les/src/wale_test_access.hpp`
- Create: `tests/unit/test_wale.cpp`
- Create: `tests/mpi/test_wale_mpi.cpp`
- Create: `tests/unit/test_wale_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: deterministic owned-first local-active global-cell IDs and owned
  count, mesh volumes,
  execution context, externally prepared nine-component lagged-velocity
  gradient and trial-density checked views.
- Produces: `WaleModel`, `WaleAttemptCoefficients`,
  `WaleCoefficientIdentity`, `WaleSummary` from §3.4.

- [ ] **Step 1: Write RED tensor oracle**

Independent test code computes:

```text
g_ij  = partial(u_i)/partial(x_j)
S_ij  = 0.5*(g_ij+g_ji)
G_ij  = g_ik*g_kj
Sd_ij = 0.5*(G_ij+G_ji) - delta_ij*G_kk/3

nu_t = (Cw*Delta)^2 * (Sd:Sd)^(3/2)
       / [(S:S)^(5/2) + (Sd:Sd)^(5/4)]
```

Cover fixed nonsymmetric tensors, orthogonal rotations and dimension scaling.

- [ ] **Step 2: Write RED exact-zero and y-cubed tests**

When both mathematical invariants are zero, `nu_t` must be bitwise `+0.0`;
there is no fixed epsilon. Run §4.6 wall sequence and require slopes
`[2.9,3.1]`.

- [ ] **Step 3: Write RED state/identity/failure tests**

Cover constant/material/ideal density inputs, negative/non-finite density,
non-finite gradient/config, metadata/input-fingerprint changes, retry metadata
with changed `dt/rho_attempt`, no committed model revision and deterministic
`WaleCoefficientIdentity`. The WALE core does not construct `u_lag`; that
flow-level contract is Task 13.

- [ ] **Step 4: Run RED**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R 'test_wale' --output-on-failure
```

- [ ] **Step 5: Implement attempt-local coefficients**

Use `Delta=cbrt(active background volume)` and
`mu_sgs=rho_attempt*nu_t`. Store coefficient workspace in move-only
attempt-local buffers; never add it to committed `FlowState`. Compute once
before momentum predictor and expose only the active-order backend-neutral
`VectorView` accessors frozen in §3.4.

- [ ] **Step 6: Keep kernels backend-neutral**

Use `ExecutionContext`; obtain input `KernelCellView` objects only inside
`with_kernel_cell_view`, and obtain output host pointers once from the
backend-neutral buffer views before the cell loop. Kernel views remain
trivially copyable and contain no virtual calls/shared ownership/checks.
Device test double proves capability rejection but is not production
registerable.

- [ ] **Step 7: Run GREEN, MPI, sanitizers and tests-off**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R 'test_wale' --output-on-failure
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_wale_mpi
done
cmake --build --preset asan -j2
ctest --test-dir build/asan -R 'test_wale' --output-on-failure
cmake --build --preset ubsan -j2
ctest --test-dir build/ubsan -R 'test_wale' --output-on-failure
cmake -S . -B build/tests-off-stage3-task12 \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/tests-off-stage3-task12 -j2
```

- [ ] **Step 8: Candidate commit**

```bash
git add les/include/hundun/les/wale.hpp \
  les/src/wale.cpp les/src/wale_test_access.hpp \
  tests/unit/test_wale.cpp tests/mpi/test_wale_mpi.cpp \
  tests/unit/test_wale_header_contract.cpp CMakeLists.txt
git commit -s -m "feat: add WALE model"
```

### Task 13: Integrate body-fitted WALE flow and close Gate 4

**Ownership:** main agent only.

**Files:**

- Modify: `flow/include/hundun/flow/stage3_flow.hpp`
- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_piso_detail.hpp`
- Modify: `finite_volume/include/hundun/finite_volume/cell_centered_fvm.hpp`
- Modify: `finite_volume/src/cell_centered_fvm.cpp`
- Create: `tests/unit/test_variable_viscosity_fvm.cpp`
- Create: `tests/mpi/test_wale_body_fitted.cpp`
- Create: `tests/numerical/test_wale_taylor_green.cpp`
- Create: `tests/mpi/test_wale_transaction.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Task 12 WALE with the topology's owned-first ordered local
  owned+ghost cell list.
- Produces: complete `none/wale` Stage 3 fixed-step/retry composition and
  canonical all-active identity layout semantics.

- [ ] **Step 1: Write RED none/wale identity tests**

Prove all background cells active, Stage 2 `BoundaryRegistry` used directly,
no surface/link/stencil/LFP/wall-force objects/counters/providers constructed,
and only WALE presence is reported.

- [ ] **Step 2: Write RED lagging transaction tests**

Capture `rho_attempt` and extrapolated `u_lag` before momentum predictor.
For startup require `u_lag=u^n`; for BDF2 require
`u_lag=u^n+(dt_attempt/dt_previous)*(u^n-u^(n-1))`. Build its nine-component
gradient with the accepted body-fitted gradient path, pass that field to
`WaleModel`, and prove the same coefficients are used by predictor, two correctors,
finalization and final residual; a second post-corrector WALE evaluation is a
failing mutation.

- [ ] **Step 3: Implement variable-viscosity deviatoric stress**

Add a `CellCenteredFvmOperators` overload that accepts
`FaceFieldView<const double> dynamic_viscosity_by_face`; retain the existing
scalar-viscosity overload and prove it is bitwise unchanged. Interpolate each
active-cell `mu_sgs` to faces once, form `mu_eff=mu+mu_sgs`, and pass the
frozen face field through deviatoric viscous residual.
Add `mu_sgs/Pr_t` and `mu_sgs/Sc_t` to h/scalar diffusion. Do not add SGS
isotropic stress outside mechanical pressure.

- [ ] **Step 4: Run fast TGV/channel GREEN**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_wale_(body_fitted|taylor_green).*fast' \
  --output-on-failure -j1
```

- [ ] **Step 5: Run full Release acceptance**

```bash
cmake --build --preset release -j2
ctest --test-dir build/release -R \
  'test_wale_(body_fitted|taylor_green).*acceptance' \
  --output-on-failure -j1
```

Enforce every §4.6 channel/TGV threshold and 1/2/4 decomposition matrix.

- [ ] **Step 6: Run failure/rollback/MPI**

Inject rank-local non-finite gradient and failed final residual. Verify lowest
rank, bitwise state/controller rollback, discarded coefficient identity and
deterministic retry recomputation.

- [ ] **Step 7: Main reviews and exact-HEAD configurations**

Review all variable-viscosity and scalar-diffusion callers. Run complete
Debug; ASan/UBSan fast cases; Release full acceptance; tests-off full build.
Confirm Stage 2 constant/material/ideal tests still use molecular viscosity
only.

- [ ] **Step 8: Commit**

```bash
git add flow/include/hundun/flow/stage3_flow.hpp \
  flow/src/stage3_flow.cpp flow/src/stage3_piso_detail.hpp \
  finite_volume/include/hundun/finite_volume/cell_centered_fvm.hpp \
  finite_volume/src/cell_centered_fvm.cpp \
  tests/unit/test_variable_viscosity_fvm.cpp \
  tests/mpi/test_wale_body_fitted.cpp \
  tests/numerical/test_wale_taylor_green.cpp \
  tests/mpi/test_wale_transaction.cpp CMakeLists.txt
git commit -s -m "feat: integrate body-fitted WALE flow"
```

Gate 4 closes only when tensor/y³, TGV, body-fitted channel, density-input
unit contracts, rollback and 1/2/4-rank evidence are accepted.

---

## Gate 5 — Variable Density and Combined IBM+WALE

### Task 14: Integrate material-density transport with LFP-GCIBM

**Ownership:** main agent only.

**Files:**

- Modify: `flow/include/hundun/flow/stage3_flow.hpp`
- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/src/stage3_piso_detail.hpp`
- Modify: `flow/include/hundun/flow/material_density_transport.hpp`
- Modify: `flow/include/hundun/flow/material_density_piso.hpp`
- Modify: `flow/src/material_density_transport.cpp`
- Modify: `flow/src/material_density_piso.cpp`
- Create: `tests/mpi/test_material_density_ibm.cpp`
- Create: `tests/mpi/test_material_density_ibm_transaction.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: accepted Stage 2 material transport/PISO, all Stage 3 IBM plans.
- Produces: material-density branch of `FixedStepStage3Flow`; no new density
  boundary model.

- [ ] **Step 1: Write RED nonzero wall-normal density-gradient tests**

Use:

```text
rho(x,y,z) = 1 + 0.10*x + 0.05*y - 0.03*z
```

on sphere/prism interfaces. Donor cell averages are positive, and analytic
normal gradient is generally nonzero. Verify `rho_wall` comes only from
fluid-side quadratic extrapolation; a homogeneous-Neumann mutation fails.

- [ ] **Step 2: Write RED conservative transport tests**

Initialize a positive nonuniform density in the closed MMS domain and advance
with divergence-free no-slip velocity and no density source. Require:

```text
wall mass flux = mathematical zero
same final FaceMassFlux for rho, rho*h and every rho*phi
active-fluid total mass relative error <= 5e-12
rho finite and positive
inactive slots bitwise +0.0
```

- [ ] **Step 3: Write RED `rho_wall`/`D_wall` retry tests**

Construct positive donors whose quadratic extrapolation becomes negative, and
another that becomes NaN through injected trial data. Both are collective
recoverable failures; clipping, epsilon replacement, zero-gradient fallback
and stencil replacement mutations must fail.

- [ ] **Step 4: Implement Stage 3 material branch**

Reuse Stage 2 conservative transport and final-flux provenance through
active-layout iteration and IBM-aware reconstruction. Compute material
`rho_wall` before each pressure correction, validate it and dynamic `D_wall`,
then freeze them for that correction. Public Stage 2 signatures remain
unchanged; header edits are limited to private detail-bridge/friend
declarations used by `FixedStepStage3Flow`.

- [ ] **Step 5: Run fast GREEN**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_material_density_ibm.*fast' --output-on-failure
```

- [ ] **Step 6: Run full 1/2/4 acceptance**

```bash
cmake --build --preset release -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_material_density_ibm acceptance
  mpiexec -n "$n" \
    build/release/test_material_density_ibm_transaction acceptance
done
```

Compare fields/final flux/fingerprints to the approved decomposition
threshold; transaction state uses bitwise equality.

- [ ] **Step 7: Main impact review**

Use codegraphf on every
`MaterialDensityTransport::advance`,
`FixedStepMaterialDensityFlow::attempt`, gradient and PISO caller. Confirm
v2 paths remain full-layout/body-fitted and Stage 3 paths never scan inactive
values for positivity.

- [ ] **Step 8: Exact-HEAD tests and commit**

Run complete Debug, ASan/UBSan fast 1/2-rank, Release 1/2/4 acceptance, then:

```bash
git add flow/include/hundun/flow/stage3_flow.hpp \
  flow/src/stage3_flow.cpp flow/src/stage3_piso_detail.hpp \
  flow/include/hundun/flow/material_density_transport.hpp \
  flow/include/hundun/flow/material_density_piso.hpp \
  flow/src/material_density_transport.cpp \
  flow/src/material_density_piso.cpp \
  tests/mpi/test_material_density_ibm.cpp \
  tests/mpi/test_material_density_ibm_transaction.cpp CMakeLists.txt
git commit -s -m "feat: couple material density to immersed walls"
```

### Task 15: Integrate ideal-gas closure with LFP-GCIBM

**Ownership:** main agent only.

**Files:**

- Modify: `flow/include/hundun/flow/stage3_flow.hpp`
- Modify: `flow/src/stage3_flow.cpp`
- Modify: `flow/include/hundun/flow/ideal_gas_closure.hpp`
- Modify: `flow/include/hundun/flow/ideal_gas_piso.hpp`
- Modify: `flow/src/ideal_gas_closure.cpp`
- Modify: `flow/src/ideal_gas_piso.cpp`
- Create: `tests/mpi/test_ideal_gas_ibm.cpp`
- Create: `tests/mpi/test_ideal_gas_ibm_transaction.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Stage 2 `IdealGasClosure`, active volumes, IBM h Neumann weights.
- Produces: ideal-gas branch of `FixedStepStage3Flow`; same
  `h=cp*T`, `rho=p0/(R*T)` closure.

- [ ] **Step 1: Write RED closed-domain active-volume tests**

Use:

```text
T(x,y,z) = 300 + 20 * E(x,y,z) * F_b(x,y,z)^2
h = cp*T
p0 = M_target*R / sum_active(V_i/T_i)
```

This has zero normal derivative at outer and immersed walls. Verify only
active volumes enter `p0`, mass target is conserved and inactive `T/rho`
zeros are not treated as physical states.

- [ ] **Step 2: Write RED open-domain fixed-p0 tests**

Use the finite sphere inlet/outlet geometry with configured positive finite
`p0`. Pressure outlet constrains only `pi`; h/T inlet authority remains
exactly one source, and redundant h/T/rho is cross-validated with the existing
tolerance.

- [ ] **Step 3: Write RED invalid-state and rollback tests**

Cover non-positive/non-finite T/rho/h, closure mismatch, dynamic-p0 failure,
IBM reconstruction failure and rank-local failure after p0 update. State,
history, controller, final flux and p0 rollback bitwise.

- [ ] **Step 4: Implement active ideal-gas closure adapter**

Reuse the accepted closure formula but supply active-cell volume iteration
and IBM-aware zero-h-gradient reconstruction. Do not add an enthalpy feedback
outer iteration, species-dependent cp/R or a second density field. Public
Stage 2 signatures remain unchanged; header edits are private
detail-bridge/friend declarations only.

- [ ] **Step 5: Run fast and acceptance matrices**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_ideal_gas_ibm.*fast' --output-on-failure
cmake --build --preset release -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_ideal_gas_ibm acceptance
  mpiexec -n "$n" build/release/test_ideal_gas_ibm_transaction acceptance
done
```

Enforce closure identities `<=1e-12`, closed mass `<=5e-12`, Stage 2 final
residuals and decomposition thresholds.

- [ ] **Step 6: Main reviews and configurations**

Review all p0/closure state and Checkpoint-facing access. Run complete Debug,
ASan/UBSan fast 1/2-rank, Release acceptance. Confirm Stage 2 open/closed
ideal-gas tests remain unchanged.

- [ ] **Step 7: Commit**

```bash
git add flow/include/hundun/flow/stage3_flow.hpp \
  flow/src/stage3_flow.cpp \
  flow/include/hundun/flow/ideal_gas_closure.hpp \
  flow/include/hundun/flow/ideal_gas_piso.hpp \
  flow/src/ideal_gas_closure.cpp \
  flow/src/ideal_gas_piso.cpp tests/mpi/test_ideal_gas_ibm.cpp \
  tests/mpi/test_ideal_gas_ibm_transaction.cpp CMakeLists.txt
git commit -s -m "feat: couple ideal gas to immersed walls"
```

### Task 16: Close the combined IBM+WALE hard gate

**Ownership:** main agent only.

**Files:**

- Modify: `flow/include/hundun/flow/stage3_flow.hpp`
- Modify: `flow/src/stage3_flow.cpp`
- Create: `tests/mpi/test_ibm_wale_wake.cpp`
- Create: `tests/mpi/test_ibm_wale_density_models.cpp`
- Create: `tests/mpi/test_stage3_retry.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: Tasks 11--15.
- Produces: all three density models with legal `LFP-GCIBM/wale` composition,
  final Gate 5 scientific verdict.

- [ ] **Step 1: Write RED combined coefficient/provenance tests**

For each density model prove:

```text
rho_attempt chosen after transport prediction/closure
u_lag chosen from committed/history and attempted dt
WALE evaluated exactly once
mu_eff frozen through predictor/two correctors/finalization/force
all conservative equations use the same final FaceMassFlux
wall traction uses the same frozen mu_eff
```

- [ ] **Step 2: Write RED retry identity tests**

Inject recoverable failure after first attempt. Retry halves `dt`, recomputes
time stencil/rho_attempt/WALE identity from unchanged committed state, discards
old force report, and commits only the successful second attempt. Eight-retry
and min-dt termination follow Stage 2 exactly.

- [ ] **Step 3: Implement combined orchestration**

Stage 3 flow order is fixed:

```text
begin trial
-> density/transport prediction and ideal-gas closure when selected
-> build u_lag from committed/history and attempted dt
-> build one nine-component velocity-gradient field through
   ImmersedReconstruction on interface rows and the accepted background
   provider elsewhere
-> compute attempt-local WALE from that gradient and rho_attempt
-> momentum predictor with IBM mu_eff
-> pressure corrector 1
-> provisional transport
-> pressure corrector 2
-> final transport from step-start state and final flux
-> final closure
-> residual/conservation/wall penetration/force
-> collective commit or rollback
```

- [ ] **Step 4: Run fast sphere wake**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_ibm_wale_.*fast|test_stage3_retry' --output-on-failure -j1
```

- [ ] **Step 5: Run full Release finite-sphere wake**

```bash
cmake --build --preset release -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_ibm_wale_wake acceptance
done
```

Use §4.6 geometry/grid/steps and all wall/residual/mass/force/WALE/MPI
thresholds.

- [ ] **Step 6: Run density-model matrix**

Constant runs the full wake. Material and ideal-gas run a shortened 20-step
acceptance on the `96x48x48` grid, still enforcing conservation, closure,
positivity, rollback and decomposition; they do not claim engineering
statistics.

- [ ] **Step 7: Main full-diff requirements and quality reviews**

Review all changes from Gate 4 accepted HEAD. Search for multiple WALE
evaluations, wrong density time layer, final-velocity coefficient refresh,
stateful model revision, non-final flux use, wall-function code and
uncommitted report persistence.

- [ ] **Step 8: Run exact-HEAD configurations**

Complete Debug; Release full Gate 5 acceptance; ASan/UBSan fast combined cases;
tests-off full build. Any product/public-header repair invalidates and repeats
complete Debug.

- [ ] **Step 9: Commit**

```bash
git add flow/include/hundun/flow/stage3_flow.hpp \
  flow/src/stage3_flow.cpp tests/mpi/test_ibm_wale_wake.cpp \
  tests/mpi/test_ibm_wale_density_models.cpp \
  tests/mpi/test_stage3_retry.cpp CMakeLists.txt
git commit -s -m "test: accept combined IBM and WALE"
```

Gate 5 closes only after constant/material/ideal-gas, retry, finite-body wake
and 1/2/4-rank evidence all pass.

---

## Gate 6 — Checkpoint, Diagnostics, Driver and Exit

### Task 17: Implement transactional Checkpoint v3

**Ownership:** bounded worker eligible.

**Files:**

- Create: `flow/include/hundun/flow/checkpoint_v3.hpp`
- Create: `flow/src/checkpoint_v3.cpp`
- Create: `flow/src/checkpoint_v3_detail.hpp`
- Create: `flow/src/checkpoint_v3_test_access.hpp`
- Modify: `flow/include/hundun/flow/flow_state.hpp`
- Create: `tests/unit/test_checkpoint_v3_header_contract.cpp`
- Create: `tests/unit/test_checkpoint_v3_protocol.cpp`
- Create: `tests/mpi/test_checkpoint_v3.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
namespace hundun::flow {

enum class IbmPresenceTag : std::uint32_t {
  absent = 0,
  lfp_gcibm_v1 = 0x00030101U
};
enum class LesPresenceTag : std::uint32_t {
  absent = 0,
  wale_v1 = 0x00030d01U
};
struct IbmStaticIdentity final {
  std::uint64_t normalized_stl_bytes{};
  std::uint64_t normalized_stl_content{};
  double length_scale_to_m{};
  config::ImmersedFluidSide fluid_side{config::ImmersedFluidSide::outside};
  std::uint64_t surface_topology{};
  std::uint64_t surface_orientation{};
  std::uint64_t surface{};
  std::uint64_t classification{};
  std::uint64_t interface_links{};
  std::uint64_t active_layout{};
  std::uint64_t triangle_ownership{};
  std::uint64_t surface_coverage{};
  std::uint64_t active_boundary{};
  std::uint64_t ghost_stencil{};
  std::uint64_t wall_quadrature{};
  std::uint64_t local_flow_pattern{};
  std::uint64_t wall_force_integrator{};
  std::uint64_t wall_force_diagnostics{};
};
struct WaleStaticIdentity final {
  std::uint64_t configuration{};
  std::uint64_t transient_field_schema{};
};
struct Stage3StaticIdentity final {
  std::uint64_t canonical_schema_v3{};
  IbmPresenceTag ibm{IbmPresenceTag::absent};
  LesPresenceTag les{LesPresenceTag::absent};
  std::optional<IbmStaticIdentity> ibm_identity;
  std::optional<WaleStaticIdentity> les_identity;
};

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext&,
    const runtime::StructuredDecomposition&,
    const mesh::MeshTopology&, const mesh::MeshGeometry&,
    const config::ImmersedFlowCaseConfig&,
    const Stage3StaticIdentity&, const FlowState&,
    const TimeControlState&, std::optional<IdealGasClosureState>,
    const std::filesystem::path&);

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext&,
    const runtime::StructuredDecomposition&,
    const mesh::MeshTopology&, const mesh::MeshGeometry&,
    const config::ImmersedFlowCaseConfig&,
    const Stage3StaticIdentity&, FlowState&,
    const std::filesystem::path&);

}  // namespace hundun::flow
```

`CheckpointV3Report/ReadResult` mirror v2 value-query style but use distinct
enums/types; v2 headers and bytes remain untouched.

- [ ] **Step 1: Write RED byte-layout tests**

Freeze:

```text
manifest.v3.bin
rank-%06u.v3.bin, where %06u is the zero-padded rank
COMPLETED published last
little-endian integer fields
binary64 FP fields
CRC-64/ECMA-182 per rank and manifest
exact EOF, no trailing bytes
```

Presence tags are canonical algorithm/version IDs. `absent` requires a
disengaged matching optional, zero section count/bytes and no fake
fingerprint; either non-absent tag requires exactly one matching identity
section. Unknown tag, tag/optional mismatch and duplicate section are
integrity failures.

- [ ] **Step 2: Write RED semantic-state tests**

Persist committed/history, dt/order/controller, final face flux, dynamic p0,
field schema, full inactive `+0.0` slots and all static fingerprints. Do not
persist BVH, raw address, cached weights, force report, `nu_t` or WALE attempt
identity.

- [ ] **Step 3: Write RED corruption/mismatch tests**

Cover missing marker/rank, CRC, size, endian marker, version, presence,
canonical schema, normalized STL bytes/content, length scale, fluid side,
surface topology/orientation, classification, interface links, active layout,
triangle ownership, coverage, active boundary, stencil, wall quadrature,
LFP, wall-force, WALE configuration/transient-schema fingerprints, ranks,
process-grid, owned-box, inactive negative zero and field schema.

- [ ] **Step 4: Write RED transaction/epoch tests**

Read entry invalidates old checked views before I/O. Every failed read leaves
field values, committed metadata, history, controller and p0 bitwise
unchanged. No stale view is dereferenced.

- [ ] **Step 5: Implement isolated v3 codec and publication**

Reuse the project CRC implementation through a shared private utility without
changing v2 serialized order. Write temporary rank files, collectively
validate, publish rank files, manifest and finally `COMPLETED`.

- [ ] **Step 6: Implement deterministic restore**

Driver rebuilds static plans before read and passes `Stage3StaticIdentity`.
Reader validates all bytes/fingerprints into temporary values, prepares field
replacement, then publishes once after collective success.

- [ ] **Step 7: Run unit/MPI/corruption GREEN**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_checkpoint_v3' --output-on-failure
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_checkpoint_v3 acceptance
done
```

- [ ] **Step 8: Run continuation and sanitizers**

For IBM-only, WALE-only and combined density models, continuous vs
write/read/continue must be bitwise equal for fields/history/next dt/order/
final flux/p0 and next recomputed IBM/WALE result. Run ASan/UBSan corruption
and 1/2-rank cases; tests-off/header consumer mandatory.

- [ ] **Step 9: Candidate commit**

```bash
git add flow/include/hundun/flow/checkpoint_v3.hpp \
  flow/src/checkpoint_v3.cpp flow/src/checkpoint_v3_detail.hpp \
  flow/src/checkpoint_v3_test_access.hpp \
  flow/include/hundun/flow/flow_state.hpp \
  tests/unit/test_checkpoint_v3_header_contract.cpp \
  tests/unit/test_checkpoint_v3_protocol.cpp \
  tests/mpi/test_checkpoint_v3.cpp CMakeLists.txt
git commit -s -m "feat: add Checkpoint v3"
```

### Task 18: Add Stage 3 structured diagnostics and exact counters

**Ownership:** bounded worker eligible.

**Files:**

- Modify: `diagnostics/include/hundun/diagnostics/structured_diagnostics.hpp`
- Create: `diagnostics/include/hundun/diagnostics/stage3_module_diagnostics.hpp`
- Create: `diagnostics/include/hundun/diagnostics/checkpoint_v3_diagnostics.hpp`
- Create: `diagnostics/src/stage3_module_diagnostics.cpp`
- Create: `diagnostics/src/checkpoint_v3_diagnostics.cpp`
- Create: `tests/unit/test_stage3_module_diagnostics.cpp`
- Create: `tests/mpi/test_stage3_diagnostics_mpi.cpp`
- Create: `tests/mpi/test_checkpoint_v3_diagnostics.cpp`
- Create: `tests/unit/test_stage3_diagnostics_header_contract.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: read-only public summaries from surface/domain/stencil/LFP/force/
  WALE/flow/checkpoint.
- Produces additive `DiagnosticModuleKind` enum values:

```cpp
immersed_surface   = 18,
ghost_stencil      = 19,
local_flow_pattern = 20,
wall_force         = 21,
les                = 22
```

and the concrete friend-free `collect_diagnostics` overload set in
`stage3_module_diagnostics.hpp`.

- [ ] **Step 1: Write RED stable-enum/header tests**

Assert every existing 0--17 underlying value, new 18--22 values and
`kDiagnosticRecordSchemaV1==1`. Re-run Stage 2 canonical-record snapshots.

- [ ] **Step 2: Write RED provider inventory/presence tests**

```text
none/wale       -> added kind 22 only
LFP/none        -> added kinds 18,19,20,21
LFP/wale        -> added kinds 18,19,20,21,22
absent module   -> no provider, instance, counter or fake zero record
```

`ActiveBoundaryLayout` remains under the existing `boundary` kind with a
stable Stage 3 instance ID; driver/checkpoint remain under their existing
kinds and report the canonical presence tags.

- [ ] **Step 3: Write RED read-only/determinism/sample tests**

Repeated local collection must be bitwise/canonical identical and perform no
collective. Explicit collective requests run 1/2/4 ranks, return same failure
class/lowest rank, obey sample budget and stable global IDs. Snapshot fields,
revisions, generations, allocation IDs, caches, controller and business
counters before/after.

- [ ] **Step 4: Implement adapters**

Record surface counts/bbox/area/closure/orientation; query ambiguity,
classification/link/active layout/coverage; active-boundary state under the
existing boundary kind; stencil/quadrature donor min/max, rank, condition,
halo and triangle coverage; LFP scale/rotation invariants, limiting-case
state, group counts and coefficient norms; PISO predictor/final wall flux,
penetration, dynamic pressure constraint and lowest failing rank;
pressure/viscous/total force, moment, area closure and three consistency
residuals with units; WALE `nu_t/mu_sgs/mu_eff` min/max/norm, zero/active
counts, config and identity; Checkpoint v3 presence/CRC/fingerprint.

- [ ] **Step 5: Add exact performance counters**

Counters:

```text
surface triangle/candidate/ray/segment queries
classification and interface-link count
stencil donor values, QR plans/rejections and maximum halo reach
Ghost constraint evaluations, Halo payload bytes/messages
replacement groups, LFP transforms and simultaneous substitutions
immersed rows, pressure wall constraints and wall quadrature points
force reductions, WALE gradient cells and WALE evaluations
Stage 3 allocations, matvecs, reductions and I/O logical bytes
```

Initialization queries must not increase during time-step apply.

- [ ] **Step 6: Run GREEN 1/2/4**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_stage3_(module_diagnostics|diagnostics_header)' \
  --output-on-failure
for n in 1 2 4; do
  mpiexec -n "$n" build/debug/test_stage3_diagnostics_mpi
  mpiexec -n "$n" build/debug/test_checkpoint_v3_diagnostics
done
```

- [ ] **Step 7: Run disabled-path and sanitizer evidence**

Diagnostics disabled/not-due must add zero allocations, collectives, field
copies or business-counter increments. Run ASan/UBSan unit and 1/2-rank;
tests-off full build/header consumer.

- [ ] **Step 8: Candidate commit**

```bash
git add diagnostics/include/hundun/diagnostics/structured_diagnostics.hpp \
  diagnostics/include/hundun/diagnostics/stage3_module_diagnostics.hpp \
  diagnostics/include/hundun/diagnostics/checkpoint_v3_diagnostics.hpp \
  diagnostics/src/stage3_module_diagnostics.cpp \
  diagnostics/src/checkpoint_v3_diagnostics.cpp \
  tests/unit/test_stage3_module_diagnostics.cpp \
  tests/unit/test_stage3_diagnostics_header_contract.cpp \
  tests/mpi/test_stage3_diagnostics_mpi.cpp \
  tests/mpi/test_checkpoint_v3_diagnostics.cpp CMakeLists.txt
git commit -s -m "feat: add Stage 3 diagnostics"
```

### Task 19: Integrate the same-executable Stage 3 driver

**Ownership:** main agent only.

**Files:**

- Create: `applications/hundun/stage3_driver.hpp`
- Create: `applications/hundun/stage3_driver.cpp`
- Modify: `applications/hundun/main.cpp`
- Modify: `applications/hundun/case_config_broadcast.cpp`
- Create: `tests/unit/test_stage3_dispatch_contract.cpp`
- Create: `tests/acceptance/stage3_flow_models.sh`
- Create: `tests/acceptance/stage3_restart_continuation.sh`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: schema v3, all Stage 3 plans/flow/checkpoint/diagnostics.
- Produces:

```cpp
int run_stage3_case(
    const CliOptions&, runtime::MpiContext&,
    const config::ImmersedFlowCaseConfig&,
    const std::filesystem::path& authoritative_case_root);
```

- [ ] **Step 1: Write RED dispatch compatibility tests**

v1 and v2 CLI/version/validate/print-resolved/run stdout snapshots must remain
byte-for-byte. v3 validate/print-resolved use v3 canonical JSON. Normal v3
run dispatches exactly once to `run_stage3_case`.

- [ ] **Step 2: Write RED three-combination driver cases**

Generate in build tree:

```text
none/wale
LFP-GCIBM/none
LFP-GCIBM/wale
```

for constant/material/ideal-gas where permitted, uniform/warped, inside/outside
and open/closed boundaries.

- [ ] **Step 3: Implement construction order**

```text
decomposition/topology/geometry/boundary
-> optional surface/query/domain
-> optional Ghost/Wall/LFP plans
-> field registry/state and canonical inactive zeros
-> optional WALE
-> solvers and FixedStepStage3Flow
-> optional Checkpoint v3 restore
-> retry loop / diagnostics / checkpoint writes
```

If IBM absent, do not construct IBM objects. If LES absent, do not construct
WALE objects. Validation completes before any time-step or output write.

- [ ] **Step 4: Implement static wall-motion test double isolation**

Only tests may construct a nonzero wall-kinematics test double behind
`HUNDUN_APPLICATION_ENABLE_TEST_ACCESS`; production schema always yields
zero velocity. Tests-off symbol scan must find no production registration.

- [ ] **Step 5: Run driver GREEN**

```bash
cmake --build --preset debug -j2
ctest --test-dir build/debug -R \
  'test_stage3_dispatch_contract|test_stage3_flow_models' \
  --output-on-failure -j1
```

- [ ] **Step 6: Run 1/2/4 restart continuation**

```bash
for n in 1 2 4; do
  bash tests/acceptance/stage3_flow_models.sh \
    build/debug/hundun "$(command -v mpiexec)" "$n"
  bash tests/acceptance/stage3_restart_continuation.sh \
    build/debug/hundun "$(command -v mpiexec)" "$n"
done
```

- [ ] **Step 7: Main complete impact review**

Use codegraphf on `main`, all three drivers, loaders/broadcast, FlowState,
checkpoint and diagnostics. Confirm Stage 1/2 paths do not instantiate or link
through Stage 3 behavior except additive libraries.

- [ ] **Step 8: Exact-HEAD configurations**

Complete Debug; Release 1/2/4 driver; ASan/UBSan fast driver/restart;
tests-off full build; standalone header; `nm` confirms no test access symbol;
`ldd` confirms only approved system/MPI dependencies.

- [ ] **Step 9: Commit**

```bash
git add applications/hundun/stage3_driver.hpp \
  applications/hundun/stage3_driver.cpp \
  applications/hundun/main.cpp \
  applications/hundun/case_config_broadcast.cpp \
  tests/unit/test_stage3_dispatch_contract.cpp \
  tests/acceptance/stage3_flow_models.sh \
  tests/acceptance/stage3_restart_continuation.sh CMakeLists.txt
git commit -s -m "feat: add Stage 3 flow driver"
```

### Task 20: Add performance evidence, capability ledger and Stage 3 gate

**Ownership:** main agent only.

**Files:**

- Create: `docs/numerics/stage3-capability-ledger.md`
- Create: `tests/support/stage3_case_generator.cpp`
- Create: `tests/support/stage3_performance_evidence.cpp`
- Create: `tests/mpi/test_stage3_performance.cpp`
- Create: `tests/acceptance/stage3_acceptance.sh`
- Create: `tests/cmake/task20_stage3_acceptance_contract.cmake`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes: all accepted Stage 3 capabilities/tests.
- Produces: exact registered Stage 3 acceptance inventory, performance JSON
  artifacts and final capability traceability.

- [ ] **Step 1: Freeze capability rows**

Every row has exactly one disposition:

```text
implemented-and-accepted
explicitly-deferred
out-of-scope
```

Rows cover schema, surface/query, classification/layout, Ghost/Wall plans,
LFP, reconstruction/operator, PISO, force, WALE, material, ideal gas,
combined flow, Checkpoint, diagnostics, driver and performance. Deferred/
excluded rows cover moving/multipart/cut-cell, wall functions/thermal walls,
production GPU/GPU-aware MPI, dynamic SGS, chemistry/species, vendor solvers
and rank-changing restore.

- [ ] **Step 2: Write acceptance inventory contract**

`stage3_acceptance.sh` has a list-only mode used only by its contract test and
a normal mode that runs the exact inventory serially. It validates source
root, Debug cache, tests enabled, build HEAD, ledger presence and every test
executable before execution.

- [ ] **Step 3: Add portable exact-counter cases**

Run:

```text
ibm_strong:
  global 48^3 sphere, ranks 1/2/4

wale_weak:
  48^3 cells/rank periodic TGV, ranks 1/2/4

ibm_wale_reference:
  global 96x64x64 finite sphere, ranks 1/2/4

warmup = 3 steps
measured = 10 steps
repetitions = 3
```

Exact bytes/messages/queries/replacements/quadrature/evaluations/collectives/
matvec/I/O counters are hard gates. Wall-clock, RSS, bandwidth and throughput
must be positive finite and are stored with compatibility metadata but have no
portable pass threshold.

- [ ] **Step 4: Add canonical performance artifacts**

Bind artifacts to source HEAD, build/test binary SHA, compiler, MPI, process
grid, geometry/config fingerprints and presence tags. Incompatible metadata
is `incomparable`, not pass/fail.

- [ ] **Step 5: Run Stage 3 acceptance inventory**

```bash
cmake --build --preset debug -j2
bash tests/acceptance/stage3_acceptance.sh
```

The inventory includes all hard-gate 1/2/4 MPI cases, all nine full laminar
order sequences, the §4.5 engineering/closed-transient matrix, WALE
acceptance, density integrations, Checkpoint, diagnostics and driver, but
excludes its own contract test.

- [ ] **Step 6: Run Release performance matrix**

```bash
cmake --build --preset release -j2
for n in 1 2 4; do
  mpiexec -n "$n" build/release/test_stage3_performance acceptance
done
```

- [ ] **Step 7: Main reviews and exact-HEAD Debug**

Check bidirectional capability traceability, exact test names/ranks/selectors,
counter formulas, no skipped test, no fast case presented as acceptance and
no unapproved capability. Run complete Debug after the final tracked change.

- [ ] **Step 8: Commit**

```bash
git add docs/numerics/stage3-capability-ledger.md \
  tests/support/stage3_case_generator.cpp \
  tests/support/stage3_performance_evidence.cpp \
  tests/mpi/test_stage3_performance.cpp \
  tests/acceptance/stage3_acceptance.sh \
  tests/cmake/task20_stage3_acceptance_contract.cmake CMakeLists.txt
git commit -s -m "test: add Stage 3 acceptance evidence"
```

### Task 21: Final Stage 3 acceptance

**Ownership:** main agent only. No implementation worker or reviewer worker.

**Files:**

- Read: entire accepted diff from the Stage 3 implementation start to
  candidate HEAD.
- Update only if evidence finds a real defect: the directly affected product,
  tests or `docs/numerics/stage3-capability-ledger.md`.
- Record ignored coordinator evidence under:
  `.superpowers/sdd/2026-*-stage3-final-acceptance.md`.

**Interfaces:**

- Consumes: accepted Tasks 1--20.
- Produces: final accepted HEAD/report SHA and a stop at the Stage 3 boundary.

- [ ] **Step 1: Audit history and workspace**

Verify every task parent chain, subject, DCO, allowed files and accepted
report. Confirm the unrelated
`docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md` remains
untouched/untracked unless the user separately changes its status.

- [ ] **Step 2: Fresh Debug configure/build/test**

```bash
cmake --preset debug
cmake --build --preset debug -j2
ctest --preset debug --output-on-failure
```

- [ ] **Step 3: Fresh Release configure/build/test**

```bash
cmake --preset release
cmake --build --preset release -j2
ctest --preset release --output-on-failure
```

- [ ] **Step 4: Fresh ASan configure/build/test**

```bash
cmake --preset asan
cmake --build --preset asan -j2
ctest --preset asan --output-on-failure
```

- [ ] **Step 5: Fresh UBSan configure/build/test**

```bash
cmake --preset ubsan
cmake --build --preset ubsan -j2
ctest --preset ubsan --output-on-failure
```

- [ ] **Step 6: Run all stage acceptance scripts**

```bash
bash tests/acceptance/stage1_acceptance.sh
bash tests/acceptance/stage2_acceptance.sh
bash tests/acceptance/stage3_acceptance.sh
```

- [ ] **Step 7: Run tests-off, offline and linkage checks**

```bash
cmake -S . -B build/tests-off-final \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/tests-off-final -j2
cmake -S . -B build/offline-final \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF \
  -DCMAKE_DISABLE_FIND_PACKAGE_Git=ON
cmake --build build/offline-final -j2
if nm -C build/tests-off-final/hundun |
    rg -q 'TestAccess|test seam|mock device'; then
  echo "tests-off executable contains a forbidden test symbol" >&2
  exit 1
fi
ldd build/tests-off-final/hundun
```

`ldd` must contain only approved system C/C++/MPI/thread/dl dependencies and no
Python/vendor solver/GPU runtime.

- [ ] **Step 8: Re-run policy/provenance and exact inventory**

Run all source-policy, standalone-header, provenance, DCO, test-registration,
diagnostic-enum and capability-ledger contracts. Confirm no skipped tests and
no public dependency on the user paper PDF path.

- [ ] **Step 9: Audit runtime/worker/process state**

Confirm every worker completed, no build/test/MPI process remains, no research
process was inspected or altered, no publication/push occurred and workspace
contains no unapproved tracked/untracked Stage 3 output.

- [ ] **Step 10: Write final acceptance report and stop**

Report:

```text
accepted HEAD and parent
all task accepted HEADs
spec/plan/report SHA-256
Debug/Release/ASan/UBSan counts and failures=0
Stage 1/2/3 gate results
1/2/4-rank numerical results and observed orders
Checkpoint continuation result
diagnostic provider coverage
performance exact counters / informative timing
DCO/policy/provenance/tests-off/linkage status
worker/process/publication audit
capability dispositions
```

Do not enter Stage 4 and do not add a commit solely to claim success. If the
ledger needed a genuine correction, commit it with DCO and repeat every
invalidated exact-HEAD gate first.

---

## 5. Spec-to-task Traceability

| Approved design sections | Implementation/tasks | Positive and failure evidence | Final acceptance |
|---|---|---|---|
| 1--3 decisions, scope and non-goals | Global constraints; Tasks 1, 20, 21 | source policy, schema-combination tests, capability disposition | no Stage 4 capability or claim |
| 4 architecture/component boundary | Tasks 2--10, 12, 19 | standalone headers, tests-off, codegraphf impact review | target/linkage audit |
| 5 schema v3/compatibility | Tasks 1, 19 | canonical round-trip, exact JSON Pointer failures, typed MPI broadcast | v1/v2 byte snapshots and same executable |
| 6.1--6.2 surface/query | Task 2 | ASCII/binary positive fixtures; topology/query/collective failures | deterministic 1/2/4 fingerprints |
| 6.3--6.4 classification/active boundaries | Task 3 | inside/outside, coverage, patch/nullspace, small-body and multi-hit failures | 1/2/4 and two process-grid evidence |
| 7 Ghost quadratic constraints | Tasks 4, 5 | cell-average reproduction, QR threshold, donor/Halo and no-fallback failures | uniform/warped 1/2/4 plans |
| 8 Local Flow Pattern | Task 6 | independent paper row, singleton/joint groups and all mutation oracles | stable 1/2/3-link results |
| 9 residual/shared-flux/storage | Tasks 7, 8 | polynomial reconstruction, inactive-payload oracle, no-double-count/shared-face tests | allocation-free 1/2/4 operator |
| 10 PISO pressure Ghost | Task 9 | dynamic `rho_wall/D_wall`, revision, nullspace, exactly-two-corrector rollback | checkerboard/residual/MPI |
| 11 constant/material/ideal-gas paths | Tasks 9, 14, 15 | conservation, positivity, closure, wall-flux and contradiction failures | density-model 1/2/4 matrix |
| 12 traction/force | Tasks 10, 11 | analytic pressure/viscous/moment, ownership and transaction neutrality | separate second-order force/consistency rows |
| 13 WALE math/lagging | Tasks 12, 13, 16 | tensor oracle, exact zero, y-cubed, identity, lagging and retry mutations | TGV/channel/finite-body wake |
| 14 transaction/collective failure | Tasks 9, 11, 13--16 | authoritative bitwise helper and every injected phase failure | lowest-rank, retry and no-third-corrector evidence |
| 15 Checkpoint v3 | Task 17 | byte protocol, presence/mismatch/corruption and failed-read transaction | bitwise same-layout continuation |
| 16 diagnostics/performance | Tasks 18, 20 | enum/provider/presence/sample/read-only/counter tests | canonical 1/2/4 artifacts |
| 17 backend neutrality | Tasks 7, 12, 20, 21 | kernel-view/header/test-double/tests-off rejection | CPU-reference claim only |
| 18 scientific acceptance | Tasks 11, 13--16 | full MMS order, laminar engineering, LES and density matrices | all frozen thresholds and process grids |
| 19 six hard gates | Gate boundaries and Tasks 3, 6, 11, 13, 16, 21 | each gate remains closed until its whole-task verdict | matrix in §6 |
| 20 coordinator/TDD | §1 and every task brief | RED/GREEN, main dual review, exact-HEAD evidence identity | accepted-head/report ledger |
| 21 independence/public sources | Tasks 1, 20, 21 | path/string policy, provenance and workspace audit | no private-source access/publication/push |
| 22 completion/claim | Tasks 20, 21 | capability ledger and complete final commands | claim text in §7 and stop |

This table is bidirectional: every row must appear in the capability ledger,
and every new product symbol must point back to one of these approved rows
and to a named test. A row without both implementation and evidence is not
`implemented-and-accepted`.

---

## 6. Hard-gate Exit Matrix

| Gate | Required final-candidate evidence | Fast cases allowed only for iteration | Acceptance evidence |
|---|---|---|---|
| 1 | schema/header/source-policy complete Debug; surface/classification Release 1/2/4; parser ASan/UBSan; tests-off | small tetra/cube and `12^3` classification | strict STL/query/coverage/inside/outside/active-boundary at 1/2/4 |
| 2 | reconstruction/LFP complete Debug; 1/2/4 plan identity; ASan/UBSan; tests-off | 14-donor single-link polynomial | uniform/warped cell-average quadratic, remote donors, 1/2/3-link replacement groups |
| 3 | complete Debug; Release full `12/24/48` sequences; rollback 1/2/4; focused sanitizers; tests-off | `12/24` sphere and reduced engineering/decomposition cases | six shape/mapping plus two translated-sphere and one inside-cavity sequences with every velocity/pressure/penetration/force/consistency order `>=1.8`; §4.5 cylinder/sphere/transient baselines |
| 4 | complete Debug; WALE unit/MPI/sanitizers; Release TGV/channel; tests-off | `32^3` TGV, `32x33x16` channel | `64^3` TGV, `64x65x32` channel, 1/2/4, y³ and transaction |
| 5 | complete Debug; focused sanitizer; Release density/combined cases | 10-step `96x48x48` wake | 50-step `192x96x96` constant wake plus material/ideal integration and retry |
| 6 | complete Debug; I/O/diagnostics sanitizer; Release performance; tests-off/linkage | fast corruption and diagnostic samples | Checkpoint bitwise continuation, canonical providers, exact counters, driver and Stage 1/2/3 full regressions |

---

## 7. Final Capability Claim

Stage 3 完成后只能声明：

```text
HUNDUN-FLOW provides, on the CPU-reference backend:
- schema-v3 additive composition;
- one static closed connected STL;
- a sharp-interface Local-Flow-Pattern Ghost-Cell IBM that passes the
  approved second-order velocity, pressure, penetration and force contracts;
- pressure/viscous/total real-surface force integration;
- a WALE LES baseline passing the approved tensor, near-wall, body-fitted
  and finite-body composition gates;
- constant/material/ideal-gas density paths with collective transaction,
  same-layout Checkpoint v3, diagnostics and 1/2/4-rank evidence.
```

不得声明：

```text
cut-cell conservation
moving or multi-part IBM
general wall functions or thermochemical walls
DNS or universal experimental agreement
production GPU/GPU-aware MPI/mixed precision
chemistry/species/TPDF-TCR/spray/particles
rank-changing restart
Stage 4 readiness beyond the frozen backend-neutral interfaces
```

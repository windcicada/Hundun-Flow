# FGMRES 递推恢复原因观测

开发基线 `453eedfc223d4b72e435c787a88c77820a0e9fbd`。
本切片只补公开 FGMRES 调用的可选观测；不改变算法、容差、restart、压力方向、
方法历史或 checkpoint。ProductDriver/runner 尚未接线，不能用本报告给 Re3900
的恢复原因分配占比，也没有新的加速结论。两相继续暂挂。

## 定位依据与假设

[上一轮 MG 数据](2026-09-08-runner-mg-profile.md)的 9502–9510 窗口，C2-r1
有 291 次迭代/M、428 次 A，以及 58 次 `norm_breakdown_restarts`。
源码核查确认旧计数至少有两个来源，不能把 58 称为失败次数。

| 待核实假设 | 区分它所需的观测 | 当前结论 |
|---|---|---|
| 后续 Arnoldi 列的单归约范数递推为负，丢弃已经花费 A/M 的当前方向，增加工作 | 负范数事件、丢弃列、该分支重建真实残差的 A、真正继续的恢复次数 | 路径存在，已有受控公开接口回归；目标 Re3900 占比待测 |
| 显式 Arnoldi 残差为零，但真实残差尚未达标，需要 happy restart | 与负范数恢复分开的 happy restart | 无故障注入的单位算子回归实际触发；目标占比待测 |
| 常规 restart 或周期/提前真实残差审计贡献额外 A | 初始、Arnoldi、unsafe、周期内候选、周期结束五类 A | 外部算子调用计数可独立对账；目标占比待测 |

CodeGraphF 当前索引未命中内部 norm helper，本次以 `rg` 和源码继续核查。
只读确认 `unsafe_recurrence && column != 0` 保留先前有效子空间，并没有把当前列
加入三角回代；column=0 的负范数及零递推范数进入已有显式重正交化。
这说明存在潜在的额外工作来源，不证明该保护可以删除或方向可以直接复用。

## 接口与最小源码范围

`include/hundun/v04_linear.hpp` 为 `solve_fgmres()` 增加末尾可选借用指针
`FgmresRecoveryObservation*`，默认 `nullptr`；`LinearSolveInvocation` 和
`LinearSolveResult` 的载荷不增加。当前 Clang/x86-64 的独立 sink 为 96 字节。
只有调用者提供 sink 才写计数；它不进入字段、arena、持久证书、Restart 或 Evidence。
调用者保证 sink 在本次调用期间有效，且不与 solver 字段/workspace 重叠。

- 每次调用先清零，预检查通过后 `available=true`。它表示进入了有效求解，
  不是成功标志；失败时只报告已经执行的前缀，仍由返回 Status 决定是否接受解。
- sink 是否存在及地址允许 rank-local，不进入 collective fingerprint，
  不引入计时器、动态分配或新的通信边界。
- `unsafe_norms` 统计负递推范数事件；`discarded_columns` 统计未使用当前方向的事件；
  `explicit_reorthogonalizations` 统计显式重正交化尝试，不要求尝试最终成功。
- `unsafe_restarts + happy_restarts` 严格对应旧 `norm_breakdown_restarts`。
  收敛返回、无进展失败、迭代耗尽不算继续的 restart；普通长度 restart 单列。
- 五类 `*_applies` 统计 A 的尝试，包括失败调用。它们与返回结果中普通 A 总数相等；
  recycle projection 的 A 原本单独记账，不能再次计入这个和。

`src/solver_krylov.cpp` 的差异只在上述既有事件处更新整数计数。
没有变更浮点表达式、reduction 参数、workspace 访问/修订或异常返回条件。
现有迭代上限为 uint32，单次调用事件数远小于 uint64 上限；无需为了观测改变求解
失败路径。默认关闭仍有参数/分支检查，未测产品开销，不能宣称成本绝对为零。

## 实际开发回归

入口沿用已确认的 `solve_fgmres()` 公共数值接口。新增回归没有调用 norm 私有函数：
第一组沿用已有故障注入驱动真实求解，检查独立三对角残差及外部 A 回调；
第二组使用正常/失败公开算子、真实 ReductionEngine 和 caller solution。

| 验证 | 实际结果 |
|---|---|
| 声明 + 空观测的 RED | `v04_solver_krylov_mpi_1` 在原因/A 归因断言失败，2.33 s；不是数值不收敛 |
| 加计数后的 GREEN | 1/2/4 ranks，3/3，0.98 s |
| 生命周期扩展 | 1/2/4 ranks，3/3，0.97 s；第一列显式恢复、后续丢列恢复、普通长度 restart、周期内/周期末 A 审计 |
| 观测等价及资源 | 所有 rank、仅 root、仅末 rank 开启：caller solution 逐值相同，非计时结果/外部 A/M/归约工作完全相同；观测开启/关闭的热路径 C++ 分配计数均为 0 |
| 失败与复用 | 非零 rank A 失败全体返回，已接受 caller solution 不变；只报告已尝试的初始 A 和首个 Arnoldi A；零 RHS 和单 rank 无效控制入口清空旧观测 |
| 非注入 happy restart | 非均匀 RHS、A=I、极小公开绝对容差：1/2/4 ranks 均为 2 iterations、1 happy restart、0 unsafe restart、5 A；仍由独立真实残差或既有失败回退合同验收 |
| ASan + UBSan | Clang Debug，1/2/4 ranks，3/3，1.84 s；sink 实测 96 字节 |

原始 stdout 见 [数据目录](data/2026-09-08-fgmres-recovery-observation/)。
Sanitizer 使用 `detect_leaks=0` 避免系统 MPI 保留对象干扰，这不是全进程泄漏/RSS 验收。
开发 Ninja 仍提示既有日志末尾不完整并重建；开发产物不作为干净发布身份。

```sh
cmake --build build-review-fixes --target v04_solver_krylov_mpi_test -j4
ctest --test-dir build-review-fixes -R '^v04_solver_krylov_mpi_[124]$' --verbose -j1
cmake --build build-review-sanitize-clang --target v04_solver_krylov_mpi_test -j4
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ctest --test-dir build-review-sanitize-clang -R '^v04_solver_krylov_mpi_[124]$' --verbose -j1
```

共同环境：`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`、
`LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu`。
所有 MPI 测试串行，原长测 128 ranks 保持 SIGSTOP，没有新 128-rank 窗口或长测恢复。

## 下一边界

先对本切片 DCO 提交做独立干净构建及相关 Krylov/PISO 调用方验收，再接入
ProductDriver 的逐 loop 出口和显式 runner 观测。后续接线必须保留失败尝试、
完整 rank/step/attempt/sweep 身份及 64-loop 溢出拒绝，默认关闭且无热路径分配。
未启用不能伪造为有效零计数；旧 Evidence/CSV 不回填原因字段。

取得修正版目标窗口的实际原因计数后，再决定是否对某个恢复分支做唯一最小实验。
本报告没有批准延长 restart、提前 spatial、放宽容差或移除保护。
COAST 同工况对标、基础流动替代和后续燃烧仍由[总目标台账](../plans/2026-09-08-coast-replacement-and-stages.md)跟踪。

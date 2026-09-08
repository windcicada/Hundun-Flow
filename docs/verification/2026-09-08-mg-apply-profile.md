# MG 内部成本观测

开发基线 `da56300a9040211b2b50851470b1f2f8c984fb2d`。
本切片进行中，尚未接入生产 runner，也没有新的逐层生产测量或加速结论。
依据是 [9501–9510 单轮观测](2026-09-08-linear-criterion-observation.md)：
C2 r1 工作更多，但单次 M 成本与其他 diagonal 轮接近。
要区分细层平滑、粗层求解、层间传递及通信，先补成本来源。

## 接口与计时口径

`NativeCartesianMgPlan::set_apply_profiling()` 在静止调用边界启停并重置本地
观测 epoch；`apply_profile()` 返回 plan-owned 借用快照。默认关闭，允许各 rank
开关不同。数据固定为至多 32 层，不写入 checkpoint、方法历史或数值证书。
成功重新 compile 替换计划、move 或销毁后，旧借用引用失效。

- `attempts/successes/failures` 区分全部 apply 尝试与成功提交的 correction。
  `prepare_batch` 不算 apply；prepared 入口检查和退出归入整个 apply。
- pre/post smooth、独立 residual、restriction、prolongation、terminal 六项
  不包含递归子调用，时间互斥。pre smooth 已保留的 defect 不再重复计 residual。
  初始化、最终 finest residual/缩放/发布和入口检查留在 apply 总时间中。
- Halo 等待/控制、ReductionEngine 和直接 MPI 是嵌套子项，不能与六项或
  apply 总时间再相加。借用 coarse halo 的通信归入该 halo 的层。
- 计数饱和或无法解释的通信身份使观测 `complete=false`，不改变求解 Status。
  有效服务上的数值失败可以有完整成本记录；不能把失败尝试全部丢弃。
- 关闭时不创建新计时器或采通信快照；固定 profile 本身仍占少量存储。
  不新增全局通信检查，不改变 prepared epoch 的失败汇总与发布顺序。

## 已实际运行的开发回归

公开入口为 NativeCartesianMgPlan，沿用已确认的公开数值接口测试范围。
所有 MPI 测试由根代理串行启动，原长测 128 ranks 保持 SIGSTOP。

| 切片 | RED | GREEN |
|---|---|---|
| 直接 V-cycle | 仅有声明/空观测时，新观测断言失败，1 rank，2.41 s | 1 rank，0.42 s；阶段次数、互斥计时、校正场/初终残差逐位不变、C++ 分配数 0 |
| prepared/失败 | 尚未包住 prepared 入口时，成功及失败累计断言失败，2 ranks，2.44 s | 2 ranks，0.41 s；rank 0 单独开启、非零 rank 非法 ticket、既有通信故障与恢复；原通信次数和输出回退保持 |
| 通信细项 | 空通信计数使逐层 Halo/ReductionEngine 对账及 replicated MPI 断言失败，2 ranks，2.46 s | 1/2/4 ranks，3/3，1.30 s；逐层差值与公共计数器精确相等，通信作为嵌套时间 |

原始 stdout 位于 [本切片数据目录](data/2026-09-08-mg-apply-profile/)。
开发构建曾提示 Ninja 日志末尾不完整并恢复；后续接受必须来自干净构建，
这些增量构建不充当发布身份。通信实现初次编译把局部参数 `level_index` 误写为
`level`，编译器拒绝后已纠正；这不是运行期数值失败。
扩大验收后的 Release 组 16/16 通过（7.67 s）；覆盖 2–7 层 V/F 的独立次数、
4-rank 延拓直接 MPI、prepared 计时排除 prepare_batch、生命周期和调用方指针表、
借用对象失效、数值失败、MG 更新合同/奇数分区及 Krylov 隔离。
NaN 数值失败和无效输入在有效服务上仍有 complete 观测；服务换绑失败则为 incomplete。
Move 的测试要求曾被委派语句混淆，核对 ownership 后明确为转移 epoch，
不是 move 后清空；未将这次说明澄清记录为数值实现缺陷。

Clang Debug ASan+UBSan 的 reuse 与 MG MPI 1/2/4 共 4/4 通过（15.05 s）。
`detect_leaks=0` 用于排除系统 MPI 保留对象，不能据此宣称完整泄漏检查或 RSS 硬预算通过。
扩大测试保留原科学工作量，没有为了观测改变平滑次数、残差门槛或通信方案。

开发验收命令（该 checkout 的两个既有构建目录）：

```sh
ctest --test-dir build-review-fixes -R '^(v04_solver_mg_.*|v04_solver_krylov_mpi_[124])$' --output-on-failure -j1
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ctest --test-dir build-review-sanitize-clang -R '^v04_solver_mg_(reuse|mpi_[124])$' --output-on-failure -j1
```

共同环境：`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`，
`LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu`。
计数饱和与重复 Halo 冷归属防护经过静态核查，未借私有变异接口伪造饱和回归；
不将这些分支记为已实际触发。库级测试不证明生产规模额外开销为零。

通信采样只接受冻结 instance identity；不要求失败后的 Halo 仍为 ready，
因此同身份 poisoned 对象可以保留失败成本。计数倒退或输入已饱和时标记不完整。
setter 在冷端排除重复 Halo 归属，既不重复计费，也不因此改变数值求解合同。

## 后续接线和接受条件

Native 层验收后，再将 compact 每-loop 摘要接到现有 CaptureSolve 出口；
完整逐层快照只在独立 step sidecar 保存，避免把 32 层数组复制进 64-loop 载荷。
runner 使用显式 opt-in，并在冷入口冻结输出开关；它不同于 Native 允许的
rank-local 开关。覆盖失败步、缺 rank/step/level、重复、首错及 flush/close。
observer 必须能把完整 step 层账与 loop 摘要对账，未启用不冒充有效零成本。

完成公开 1/2/4-rank、prepared/direct、V/F、通信、更新/重建、失败与资源回归后，
从干净候选运行同起点的一轮 Re3900 窗口。只有实际热点和数值不变证据成立后，
才选择一个最小性能实验。基础流动替代验收、COAST 对照和后续燃烧仍属
[持续目标](../plans/2026-09-08-coast-replacement-and-stages.md)，没有缩减为观测功能验收。
两相 Stage 6 按用户 2026-09-08 最新指示暂挂，不在当前接入范围内。

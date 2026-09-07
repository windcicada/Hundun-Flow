# 真实线性停止阈值：观测回归

基线 `86542bb96678ec844ae5ac95d8f6391993da239e`。本切片只增加观测，
不改变压力—焓方法、物性、初猜保护、收敛控制或 checkpoint 合同。
本文件记录局部验收；生产规模窗口及 COAST 比较仍需独立证据。

## 缺口与最小修改

旧 loop CSV 的 `linear_initial` 是 `||b-A*x0||`，不能据此计算实际
停止阈值 `max(atol, rtol*||b||)`。热初猜及已有初猜保护还可能改变这个初始残差。
已有 `LinearSolveResult::convergence_limit` 属于附加物理审计，不覆盖它。

- `solver_krylov.cpp::convergence_limit()` 复用三种求解器已经算出的全局 RHS 范数，
  在零 RHS 快速路径前保存 canonical criterion；没有新增归约、场遍历或动态分配。
- `LinearSolveResult` 新增可用标志、RHS 范数、atol、rtol、真实残差上限五个定长字段。
  preflight、RHS 范数失败或阈值溢出时保持不可用及全零；有效零 RHS 可有零阈值。
- runner 使用 17 位精度写入同名 CSV 观测列，`RUN.meta` 的 observation schema 升为 4。
  observer 保留 V3 输入支持，明确旧行的 criterion 不可用；V4 要求完整字段、
  有限非负数、公式一致、独立来源/步/rank 覆盖与跨 rank 一致性。
- 公共求解器允许 atol 或 rtol 单独为零，以及 rtol≥1；观测沿用该合同。
  双零容差仍非法。真实零阈值必须精确为零，不能用舍入余量接受非零值。
- 失败求解也能有合法停止阈值；observer 不以最终残差高于该阈值为理由丢弃故障证据。

原 Evidence V8 显式编码和哈希不包含这些新增字段；没有对结构体原始字节序列化。
它们不改变存储状态、速率、面通量或时间语义，因此不改方法历史签名。
64-loop 固定容量仍保留，溢出仍拒绝完整归因；本切片没有解决分段输出和总峰值预算。

## RED → GREEN 证据

公开入口为 PCG / FGMRES / BiCGStab、真实 runner CLI、observer CLI。

| 验证 | RED / 修正 | 最终结果 |
|---|---|---|
| RHS 与初猜分离 | 仅加声明与回归时，旧求解器不发布 criterion，新断言失败 | 67 单元独立 `A=I,b=2`；零/暖/精确初猜、相对项/绝对下限及零 RHS，1/2/4 ranks 通过 |
| 失败与公共范围 | 扩展合法容差、溢出、preflight、单 rank NaN RHS 和调用方输出不变 | 三算法均通过；附加物理 audit 拒绝后仍保留独立 canonical criterion；guard 统计无新增 C++ 分配 |
| 测试夹具纠正 | 初次 GREEN 中 BiCGStab 被夹具错误的 `fixed_spd` 声明拒绝 | 改用其既有 `fixed_general` 合同；不是生产求解器缺陷 |
| 实际 CSV | 旧 runner 在真实 2-rank CLI 的 schema=4 断言处失败 | 非零统计样本：同方法→方法恢复→新 epoch 精确续算，1/2/4 ranks 通过，源附件哈希保持 |
| observer 边界 | 初版校验错拒合法纯绝对容差；舍入余量曾误接受 expected=0、limit=5e-324 | 两个局部 RED→GREEN 后，88 项 CLI 检查通过；旧 V3、V4、非法阈值、跨 rank、partial 与流式路径覆盖 |

Release 最终选定组 7/7 通过（40.57 s），随后只重查最终 observer 版本，88 checks 通过。
Clang Debug ASan+UBSan 的 Krylov 1/2/4-rank 组 3/3 通过（1.82 s）。
关闭 LeakSanitizer 以避免系统 MPI 保留对象干扰；这不是进程 RSS 或完整内存预算验收。
observer 的 8/2048 步 details-output 测试 Python traced peak 为 2,806,838 / 3,839,510 bytes，
仅表示该合成流式回归的分配规模；默认保留全部 loop 对象的汇总模式不宣称常量内存。

命令（工作目录为本开发 checkout）：

```sh
cmake --build build-review-fixes --target v04_solver_krylov_mpi_test v04_thin_domain_runner -j4
ctest --test-dir build-review-fixes -R '^(v04_solver_krylov_mpi_[124]|v04_runner_statistics_epoch_mpi_[124]|v04_solver_observation_cli)$' --output-on-failure -j1
ctest --test-dir build-review-fixes -R '^v04_solver_observation_cli$' -V
cmake --build build-review-sanitize-clang --target v04_solver_krylov_mpi_test -j4
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ctest --test-dir build-review-sanitize-clang -R '^v04_solver_krylov_mpi_[124]$' --output-on-failure -j1
```

MPI 测试使用 `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`，
`LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu`。
原 128-rank 长测全体 SIGSTOP，所有 MPI 检查串行，未恢复或替换其程序。

原始日志见 [本切片数据目录](data/2026-09-08-linear-criterion/)。两项 observer
开发期 RED 仅保存在代理工具 stdout 中；最终 88-check CLI 日志已独立落盘，
不把不存在的 RED 文件记为完整归档。

## 后续验证边界

本地回归通过不等于新性能排名。接下来从只读 checkpoint，以同方法精确续算，
在新目录运行一次有界窗口，验证 V4 完整来源和数量，并记录实际阈值与 A/M 成本。
当前 C2 refinement 已有比零初猜差则回退的保护；不能重复实现该功能，
也不能从覆盖后的 `linear_initial` 反推原始初猜质量。
数值优化需由新测量定位后单独实验、单独提交；本提交不包含速度提升结论。

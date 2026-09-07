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

## 干净候选与单轮生产观测

观测提交 `16fd8f085e845628a545dc7d4d371b3064e9ef38`，
tree `8a53030af5f42f63dff6dac7fd3ce134aaa8e0c9`。
独立干净 checkout `hundun-flow-criterion-accept-20260908` 使用与原冻结程序相同的
Clang 15/libc++、Release、`-march=znver3 -mno-fma -ffp-contract=off`，
7/7 选定回归通过（43.99 s）。后续文档提交不是新的已测二进制身份。

- runner SHA-256：`8cfc9fe3a613a84e677fe536a353779040b42bf525650ae10373616e23d33545`。
- build manifest SHA-256：`d235d3a3974eec4fd88ce62a03ff24ff3ed6fbc5088209c80652dae8b0140f3f`。
- 只读源为原长测 `generation-9500-139424060480232`，manifest
  `ed4600dc5068bb3118acd14aeefffcfde527a47383c23c3c36fee709f47aba78`。
- 新目录 `pilot-criterion-9500-9510-20260908`，128 ranks，10 步，仅一轮；
  原网格、变物性、固定 dt、门槛、refinement 及持久化设置未改。
  没有 method-recovery；来源/目标方法签名均为 `12213963202598979269`。
- 10 步 BDF2，无重试/恢复/降阶；统计 epoch 保留 7000、采样起点 17001、样本数 0。
  `conservation.csv` 的积分账本按既有公共合同另从重启点 9500 开始，不与统计 epoch 混同。
  复用旧检查器时一度误将两种 epoch 都写为 7000，纠正检查器后通过，生产实现未改。
- 全窗口连续性/能量归一化最大值分别为 `2.61757e-7 / 8.80580e-7`；
  流体温度范围 `299.97904–300.02367 K`，固体占位区域极值及极值位置没有漂移。
  本观测改动未重复全固体 checkpoint MPI 探针，不把区域极值检查写成全载荷逐单元比较。
- V4 observer 验证 10 步、70 logical loops、每 loop 128 ranks，来源和覆盖完整；
  原 Evidence 的 runtime 校验通过。末次 flush/close 成功，程序正常退出。
- 9510 Visit 的 128 块引用集合、LittleEndian 二进制长度、有限字段、单调坐标、
  无重叠分区及总计 6,070,272 单元通过检查。末次 checkpoint 有 128 个非空 rank 文件、
  manifest、完成标记、accumulator/statistics；代次 `generation-9510-146520080763070`。
  manifest SHA-256 为 `f8bb868eb6c006b06a5b1e1f6ff9718071e0a65b5762ff65e139628806e88ccc`。
- 133 个源 checkpoint/统计文件及冻结输入/程序哈希复核不变。原长测仍为全体 SIGSTOP，
  短测期间仅候选的 128 ranks 活跃，未替换原长测或从较早的 pilot 点重启它。

### 新成本数据及口径

进程总墙钟 111.45 s，含启动、恢复与末次输出；不能直接除以十当作纯求解成本。
重启会按既有合同清除非持久的 warm authority，首步单列；
旧长测的相同步号因此不是相同求解初猜基线。本轮没有算法 A/B 或 COAST 速度比较。

| 窗口 | mean(max-rank full) | mean(max-rank advance) |
|---|---:|---:|
| 首步 9501 | 9.774696 s | 9.762793 s |
| 9502–9510，9 步 | 10.845357 s | 9.562760 s |

9510 含最终输出，full=20.903768 s、advance=9.452607 s；该步 rank-mean
Visit/checkpoint 分别 2.159710/9.279872 s。不能将短窗的一次末次输出频率外推为长测均摊成本。

后 9 步每步 rank-mean：Krylov 4.014488 s，内部 A/M 为 1.389056/2.015330 s；
候选装配 2.278196 s；MG refill/copy 为 0.075388/0.011184 s，copy 包含于 refill；
Schur prepare 0.227973 s，最终动量 0.201856 s。嵌套子项不得重复加到总墙钟。

| 求解类别（后 9 步） | loops | 平均迭代 | solve ms/loop | A/M ms/apply | breakdown/loop |
|---|---:|---:|---:|---:|---:|
| C1 spatial | 9 | 22.667 | 944.296 | 17.263 / 11.394 | 5.000 |
| C2 r0 diagonal | 9 | 20.111 | 394.856 | 3.534 / 11.490 | 4.667 |
| r1 diagonal | 9 | 32.333 | 650.022 | 3.559 / 11.313 | 6.444 |
| r2 diagonal | 9 | 29.111 | 561.479 | 3.554 / 11.060 | 6.333 |
| r3 diagonal | 9 | 26.889 | 534.483 | 3.567 / 11.287 | 5.667 |
| r4 diagonal | 9 | 27.556 | 535.345 | 3.564 / 11.122 | 6.111 |
| r5 diagonal | 8 | 20.500 | 401.784 | 3.566 / 11.345 | 4.375 |
| r6 diagonal | 1 | 17.000 | 331.797 | 3.594 / 11.398 | 4.000 |

63 logical loops、1609 迭代、2358 A/1609 M；没有 C2 spatial 样本。
现有 Evidence JSON 与 CSV 按 step/corrector/refinement ordinal 对应，
工作次数和初末残差精确匹配。breakdown 为既有递推恢复计数，不是物理发散次数。

实际停止阈值带来了两点此前不能从 CSV 直接证实的结论：

- r1 的 rtol=`1e-6`，`||b||=0.014628–0.015590`，门槛为 `1.4628e-8–1.5590e-8`。
- r2–r6 的全部 36 loops 已由 atol=`1e-8` 主导，单看 rtol 大小会误判实际松紧。

C2 r0 的连续性/能量平均收缩比为 `0.009466 / 1.094405`，r1 为
`1.187417 / 0.161684`；分量不逐次同时单调收缩，不据此单独判为实现缺陷。
所有 refinement 的记录初残都等于 RHS 范数；因为已有保护会覆盖初残，
这不能证明保护触发过多少次。

### 数据身份与后处理限制

外部完整审计目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/criterion-observation-20260908`。
仓库内保留[紧凑验收收据](data/2026-09-08-linear-criterion/pilot-acceptance.json)。
`observation.json` / `loops.jsonl` 来自本候选的正式 CLI；VTK/区域检查使用封存的
`check-pilot-output.py`，成本复算使用事后整理并实际重算的 `summarize-cost.py`。

首份一次性 `cost-summary.json` 经 JS 中转时舍入了 52 个 uint64 身份字段，
不是求解器或原始 Evidence 的问题。原文件保留只读，不作为身份权威；
新的 `cost-summary-reproduced.json` 由原生 Python JSON 写出，记录逐项纠正清单。
全部成本、计数及其他既有字段一致，source history signature 精确保留；
复算脚本 hash `238e0d491c8ed985c4591900d04ee04bb9f185e4426b84ebfafa492204729f99`，
复算产物 hash `3ce62a7a7ea7dc361069e3dad0e68cd076d91e87f5b92f5140a3b48a1feb9935`。

## 下一切片与验收边界

本地回归和单轮观测通过不等于新性能排名。下一切片只增加默认关闭的
M 内部 MG 层级归因，区分平滑、残差、传递、粗解和通信，再据新证据选择一个优化。
当前 C2 refinement 已有比零初猜差则回退的保护；不能重复实现该功能，
也不能从覆盖后的 `linear_initial` 反推原始初猜质量。
数值优化需由新测量定位后单独实验、单独提交；本提交不包含速度提升结论。

# Re3900 审阅核对、首批修复与后续调试

日期：2026-09-05。审阅基线：`62ec29b359d78d6b3bff9af00bf48a930a0d6673`。

本批修复输出阶段的 MPI 异常一致性、Visit 大载荷双缓冲，以及六处标量目录打包错配。没有修改时间步、物性、压力求解算法、SIMPLE refinement 上限或物理残差容限。这些修复不构成当前 Re3900 长测通过的证据。

## 1. 当前长测的停止原因

核对的实际文件为：

```text
/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/case/
/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/long-simple-35000-iohotfix-62ec29b-from500.launch.log
```

该配置为固定时间步 `dt = 1.3808912271980336e-5 s`，`transported_scalars = []`。物性、COAST 原生网格、456 × 256 × 52 网格、128 rank 和 Re3900 参数保持不变。因此，自适应时间尺度输入缺陷和混排标量错配都不能解释这次长测的停止。

服务 `hundun-long-re3900-iohotfix-20260904.service` 的现场状态为 `failed / exit-code / 7`。日志给出的直接证据是：

| 项目 | 记录 |
| --- | --- |
| 失败目标步 | 703 |
| Status / failed_stage | `6/10210` / `54` |
| 尝试数 | 1 |
| C2 refinement 次数 | 12 |
| 最后选中步长系数 | 2 |
| 最后连续性残差 | `4.90634e-7` |
| 最后能量残差 | `1.22875e-6`，高于 `1e-6` |
| 本次失败尝试最大 rank 耗时 | `14.722812894 s` |
| stage 50 最大 rank 耗时 | `12.645943496 s` |

这是压力—能量终止检查未通过的受控退出。不能据此认定 OOM、段错误或输出失败；也不能仅用最后一个 refinement 样本确定收敛慢的更深层原因。不同阶段的最大 rank 时间可能来自不同 rank，不能当作同一 rank 的严格可加耗时分解。

## 2. 本批修改与验证

### 标量映射

`core_product_freeze.cpp` 的冷启动依赖重建、预测后交换、C2 IBM donor 刷新、最终非对流速率及其 IBM donor 路径，统一按原始目录打包。组分和被动标量仍保留各自的紧凑工作数组；没有增加新的运行时映射容器。低阶 predictor slow-path 使用自己的组分优先协议，未改动其打包顺序。

修复前，实际产品路径对 `[tracer, sp0, sp1]` 和 `[sp0, tracer, sp1]` 出现 `2/736`（stage 40）或 IBM 初始化 `2/14504`。修复后，三种目录排列在普通开边界和 IBM 算例中均接受两步；在各自的同一分区下，所有输出字段按名称比较完全一致，FieldId 在推进前后保持不变。测试分别覆盖 1、2、4 rank；这不是长时间湍流逐点跨分区一致性测试。

### 输出异常与 MPI

Visit、Evidence、Screen、Monitor 和共享目录创建采用“本地工作及异常转换 → 全体 Status 汇总 → 下一阶段”的顺序。Evidence 编码成功后才进入全局哈希比较；正常路径保留全 rank 编码与一致性检查。字符串流开启失败异常，避免内存不足后发布截断文本。

故障注入只存在于测试可执行文件的分配器替身，不增加生产程序环境开关。对每条正常输出路径测得的分配位置逐一注入一次失败，分别选择最后一个 rank 和 rank 0。验证全体 Status 一致且为失败、不发布 `.visit` 索引或日志记录，并能随后成功执行正常写入。

修复前，2 rank 在首个注入点出现 `MPI_Allreduce / MPI_ERR_OTHER`，退出码 16。修复后，四种 writer 的 2、4 rank 测试全部通过。部分标准库字符串流将 `bad_alloc` 转成 I/O 异常，现有身份校验也可能返回 `invalid_plan`；本批保证一致失败，尚未统一这些错误的精细分类。

这只覆盖 writer 内部的本地异常边界。调用者在进入 writer 前的路径/文本构造、产品编译中的 arena 分配，以及其他 MPI 阶段，仍需各自的故障注入验收，不能宣称整个应用已具备完整的内存不足恢复能力。

### Visit 内存

先计算偏移和载荷长度，再一次分配完整文件缓冲区，直接编码字段和坐标。仍显式输出 LittleEndian；删除了完整 appended 副本和逐字节 `push_back` 扩容。格式化流在大载荷分配前销毁，文件缓冲区在生成索引前释放。大缓冲区分配前的容量检查包含同时存活的 prefix 和偏移数组，不再只检查最终文件长度。

下表是原有小型 I/O fixture 的新增 C++ 分配峰值，不是 Re3900 作业的 RSS，也不是端到端性能测量。两边输出文件逐字节一致。

| 输出 | 文件字节数 | 修复前新增分配峰值 / B | 修复后新增分配峰值 / B |
| --- | ---: | ---: | ---: |
| 均匀网格 VTI | 4,387 | 9,779 | 4,915 |
| 拉伸网格 VTR | 42,997 | 110,069 | 43,765 |

同一测试构建中的文件 FNV-1a 校验值分别为 `16284654938302156485`、`4264152898443382635`。另测最终文件大小恰好等于 staging 预算时应拒绝额外元数据开销，以及含 NaN 的字段不得发布输出。格式化元数据、分配器开销和求解器常驻内存的完整进程预算仍属后续工作；本批未实现跨输出步复用 staging。

### 可运行检查

使用现有 Release / Clang / libc++ 检查构建，没有启用 fast-math。以下 14 项检查单轮通过，总墙钟 13.34 s；其中包括 3 个标量排列 MPI 测试、8 个输出异常 MPI 测试、普通 I/O 检查、产品 freeze 检查和 thin runner 自检。

```bash
cd /home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity
cmake --build build-re3900-fix-test --target \
  v04_io_product_path_test v04_core_product_freeze_mpi_test \
  v04_thin_domain_runner v04_core_product_freeze_test -j2
LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu \
ctest --test-dir build-re3900-fix-test --output-on-failure -j4 \
  -R 'v04_(scalar_catalog|io_.*allocation_failure|io_product_path$|core_product_freeze$|thin_domain_runner_self_test)'
```

生产库和测试库均构建成功，`git diff --check` 通过。本批未跑全仓库测试、ASan/UBSan 全量构建或新的 Re3900 长测。

## 3. 后续仍按两部分推进

| 部分 | 下一项工作 | 验收依据 |
| --- | --- | --- |
| 直接影响当前 Re3900 | 定位第 703 步 C2 能量收敛不足：在临近故障的检查点重放中记录完整 C1/C2/refinement 轨迹、真实线性残差、候选数、Jacobian scope 与独立物理残差 | 区分方向近似不足、线性解精度不足及残差/线性化不一致；不以增大迭代上限或放宽容限代替诊断 |
| 直接影响当前 Re3900 | 对确认的根因做局部修复；从保留的 checkpoint 进行一轮有目标的跨故障时刻检查 | 在原网格、物性、边界、固定 dt 和 `1e-6` 物理容限下跨过原停止点，输出/checkpoint 成功，随后才据此恢复长测 |
| 应用、记录与通用正确性 | 在每次提交后立即更新已接受步数/时间；区分 numerical、Visit、Restart、resource、Evidence 失败 | 输出失败时仍报告已经完成的步号，不把文件失败写成数值求解失败 |
| 应用、记录与通用正确性 | 给普通自适应 CLI 接入最近接受状态的离散时间尺度；CFL 缩放只做一次，按显式项选择扩散限制，保留事后 CFL 审计 | 静止流、速度变化、非均匀网格、一次拒绝后的增长冷却及跨 rank 测试；固定 dt 路径保持原合同 |
| 应用、记录与通用正确性 | 修复 `valid_boundary()` 分配与 `noexcept` 的冲突；逐阶段补产品编译/应用调用层异常汇总 | 单 rank 注入失败，全体一致退出，无 collective/free/finalize 挂起 |
| 记录与性能分析 | 普通 app 复用节点通信器；补计算、输出、完整步墙钟和 halo 控制归约计数；明确采样和持久化周期 | 总墙钟能与分项对账；不把 `--output-interval 0 --restart-interval 0` 写成完全无 I/O 基线 |
| 后续内核优化 | 按新测量再选择 halo 控制通信、multidot、MG 或 workspace 生命周期复用 | 同物理推进时间比较 HUNDUN PISO、SIMPLE 和普通可压缩 COAST；每候选单轮，结果与未完成状态写入测试表 |

完整内存预算、明确初始状态、BDF1/BDF2 实际使用计数、压力—能量方向导数测试、发布策略/运行证据合同仍未完成。这里的计划不是已实现功能或通过的测试结论。

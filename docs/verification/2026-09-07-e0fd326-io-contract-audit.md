# e0fd326 之后的 I/O 与入口合同核查

基线：`e0fd326d3ed6cccd52f812c113712f34a08cec95`，实际工作区
`/home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity`。
原工作区在本轮开始时干净。未改物性、网格、时间步、残差门槛、
refinement 容量或 checkpoint 发布/持久化协议。

## 1. 运行器日志结束状态：实际复现并修正

`tools/v04_thin_domain_runner.cpp::run()` 的原结束路径没有检查四个
root 流的 close 状态；performance 和逐 rank 的 solver 日志依赖析构关闭。
真实两 rank CLI 在持久化 step 1 之后，向 root `force.csv` 的 libc
`fclose` 注入 EIO：原程序两个 rank 均返回 0，仍打印 `COMPLETED`。
这是实际注入复现，不是一次真实磁盘故障记录。

最小修改：`complete_logs()` 显式处理所有实际打开的流，逐流 flush/close，
保存首个本地错误，再选出失败 rank 并广播错误。disabled 流不参与。
运行主体提前返回也先经过相同清理边界；原求解/输出失败退出码优先，
不因清理失败变成成功，也不改写已接受物理状态。只有全体关闭成功，
才完成日志只读化与目录同步并打印 `COMPLETED`。

“日志 flush/close 成功”只表示缓冲输出和关闭成功，不声称普通 CSV
已经 fsync 持久化。checkpoint 的文件同步、目录同步和 current 发布合同不变。

回归入口：`v04_runner_log_completion_mpi_2`。39 次独立 CLI 检查通过：
普通/观测开关、六类文件的 EIO/ENOSPC、末次显式刷新、关闭、
先 EIO 后 ENOSPC 的首错保留；含 rank 1 的 solver 文件。
31 个注入失败用例全 rank 返回 6，不打印完成；8 个成功/未启用用例正常完成。
每个用例的持久化 rank checkpoint 字节 SHA-256 与正常基线完全相同。

实际证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/module-io-audit-20260907/`。
`close-red-checkout/results.json`、`force-close-eio.log` 记录 RED；
`close-green/results.json` 和每次 CLI 日志记录 GREEN。

## 2. Restart 读取上限：实际复现并修正

原 `read_file()` 在核对清单长度前按 `fstat.st_size` resize。
四 rank 回归将独立 rank 3 文件扩成 1 TiB 稀疏文件；测试分配器在
16 MiB 以上直接拒绝并计数，避免实际耗尽内存。旧代码返回
`allocation_failure/10301`，rank 3 记录一次超大分配尝试。
修正后全 rank 返回 `io_failure/10306`，超大分配尝试为零。

`io_restart.cpp` 的修改边界：

- `read_file` 在 resize 前检查普通文件、可表示长度、期望长度和容量。
  非阻塞打开使 FIFO 能在 fstat 阶段拒绝；读满后再检查 EOF，拒绝读期间追加。
  Writer 回读也使用清单期望长度和已预留的最大 rank 容量。
- `current` 在读取前限制为 256 B。manifest 有独立字节和 source-rank 数上限；
  解析 rank 目录前核对 40 B/record 与文件剩余长度，不能由恶意计数触发大分配。
  现有 64 个状态/64 个速率目录上限保留。
- `plan_read_bulk` 在 rank 读取和新 image 分配前计算独立 bulk 预算。
  包含旧 out 的实际容量、清单、原始 rank 字节、解析块、新 image 和覆盖数组。
  原始 manifest 完成完整性检查、解析和 SHA-256 后释放；恢复遍历重新核对
  rank 文件长度/哈希/公共头，与完整来源扫描保持一致。
- `RestartReadLimits` 是显式、可配置、rank-local 容量，任一 rank 超限时统一失败。
  默认 bulk=1 GiB、manifest=8 MiB、来源 ranks=65536；零值非法。
  `RestartReadReport` 区分 retained/new image 和 peak bulk 上界。

实际四 rank `v04_io_restart_mpi_test` 全套通过，日志为证据目录下
`restart-green.log`：V1/V2/V3、1→2/2→4/4→1/4→4 恢复、原有持久化和保留代次
测试、单 rank 系统调用故障，新增稀疏超大 current/manifest/rank、错误清单长度、
恶意计数、FIFO、空 rank 文件和 directory-as-current；全部在大分配前拒绝。
Writer 回读的超大 fstat 长度注入也被拒绝，current 仍指向前一个持久化步骤。
预算等于报告上界时成功；仅 rank 3 减少 1 B 时全体拒绝，旧 out 的完整字段/通量
载荷与存储地址保持不变。被修改的是独立测试副本，原测试 checkpoint 仍可读取。

预算边界必须说明：这是读取函数的 bulk 分配规划上界，不是进程 RSS 或全产品
硬预算；不包括小路径字符串、分配器/MPI 内部资源、ProductDriver 恢复物性暂存、
运行 arena 外数组和 halo。这些阶段仍须在后续总峰值模型中按同时存活关系加入，
不能把这里的 reader 预算当成完整运行内存验收。

## 3. MPI 公共接口：合同核查与防御性修正

`NativeCartesianMgPlan::update_coefficients()` 的可选 `counters` 是本地统计接收器，
不属于算子/层级的全局身份。原代码确实在 `counters != nullptr && replicated_coarse`
内调用 consensus，可能使不同 rank 的 collective 顺序错位。
现在只让指针控制本地计数运算；replicated coarse 分支中的结果汇总由所有 rank 参加。
原有生命周期合同、系数检查、层级发布及 prepared epoch 汇总均保留。
此项旧路径属于静态确认的通信风险；未故意提交会挂起的旧 MPI 作业，也未用 timeout
充当复现证据。

公共 MG 回归在真实 replicated coarse 路径中分别让 rank 0、最后一个 rank 传空指针；
所有 rank 成功刷新且仍能 apply，内部刷新次数正确、外部可选计数仅写入非空接收器，
层级持久地址不变。2/4 ranks 两项通过，见 `mg-contracts-green.log`。
最初测试夹具走了各向异性线松弛，路径断言正确报失败；改用等距夹具后确认进入
replicated coarse，并非放宽路径断言使测试通过。

普通 `ApplicationService::run()` 现在在冷入口比较七个控制量：steps、两个输出周期、
是否 restart、两种恢复策略、是否显式初态。有效但不同的控制量一致返回
`invalid_case/10506`，不会先编译产品、生成输出或推进。
路径需要指向同一逻辑共享资源，但挂载路径拼写、内存地址、局部 patch/workspace
不做逐位一致要求；`LocalTimeLimits` 保持 rank-local 候选尺度。

`app_control_contract_mpi_test --missing-case` 用共同缺失输入安全确认旧入口没有先做
控制检查，七种情况均进入后续输入读取而非目标 cold rejection；避免故意运行错序循环。
修正后以真实合法 case 测七类控制差异，2/4 ranks 均在 input 阶段一致拒绝，
accepted_steps/attempts=0 且没有输出目录；不同本地时间尺度的一步运行正常接受。
见 `app-contract-red.log`、`control-contracts-green.log` 中应用两项的记录。

## 独立运行与构建观察

本轮开始前，冻结服务 `hundun-re3900-module-reviewed-20260907.service`
已于 15:42:58 CST 退出，状态 failed/exit-code 7；尝试第 7232 步时
返回 `status=5/804 stage=15`，已接受步为 7231，最新持久化 checkpoint 为 7000。
本轮未停止、替换或重启该程序。此退出未经过正常日志完成分支，
不能据此把它归因于上述 close 缺口；其数值根因尚未确定。

从基线 `git archive` 解压的无 `.git` 源码可独立 configure、构建
`hundun`、`v04_thin_domain_runner` 和 `v04_app_case_test`。
但真实 runner 初始执行返回 `candidate_identity_status=2/10505`：
当前运行证据身份要求有效 Git commit/tree，而归档构建记录为 unavailable。
本轮未伪造 Git 身份，也未把旧冻结程序的验收标签移给此构建。
这是与旧入口清理无关的归档运行限制，需要单独定义归档来源身份合同。

## 最终干净 checkout 回归

验证的代码提交为 `33bbc6d95b6f47f0d66aa5815dac32da3c659224`。
从该提交创建新的 detached `clean-checkout`，使用全新 `clean-build`；
两者均在本轮证据目录下，没有复用 `build-review-fixes` 对象文件。
configure 与生产/测试目标构建通过，见 `clean-configure.log`、`clean-build.log`。

工具链：CMake 3.31.12、Ninja、Clang 15.0.6 + libc++、C 编译器 GCC 7.5.0、
Open MPI 2.1.1（CMake 检测的 MPI 接口版本 3.1）。Release、tests=ON、
HYPRE/ASan/UBSan=OFF；C++ 额外参数为：

```text
-stdlib=libc++ -isystem /home/wyf/.local/opt/hundun-toolchain/clang/include/c++/v1 -march=znver3 -mno-fma -ffp-contract=off
```

未使用 fast-math。构建最多 4 个作业；MPI 回归串行、每次最多 4 ranks，
使用现有独占锁。冻结长测在开始前已失败，不与运行中的长测争抢算力。

`clean-ctest.log`：26/26 通过，实际耗时 444.10 s，覆盖：

- 本轮日志关闭 CLI（内部 39 个正常/故障场景）、读取上限与原 Restart 全套、
  Application/MG 2/4-rank 合同、有符号标量 1/2/4 ranks；
- V3 observer、性能对账、runner self-test、target identity；
- 普通应用同方法→方法恢复→再次精确续算 1/2/4 ranks、统计新 epoch 1/2/4 ranks；
- checkpoint 和启用 Visit 的真实 CLI 分配失败、应用初始化/驱动/输入解析、
  跨 rank 重启续算及比较器。

这是针对性验收，不是仓库全部 CTest 或本轮 sanitizer 验收。
39 个日志场景已包含在 26 项中的一项，不把两者相加当成独立测试总数。

README 最小算例按显式初态 `101325,300,0.1,0,0`、10 步、两个输出周期均为 10
实际执行：`--version`、`validate --dry-plan`、`run` 均 exit 0，输出 Visit 和 Restart，
见 `clean-version.log`、`clean-minimal-validate.log`、`clean-minimal-run.log` 和
`clean-minimal-run/`。这仅是文档/入口检查，不是 Re3900 或科学精度证明。

## 交付边界与提交

| 类别 | 提交 | 内容 |
| --- | --- | --- |
| 正确性 | `ce9580a` | 运行器日志显式结束与全 rank 完成状态 |
| 正确性 | `e70fc1c` | Reader/Writer 回读分配前检查与 reader bulk 预算 |
| 接口防御 | `4b0dab6` | rank-local MG counters；应用冷入口控制合同 |
| 独立实验/观测 | `33bbc6d` | 有符号标量、固定几何实验、已有 C2 数据分析 |
| 诊断正确性 | `9a2e128` | 按 StatusCode 和 detail 一起选择初态冲突提示 |

生产数值算子没有随实验改变。第 4–7 项的量测、方法设计和剩余限制见
[独立后续报告](experiments/2026-09-07-c2-scalar-followup.md)。
新增公共读取容量的说明见 [Restart API](../api/restart-schema.md)。
本轮没有推送 GitHub，也没有重新提交长测。

## 归档检查发现的诊断缺口与最终补验

在 `33bbc6d` 的全新源码归档中，三个目标 `hundun`、`v04_thin_domain_runner`、
`v04_product_scalar_contract_experiment` 构建通过；最小输入 validate 通过。
普通 CLI 的 run 在 `runtime_identity` 阶段返回 invalid_plan/10505，
committed_step=0、attempts=0，却同时提示“conflicting boundary-derived initial state”。
当时命令已经给出合法 `--initial-state`，故提示确实误导。

根因是 `app_main.cpp::finish()` 只看 detail=10505；初态冲突使用
invalid_case/10505，而 `app_identity.cpp` 的身份失败使用 invalid_plan/10505。
最小 diff 只给提示条件加上 `status.code == StatusCode::invalid_case`。
没有更改错误编号、身份接受条件或任何数值方法。

新增公共 CLI 回归 `app_failure_hint_cli_test.py`：真实冲突初态仍应给出提示，
真实无 Git 身份的归档程序不应误提示初态。`failure-hint-red.log` 记录第二个断言
在旧代码失败；修正后 `failure-hint-green.log` 两个场景都通过。

最终代码提交为 `9a2e1280dc3c3f9ce71fc575e0e11e97108b0f12`。
重新创建干净 `release-checkout` 和无 `.git` 的 `release-source`，各用独立新 build；
两者均 configure 并构建上述三个目标成功。所有日志仍在本轮证据目录，
以 `release-checkout-*`、`release-source-*` 命名。这些目录名称不代表正式发布或科学验收。

| 对象 | 实际检查 | 结果 |
| --- | --- | --- |
| 新 checkout | hint、普通重启链 2 ranks、runner self-test、target identity、有符号标量 | 5/5，25.02 s |
| 新 checkout | README 最小算例，显式初态，10 步与 Visit/Restart | exit 0，`release-minimal-run.log` |
| 新归档 | 有符号标量公共 ProductDriver 路径 | 1/1，2.20 s；`release-source-ctest.log` |
| 新归档 | 最小输入 `validate --dry-plan` | exit 0，`release-source-validate.log` |
| 新归档 | 普通 CLI 身份失败诊断 | 严格拒绝、零步/零尝试、没有错误初态提示 |
| 新归档 | 圆柱 runner 一步请求 | `candidate_identity_status=2/10505`，exit 5，不打印完成 |

第 1–3 项的生产源文件在 `33bbc6d` 到 `9a2e128` 之间完全未变；
26 项主回归与后续 5 项补验按各自提交记录，不把重叠测试或归档拒绝当作额外成功样本。
后续提交只更新本文、读取 API 和 README/构建说明，没有重建或替换冻结程序。

未解决的是**归档正式运行的来源身份合同**，不是这个提示条件。
需要明确随包来源元数据和实际内容校验规则后再接入；本轮没有用占位 SHA、
移除身份校验或照搬旧验收标签。完整产品峰值、MG 层级归因、曲面标量方法、
同物理时间动量 dt 细化及长测第 7232 步的数值根因仍未闭合。

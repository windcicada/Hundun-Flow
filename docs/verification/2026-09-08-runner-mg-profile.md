# Runner 的 MG 层级成本观测

基于已验收的 [ProductDriver 接线](2026-09-08-product-mg-profile.md)，
本切片只增加可选观测和读取器校验；没有更改方程、物性、dt、残差门槛、
refinement 容量、方法历史签名或 checkpoint 格式。Re3900 生产测量尚未执行。

## 接口和计数含义

runner 增加 `--observe-mg-cost`，必须同时给出 `--observe-performance`。
默认关闭。两个开关在现有冷入口共识中核对，不进入科学/统计 fingerprint；
底层 ProductDriver 仍允许 rank-local 观测开关，runner 的文件与通信协议要求一致。
重复参数、缺少父开关、跨 rank 开关不一致均在推进和创建运行目录前拒绝。

| 文件 / 函数 | 责任 |
|---|---|
| `tools/v04_thin_domain_runner.cpp::write_mg_layout()` | 从 ProductDriver 冻结几何读取每 rank 的层数、全局/局部尺寸、coarsening、line mask；一次冷 Gather 后写独立 `mg-layout.meta`，其 SHA-256 绑定到 `RUN.meta` |
| `write_mg_loop()` / `write_mg_step()` | 每-loop 紧凑差值沿用 corrector/refinement/attempt/sweep/dt；`mg-rank-R.csv` 每步一个 total 行（level=-1）及全部冻结层的行。total 和 level 的字段互斥，不重复计总成本 |
| `complete_logs()` | 新流与既有日志一起显式 flush/close、保留首错、汇总全 rank 结果；普通日志关闭不是 checkpoint fsync 的替代品 |
| `tools/v04_solver_observe.py` | schema 5 检查独立布局、全 rank/步/层覆盖、来源哈希、逐-loop与逐-step计数/时间；schema 3/4 仍兼容 |

MG 关闭时仍输出 schema 4，原 performance/solver CSV 列保持不变；RUN.meta
增加关闭标志。开启时 schema 5 增加紧凑 MG 列和 sidecar，各行均绑定 RUN.meta 哈希。
失败尝试沿原有输出边界记录；写日志失败不改变已接受载荷，也不覆盖数值失败退出码。

六个非递归阶段（pre/post smoothing、residual、restriction、prolongation、terminal）
是 apply 总时间的互斥子集；halo、reduction、direct MPI 是嵌套成本，不能再加一次。
direct MPI 包括 coarse exchange，不能全归到 terminal solve。
MG profile 不包含 Fresh projection，也不包含已有单列的 refill/copy。

层级 sidecar 是整步汇总，不提供任意单个 loop 的完整逐层分布。每-loop 摘要仍可区分
C1/C2、refinement、retry 和组成 sweep，并额外保留 finest pre/post smoothing。
无法对账、未初始化、饱和、缺层、缺 rank、截断尾步一律不成为完整归因。
显式 `--allow-partial` 只保留已验证的完整步，complete=false；它不宣称数值失败已修好。

## 存储与使用边界

冷布局使用每 rank 257 个 int32 的固定 wire；128 ranks 时 root Gather 缓冲
为 131584 字节，另有文本编码存储。分配/路径构造均在纯本地阶段内，成功汇总后才通信。
不是全产品峰值预算。推进前后各复制定长累计 profile，关闭时不执行复制或逐层差值；
没有新增全场数组，也不把32层数组放进每个 loop。固定栈容量仍存在。
本轮没有宣称全 runner 的零 C++ 分配或零进程 RSS 增量。

新增 CSV 的序列化和刷新仍落在既有步计时范围（目前随 loop 输出计入 observables），
并非没有 I/O 成本。生产比较仍需完整进程墙钟，且两侧实验使用一致的观测配置。
完整日志按步流式读取、1 MiB 块增量哈希；新增汇总容量随 rank×level 而非步数增长。
使用 `--details-output` 避免默认 loops JSON 随窗口累积；64-loop 溢出仍拒绝完整归因，
本切片没有实现日志分段，也没有减少物理迭代来适配容量。

## 已执行回归

MPI 全部串行调度；原128-rank长测保持 SIGSTOP，旧运行目录、checkpoint 和程序未改。
测试入口使用已确认的真实 CLI 和公开 observer，未通过私有求解函数造状态。

- CLI RED：原 runner 对新参数返回2，未进入求解，见
  [输出](data/2026-09-08-runner-mg-profile/red-cli-on.log)。首次实现后1/2/4 ranks
  3/3通过（10.05 s），检查独立布局、逐层/逐轮对账及开关前后 checkpoint 全 rank 字节一致。
- 扩展真实 CLI 与统计恢复链：6/6通过，61.91 s，见
  [逐测试日志](data/2026-09-08-runner-mg-profile/release-cli-epoch-LastTest.log)。
  16³静止周期 IBM 夹具各 rank 的两步窗口均实际执行6次 MG apply；不是全零计数自证。
  包括新 observer 读取、开关缺失/重复及 MPMD 跨 rank 不一致的冷入口拒绝。
- 新 observer 合成账本先 RED（不支持 schema 5），后 GREEN。33项拒绝检查覆盖
  缺整 rank、缺/重层、非法/饱和计数、来源错配、重新绑定后的非法布局、截断尾步、
  逐层账目和 finest 错配、partial 只留前一完整步、失败输出清理及源文件不变。
  合成 ns 是算例期望值，不是 CFD 性能实测。损坏 metadata 的 partial 测试还发现
  新版初稿 `schema` 未初始化异常，已 RED→GREEN 修正。旧 observer 88项检查通过。
- libc I/O 故障：14次一小步调用，含4次正常/关闭观测对照和10次失败。
  root/非零 rank 的 MG 流分别覆盖 flush/close EIO、ENOSPC 和双错误首错保持。
  所有注入都实际触发，所有 rank 返回6，没有 COMPLETED；已持久化 checkpoint 与对照
  逐字节一致。见 [退出记录](data/2026-09-08-runner-mg-profile/log-failure-results.json)、
  [root](data/2026-09-08-runner-mg-profile/log-first-error-root.log)、
  [非零 rank](data/2026-09-08-runner-mg-profile/log-first-error-nonroot.log)。
- Clang Debug ASan+UBSan 真实 CLI 1/2/4 ranks：3/3，114.03 s，见
  [日志](data/2026-09-08-runner-mg-profile/asan-cli-LastTest.log)。`detect_leaks=0`；
  不宣称 LeakSanitizer 或总 RSS 验收。

CTest/故障 stderr 源日志已归档；退出记录 JSON 仅补仓库末尾换行。
原 JSON SHA-256 为 `01081710965ff5795f10b8778a2da7ded5114b61f1338a41c75f7bef6cf4dcbd`。
以上不是全仓库测试通过，也不是 COAST 替代或性能优化已完成。

## 下一验收

观测源码分项 DCO 提交后，从独立干净 checkout 构建 runner 并运行上述相关 CLI、
observer 和 I/O 测试。确认身份与输出后，以单轮 Re3900 窗口定位 M 的细层平滑、
粗层求解和通信成本，再只选一个有数据依据的最小优化。
两相继续暂挂，当前没有接入任何燃烧或两相源码。

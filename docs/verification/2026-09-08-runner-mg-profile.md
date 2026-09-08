# Runner 的 MG 层级成本观测

基于已验收的 [ProductDriver 接线](2026-09-08-product-mg-profile.md)，
本切片只增加可选观测和读取器校验；没有更改方程、物性、dt、残差门槛、
refinement 容量、方法历史签名或 checkpoint 格式。独立干净验收13/13和单轮
Re3900 9500→9510观测已完成；下文区分测得成本与尚未实施的优化。

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

## 干净构建与生产窗口验收

DCO源码提交 `2ea65b6`，干净验收HEAD `461d7631aa251969ddfc4e639bcadb74c640070c`，
tree `e1e46205ca609ad32e35875bb30097c1b787b008`。
独立checkout为 `/home/wyf/code_dev/.worktrees/hundun-flow-mg-runner-accept-20260908`；
新Ninja Release目录 `build-accept`，Clang/libc++，`-march=znver3 -mno-fma -ffp-contract=off`，
ASan/UBSan/Hypre关闭，测试开启，构建 `-j4`。未沿用开发目录的Ninja恢复日志。

干净CTest **13/13通过，102.90 s**，范围为旧/新observer、runner self-test、
MG CLI与统计epoch链各1/2/4 ranks、MG日志关闭故障2 ranks、ProductDriver retry各1/2/4 ranks。
见[原始日志](data/2026-09-08-runner-mg-profile/clean-acceptance-LastTest.log)和
[构建manifest](data/2026-09-08-runner-mg-profile/clean-runner-build-manifest.txt)。
runner SHA-256 `50ae573d1b8a39e606310e4a7dcfd803cbed96decb772a7308f459ed847f711a`，
manifest SHA-256 `8d4ba4897e575c185989720fcb76ef8d644510a269ff12fae58a67829157bfa2`。
manifest中的head/tree是带既有前缀的SHA-256，不是明文Git SHA；已独立重算匹配。

生产窗口仅运行一次，目录为
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/pilot-mg-9500-9510-20260908`。
使用原物理输入、原9500 checkpoint，128 ranks，显式performance/MG观测，精确续算10步。
方法签名仍为 `12213963202598979269`；未执行方法恢复或清空统计。
完整进程墙钟 **110.86 s** 包含恢复和末次Visit/checkpoint，不是十倍纯advance。

- 10步BDF2、每步1次attempt、无retry；V5 observer验证全部128 ranks、10步、70 logical loops、6层。
- 连续性最大 `2.617570949848936e-7`，能量最大 `8.805795074746012e-7`；EOS、质量目标和gauge报告为0。
  流体温度 `299.9790418656781–300.023665045917 K`；固体占位区域极值和位置没有变化。
- 末次128份Visit覆盖6070272个单元，字段/有限性/坐标/范围通过；runtime validator返回0；日志正常关闭并打印COMPLETED。
- 新generation `generation-9510-158606054634033`。与此前criterion观测窗口的同一步比较，
  **128份rank checkpoint、manifest、128份Visit、statistics及accumulator均逐字节一致**。
  complete附件仅generation名称不同，两侧均正确指向各自清单；不是全目录逐字节相同。
- 133份原checkpoint/统计来源文件及冻结程序/输入哈希前后匹配。原长测128 ranks仍SIGSTOP，
  其已接受状态更前进，不能拿pilot 9510替换原进程。

精确来源、命令、字节比较和产物哈希见[验收收据](data/2026-09-08-runner-mg-profile/PILOT_ACCEPTED.json)。
读取核查脚本和成本汇总脚本一并归档。原始CSV与逐loop明细保留在上述运行目录和同级
`mg-observation-20260908`；`FROZEN.sha256`未改写，后处理产物由独立`ANALYSIS.sha256`绑定。
这是观测接线的验收，不是新的长期稳定性、全产品硬内存预算或COAST替代验收。

## 新版成本排序与下一最小实验

下表取9502–9510的9步，每步先对128 ranks取均值再平均；9501为恢复首步单独保留。
这9步不是“统计稳态”。纯advance的逐步max-rank均值 **9.526839 s/步**；
全部实际阈值、收缩率、C1/C2/refinement和6层分布见[数值表](data/2026-09-08-runner-mg-profile/cost-window.json)。

| 成本 | s/步 | 解释 |
|---|---:|---|
| Krylov | 4.045411 | 包含A、M及正交化，不与其子项相加 |
| A / M wrapper | 1.395209 / 2.054118 | M内部Native MG为2.053430 |
| 候选装配 | 2.271481 | 独立于上述Krylov阶段 |
| MG pre / post smoothing | 0.703746 / 0.350072 | 合计占Native MG的51.32%；最细层两项合计0.776889 |
| MG restriction / prolongation | 0.306723 / 0.104016 | 互斥阶段 |
| MG terminal | 0.200820 | 9.78%；不是当前最大项 |
| MG六阶段之外 | 0.388052 | 包含入口/出口、初始化、最终projection等，尚未细分，不能全称为清零成本 |
| MG halo wait / direct MPI / reduction | 0.353292 / 0.117029 / 0.065286 | 嵌套成本，不追加到六阶段总计；direct MPI并非全在terminal |
| MG refill / 其中copy | 0.074407 / 0.011250 | copy已包含在refill，不直接swap存储 |
| 最终动量 / terminal metrics / 边界账本 | 0.202341 / 0.005116 / 0.033148 | 单列新增数值审计成本，不取消检查 |

独立residual阶段为0，因为最终缺陷已在pre smoothing内装配；prepared halo的control列为0
不代表MG没有其他MPI一致性检查。实际分区为16×8×1；全局层尺寸为
456×256×52、228×128×26、114×64×13、57×32×7、29×16×4、29×8×2。
所有层line mask为0，生产policy是Chebyshev pre=1/post=2、F-cycle，不是公共结构默认的red/black。

C2-r1平均32.33次迭代，而r0为20.11；两者M/apply分别11.405和11.737 ms。
该证据支持先看工作次数与平滑单次成本，**不支持**提前spatial或放宽阈值。
窗口仍无C2 spatial样本。未与COAST做新对照，不从110.86 s与旧墙钟差值声称加速。

本观测之后选择的实验限定为Chebyshev单stage且保留最终缺陷的路径：检查中间direction写入
是否在读取前被最终缺陷覆盖，尝试只消去这一项无消费者的写入。
保留FP求值/有限性检查、copyback、view地址、halo顺序和全部残差门槛。
该[实验](2026-09-08-mg-single-stage-store.md)随后通过公开回归和干净构建，
但单轮总墙钟没有收益，已撤回算法改动，仅保留回归和证据；没有将逻辑写入数量当成硬件流量。
这不是已复现数值缺陷，也没有耗时收益结论。现有streamed stencil、multidot不重复实现。
两相继续暂挂，当前没有接入任何燃烧或两相源码。

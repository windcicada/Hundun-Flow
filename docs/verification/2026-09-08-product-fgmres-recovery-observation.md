# ProductDriver 与 runner 的 FGMRES 恢复原因观测

开发基线 `ab1fcad5d03ebd8aed4b2818d6f338b98e2f3f8e`。
本切片把[已验收的 Native sink](2026-09-08-fgmres-recovery-observation.md)
接入现有有界逐 loop 报告和真实 CLI；不改变数值方法、门槛、restart 长度、
refinement 容量或持久化合同。两相继续暂挂。

## 接线与合同

| 位置 | 最小变更与责任 |
|---|---|
| `v04_flow.hpp` / `solver_piso.cpp::PisoPressureSolveEpoch::solve_prepared` | 末尾可选借用指针；入口清空，FGMRES 时向下传递，BiCGStab 保持 unavailable；早期拒绝不能遗留旧观测 |
| `v04_app.hpp` / `core_product_freeze.cpp` | `set_pressure_recovery_observation(bool)` 在无活动 proposal 时开关；允许 rank-local，不新增通信。现有 CaptureSolve 收集每次求解的成功或失败前缀，现有 attempt/sweep 复制路径保留它 |
| `v04_thin_domain_runner.cpp` | 显式 `--observe-fgmres-recovery`，必须同时启用 performance；开关在冷入口全 rank 一致验证。可独立于 MG 打开，无新输出文件 |
| `v04_solver_observe.py` | V6 读取新增的 13 列；五类 A 与总 A 对账、两类恢复与旧 restart 对账，并检查事件界限、rank 一致性、来源、全步覆盖与溢出 |

V6 元数据冻结压力算法和两个可选观测开关；MG 开启时仍执行完整的层级对账。
仅 MG 保持 V5，仅 performance 保持 V4，旧 V3–V5 回归保留。
新增列包括 availability、旧 norm restart 和 11 个 Native 原因/工作计数；
只在开启时写入，不给旧 CSV 或 Evidence 回填原因。非 FGMRES 即使有普通 A/M
工作，恢复摘要仍是 `available=false, counts_by_rank=null`，不能伪装为零恢复样本。

五类 A 不含原本独立计数的 recycle projection A。负递推范数、丢列、显式
重正交化、unsafe/happy/长度 restart 含义沿用 Native 合同；恢复计数不是失败数。
64-loop 容量不增加；出现 dropped loop 时仍拒绝完整归因，partial 只保留已完整
校验的步。不能为了日志容量减少求解器物理迭代。

## 存储与关闭路径

当前 Clang/x86-64 实测 sink 为 96 字节，solve observation 为 360 字节。
每个 14-slot attempt 报告因此增加 1344 字节，每个 64-slot step 报告增加
6144 字节；还存在现有报告副本和局部 CaptureSolve 的定长成本。
这些数字是每个报告对象的新增载荷，不是进程峰值模型。
关闭时仍有固定存储、初始化/复制及指针分支，不能称为绝对零开销。
没有新增全场数组、halo、计时器或热路径动态分配点。Native 热路径分配为零的
实测见前一报告；本切片尚未以新 Re3900 pilot 量化整个产品开销。

现有 solver stream 的阶段保护和显式 flush/close 汇总也覆盖新列。
真实 CLI 的 root/非零 rank 末次 flush/close EIO、ENOSPC 及先 flush 后 close
双故障检查保持首错，一致返回失败、不打印 COMPLETED；durable checkpoint
字节保持不变。关闭成功不替代 fsync，checkpoint 合同未降级。

## 实际开发验证

测试沿用用户已确认的 ProductDriver、公开数值接口、真实 CLI 和 observer CLI。
原始日志在[数据目录](data/2026-09-08-product-fgmres-recovery/)。

| 组别 | 实际证据 |
|---|---|
| Product RED | 仅声明/空接线时 1-rank 测试因缺失恢复记录失败，3.18 s；retry/control 已接受载荷仍相同，不是数值失败 |
| Product GREEN | 接线后 1/2/4 ranks 3/3，3.14 s；仅 root 或末 rank 开启、真实 dt retry/耗尽、首末尝试、下一 BDF2 步和已接受载荷回退 |
| runner RED | 对旧干净 `74e262d` runner 传新开关返回 usage/2，未运行新观测；不是原 CFD 故障 |
| runner GREEN | 2 ranks 的 MG+恢复、仅恢复两种路径均实际执行总计 12 次 Arnoldi；与关闭时 checkpoint 字节一致，observer complete，11.66 s |
| 扩展 Release | Product 1/2/4、runner 1/2/4、2-rank 日志关闭及 V6 reader，共 8/8，74.12 s |
| 旧 reader 兼容 | V3/4 的 88 checks、V5 的 33 checks，共 2/2，5.77 s |
| 新 reader | 28 项拒绝检查，另验有/无 MG、非 FGMRES、partial 整步、details-output 等价及原输入只读；这是合成账目，不是 CFD 耗时证据 |
| ASan/UBSan | Product 1/2/4 和 runner 2 ranks，共 4/4，181.74 s；之后包含真实 BiCGStab 的完整 runner 扩展再通过 1/1，88.59 s |
| 真实非 FGMRES 合同 | 补正输入后的 2-rank CLI 通过，16.99 s；BiCGStab 实际 A+M=40，恢复 unavailable，开关前后 checkpoint 相同 |

两个测试输入问题单独记录，不修改生产算法来迎合夹具：

- V6 截断测试起初只删 CRLF 的 LF，剩余 CR 在 Python 通用换行模式下仍是合法
  终止符。改成删掉整个行终止符，生产 reader 不变。
- 真实 BiCGStab 扩展最初继承 FGMRES `krylov_restart=12`，在关闭观测的基线
  就被合法性校验拒绝（CLI exit 4）。`valid_solver()` 的既有合同要求非 FGMRES
  为 0；夹具改成 0，MG scaling 仍由公开 case parser 选择 unit-linear。
  这一路输入拒绝的 runner 首次返回没有打印具体 status，属于另行登记的诊断缺口，
  本次不借观测接线改动输入控制路径。

环境：Clang 15/libc++，Release 未启用 sanitizer/fast-math；Debug 使用 ASan+UBSan，
`detect_leaks=0`，不是 LSan/RSS 验收。MPI 全部串行，最多 4 ranks；原 128-rank
长测保持 SIGSTOP，原程序、进程状态和 checkpoint 不替换。CodeGraphF 已同步接线源码。

## 干净候选验收

本地 DCO 提交 `7c03a54a9f1061d1de0e0922ed38a3bcc8f2412f`，tree
`9d502ca0023b10704aa721d94cbbb5a33e7b6db6`。独立 checkout
`/home/wyf/code_dev/.worktrees/hundun-flow-product-recovery-accept-20260908`
构建与测试前后 Git 干净；全新 Ninja Release，Clang 15.0.6/libc++，
`-march=znver3 -mno-fma -ffp-contract=off`，HYPRE/sanitizer 关闭。
生产 `hundun`、runner 和所需测试目标均构建成功，没有开发 Ninja 日志恢复告警。

现有 Krylov/PISO、Product retry/方法恢复、普通应用恢复、runner 统计 epoch、
新开关/关闭及 V3–V6 reader 共 **31/31** 通过，191.10 s。真实 FGMRES 与
BiCGStab producer 均覆盖 1/2/4 ranks。源码/说明的 diff whitespace 检查通过；
原始工具/CTest 日志的尾部空格和换行按原样保留，不将原始日志检查误报为源码检查。

- runner SHA-256：`345cca4802286a7ae1b2bd39b7c25afceeedb8a13fc60b5fdf71f80d92f0f67a`。
- manifest SHA-256：`c352d1c4ff748b5ed33d3d068115d8de54f3b8e61e5efe70cfb3ffb10a3b2e21`。
- core content SHA-256：`f5fd1d2f609024154c5c38718690549fdbf7cecceef9577b3c257158fc66ea4a`。
- manifest 两项 clean=true，带前缀的 head/tree 哈希独立重算相同；归档 manifest
  与构建原件 SHA-256 一致。[本地验收记录](data/2026-09-08-product-fgmres-recovery/LOCAL_ACCEPTED.json)
  限于该观测接线，不是完整产品验收或压力算法变更。

## 单轮目标窗口

干净验收后，独立目录 `pilot-fgmres-9500-9510-20260908` 的新原因观测配置
已完成 9500→9510、128 ranks，**仅一轮，总进程墙钟 113.58 s**。
启动前后 133 个源 checkpoint/附件文件 SHA-256 全部保持；冻结的新二进制、
manifest、reader、运行脚本和原物理配置也通过哈希回查。
这是为补缺失原因计数，不重跑已经结束的相同观测配置。

实际验收：10 步全部 BDF2、每步 attempt=1，无 retry/方法恢复；70 个 logical loops
和 128 ranks 的 V6 来源/步/loop/六层 MG 账目完整。最大连续性/能量残差分别为
`2.617570949848936e-7` / `8.805795074746012e-7`，均满足原门槛。
流体 T 范围 `299.9790418656781–300.023665045917 K`，固体占位极值不漂移。
本次只验区域极值，不把它写成新的 source-to-final 固体全载荷检验。

128 份 checkpoint、manifest、128 份 Visit、accumulator 和 statistics 与之前
`pilot-mg-9500-9510-20260908` **逐字节一致**；旧 CSV 全部非计时/非来源列也相同。
完成声明只允许 generation 名不同，各自指向通过验证的 generation。
runtime validator 返回 0；沿用 epoch=7000、采样起点=17001，样本仍为 0。
新 generation 为 `generation-9510-168455369402501`，只属于本 pilot。
原 128 个 SIGSTOP ranks、较前进的内存状态和源 9500 checkpoint 均未替换。

结果见[目标窗口验收](data/2026-09-08-product-fgmres-recovery/PILOT_ACCEPTED.json)、
[物理/输出检查](data/2026-09-08-product-fgmres-recovery/pilot-output-verified.json)与
[紧凑原因归因](data/2026-09-08-product-fgmres-recovery/RECOVERY_ATTRIBUTION.json)。
完整原始 per-rank CSV、`observation.json`、`loops.jsonl` 和 `cost-summary.json`
保留在基准目录 `trial-D0p02-zpi2-52/fgmres-observation-20260908`，收据记录其哈希。
归档脚本与基准目录原件字节一致，仅用于该冻结目录，不宣称可从任意目录直接重放。

## 新测量与下一项最小实验

首步恢复暖态成本单列，以下为 **9502–9510 后九步**，不称为统计稳态。
一次 logical loop 的计数只算一次，不乘 128；时间为 rank 平均，advance 另取每步
max-rank 后平均。A/M、MG 和通信是嵌套成本，不能将它们全部相加。

| 分组 | loops | 迭代/M | A | 丢列 / unsafe restart | happy / 长度 restart | 单次 M 均值 |
|---|---:|---:|---:|---:|---:|---:|
| C1 spatial | 9 | 204 | 300 | 45 / 45 | 0 / 0 | 11.591 ms |
| C2-r0 diagonal | 9 | 181 | 275 | 42 / 42 | 0 / 0 | 11.689 ms |
| C2-r1 diagonal | 9 | 291 | 428 | 58 / 58 | 0 / 0 | 11.486 ms |
| 全部 | 63 | 1609 | 2358 | 347 / 347 | 0 / 0 | 见各 loop 组，不混淆方向 |

此前无法区分的 347 次恢复，现已确认全部进入后续列的 unsafe 分支，
不是 happy restart，更不是 347 次时间步失败。347 个已付出 A/M 的方向被丢弃，
占 Arnoldi 次数的 21.57%；另有 347 次 unsafe 真实残差 A，占总 A 的 14.72%。
负递推范数共 373 次；显式重正交化共 42 次，包含 column=0 的 26 个负事件和
16 个零范数事件。计数对应程序分支，不代表流场出现负物性或负真实范数。

C2-r1 的 428 次 A 完整分为 `9 初始 + 291 Arnoldi + 58 unsafe 残差 + 70 周期内残差`，
cycle-end A 为 0。r1 单次 M 没有比 r0 更贵，额外工作首先在次数上；其真实线性
阈值 `1.4628e-8–1.5590e-8` 也明显小于 r0 的约 `1.5e-5`，这是现有控制，不能为
减少次数而放宽。r0 的候选后能量分量均值比为 1.0944，r1 为 0.1617；中间分量
增加由现有联合筛选/最终审计约束，不能仅凭这两个比值宣布算法错误。

本窗口 max-rank advance 均值 9.757750 s；rank-mean Krylov/M/A/候选分别为
4.062262 / 2.046287 / 1.397141 / 2.365309 s/步。最终动量装配 0.205999 s/步，
账本 0.034380 s/步。MG pre/post 平滑占 MG apply 的 50.71%，仍不支持优先重做
已撤回的 MG copy/store 实验。113.58 s 与旧配置的 110.86 s 是两次不同观测配置
的单轮记录，未隔离系统波动和观测成本，不宣称提速或精确的观测开销百分比。

下一项只做一个局部假设验证：**后续列 unsafe 时，能否复用已存在的显式重正交化，
保留已经算出的方向，而不立即丢列重启**。这是待试验方向，尚未实现/验收。
先用公开 FGMRES 接口的近相关、非对称小算子比较独立真实残差与 A/M/归约工作，
再验受控 unsafe 列、首/非零 rank 失败、1/2/4 ranks、输出不变和工作区分配合同。
非有限性、MPI 顺序、最终真实/物理审计、迭代及 restart 上限保持；显式恢复本身
需要额外归约，不能按 21.57% 丢列比例预报速度收益。
局部证明通过后才独立提交实验，进入同物理窗口单轮比较；无总成本收益则撤回。

这些原因与五类 A 计数用于选择下一项最小实验，不是逐类计时，
也不是可直接消除的成本。当前不提前 spatial、
不放宽阈值、不增加 restart，不给出加速百分比或 COAST 替代完成结论。
本窗口为零标量 Re3900，不把这组成本外推到通用组分或燃烧路径。
Evidence/checkpoint 格式与方法签名不变；完整 COAST 对标和燃烧接入由
[目标台账](../plans/2026-09-08-coast-replacement-and-stages.md)继续跟踪。

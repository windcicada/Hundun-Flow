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
| ASan/UBSan | Product 1/2/4 和 runner 2 ranks，共 4/4，181.74 s；该次 runner 尚未加入后续真实 BiCGStab 扩展 |
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

## 下一验收边界

提交观测代码后用独立干净 checkout 构建，验收现有 Krylov/PISO、Product retry、
重启方法链、runner/关闭和各代 reader。通过之前不把开发二进制当新冻结程序。
之后只启动一次带新原因观测的 9500→9510、128-rank 诊断配置；这是为补缺失原因
计数，不重跑已经结束的相同观测配置。前一步恢复暖态成本单列，后九步用于局部归因。

取得真实原因及五类 A 的成本分配后才选择下一项最小实验。当前不提前 spatial、
不放宽阈值、不增加 restart，不给出加速百分比或 COAST 替代完成结论。
Evidence/checkpoint 格式与方法签名不变；完整 COAST 对标和燃烧接入由
[目标台账](../plans/2026-09-08-coast-replacement-and-stages.md)继续跟踪。

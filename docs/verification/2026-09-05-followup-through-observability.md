# 重启、receipt与性能观测结果

起点：`0b0eff11c8365e570725d6220d031c84a6aef4e0`，
分支`codex/v04-restart-receipt-observability`。用户要求推进到剩余回归和性能观测，
然后给出优化方案；不实施内核优化、不推送、不长测、不修改原始checkpoint或历史凭据。
沿用上一轮Release/独立ASan+UBSan工具链。已有AGENTS.md和.codegraphf不纳入提交。

本轮已完成上述1–3项。Release全量277/277、ASan+UBSan重点22/22通过；
目标Re3900从旧checkpoint的743步恢复到746步，3步均一次接受并保持BDF2。
这是一次带观测的短诊断，不是长期稳定性、实验统计量或超越COAST的证明。

## 1. 重启兼容

- 默认严格读取原69d8eee小算例：10307，restart_load，0次推进。
- 最初只重建schema/product尾部的兼容实现能通过合成测试，但原始文件仍失败。
  独立文件重放发现BoundaryCompiler的semantic_hash也包含registry fingerprint；
  现复用同一哈希实现重建旧边界身份，不分配第二套边界/求解工作区。
- 兼容必须显式选择mg-bundle-ghost-v1；只接受完整格式2历史，默认行为不变。
- 公共合成测试验证旧身份拒绝/显式接受、错误身份拒绝、全部场与速率历史、面通量、
  时间控制、压力参考保持一致，以及新checkpoint使用当前身份并通过严格加载。
- 1/2/4 ranks：3/3通过。实际原二进制checkpoint在2 ranks从步1续到步2，
  随后新checkpoint默认严格路径4 ranks续到步3。
- thin runner的完成标记和accumulator仍核对来源身份，输出标明来源和迁移，
  新文件使用当前身份。目标Re3900的743步原文件也已成功迁移，见第4节。

日志前缀`/tmp/hundun-followup-`；原始结果和哈希收录于第6节的证据包。

## 2. receipt

用户确认通过receipt-validate CLI验证。新增13项CLI检查，先确认旧程序缺失
历史模式，再实现显式`--historical`、根目录映射及内容寻址附件库。
损坏、缺失、重复映射、不完整却要求complete均拒绝；默认严格绑定当前工具不变。

- R4原文件不改写，SHA-256仍为
  `1dfdeb16b92f4a29a5dc2a0c4daa146499e8decb1af5df197b652203e36b2e3f`。
- 两份历史工具从59985db精确恢复为`.snapshot`附件，其哈希与R4一致。
  历史附件不执行、不补水印；当前生产工具水印保留。
- 原R4的14项附件通过历史校验，当前校验器与原提取器身份分别报告。
- 新生成[R5凭据](v0.4-literature-data-receipt-r5-partial.json)，绑定17项附件，
  包括父R4、历史提取工具、当前工具；默认严格校验通过。
  这次没有重新提取实验数据，也没有改变参考数值或完整性判定。
- R4/R5均为`complete=false`；R5的`--require-complete`按预期退出2。
  Norberg未完成项没有因路径修复而消失。
- [复现命令及历史快照说明](literature-r4-artifacts/README.md)。
  原始PDF/数字化附件仍是外部依赖；不将缺少附件的环境伪报为通过。

## 3. 剩余回归和观测

已补充并验证：

- `initialize_restart`真实new分配失败扫描：Cartesian/IBM、1/2/4 ranks，6/6。
  每次失败后在同一driver上重试，核对历史并检查跟踪的C++/通信资源归零。
- 新增`RestartReader::load`扫描时，原程序rank 0首次分配失败触发SIGABRT。
  根因包括read_current_name的noexcept、广播前本地分配未汇总、以及读入/重分区
  构造的外层catch。现按本地阶段汇总错误，保留原输出对象，并区分分配与文件错误。
  完整扫描Cartesian/IBM、1/2/4 ranks，6/6；2-rank示例扫描108/100个分配点。
- 已建立时间历史的非均匀重启夹具：外推alpha=2被拒绝，alpha=1梯级候选被接受，
  整步一次完成；1/2/4 ranks均覆盖。与原候选计数及应用计时检查合计7/7。
- POSIX边界故障注入确认open/write/fsync的单次EINTR原先会终止输出。
  现重试EINTR、完成短写，但不重试close；ENOSPC等首个错误不会被close覆盖。
  Visit/Screen/Monitor/Evidence提供固定容量的errno、操作、路径、失败rank报告，
  即使仅rank 0请求报告也保持同一MPI流程。实际目录冲突在ApplicationService中
  报告为已接受步之后的Visit失败。相关新增及原分配回归24/24。
- ApplicationRunReport/CLI区分完整本地运行时间与原advance计时；计时包含初始化、
  输出、资源、Evidence和清理。数值候选细分8个不重叠的包含式阶段，Schur求解分为
  准备、Krylov、关闭/焓恢复。时间只在阶段边界采样，没有逐单元计时或额外候选MPI。
- thin runner可选`--observe-performance`保存各rank的整步及子阶段数据。
  其计时不包含写这份观测记录本身的gather/write，不能据此声称完全无I/O。
- 可选PMPI动态附件只在advance区间统计9种阻塞collective及调用点，不增加collective。
  独立2-rank检查确认区间内1次Allreduce、2次Bcast、1次Barrier被准确计数，区间外忽略。

尚不宣称RestartWriter所有POSIX错误上下文/分配路径均已完成审计。

最终验证使用同一套已格式化的C/C++源码：

| 验证 | 结果 | 范围和限制 |
| --- | --- | --- |
| Release全量CTest | 277/277，708.22 s | 包括本轮新增测试和现有数值、MPI、receipt测试 |
| Debug ASan+UBSan重点CTest | 22/22，593.20 s | 重启、真实分配失败、输出、候选和应用路径；关闭LeakSanitizer，不据此宣称无泄漏 |
| Re3900真实旧checkpoint短诊断 | 单轮3步，退出0 | 128 ranks，最终Visit和checkpoint均完成 |
| 后处理工具 | 自检和原始数据对账通过 | 不增加物理试算轮数；计时按同一rank的完整步分解 |

Release使用Clang 15.0.6/libc++、OpenMPI 2.1.1，`-O3 -DNDEBUG -march=znver3
-mno-fma -ffp-contract=off`，没有开启fast-math。完整CMake缓存、二进制哈希和
测试详细日志随证据包保存；ASan/UBSan与性能运行分开构建、分开执行。

## 4. Re3900单轮诊断

### 4.1 算例与身份

本轮只运行Hundun SIMPLE，不更改算例输入、数值验收要求或原始checkpoint。

| 项目 | 本轮配置 |
| --- | --- |
| 圆柱、来流 | D = 0.02 m，Uc = 2.89668 m/s，Re = 3900 |
| 计算域 | 20D × 10D × (π/2)D，即0.4 × 0.2 × 0.031415926535897934 m |
| 坐标范围 | x ∈ [-0.1, 0.3]，y ∈ [-0.1, 0.1]，z ∈ [0, 0.031415926535897934] m |
| 网格 | 456 × 256 × 52，共6,070,272单元；沿用原COAST格式坐标文件 |
| 并行 | 128 ranks，16 × 8 × 1，`--bind-to core --map-by core` |
| 算法 | SIMPLE，2个压力校正，BDF2，固定dt = 1.3808912271980336e-5 s |
| 物性与模型 | 原COAST native air输运、NASA7热力学；Vreman_wall_function，adaptive-order IBM |
| 窗口 | 743 → 746，t = 0.010260021818081574 → 0.010301448554897516 s |
| 输出 | 配置周期500步；runner在本次结束时强制写最后一步Visit/checkpoint |

算例根目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52`。
输入为该目录的`case/`和`statistics-long-20plus50D.d`。
没有用常cp、常分子黏度或常导热系数替换COAST物性。

来源为`audit-followup-734to743/Restart/generation-743-415824991647054/manifest.bin`，
SHA-256为`de7eed33c99cf9450a5c119c1fad32c3e747a1597727eb097ee288aae3f1ab0e`，
该manifest的哈希在运行后复核未变。显式选择`mg-bundle-ghost-v1`后，来源plan/schema为
`1946063505981049456` / `16950592525103096368`，新product为
`8264025224147819262`。`restart_requires_recovery=0`，未丢失历史或退回启动步。
该兼容开关只处理已确认的旧MG存储身份，不允许绕过其他配置、网格或耦合路径的身份检查。

运行二进制来自提交`9fe82f39ef8672c30123b88c307a22ae592d7a3f`，
thin runner SHA-256为
`57e9754f72a26f907952a19425a2710d27b11f3284fcf1a6b7a106d97e1cb598`；
PMPI附件SHA-256为
`1f9629904b79fa6481c8a5271a3c21193383124ae379844c7774dab29da18c87`。
后续文档/后处理提交不会被标作这次二进制的构建来源。

新结果分别写在上述算例根目录下：

- `audit-observability-743to746-20260905/`：运行身份、384行rank-step性能记录、
  Evidence、health/force/probe、最终场和checkpoint。
- `audit-observability-743to746-20260905-mpi/`：128份原始MPI调用点记录。

### 4.2 完整耗时

| 步号 | 全rank最大advance / s | 全rank最大完整步 / s | 完整步关键rank |
| --- | ---: | ---: | ---: |
| 744 | 17.395148 | 17.406066 | 3 |
| 745 | 15.738122 | 15.748553 | 44 |
| 746 | 17.031529 | 28.285664 | 0 |

3步最大advance的算术平均为**16.721600 s/步**。这里只对同一轮中的3个步取均值，
没有重复三轮取中位数。完整步总和61.440283 s；启动到MPI退出的`/usr/bin/time`
墙钟为71.81 s。两者约10.37 s差额包括加载、启动/结束以及观测自身开销，
不能全部归给初始化，也不能把advance当作程序总耗时。

746步关键rank的分解：advance 17.031523 s，观测量0.010011 s，
Visit **2.099422 s**，Evidence/资源0.000545 s，CSV 0.000071 s，
checkpoint **9.144091 s**。这些是同一rank的六个不重叠区间，
不是把各rank各阶段的最大值相加。完整步计时尚不含写本份性能记录的gather/write；
总进程墙钟包含这些开销。

本次结束强制输出一次，不能把20.480094 s的完整步均值当成500步输出周期下的
常态吞吐量。若输出成本相同，Visit+checkpoint按500步摊销约0.0225 s/步，
这只是周期换算，尚非长时实测。

### 4.3 模块耗时

以下时间为128 ranks × 3步的本地算术平均；占比以同口径的advance均值
16.721586 s为分母。子阶段是包含式计时，例如边界阶段同时包含halo、边界处理
和派生输运准备，不能把它直接称为纯通信或纯物性耗时。

| 范围 | 平均 / s每rank每步 | 占advance |
| --- | ---: | ---: |
| Krylov求解 | 8.650889 | 51.74% |
| 候选评估合计 | 3.808093 | 22.77% |
| Schur/生命周期准备 | 0.680251 | 4.07% |
| 关闭/焓恢复 | 0.029047 | 0.17% |
| 候选内部：残差装配/审计 | 1.041883 | 6.23% |
| 候选内部：状态halo/边界/派生输运 | 0.893727 | 5.34% |
| 候选内部：通量 | 0.735916 | 4.40% |
| 候选内部：物性/精确准备 | 0.395299 | 2.36% |
| 候选内部：证书 | 0.305439 | 1.83% |
| 候选内部：最终等价性/哈希 | 0.256481 | 1.53% |
| 候选内部：速度/校正halo | 0.178541 | 1.07% |
| 候选内部：初始标量复制/revision | 0.000766 | <0.01% |

“候选内部”各行已包含在候选合计中，不能再次相加。压力—能量阶段stage 50
本地均值约14.547 s，占advance约87%；上表的求解/候选也包含于这个阶段。

| 步号 | 基线评估 | 外推评估 | 梯级评估 | 不完整评估 | 被拒外推 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 744 | 12 | 10 | 4 | 0 | 2 |
| 745 | 11 | 8 | 3 | 0 | 0 |
| 746 | 12 | 9 | 4 | 0 | 1 |

这些数值来自真实候选路径，所有rank一致。细分报告只保留最近一次数值尝试；
本窗口每步仅一次尝试，因此没有遗漏较早重试的耗时。
health还记录每步2次压力主求解，refinement分别10、9、10次，线性迭代分别
321、314、331次。下一轮需要同时解释这些额外校正的必要性和单次迭代成本，
不能只把配置里的“2个压力校正”当成实际求解工作量，也不能直接删减refinement。

### 4.4 通信对账

PMPI只观测`driver.advance`内的Allreduce、Bcast、Reduce、Barrier、Gather、
Allgather、Allgatherv、Alltoall、Alltoallv。384个rank-step覆盖完整，
记录丢弃为0，调用点均可映射到二进制符号。它不覆盖全部MPI类型或点对点halo等待。

| 步号 | 每rank观测次数 | 现有汇总次数 | 差额 | 每rank阻塞MPI平均 / s |
| --- | ---: | ---: | ---: | ---: |
| 744 | 9545 | 8415 | 1130 | 3.410879 |
| 745 | 8995 | 8021 | 974 | 2.902052 |
| 746 | 9674 | 8579 | 1095 | 3.035486 |

按实际调用者归类后的每rank每步均值如下。纳入现有计数的五类累计数与报告精确一致，
未把未归属调用塞进Reduction/MG等无关计数。

| 调用归属 | 平均次数 | 平均阻塞时间 / s | 已计入汇总 |
| --- | ---: | ---: | --- |
| ReductionEngine（含packet、contract、consensus） | 4387.00 | 2.197396 | 是 |
| 普通结构化halo控制 | 1972.00 | 0.516442 | 是 |
| MG周期及构建 | 1872.00 | 0.183283 | 是 |
| IBM donor控制 | 106.33 | 0.001784 | 是 |
| 物性predictor | 1.00 | 0.010374 | 是 |
| Coupler/BoundaryFinalizer直接调用 | 961.33 | 0.165494 | 否 |
| 其他直接调用（时间控制、事务、静态边界等） | 105.00 | 0.041366 | 否 |

九种阻塞MPI的平均3.116139 s约占advance的18.6%，**已经包含在求解/候选时间内**，
不是另加的耗时。PMPI耗时不含附件自身的调用点解析成本，而advance包含它。
归类依据是函数所属模块，不是看到参数类型含`ReductionEngine&`就归给Reduction。

漏计部分平均1066.33次但仅约0.207 s；普通halo控制约0.516 s。
因此“先删掉约一千次错误检查就能明显追平COAST”没有本轮实测依据。

### 4.5 数值和内存预判

| 步号 | 尝试次数 / 实际阶数 | 连续性残差 | 能量残差 | 已提交流出CFL |
| --- | --- | ---: | ---: | ---: |
| 744 | 1 / BDF2 | 8.456962e-7 | 2.033915e-7 | 0.220111 |
| 745 | 1 / BDF2 | 9.343137e-7 | 8.233027e-7 | 0.220689 |
| 746 | 1 / BDF2 | 8.996885e-7 | 2.156118e-7 | 0.220651 |

连续性/能量阈值仍为1e-6，CFL阈值0.8；EOS残差为0。
温度约299.9826–300.0383 K，密度约1.238978–1.239248 kg/m³。
没有重试、重启恢复降阶或时间阶数回退。结论限于：这次旧checkpoint恢复后
3步的状态与终端审计正常，最终输出成功；不等同于空间/时间精度或实验吻合验证。

三步的最大rank RSS高水位均为226.3828 MiB。节点内各rank高水位求和从
28,583,550,976增至28,584,632,320 bytes，约26.62 GiB。
`ru_maxrss`是进程历史高水位，各进程峰值不一定同时发生，不能称为同步节点峰值。
三步不足以证明RSS不再增长或不存在泄漏。当前t* = tUc/D约1.492，仍在发展段，
`samples=0`、`statistics_eligible=0`，不产生正式统计/实验数据验收。

## 5. 下一轮具体优化方案（尚未实施）

### 5.1 先做Krylov内部定位，再选择一个内核改动

Krylov占约52%，是第一优先级。下一轮在现有求解入口细分算子应用、MG预条件、
正交化、多点积、真实残差重算及邻居等待；同时保留iterations、operator/preconditioner
applies、breakdown/restart和物理验收拒绝原因。每一层采用可对账的区间，
先区分迭代次数偏多还是每次迭代偏贵，不能仅凭“Krylov耗时高”指定根因。

若确认正交化的反复全场扫描是热点，再实施分块multidot，保持现有补偿求和、
有限性检查和批量归约；若主要成本来自MG或算子，则先修相应实测热点。
每次只改变一个内核，用真残差、守恒、1/2/4 ranks及单轮相同物理窗口验收。
不预先承诺multidot提速，不换PCG、不盲目增大restart，不使用fast-math。

### 5.2 减少候选重复工作，不放宽验收

候选评估约23%，其中残差装配、边界/派生输运和通量合计约2.67 s。
优先检查这些阶段重复的网格/面遍历、静态索引构造和不随候选变化的数据。
可缓存的首先是已编译静态几何和字段槽位映射，不是未经验证的物理结果。

代码已有alpha=0物性复用和已接受候选复用，不能重复把它们列为新收益。
若进一步跨校正复用基线，必须绑定完整状态、通量、dt、边界、IBM、方程和权限身份，
保留“不更新revision却污染字段”的拒绝回归；失败不能影响回滚和MPI一致性。
不减少必要的连续性/能量审计，不通过常物性改变COAST对照条件。

### 5.3 通信和I/O按实测成本排序

- 通信：先给961.33次耦合/边界直接调用与105次其他调用补准确归属。
  仅在阶段边界可证明安全时合并consensus；单rank失败必须全体一致退出。
  计数下降不是验收标准，还要测总墙钟和等待时间。
- 输出：Visit约2.10 s、checkpoint约9.14 s。先补RestartWriter尚未覆盖的真实
  分配/POSIX故障，再评估可复用staging或分块写出。保留checkpoint持久化和
  完成标记语义，不用取消fsync来制造比较优势。
- 内存：先区分常驻数组、arena外数组、MPI/halo和I/O峰值，避免重复计数。
  阶段工作区别名复用排在上述工作之后，必须同时检查view/revision、证书、回滚
  和未完成通信。短窗口高水位不作为“内存问题已彻底解决”的证据。

仅用于排列投入优先级的算术估算：若Krylov本身快20%，按当前分布约节省
1.73 s/步（advance的10.3%）；若候选本身快20%，约节省0.76 s/步（4.6%）。
这不是已实现或承诺的加速比。

### 5.4 保留统一性能比较列表

后续沿用此表更新，不用旧的不等配置数据补空格。每配置仅运行一轮，以推进相同
物理时间的总成本为主，同时记录完整墙钟、advance、迭代/候选次数和RSS。

| 路径 | 本轮状态 | 可用时间 | 比较资格 |
| --- | --- | --- | --- |
| Hundun SIMPLE | 743→746，3步完成 | advance平均16.7216 s/步，完整窗口71.81 s | 带PMPI诊断；不是正式无探针性能基线 |
| Hundun PISO | 本轮未运行 | — | 待同条件单轮测试 |
| COAST普通可压缩 | 本轮未运行 | — | 待同条件单轮测试，保留在比较列表 |

正式比较应固定同一网格、D/U/Re、输运/热力学、边界、rank绑定、dt和物理时间窗口，
分别保留各程序的数值正确性要求。诊断探针与计时基线分开标记，不能将带探针的
本轮结果直接与无探针COAST作加速比。SIMPLE/PISO的配置身份不能互相伪装；
需从各自合法且等价的初态/时间历史开始，或使用另经验证的转换方案。

## 6. 证据与复核

[机器结果与构建身份](2026-09-05-followup-evidence/results.json)、
[原始文件清单](2026-09-05-followup-evidence/manifest.json)、
[SHA-256校验表](2026-09-05-followup-evidence/SHA256SUMS)。
包内保留最终Release/Sanitizer详细日志、初次复现失败和中间调试日志、源差异、
实际运行CSV/Evidence、128份PMPI记录以及符号级解析结果。日志文件名中的green
不作为通过依据：例如早期`restart-legacy-green.log`实际仍失败，最终通过的是
`restart-legacy-green2.log`及真实目标运行日志。

包内79份日志，SHA256SUMS覆盖220个文件；校验表和218份压缩载荷的原始/存储哈希
均已复核。另以原始TSV/CSV独立核对384条rank-step的MPI次数、耗时、差额和
完整步关键rank分解，均一致；文档引用的本地文件存在。

大体积Visit/checkpoint留在新运行目录，不提交到Git；历史PDF等外部数据不复制进包。
打包器中的`historical_receipt_read_only_audit`只比较当前路径，不是历史校验器，
历史工具与当前工具的哈希本来就可以不同；是否通过历史校验以绑定的
`receipt-r4.log`及显式CLI结果为准。

在仓库根目录复核包完整性：

```sh
(cd docs/verification/2026-09-05-followup-evidence && sha256sum -c SHA256SUMS)
```

真实运行完整命令保存在包内`logs/hundun-followup-re3900-time.log.gz`。
只重新解析已有运行、不启动CFD：

```sh
case_root=/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52
python3 tools/v04_performance_observe.py \
  --performance "$case_root/audit-observability-743to746-20260905/performance.csv" \
  --binary build-re3900-fix-test/versions/v0.4/v04_thin_domain_runner \
  --mpi-directory "$case_root/audit-observability-743to746-20260905-mpi"
```

解析要求使用上述SHA-256对应的二进制，以保持调用点符号解释一致。
本轮停止于测试、观测和优化方案；没有推送main、开启长测或实施上述性能内核改动。

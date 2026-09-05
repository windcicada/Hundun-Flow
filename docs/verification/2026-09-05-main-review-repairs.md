# main 审阅清单：修复、观测与剩余工作

本轮输入为用户附件《工作基线与执行原则》。以下“本轮”只指从
`69d8eee715c9565ccc033ea7fb95d2bf3b95aae6` 开始的这次工作，不沿用历史测试作为当前通过证据。
已实施第一批的核心正确性/存储修复，以及第二批的部分观测修正。
**完整分类内存预算、全部通信归类、细粒度墙钟和第三批物理设计尚未全部实现。**
没有进行大规模算法重写、新的 Re3900 长期统计或 COAST 正式速度对照。

## A. 基线、范围与环境

| 项目 | 本轮记录 |
| --- | --- |
| 实际 origin/main | fetch 得到 `69d8eee715c9565ccc033ea7fb95d2bf3b95aae6`，与附件参考相同 |
| 活动工作区 | `/home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity`；原 `/home/wyf/code_dev/hundun-flow` 路径已不存在 |
| 分支 | `codex/v04-mpi-memory-review`；本地提交或差异身份见证据包，不推送 |
| 初始脏文件 | 用户已有 `AGENTS.md` 修改和未跟踪 `.codegraphf/`，原样保留，不纳入修复提交 |
| Release | Clang C++ 15.0.6 / libc++；C 为 GCC 7.5.0；Open MPI 2.1.1；tests=ON，ASan/UBSan=OFF |
| Release 选项 | `-stdlib=libc++ -march=znver3 -mno-fma -ffp-contract=off -O3 -DNDEBUG`，无 fast-math |
| 检查构建 | Clang C/C++ 15、libc++、Debug、ASan+UBSan；外部 MPI 的 LeakSanitizer 关闭，不能据此声称无泄漏 |
| 实测起点 | 重新构建后 CTest **216/217** 通过，43.47 s；一项历史 receipt 路径失效，未记为通过 |
| 不变范围 | Re3900 网格、物性、边界、固定 dt、迭代容量、物理残差门槛；历史政策、日志和 checkpoint |

最终受测源码的结果：

| 验证范围 | 结果 | 进程退出码 / 实际总时长 |
| --- | --- | --- |
| Release 全部 CTest | **246/247 通过**；包含18组真实分配失败全扫描，覆盖普通/IBM产品的compile/create/initialize及1/2/4 ranks | CTest 8 / 445.67 s |
| 最终改动相关 ASan/UBSan 回归 | **14/14 通过**，覆盖公共报告、MG容量、内存阶段观测、候选记录和halo | CTest 0 / 112.03 s |

Release 唯一失败为 `v04_literature_partial_receipt`，基线也失败，并非本轮新增数值回归。
旧绝对路径映射到活动工作区后，14个条目均找到文件，12个哈希相符，两个提取脚本与历史哈希不符。
本轮保留失败和历史凭据，不修改旧哈希使其通过。早一轮较广sanitizer筛选为26/27，
唯一MG收敛超时随后单独运行通过；这不等同于247项全部执行了sanitizer。

探索使用活动工作区的 CodeGraphF CLI，状态检查/同步后再用源码核对。
没有调用 ponytail、创建 release/tag、修改全局运行环境或提交长期统计任务。
使用已确认的 ApplicationService、ProductCompiler/CompiledCasePlan、ProductDriver 和 MPI writer 公共接口。
报告结构按 codebase-design 的接口分工，把最终结果、尝试历史和停止原因分成独立含义。

## B. 逐项处理记录

“复现”表示真实路径的失败或分配观测；“控制流确认”不等于覆盖所有输入组合。

| 位置 → 原始疑点 | 核查/复现证据 → 影响 | 修改 → 回归与剩余风险 |
| --- | --- | --- |
| `ProductCompiler::compile` 本地分配后进入 collective | 单 rank 的对象、字段目录、图/arena、面数组、halo 描述分配失败，超时或 MPI 错序 | 分离本地与 MPI 阶段；本地 catch 转 Status，在同一边界按最低失败 rank 汇总再继续/退出，不发布半成品。普通/IBM × 1/2/4 全分配扫描通过；不宣称任意 MPI 故障可恢复 |
| `compile_graph`、`core_execution_graph::allocate_workspace` | 分配失败穿过错误的 noexcept，abort | 移除私有分配函数不成立的 noexcept，交给外层本地异常转换；随后先汇总再进入 collective |
| `FieldRegistry` 复制 | 当前 Clang15/libc++ 上，schema/BoundaryCompiler 定点复制失败后留下一个跟踪分配 | 先完整构造容器成员，再复制赋值，让成员析构负责部分分配。全扫描结束时跟踪 C++ 对象/通信器/持久请求归零；不是对所有标准库版本的结论 |
| `physics_input` 指纹及 thermo/transport compile | bad_alloc 被吞成指纹0，再报 invalid_plan/801 | 增加返回 Status 的内部指纹接口，原哈希算法不变；传播 allocation_failure。其他旧标量指纹调用者仍需审计 |
| `ProductDriver::create` / `FinalForceCache::bind` | 对象分配失败，其他 rank 进入 Comm_dup，挂起 | 创建、事务和本地绑定分阶段汇总；ForceCache 在 Comm_dup 前一致检查对象分配 |
| `IbmPressureOperator::bind_internal` | solid_cells 扩容失败穿过 noexcept，abort | 内部捕获；候选 RAII 所有权，完全成功后发布。IBM create 全分配扫描通过 |
| `visit_domain_clipped_quadrature` | donor 收集回调 bad_alloc 穿过模板 noexcept，abort | 移除私有遍历模板错误的 noexcept，原 SurfaceQuadratureCompiler 本地 catch/汇总接管 |
| `RemoteDonorExchangePlan::analyze` | needs/目录分配失败，其余 rank 进入 Alltoall/v，MPI_ERR_OTHER | 请求目录、通信计数/缓冲、候选发布分别本地构造并汇总；MPI 调用在阶段之外 |
| Fresh 初始化投影构造 | 对象、union-find、标签和交换缓冲分配失败后挂起 | 分阶段分配并汇总，通信前预留表面缓冲。初版按所有方向预留破坏了表面积容量合同，整库回归发现后改为仅考虑实际通信方向；保留原断言，1/2/4 和 sanitizer 通过 |
| `NativeCartesianMgPlan::compile` | IBM initialize 第3个分配位置失败，build_levels 跳过 rendezvous，其他 rank 进入粗层通信 | 层级本地构造后立即 catch+consensus；本地层级认证后也先汇总再进入 replicated operator。IBM initialize 全扫描通过 |
| `register_fields` 的 mg_arena | 内部已含 ghost/对齐，外层线性 bundle 再加三维 ghost，真实 owner 多分配 | 仅把外层 ghost 改为0；内部每层 ghost=1、stride/offset/halo 不变。实际 owner 字节与全 level/slot 区间回归通过 |
| App / Driver 公共报告 | 复用残留旧结果；早期 phase 不准；重试历史可能被误读为最终失败 | 公共入口初始化报告，补早期 phase；共享 StepCompletionReport 分开 outcome、first/last failure、stop_reason；保留并说明旧 Driver failure 语义。输出失败仍保留已提交步和时间 |
| app_init_case 旧断言 | 要求早期失败保留旧 product，与本轮要求冲突 | 仅替换该合同，检查清零、input phase、无虚构尝试、Status一致；CaseValidationReport 的独立保留旧结果合同不变，不删除测试或放宽数值断言 |
| 候选 loop 工作记录 | 实际14次基线、322次梯级评估，未接线时报告为0；其中280次中途失败 | 新增 baseline/extrapolation/ladder/incomplete/拒绝计数及本地累计评估时间；每个成功 refinement 和最后失败 loop 保留成本，无新增逐候选 collective |
| 外推后回退写 candidates[0] | 控制流确认梯级重新写 slot0，原 sample_count 不等于总评估数 | 独立保存外推 alpha、样本、评估/选择状态；保留并说明 legacy slots 的语义。已测接受外推和大量中途失败；尚未动态覆盖“外推拒绝后梯级接受”的特定序列 |
| App / thin runner collective 子合计 | PMPI 确认普通 halo 每组 begin/finish 为4次 Allreduce，原合计遗漏；donor preflight 同样遗漏 | 汇入 halo 控制计数，新增 donor 控制计数；两个入口同步更新。仍有未归类直接调用，明确是 tracked subtotal，不冒称完整 MPI 总数 |

### MG 容量与重启身份

16³ 公共产品夹具，2 ranks，每 rank：

| 量 | 修复前 | 修复后 |
| --- | ---: | ---: |
| 内部 MG bundle / doubles | 26,880 | 26,880 |
| 外层 span / doubles | 242,064 | 26,880 |
| 包含该字段的真实 arena owner / bytes | 8,401,920 | 6,680,448 |

该 owner 少申请1,721,472 bytes，约20.5%。这是一个分配对象的容量收益，不是整个进程 RSS 降幅或求解提速。
回归逐一核对每层4个工作槽的低/高 ghost 地址、不重叠区间，并执行初始化和推进。

保留的 main 原二进制 SHA-256：
`067f928d9919832f8afb2c707361e69064e8f416a267d7bbfbde05d694fe9e0d`。
它生成的小算例旧 checkpoint 被新程序在 restart_load 阶段拒绝：invalid_plan/10307，0次推进。
新程序自身 checkpoint 可继续到下一步；现有分段/连续、BDF历史和重启回归另有覆盖。
**没有迁移/覆盖用户 Re3900 checkpoint，没有伪造旧 schema/fingerprint。**

### 已核查但不是本轮新修复

- main 的对流限制已取 accepted_n 密度与 committed final flux，固定 dt 直接跳过扫描，不能重复声称仍完全使用默认对流常数。
- α=0 的生产闭合复用、选中候选 scratch/证书直接复用、接受后退出、stationary/no-op 接受均已存在，本轮没有重复实现。
- FGMRES 已批量归约，add_scaled_basis 已融合更新；逐 basis 的 local_dot 仍是实验对象。
- MG refill 与 active/inactive 操作真实存在，hierarchy_rebuilds 少不代表数值更新便宜。
- Visit 单缓冲、输出阶段异常汇总、节点内通信器复用已存在，不能把更早快照的问题当成本 main 的新缺陷。

## C. 测试与性能证据

从活动仓库根目录执行，工具链及选项以证据包的构建缓存为准：

```sh
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu
cmake --build build-re3900-fix-test -j 6
ctest --test-dir build-re3900-fix-test --output-on-failure -j 16

# 0/-1 为起点/自动测定终点；可换 create/initialize 与 cartesian。
timeout -k 3s 1800s mpirun -n 2 \
  build-re3900-fix-test/versions/v0.4/tests/v04_product_allocation_failure_mpi_test \
  0 -1 compile immersed
ctest --test-dir build-re3900-fix-test -V -R '^v04_product_memory_profile_'
ctest --test-dir build-re3900-fix-test -V -R '^v04_product_pressure_energy_retry_mpi_'
ctest --test-dir build-re3900-fix-test -V -R '^v04_parallel_halo_mpi_'

ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
ctest --test-dir build-audit-sanitize --output-on-failure -R '^v04_product_mg_capacity_mpi_'
```

所有扫描有超时。早期一次外层 timeout 未及时回收 mpirun，之后使用带 kill-after 的命令，仅清理本任务创建的进程。
检查实际注入最低 rank、返回 code/detail、对象/原计划身份、析构后的 C++ 分配、Comm_dup/split_type 及持久请求。
**不是 MPI 内部 malloc 全审计，也不模拟 rank 死亡恢复。**

| 日志族（包内文件名前缀均为 hundun-review-） | 含义 |
| --- | --- |
| 69d8eee-* | 本轮 main 基线重新构建/测试 |
| mg-red / mg-green / mg-bounds-green | 容量 RED/GREEN 和逐层区间 |
| product-alloc-* / product-ibm-* / product-initialize-* | 定点超时、abort、误分类、GDB栈及后续通过记录 |
| first-batch-ctest* | 18组普通/IBM compile/create/initialize × 1/2/4 全扫描通过；初次整库五个新增失败随后修正，不隐藏 |
| candidate-work-red* / candidate-*green* | 计数未接线 RED；梯级/中途失败/外推观测 GREEN |
| memory-pmpi-* / halo-accounting-* | 唯一分配阶段峰值及 PMPI 对账，非性能跑分 |
| sanitizer-* | 首轮26/27，MG convergence超时；随后单独运行通过，无 sanitizer 报错 |
| final-* | 最终重建及测试；实际通过数、退出状态、二进制哈希以机器结果为准 |

### 通信对账

零初速小夹具，rank0 实际调用数；核对自身的调用不进入被测 advance。第二步与首步的 BDF 路径不同。

| 夹具 / ranks / 步 | PMPI实际 | 原有子合计 | 补halo | 补donor | 修正子合计 | 待归类 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cartesian 8³ / 2 / 1 | 476 | 283 | 92 | 0 | 375 | 101 |
| Cartesian 8³ / 2 / 2 | 498 | 307 | 92 | 0 | 399 | 99 |
| IBM 16³ / 2 / 1 | 597 | 380 | 92 | 17 | 489 | 108 |
| IBM 16³ / 2 / 2 | 637 | 422 | 92 | 17 | 531 | 106 |

PMPI 观测 Allreduce、Bcast、Allgather/v、Alltoall/v、Barrier、Comm_dup；不是消息字节或耗时。
首步包含冷 MG 建立，不能把首步与后续步差异当成纯热路径收益。差额仍需沿 time/coupler 直接调用拆分，没有塞入其他模块冒充对账完成。
thin runner 只同步计数逻辑并构建/自检，本轮没有执行新的目标 Re3900 性能窗口。

### 内存：已计量的边界与缺口

product_memory_profile 跟踪同时活跃的唯一 C++ 分配，分开 compile/create/initialize、两步advance、
Visit、restart write/read 和析构。cpp_phase_peak_bytes 是各阶段峰值，阶段之间不相加；
cpp_live_bytes 已含 arena，不能再加一次 arena；只读 view 不形成第二份场数据。

本轮产品的对齐分配有两个查明的所有者：arena 与 FaceFluxStorage。
后者为3份 final/history 和6份 workspace，回归验证
`observed_aligned_bytes == arena_bytes + face_flux_bytes`，不把所有对齐分配统称 arena。

| 唯一所有者 / 活跃期 | 理论及实测口径 | 剩余缺口 |
| --- | --- | --- |
| StateLayers / 全程 | ArenaLayout::total_doubles()*8，包括状态副本和 Krylov/MG bundle；MG owner 独立核对地址/字节 | 不另加只读 level view 或 Krylov 向量 |
| FaceFluxStorage / 全程 | 每副本8*(align8(nx+1)*ny*nz + align8(nx)*(ny+1)*nz + align8(nx)*ny*(nz+1))，共9份 | 与对齐分配总量对账 |
| pressure/energy 面数组、分支及 compiled 缓存 / 常驻 | vector.capacity()*sizeof(T)，已在唯一分配总量 | 尚需逐所有者标签拆分 |
| MG active/inactive coefficients 与 replicated coarse | 位于 arena 外，已在C++峰值/存活量 | 独立系数容量、refill/复制时间尚未导出 |
| halo/donor 发送接收/目录、几何/IBM | C++分配已计入，请求/通信器跟踪生命周期 | MPI内部缓冲不在C++ new观察范围 |
| Fresh / restart / I/O临时 | 每阶段重置峰值基准，同时保留常驻量，不相加互斥峰值 | Evidence编码、实际大网格输出/迁移重启仍需覆盖 |
| 运行时余量 | 单列当前statm采样和ru_maxrss历史高水位 | allocator碎片、libc/MPI malloc、栈、共享库、输入模型、观测器静态表不能用solver请求字节替代 |

本轮获得可复用的唯一分配/阶段计量，**尚未得到Re3900全类别、全节点同步峰值预算**。
ru_maxrss 求和只是各rank不同历史高水位的聚合；析构后跟踪量归零或短程高水位稳定都不能证明长期无泄漏。

### 正式性能列表

| 程序/路径 | 本轮相同物理窗口结果 |
| --- | --- |
| HUNDUN PISO | 未进行目标规模正式对照 |
| HUNDUN SIMPLE | 未进行目标规模正式对照 |
| COAST普通可压缩 | 保留比较项，本轮未运行 |

不从小型故障夹具外推Re3900吞吐量。历史stage50约86.8%不是本提交实测，也非完整墙钟占比。
后续每个配置只测一轮；先证明起始状态/BDF历史、网格、物性、边界、dt规则、物理窗口相同，
说明重启缓存重置及预热差异。旧checkpoint不兼容时不能绕过检查凑成“同一起点”。

## D. 复查材料

证据目录：[2026-09-05-main-review-evidence](2026-09-05-main-review-evidence/)。
results.json 为机器结果；manifest.json 列出原始日志名/大小/SHA-256及对应压缩文件哈希；
SHA256SUMS 校验包内文件。完整保留RED/GREEN/超时日志，不只给本机路径。
implementation.patch.gz 加基线 SHA 标识受测实现，构建缓存、命令和二进制哈希也在包内。
最终Git提交的文档/证据变更不反向改写已记录的历史二进制身份。
最终提交SHA随交付回复给出；包内HEAD字段是打包时的基线HEAD，不冒充最终提交。

## E. 下一轮具体复查项与设计方案

1. **剩余构造组合。** ProductCompiler的不同boundary/species/非均匀分解、donor bind前置失败、
   initialize_restart的全分配扫描尚未穷举。补多组分/长字符串、空局部流体域、已有out对象，
   检查失败后旧对象/已提交状态不变。通用Status还不携带完整失败rank/申请字节，需要错误上下文；
   不能把测试注入rank冒充生产报告提供的rank。

2. **通信与墙钟。** 上述夹具仍有99–108次/步未归类，按ordinary/prepared/donor和直接time/coupler调用所有者拆分，
   不删除安全检查。App当前计时含time_control+advance，thin runner仅advance，两者都未包含完整后续输出。
   增加本地资源/Evidence/场/checkpoint/整步计时，低频汇总；各阶段rank最大值之和不当成整步rank最大值，
   Evidence写完自身的时间放后续记录或独立摘要。

3. **候选与内核实验。** 尚未拆Schur准备、线性解、物性、边界/通量、残差、哈希、复制时间。
   补“拒绝外推→接受梯级”的公共夹具。solver_krylov的分块multidot保留补偿求和、有限性、真实残差/正交化恢复。
   MG active/inactive减复制须先证明view/halo/证书地址及生命周期合同。
   单次只改一个热点、单轮受控测量，不增restart、不换PCG、不启用fast-math、不放宽物理门槛。
   提前启用IBM空间后备只能用全局一致收缩率与剩余预算，保持现有容量；本轮未改变固定切换规则。

4. **I/O错误上下文。** io_output_detail的write/sync仍返回bool，EINTR及fsync/close errno可能丢失。
   保存首次syscall/errno/path/rank，再经现有阶段边界汇总，测试单rank open/write/fsync失败和短写。
   不盲目重试close，不改变checkpoint持久化及一致提交。普通日志批量刷新与checkpoint分开验证。

5. **物理时间尺度。** 保留accepted rho/committed flux的对流适配器；其余LocalTimeLimits默认1.0不是实际物理尺度。
   黏性/热/组分尺度逐一关联预测器和校正的离散贡献、时间层和显隐式路由，后置隐式不自动免除前置显式限制。
   用速度/网格/物性缩放和1/2/4一致性测试后再接入；低马赫隐式路由不无条件加声学上限。
   外部limits调用者及fixed-dt合同不变。先记录重试/BDF降阶恢复，再决定增长冷却，不猜经验参数。

6. **显式初场。** App仍遍历boundary选最后有效温度/入口速度。
   提议独立initial-state：速度、压力参考、温度或焓（同给需一致）、按字段身份的组成/被动标量。
   用同一ThermodynamicsPlan闭合rho/T/h；独立组分全零可能表示依赖组分为1，不能直接判非法。
   默认规则版本化、确定且不受无关边界改变影响，单独说明迁移。
   测边界排列/无关面温度、显式T/h一致性、多组分与重启；本轮未改生产初场。

7. **线性化与缓存。** Cartesian对冻结物性/通量/limiter分支的同一残差做方向差分与步幅扫描，
   不用全量重新闭合残差去验收冻结线性化。IBM按PressureEnergyJacobianScope验证承诺块及候选收缩，
   不冒充完整Jacobian。切换处做分支一致单侧测试和独立跨分支可接受性测试。
   缓存失效覆盖BDF/time generation、state/flux revision、物性、limiter、backflow、IBM stencil、workspace身份。

8. **发布与新政策。** 历史V1.0 policy仍要求**456×256×104 / Evidence V7**；默认writer为**V8**，
   历史当前诊断算例为**456×256×52**，不能互换。原政策/附件哈希不变。
   下次正式统计前另建版本，固定变更范围、共同物理窗口、状态/重启迁移验证、计数定义和门槛，先登记再运行。
   历史receipt首先因旧绝对仓库路径失效而失败。本轮只读地把旧仓库路径映射到活动工作区再核对，
   14个条目均找到文件，其中12个SHA相符，2个提取脚本的当前内容与历史SHA不符，未把这项测试改成通过。
   映射只用于审计，不改原receipt；恢复历史工具版本与建立新凭据需分别处理，不能直接更新旧哈希。

本轮没有可宣称速度收益的FGMRES/MG/候选算法实验，没有放宽标准。
未完成项保留为具体设计/测试工作，不作全部正确、彻底稳定或长期无泄漏结论。

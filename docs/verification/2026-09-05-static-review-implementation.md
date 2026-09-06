# 静态审阅后的实现与 Re3900 长测准备

用户已授权任务包 A–E，以及验证完成后的 Re3900 长测。基线为
`b779bff4bbff7067e691952b25c90cbc86be7795`。本文件随实施更新；未完成项
不作为验收结论。原始 checkpoint、receipt、附件及旧二进制不改写。

## 实施清单

- [x] A1：Writer、启动身份及 runner 的本地异常/collective 边界。
- [x] A2：Restart POSIX 首错、发布/持久化/清理状态。
- [x] A3：数值代次保留、严格清理范围。
- [x] B0：SIMPLE C1 计时归属。
- [x] B：按求解类型、corrector/refinement/attempt 细分成本，选择实测热点。
- [x] C：字段绑定和候选 workspace 的最小结构整理及失效合同。
- [x] D：Restart 编码峰值/预算与 cell/face 恢复预筛选。
- [x] E1：receipt 非空参考校验，保留旧工具和 R4/R5 历史身份。
- [x] E2：按目标记录 C/C++、runner 构建身份。
- [x] E3：借用快照生命周期合同。
- [x] E4：显式初场与时间尺度设计（不加入未经论证的公式）。
- [x] 公共接口重点回归、sanitizer 检查及 Re3900 短窗。
- [x] 独立运行目录启动 Re3900 长测，记录 PID、命令、二进制和输入身份。

## 不变条件

D=0.02 m，Uc=2.89668 m/s，Re=3900；456×256×52，
20D×10D×(π/2)D；保留 COAST 网格、变热力学/输运物性、边界、
固定 dt、完整 BDF 历史及数值验收。每配置一轮，不设耗时门禁，
不以更换物性、删减检查或降低持久化要求提速。不推送 GitHub。

## 本轮证据

实施开始时源码与基线一致；工作区已有 AGENTS.md 和 .codegraphf/ 改动，
保留不动。Release 和 sanitizer 使用独立构建目录，避免覆盖旧二进制。
已完成的独立回归（以下是本轮执行结果，不引用旧轮次通过数）：

| 项目 | RED | GREEN / 状态 |
| --- | --- | --- |
| Writer 真实分配失败 | 2 ranks，index=0，SIGABRT/134 | 发布/清理报告扩展后，1/2/4-rank 全扫描通过（26.20/38.06/44.13 s）；最终编码版随全套回归再验证 |
| 启动身份 | ApplicationService，2 ranks，每 rank 2 次分配 | 同入口每 rank 0 次分配，正常接受一步；哈希输入字节顺序不变 |
| 保留顺序 | keep=2 留下 99/101；keep=3 留下 98/99/101 | 4-rank Writer 测试保留 100/101 及 99/100/101 |
| POSIX | open/write/fsync EINTR 失败；正短写原本通过 | 21 组通过，含 errno 首错、零写和 current 切换后的目录 fsync EIO |
| 发布报告 | after_current_switch 被报告为默认未切换 | 已提供可见/持久化状态，清理另报 warning；后续回归继续核验 |
| receipt 空集 | CLI 错误接受 complete | 14 项 CLI 检查通过；R5 显式历史校验通过，仍 complete=false |
| runner 实际 CLI | 旧 runner（同一新 core，用于隔离外层控制流）root 第一个分配点 SIGABRT；新阶段汇总版 index=2 出现截断文件名但仍完成 | backtrace 确认 stringbuf/ostringstream 静默 badbit；启用编码流异常后 2-rank 117+36 点全扫描通过；不是放松测试期望 |
| 恢复预筛选 | 4→4 时每个 rank 解码全部来源 | V1/V2 及 1→2、2→4、4→1、4→4 恢复通过；31 组 Writer/Reader POSIX 检查、精确预算边界、全部保留代次 Reader 检查通过 |
| SIMPLE C1 计时 | 实际接受的 SIMPLE+IBM 步，C1 invoked=true，但 Krylov 阶段时间为零 | C1 路径补齐 phase(1)/phase(2)，同测试通过；新增局部细分和跨重试总成本 |

记录暂保存在独立的 `build-review-fixes/`，结束时整理可复查材料。
最初 Python 3.11 不在 PATH，以及未设置 libc++ 运行库路径导致的退出 127，
属于测试启动配置问题，不计为缺陷复现。实际 Python 为现有 3.6；MPI 测试使用：

```sh
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu:/home/wyf/.local/opt/hundun-toolchain/clang/lib
```

## E：证据与公共合同

原 R5 SHA-256 保持
`0dae4f3d4e9bd6ec244f0197bb2b5ecd0656cc2f3b359a3346b71340ee1e87cf`。
其绑定脚本保存为只读历史数据
`historical-extractors/4f0debb1fcf6ca0b062f28cd700f8f22d0d1a0d5dde7f76cc6d723cb0c34368a.snapshot`，
内容 SHA 与文件名一致。校验不执行该脚本，也不更改 R4/R5：

```sh
python3 tools/v04_literature_extract.py receipt-validate \
  docs/verification/v0.4-literature-data-receipt-r5-partial.json \
  --historical --artifact-store docs/verification/historical-extractors --json
```

快照合同已在公共声明说明：同步借用，不跨 advance/同类再次取快照/销毁保留；
move 构造转移所有权不移动底层存储，move 赋值会释放目标的旧存储。
当前未引入租约或异步复制，不能宣称 Writer 会自动识别所有过期借用。

### F10 设计边界

本任务包授权的是设计评估，不在此次固定 dt 重启算例中改变初始化或增加时间步公式。
后续应用配置应拥有显式的温度、压力参考、速度及具名独立标量，再适配既有
`DriverInitialState`；入口与初态分别指定。兼容默认规则必须具名、确定，并对
多入口歧义给出诊断，而不是继续隐式按边界覆盖。独立组分全零可代表依赖组分为一，
不作一律拒绝。

每项时间尺度需要登记离散算子、接受状态/系数来源、显式或隐式处理及其用途。
保留现有真实对流限制和接受后 CFL 审计。扩散隐式项不能一律套用显式稳定步长；
声学限制也不能仅凭流动模型名字补入。重试冷却须有 dt/残差轨迹支持才实施。
验收通过 ApplicationService/ProductDriver 覆盖初态与无关边界解耦、组分/焓闭合、
缩放、1/2/4-rank 提案、重试及 BDF 降阶/恢复。目前这些新默认行为尚未实现。

## C：本轮采用的结构边界

本轮只提取标量的类型化绑定描述，未重写 ProductDriver，也未引入新的调度器。
`ScalarBinding` 在注册时冻结 FieldId、role、role_index、catalog_slot；普通 halo 和
IBM rate donor 的打包/取回共用它。FieldId 和 schema 顺序不变。编译阶段原有
arena、ghost reach、图资源和候选存储别名检查继续负责容量与合法性。

候选工作区保留既有 `PressureEnergyCandidateArtifacts` / LoopResult 的边界：

| 权威/数据 | 合法作用域与失效条件 | 此次处理 |
| --- | --- | --- |
| rank-invariant 方程、物性、边界语义 | 换物性、边界、离散或算子类型后重编译 | 不改语义哈希，不用存储地址替代 |
| rank-local 字段 view、halo 槽位 | runtime revision、storage identity、ghost reach 必须与冻结绑定相符 | 角色/槽位共用描述；不跨工作区冒用字段身份 |
| 候选 thermo/flux/energy 证书 | state/flux revision、dt/BDF/time generation、BC/backflow、limiter 分支、IBM stencil 或 workspace epoch 改变后重建 | 复用既有证书验证和候选 replay；未建立跨候选新缓存 |
| 压力方向与热力学候选 | 下一次候选写入会失效；拒绝后清理待提交权限和 warm-start bit | 保留事务独占提交与回退，不用日志对象发布场 |
| borrowed 输出快照 | advance、同类取快照、move 赋值目标释放及销毁 | 公共接口明确寿命；不新增全场复制 |

进一步合并整个压力—能量 workspace、自动派生全部阶段清单属于后续结构工作。
目前没有证据要求为当前长测新增别名复用或改缓存失效规则。

## D：内存与恢复观测口径

单 rank 小网格相同 memory-profile fixture：Writer 的 C++ 活跃增量峰值由
208,312 B 降为 77,840 B（约 -62.6%）。这不是整进程 RSS 的同幅度下降。
Reader 在该 1-rank fixture 无预筛选收益，临时峰值基本不变；不能冒充重分区收益。
Writer 预估 V1/V2 确切长度，预留后 move；释放已写 rank 载荷，再在 root 复用一个
验证缓冲。生产 app/runner 传入单独的 V2 staging 预算，包含两层主变量、速率、
面通量以及 root gather/manifest/verification 的 bulk 数组；不重复计借用字段。
`peak_bulk_staging_bytes` 不含路径字符串、分配器元数据、MPI 和求解器常驻内存，
因此不称为完整 RSS 上限。Reader 新增逻辑读字节与解码块数，保留全档校验；
读取缓存不等于实际磁盘流量。更完整的 Reader 总内存硬预算仍需独立 owned-image 合同。

## B：单轮算例观测

独立目录 `review-kernel-observe-746to749`，128 ranks、固定原输入和 checkpoint，
3 步一次运行，67.28 s 全进程墙钟；平均 advance 17.0441 s。
每步 C1 pressure=1、双对角 Schur=7、空间 Schur=4。
各 rank-step 均值：prepare 0.6210 s、Krylov 9.1156 s、recovery 0.0311 s，
candidate 3.7376 s。Krylov 内 A apply 3.8905 s、M apply 3.6184 s，
主 Arnoldi dot/reduce/update 为 0.2613/0.5543/0.0988 s。
A/M 为含通信成本；这些分量不与 PMPI 时间相加，也不相加独立 rank 最大值。
本轮只尝试保留每个点积补偿求和顺序的 x-strip multidot，未放宽容差或更换物性。
正确性已通过逐位对照和非有限值检查；目标算例收益尚待最终单轮验证。

构建验证中的混合 GCC 7 C + Clang 15 C++ sanitizer 链接缺少旧 UBSan 符号，
属于检查构建工具链不匹配；改用同一 Clang C/C++ 的独立目录，不改 Release 物理配置。
identity fixture 最初同秒文件修改未触发 Makefile 重配置，测试显式跨越时间戳粒度后
验证 C-only/runner-only 的实际增量更新；该启动问题不计为源代码 RED。

## 最终回归进度

- Release 全套：287/287，561.43 s，包括真实分配故障扫描及 1/2/4-rank 回归。
- 同一 Clang C/C++ 的 ASan/UBSan：11/11，111.12 s，涵盖 Writer/Reader、
  重启存储迁移、压力—能量重试、Krylov 和逐位 multidot 对照。另运行一次
  2-rank SIMPLE+IBM refinement sanitizer 检查通过。
- 最后仅构建身份材料更新后，身份/runner/性能解析公共入口：9/9，10.98 s。
- 最终 runner 的固定长度路径 CLI 全扫描：2 ranks，root 131 点、另一 rank
  36 点，全部通过。此前未固定长度的夹具在 index=129 因路径跨过 small-string
  阈值而未触发分配；该记录不算程序失败，也不计作成功注入。
- sanitizer 使用 `ASAN_OPTIONS=detect_leaks=0`（MPI 运行时残留不由 LSan 验收），
  `UBSAN_OPTIONS=halt_on_error=1`。不能将这组结果称为泄漏检查通过。
- 程序仍在带有未提交源码改动的工作区；V2 manifest 记录 source content hash、
  C/C++ 配置、目标入口 SHA 和 target_source_clean=false，另由 executable SHA 绑定
  实际二进制。不把 b779bff 的 HEAD 误当成本次未提交程序的完整源码身份。

## 修改后单轮与长测提交

`review-multidot-746to749` 使用冻结程序及逐字节相同的输入，从同一 746 checkpoint
推进三步。整个进程 64.85 s，平均 advance 16.8788 s/步。主 Arnoldi dot
0.22457 s（前值 0.26127 s，约 -14.0%）；全步约 -0.97%。这是一次配对诊断，
不为 1% 总体差异建立统计显著性结论。保留逐点补偿求和与有限性检查，未调整
求解容差、物性、网格、dt 或任何接受标准。

更重要的是：两次 749 checkpoint 的 **128 个 rank 文件加 manifest 共 129 个文件
逐字节相同**，health.csv 除 `max_rank_step_ns` 外所有值相同。
三步均一次接受，BDF2 无降阶；连续性最大 9.1953e-7、能量最大 2.4982e-7，
EOS 误差 0，均满足当前 1e-6 门限。T 为 299.9828–300.0393 K，压力为
107114.31–107142.32 Pa，密度保持正值。最高 rank RSS 最后两步均 235,319,296 B；
三步窗口不是长期内存或湍流统计收敛证明。

长测独立目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/long-review-749to35000-20260905/`

启动脚本与比较清单在 `2026-09-05-review-evidence/`。从修改后最新有效的 749 checkpoint
推进 34,251 步至 35,000，128 ranks，core 绑定；每 500 步 Visit/checkpoint，checkpoint
保留两代，所有路径独立。估计纯推进约 6.69 天，实际随迭代数、系统负载及输出变化。
此处未设置耗时门禁；初始化和每一步仍执行原有数值检查及重启完整性校验。

### 启动后核验（长测尚未完成）

2026-09-05 23:38:10 +08:00 已启动上述长测；外层 time PID 为 10097，
mpirun PID 为 10099。进程检查确认其下有 128 个非 zombie rank 进程。
`RUN.meta` 确认从 749 步恢复，未迁移存储布局，也不需要重启恢复降阶；
运行时 executable SHA 与冻结清单一致。源码仍为本地未提交修改，
基线 HEAD 不能单独代表此二进制，完整源码 archive 和逐文件哈希共同留档。

截至 2026-09-05 23:48:42 +08:00，完整 Evidence 行确认第 750–786 步
连续接受，共 37 步；全部保持 BDF2，无 retry、restart recovery 或 temporal fallback。
逐行核查最终物理审计均通过：连续性最大 9.899950683e-7，能量最大
9.541220803e-7（各自阈值 1e-6），提交态对流 CFL 最大 0.244885148（阈值 0.8）。
这段平均推进耗时 16.7375 s/步，最高 rank RSS 为 235,352,064 B。
这些仅是启动段实测，不构成长时间稳定性、内存增长或统计收敛结论。

机器可读启动快照见 `2026-09-05-review-evidence/launch-status.json`；
持续运行记录以独立运行目录中的 `evidence.jsonl`、`health.csv`、
`performance.csv` 和后续 checkpoint 为准。CSV 存在缓冲，读取时不将
尚未写完的末行当成完整步，也不据此判断 Evidence 已确认的步丢失。

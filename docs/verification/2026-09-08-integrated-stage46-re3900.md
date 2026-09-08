# 主线与 Stage 4–6 模块整合及 Re3900 续算

> 状态修正：用户随后明确要求完成产品耦合接线才算合并完成。本记录只证明源码归集及单相兼容性，不满足该验收条件。新长测已再次 SIGSTOP（完整健康前缀至 10577，独立恢复源为 10510）；产品接线和耦合验收完成前禁止恢复圆柱计算。既有数值证据保留，不将其改写为耦合证据。

按用户本轮要求，将 `codex/stage4-6-coast-replacement@f3ad77b` 与已验收
主线 `25b3ea9` 合并，并以新构建接替圆柱长测。整合源码提交
`f4f04eccc4f159cf7aa55e4244cb36b7e642cb77`，两个父提交均保留；
本文和证据归档是随后文档变更，运行程序仍明确绑定该源码提交。

## 整合范围与构建

开发工作树为 `/home/wyf/code_dev/.worktrees/hundun-flow-integrated-stage4-6-re3900`，
分支 `codex/integrated-stage4-6-re3900`。原主线和供体开工时均干净；
供体保持只读，整合通过独立 worktree 完成。
共同祖先为 `86542bb96678ec844ae5ac95d8f6391993da239e`。
唯一重叠文件 `versions/v0.4/tests/CMakeLists.txt` 自动合并成功，其余文件分别
与对应分支逐文件比对一致。主线流动 driver、求解器、Restart 和 runner 保留。

[P1–P8 供体台账](../handoff/2026-09-08-later-stages-portable.md)所列普通值模块
已进入同一源码版本和 core 构建，包括 chemistry/PaSR/ESF/TCR、液滴、TAB、
事件、候选迁移及交换组合器。它们仍采用候选接口；本次没有新增产品
chemistry source admission、反应/喷雾 driver 接线或磁盘持久化格式。
圆柱继续使用原单相变物性配置。

全新 `build-integrated` Ninja Release 构建了 `hundun`、
`v04_thin_domain_runner` 和相关测试。采用与长测基线相同的 Clang 15/libc++、
`-O3 -march=znver3 -mno-fma -ffp-contract=off`；HYPRE、sanitizer 关闭。
构建 manifest 的 core/target clean 均为 true，head/tree 的具名前缀哈希独立重算一致。
Cantera 仍是单独的 GCC/libstdc++ ABI 目标；本次 Clang 产品构建关闭该选项，
未重跑 Cantera conformance，也不把供体既有 conformance 记成本轮检查。

| 检查 | 本轮结果 |
| --- | --- |
| 可移植模块及公开头 | 32/32 |
| V1 单物理 | 23/23 |
| V2 组合 | 25/25；六组 1/2/4-rank 已接受状态指纹一致 |
| 既有流动、Krylov、PISO、方法历史、runner、日志关闭和 V3–V6 observer | 31/31 |
| 构建布局与合并范围 | current-source fixture、源码空白检查通过 |

共 111 个唯一检查通过。第一轮有 6 个测试因构建目标清单漏列可执行文件而
Not Run，105 个已执行项通过；随后补建并仅执行这 6 项，全部通过。
[首轮 XML](data/2026-09-08-integrated-stage46/tests.xml)及
[补跑 XML](data/2026-09-08-integrated-stage46/tests-missing.xml)均原样保留。
这是工具清单遗漏，不是产品断言失败；没有隐藏首次 CTest 非零退出。

## 暂停、保存与同恢复点对照

基准根目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52`。
完整原始证据位于其下 `integrated-stage46-20260908`。

旧 `hundun-re3900-observed-long-20260908.service` 的全部 128 ranks 经 SIGSTOP
验证暂停，完整健康记录至 10601；最新持久点为 10500。
`generation-10500-176021483306485` 的 128 个 rank 文件、manifest、current、
complete、statistics 和 accumulator 共 133 文件独立复制至
`integrated-stage46-20260908/source-10500`，逐文件 SHA-256 一致并同步。
10500 的 128 个 Visit 文件及索引齐全。统计 epoch=7000，development=10000，
sampling_start=17001，样本数仍为零，V3 方法历史签名为 `12213963202598979269`。

确认副本、构建及回归后，旧 unit 被管理性停止；MainPID=0、inactive，旧 MPI
进程全部退出。原输出保留为中止前缀，不能当作正常流关闭或 COMPLETED 验收。
10501–10601 的原内存尾段由新作业从最近持久点重算，不伪装为从 10601 精确恢复。

| 单轮配置 | 结果 |
| --- | --- |
| 新版本 `pilot-integrated-10500-10510-20260908` | 128 ranks，10500→10510，106.13 s，COMPLETED |
| 旧冻结程序 `control-mainline-10500-10510-20260908` | 同一个只读 10500 源，128 ranks，111.18 s，COMPLETED |
| 新版本健康 | 10 步均 BDF2、attempt=1、无 retry/方法恢复；最大连续性 `2.5462782876084995e-7`、能量 `9.728335760201628e-7`，均满足原 `1e-6` 门槛 |
| 持久载荷对照 | 260 文件逐字节一致：128 rank、manifest、statistics、accumulator、128 VTR、Visit 索引；complete 仅 generation 名不同 |
| CSV 对照 | 10 行 health、10 行 force、30 行 probe、10 行 conservation 的物理列及稳定 revision 相同；计时单列 |
| 观测与证据 | 新旧 V6 observer 均 complete，10 步、128 ranks；新版本 63 loops。两组 runtime validator 返回 0 |
| VTK | 实际 appended binary blocks、分区覆盖/无重叠、6,070,272 cells、有限值与坐标递增检查通过 |

每个配置只执行一次，串行使用互斥计算锁。两个计时不是统计性能验收，不据此
宣称加速。原连续轨迹与新重启窗口有微小数值差异；同恢复点对照则载荷完全一致。
现有恢复实现会清除瞬态压力修正 warm-start 权限（`core_product_freeze.cpp`），
因此连续运行与重启运行不能未经对照就假定逐位相同。

比较器首次因 force 的 `certificate_terminal_state`、`certificate_state` 不同而
拒绝。源码确认：terminal audit 包含 rank-local 字段绑定，绑定哈希包含 base
地址；force state 再混入 terminal audit（`solver_piso.cpp` 的
`mix_complete_view_identity`、`semantic_composition_rank_local_binding` 和
terminal audit；`solver_ibm_force.cpp::force_state_revision`）。这两列是进程内
权限标识，不是跨进程物理校验值。比较器仅明确排除这两列及计时，保留双方原值；
其余力证书列、实际力、探针和物理输出仍逐值比较。

辅助的“旧连续轨迹”probe 比较起初只以 step 配对，漏了 station；已修正为
(step, station)。完整记录中的 `PROBE_COMPARISON_CORRECTED.json` 明确取代
旧辅助 JSON 的 probe 部分。最终同恢复点验收使用完整有序行，不依赖该辅助结果。
固体区域极值和 VTK 有限性已检查；这不是新增 source-to-final 固体全载荷专项验收。

## 新版本长测现场

短测 10510 的 133 文件再次独立保存为 `source-10510` 并逐文件校验。
随后启动 `hundun-re3900-integrated-long-20260908.service`，MainPID=165884；
新目录为 `long-integrated-35000-20260908`。新作业从 10510 精确续算 24490 步，
目标 35000，128 ranks；原程序已退出，没有同时运行两组 128 ranks。
冻结的程序、manifest、输入及工具保存在 `integrated-stage46-20260908/frozen`。

启动只读快照完整健康记录为 10511–10522：12 步均 BDF2、无重试，最大连续性
`2.8302423439211935e-7`、能量 `9.296890369296903e-7`。
公开 Reader 实际读取 10510 源，V3 source/target 签名一致，
policy=`require_compatible`，无 method recovery，统计 epoch 沿用。
Runtime evidence 验证通过；observer 完整验证 13 步，唯一 partial issue 为
下一步 10524 尚未完整，保留 RUN.meta 的真实 35000 目标。
程序和源文件哈希回查未变。128 ranks 的 RSS 求和快照为 27,325,004 KiB，
只是当时占用，不是泄漏趋势验收。

下一次 checkpoint/Visit 在 11000。后续检查应核对新的两代恢复点、统计附件、
RSS 趋势及残差/工作量；继续运行时只读观察，不替换冻结程序。
目前通过的是合并、回归、同恢复点兼容性和长测启动验收，尚未到 35000 步，
不宣称完整燃烧/喷雾耦合、科学统计充分发展或 COAST 替代验收完成。

紧凑收据和实际脚本位于
[data/2026-09-08-integrated-stage46](data/2026-09-08-integrated-stage46/)。
脚本绑定上述冻结基准目录，不是任意目录可重放的启动器；launch 脚本拒绝覆盖
已存在输出。完整构建/测试日志、原始运行 CSV 和观察报告保留在基准证据目录。

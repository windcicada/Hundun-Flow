# ProductDriver 的 MG 成本观测接线

基线 `95eb7d9a522098b618d6acd46175b01e51e5716b`；继承已验收的
[Native MG profile](2026-09-08-mg-apply-profile.md)。本切片不改变数值方法，
尚未接入 runner CSV/逐层 sidecar，未进行新的生产性能测量。

## 接口及容量

`ProductDriver` 增加 rank-local 冷配置 `set_pressure_mg_profiling(bool)`、
借用累计记录 `pressure_mg_profile()` 和冻结层几何 `pressure_mg_level()`。
配置允许各 rank 不同；没有新 collective，不参与科学 fingerprint。
active proposal 时拒绝配置，每次成功 setter 重置观测 epoch，不修改方法状态。

层几何来自 ProductCompiler 已冻结的 `MgWorkspaceRequirements`，create 后即可读取，
不触发 lazy MG compile。尚未编译时分别报告 requested enabled、initialized=false、
cumulative=nullptr。无效对象/index 返回 invalid_plan，层输出不写；move 转移 owner。
view 中的布尔/层数是读取时的快照，不会自动随以后初始化刷新。累计值需在 advance、
reset 或初始化之前复制，不能把借用引用当成历史样本。

每个非嵌套 `solve_pressure_energy` 入口复制一份定长累计 baseline，退出时折叠差值，
沿既有 corrector/refinement、attempt、composition sweep 和实际 dt 记录保存。
不因 dt retry 清空累计记录，失败尝试亦计入。完整32层数组不放进64-loop报告。

`MgSolveProfile` 为200字节，`MgApplyProfile` 为4664字节（本次 Clang/libc++ 构建实际 sizeof）。
新增一份 Impl baseline；64-loop数组的新增字段载荷12800字节，另有既有 attempt/report
对象中的摘要副本。没有新增全场数组。关闭时跳过 baseline copy 和逐层差值计算；
固定对象/摘要的清零与报告拷贝成本仍存在，lazy compile 后也有一次 profile 初始化。
这些数值不是完整产品峰值/RSS预算，也不证明生产额外开销为零。

紧凑摘要保存 apply尝试/成功/失败、总时间、六项非递归阶段、通信嵌套成本和 finest
pre/post smoothing。它不保存任意单个 loop 的完整逐层分布；计划中的每-step sidecar
只提供整步各层成本。Fresh projection 和 MG refill/copy 不属于此累计对象，
后者继续沿用已有独立计时。

差值检查 enabled/层数身份、倒退、UINT64_MAX和汇总溢出；无法解释时 complete=false，
不改变求解 Status。默认关闭路径无新分配表达式；本切片未新增全 ProductDriver 的
C++ malloc拦截实验。Native apply 的零热分配证据仍仅按原报告范围引用。

## 实际回归

使用用户此前确认的 ProductDriver 公共接口，在现有真实重试夹具上增量验证。
MPI 全部由主任务串行调度，原128-rank长测保持SIGSTOP，未替换冻结程序。

- RED：接口声明/空实现时2 ranks失败（2.04 s），MG=0；原物理解逐位比较通过。
- 首次 GREEN：1/2/4 ranks，3/3，3.14 s。rank0单独开启；24项摘要逐项与
  Native累计差值精确对账，apply次数与公开线性报告一致。原失败回滚及后续BDF2
  物理解逐位相同。各分区初窗/累计apply为70/85、90/108、94/114，不跨分区强求相等。
- 契约扩展复用已有重试耗尽窗口，最后一个rank单独开启。两个拒绝尝试均保留，
  move、reset、disable、invalid/default接口后检查已接受完整Restart载荷及flux不变。
  不增加推进次数；这是已有实现的合同覆盖，不伪造第二次RED。
- Clang Debug ASan+UBSan retry 1/2/4：3/3，123.65 s，包含契约扩展。
  `detect_leaks=0`，不宣称LeakSanitizer或进程RSS硬预算验收。
- 扩大相关组首次12/15：MG容量、存储恢复、core产品冻结1/2/4及有符号/变物性/
  重启标量2rank通过。失败的三个Fresh夹具另行确认为
  [旧固体内部通量预期](2026-09-08-restart-inactive-flux-fixture.md)，测试修正与观测分开提交。

原始stdout及逐测试日志见[数据目录](data/2026-09-08-product-mg-profile/)。
最终Release的retry1/2/4与Fresh串行/MPI1/2/4共7/7通过（6.21 s），覆盖完整生命周期。
Fresh组另有ASan+UBSan 4/4通过（59.73 s）。连同前述相关12项，本轮影响范围内
19个不同Release测试均通过；不是全仓库测试通过的声明。
Release增量目录仍有历史Ninja末尾恢复提示，不能作为干净发布身份；本轮生产runner
已成功构建，但未运行新生产窗口。

## 独立干净验收

在 detached `d199f9681971e89b3fe1872f95f4b92f2d477b81` 的独立 checkout
`hundun-flow-mg-driver-accept-20260908` 重新 configure/build，构建前后 Git 均干净。
树为 `3068596aa0e4e5798a3b348f3b070211cc2d9cfe`。Clang 15/libc++、Release、
`-march=znver3 -mno-fma -ffp-contract=off`，测试开启，ASan/UBSan/HYPRE 关闭。
生产 runner 及相关测试 target 均构建成功；没有增量目录的 Ninja 恢复提示。

同一相关组 **19/19 通过，21.85 s**，MPI 串行调度。
[stdout](data/2026-09-08-product-mg-profile/clean-acceptance-19.txt)、
[逐测试原始日志](data/2026-09-08-product-mg-profile/clean-acceptance-LastTest.log) 和
[runner manifest](data/2026-09-08-product-mg-profile/clean-acceptance-runner-manifest.txt)
已保存。这是 ProductDriver 观测接线及 Fresh 夹具合同的验收，不包含后续 runner MG CSV
实现，不把验收标签自动转移给新二进制，也不是全仓库或生产性能验收。

- runner SHA-256：`56e9b662a849f38fce25f6dc8dab79dd1c4ecca5c0e3da688c9622a6cb1f942d`
- manifest SHA-256：`689ce95fe91ad985bc8f7825418e2f9cc2cf0532a4a35fe97cd76db484242fb7`
- core content SHA-256：`5a7e8c45b71f962d71f867b6f030ed9052c8785c2f816f7a887725a84e89104c`

manifest 中的 head/tree 是带 `hundun-git-head-v1:` / `hundun-git-tree-v1:` 前缀的
SHA-256，不是原始 Git SHA；独立重算与 manifest 一致。source clean 两项均为 true。

## 下一步边界

完成runner opt-in、紧凑loop字段及每step层级sidecar；冷阶段从冻结几何形成独立
rank/level清单和来源哈希，不能让两份不完整CSV互证完整。包含I/O关闭/失败汇总、
计数饱和/缺层/缺rank/partial校验后，再用一轮同窗口Re3900测量选择最小优化。
不提前切换spatial、不调整阈值，也不将此观测接线称为已经提速。
两相按用户最新要求暂挂，本轮未接入任何燃烧或两相源码。
